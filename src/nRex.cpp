/*
============================================================================
nRex
============================================================================
Copyright (C) 2016 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include "htslib/vcf.h"


// Config arguments
struct Config {
  bool selTranscripts;
  bool clinvarPathogenic;
  uint32_t maxCarrier;
  std::string severity;
  boost::filesystem::path vcf;
  boost::filesystem::path transcripts;
  boost::filesystem::path outfile;
};

  
template<typename TMap>
inline bool
getCSQ(std::string const& header, TMap& cmap) {
  std::string delimiters("\n");
  typedef std::vector<std::string> TStrParts;
  TStrParts lines;
  boost::split(lines, header, boost::is_any_of(delimiters));
  TStrParts::const_iterator itH = lines.begin();
  TStrParts::const_iterator itHEnd = lines.end();
  for(;itH!=itHEnd; ++itH) {
    if (itH->find("##INFO=<ID=CSQ,")==0) {
      std::string delim(",");
      TStrParts keyval;
      boost::split(keyval, *itH, boost::is_any_of(delim));
      TStrParts::const_iterator itKV = keyval.begin();
      TStrParts::const_iterator itKVEnd = keyval.end();
      for(;itKV != itKVEnd; ++itKV) {
	size_t sp = itKV->find("=");
	if (sp != std::string::npos) {
	  std::string field = itKV->substr(0, sp);
	  if (field == "Description") {
	    std::string desc = itKV->substr(sp+1, itKV->size() - sp - 2);
	    size_t colon = desc.find(":");
	    if (colon != std::string::npos) {
	      std::string format = desc.substr(colon+2);
	      TStrParts columns;
	      boost::split(columns, format, boost::is_any_of(std::string("|")));
	      TStrParts::const_iterator itC = columns.begin();
	      TStrParts::const_iterator itCEnd = columns.end();
	      int32_t i = 0;
	      for(;itC != itCEnd; ++itC, ++i) {
		cmap.insert(std::make_pair(*itC, i));
		//std::cerr << i << ": " << *itC << std::endl;
	      }
	    }
	  }
	}
      }
    }
  }
  return true;
}


template<typename TStrParts, typename TMap>
float
popMax(TStrParts const& vals, TMap const& cmap) {
  float popmax = 0;
  std::string tmp[] = {"gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "gnomADg_AF_AFR", "gnomADg_AF_AMR", "gnomADg_AF_ASJ", "gnomADg_AF_EAS", "gnomADg_AF_NFE", "gnomADg_AF_OTH"};
  typedef std::set<std::string> TAfSet;
  TAfSet afs(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
  for(typename TAfSet::const_iterator afIt = afs.begin(); afIt != afs.end(); ++afIt) {
    if (cmap.find(*afIt) == cmap.end()) std::cerr << *afIt << " does not exist!" << std::endl;
    std::string afstr(vals[cmap.find(*afIt)->second]);
    TStrParts afparts;
    boost::split(afparts, afstr, boost::is_any_of(std::string("&")));
    for(typename TStrParts::const_iterator afP = afparts.begin(); afP != afparts.end(); ++afP) {
      size_t sp = afP->find(":");
      std::string freqstr;
      if (sp != std::string::npos) {
	TStrParts freq;
	boost::split(freq, *afP, boost::is_any_of(std::string(":")));
	if (freq.size()==2) freqstr = freq[1];
      } else freqstr = *afP;
      float afval = 0;
      if (freqstr.size()) {
	if (freqstr != ".") {
	  afval = boost::lexical_cast<float>(freqstr);
	  if (afval > popmax) popmax = afval;
	  //std::cerr << *afIt << "\t" << afval << "\t" << afstr << "\t" << popmax << std::endl;
	}
      }
    }
  }
  return popmax;
}

template<typename TStrParts, typename TMap, typename TConsequence>
bool
candidateVar(TStrParts const& vals, TMap const& cmap, TConsequence const& consequence_filter) {
  bool pass = false;
  TStrParts consvals;
  boost::split(consvals, vals[cmap.find("Consequence")->second], boost::is_any_of(std::string("&")));
  for(typename TStrParts::iterator consIt = consvals.begin(); consIt != consvals.end(); ++consIt) {
    if (consequence_filter.find(*consIt) != consequence_filter.end()) {
      pass = true;
      break;
    }
  }
  return pass;
}


int main(int argc, char **argv) {
  Config c;
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("max-carrier,n", boost::program_options::value<uint32_t>(&c.maxCarrier)->default_value(10), "max. carrier sample to list in the output")
    ("severity,s", boost::program_options::value<std::string>(&c.severity)->default_value("missense"), "max. severity to report [ptv, missense, splice, all]")
    ("transcripts,t", boost::program_options::value<boost::filesystem::path>(&c.transcripts), "list of selected transcripts (if not provided canonical transcripts are used)")
    ("pathogenic,p", "Always include ClinVar pathogenic variants")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("variants.tsv"), "output file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcf), "input VCF file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << std::endl;
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <variants.vcf.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  }

  // ClinVar Pathogenic
  if (vm.count("pathogenic")) c.clinvarPathogenic = true;
  else c.clinvarPathogenic = false;
  
  // Load selected transcripts
  std::set<std::string> seltrans;
  c.selTranscripts = false;
  if (vm.count("transcripts")) {
    if (!(boost::filesystem::exists(c.transcripts) && boost::filesystem::is_regular_file(c.transcripts) && boost::filesystem::file_size(c.transcripts))) {
      c.selTranscripts = true;
      std::string line;
      std::ifstream ifile(c.transcripts.string().c_str());
      if (ifile.is_open()) {
	while (getline(ifile, line)) seltrans.insert(line);
	ifile.close();
      }
    }
  }

  // Limit output up to this severity (maxSOterm)
  std::string maxSOterm = "intergenic_variant";
  if (c.severity == "ptv") maxSOterm = "transcript_amplification";
  else if (c.severity == "missense") maxSOterm = "missense_variant";
  else if (c.severity == "splice") maxSOterm = "splice_region_variant";
  std::string consSO[] = {"transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","splice_region_variant","incomplete_terminal_codon_variant","stop_retained_variant","synonymous_variant","coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_ablation","regulatory_region_amplification","feature_elongation","regulatory_region_variant","feature_truncation","intergenic_variant"};  
  std::vector<std::string> consequence_all(consSO, consSO + sizeof(consSO) / sizeof(consSO[0]));
  std::set<std::string> consequence_filter;
  for(uint32_t i = 0; i < consequence_all.size(); ++i) {
    consequence_filter.insert(consequence_all[i]);
    if (consequence_all[i] == maxSOterm) break;
  }
  
  // Load bcf file
  htsFile* ifile = bcf_open(c.vcf.string().c_str(), "r");
  if (ifile == NULL) {
    std::cerr << "Fail to load " << c.vcf.string().c_str() << std::endl;
    return 1;
  }
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  char* htxt = NULL;
  int32_t hlen;
  htxt = bcf_hdr_fmt_text(hdr, 1, &hlen);
  typedef std::map<std::string, int32_t> TColumnMap;
  TColumnMap cmap;
  if (!getCSQ(std::string(htxt), cmap)) {
    std::cerr << "Fail to parse BCF header." << std::endl;
    return 1;
  }
  
  // Read SV information
  int ngt = 0;
  int32_t* gt = NULL;
  int32_t ncsq = 0;
  char* csq = NULL;

  std::cout << "chr\tpos\tref\talt\texisting_variation\tsingleton\tgt\tsymbol\texon\tstrand\t";
  std::cout << "biotype\tconsequence\tclin_sig\tpopmax\thgvsc\thgvsp\timpact\t";
  std::cout << "polyphen\tsift\tLoFtool\tMaxEntScan(Ref,Alt,Diff)\tcanonical\t";
  std::cout << "carrier\taf\tmissingrate" << std::endl;

  // Parse VCF records
  bcf1_t* rec = bcf_init();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_INFO);

    // Only bi-allelic
    if (rec->n_allele == 2) {
      // Parse VEP
      typedef std::vector<std::string> TStrParts;
      std::string vep;
      if (bcf_get_info_string(hdr, rec, "CSQ", &csq, &ncsq) > 0) vep = std::string(csq);

      // Iterate all transcripts
      TStrParts mgenes;
      boost::split(mgenes, vep, boost::is_any_of(std::string(",")));
      for(TStrParts::const_iterator mgIt = mgenes.begin(); mgIt != mgenes.end(); ++mgIt) {
	TStrParts vals;
	boost::split(vals, *mgIt, boost::is_any_of(std::string("|")));

	// Find transcript
	std::string transcript("NA");
	if (vals[cmap.find("Feature")->second].size()) transcript = vals[cmap.find("Feature")->second];
	std::string canonical("NA");
	if (vals[cmap.find("CANONICAL")->second].size()) canonical = vals[cmap.find("CANONICAL")->second];
	  
	// Is this a selected transcript?
	if (c.selTranscripts) {
	  if (seltrans.find(transcript) == seltrans.end()) continue;
	} else {
	  if (canonical != "YES") continue;
	}

	// Is this a candidate consequence?
	std::string clinsig("unknown");
	if (vals[cmap.find("CLIN_SIG")->second].size()) clinsig = vals[cmap.find("CLIN_SIG")->second];
	bool pathogenic = false;
	if ((c.clinvarPathogenic) && (clinsig.find("pathogenic") != std::string::npos)) pathogenic = true;
	if ((!candidateVar(vals, cmap, consequence_filter)) && (!pathogenic)) continue;

	// Check pop max
	float popmax = popMax(vals, cmap);
	if ((popmax >= 0.01) && (!pathogenic)) continue;
	  
	// Debug
	//for (typename TColumnMap::const_iterator cIt = cmap.begin(); cIt != cmap.end(); ++cIt)
	//if (cIt->second < (int32_t) vals.size()) std::cerr << cIt->first << ',' << vals[cIt->second] << std::endl;
	//std::cerr << std::endl;

	// Get allele frequency
	bcf_unpack(rec, BCF_UN_ALL);
	bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
	int32_t ac[2];
	ac[0] = 0;
	ac[1] = 0;
	std::set<std::string> carrier;
	int32_t uncalled = 0;
	int32_t singlecarrier = 0;
	for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	  if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	    int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	    ++ac[bcf_gt_allele(gt[i*2])];
	    ++ac[bcf_gt_allele(gt[i*2 + 1])];
	    if (gt_type != 0) {
	      carrier.insert(hdr->samples[i]);
	      singlecarrier = gt_type;
	  }
	  }
	}
	if (carrier.size() >= 1) {
	  // Compute GT stats
	  float af = (float) ac[1] / (float) (ac[0] + ac[1]);
	  float missingRate = (float) uncalled / (float) bcf_hdr_nsamples(hdr);
	  std::string gtstr = "NA";
	  if (carrier.size() == 1) {
	    if (singlecarrier == 1) gtstr = "het";
	    else if (singlecarrier == 2) gtstr = "hom";
	  }
	  
	  // Get other fields for output
	  std::string polyphen("unknown");
	  if (vals[cmap.find("PolyPhen")->second].size()) polyphen = vals[cmap.find("PolyPhen")->second];
	  std::string sift("unknown");
	  if (vals[cmap.find("SIFT")->second].size()) sift = vals[cmap.find("SIFT")->second];
	  std::string loftool("NA");
	  if ((cmap.find("LoFtool") != cmap.end()) && (vals[cmap.find("LoFtool")->second].size())) loftool = vals[cmap.find("LoFtool")->second];
	  std::string mescan("NA");
	  if ((cmap.find("MaxEntScan_diff") != cmap.end()) && (vals[cmap.find("MaxEntScan_diff")->second].size())) {
	    mescan = vals[cmap.find("MaxEntScan_ref")->second] + "," + vals[cmap.find("MaxEntScan_alt")->second] + "," + vals[cmap.find("MaxEntScan_diff")->second];
	  }
	  std::string hgvsp("NA");
	  if (vals[cmap.find("HGVSp")->second].size()) hgvsp = vals[cmap.find("HGVSp")->second];
	  std::string hgvsc("NA");
	  if (vals[cmap.find("HGVSc")->second].size()) hgvsc = vals[cmap.find("HGVSc")->second];
	  std::string impact("NA");
	  if (vals[cmap.find("IMPACT")->second].size()) impact = vals[cmap.find("IMPACT")->second];
	  std::string exvar("NA");
	  if (vals[cmap.find("Existing_variation")->second].size()) exvar = vals[cmap.find("Existing_variation")->second];
	  std::string strand("NA");
	  if (vals[cmap.find("STRAND")->second].size()) strand = vals[cmap.find("STRAND")->second];
	  std::string exon("NA");
	  if (vals[cmap.find("EXON")->second].size()) {
	    exon = vals[cmap.find("EXON")->second];
	    std::replace(exon.begin(), exon.end(), '/', ',');
	  }
	  std::string biotyp("NA");
	  if (vals[cmap.find("BIOTYPE")->second].size()) biotyp = vals[cmap.find("BIOTYPE")->second];
	  std::string symb("NA");
	  if (vals[cmap.find("SYMBOL")->second].size()) symb = vals[cmap.find("SYMBOL")->second];
	  std::string cons("NA");
	  if (vals[cmap.find("Consequence")->second].size()) cons = vals[cmap.find("Consequence")->second];
	  
	  // Build carrier string
	  std::string carrierstr;
	  if (carrier.empty()) carrierstr = "NA";
	  else {
	    std::set<std::string>::iterator iter = carrier.begin();
	    std::set<std::string>::iterator iterEnd = carrier.begin();
	    carrierstr = *iter;
	    if (c.maxCarrier < carrier.size()) std::advance(iterEnd, c.maxCarrier);
	    else iterEnd = carrier.end();
	    for(++iter; iter != iterEnd; ++iter) carrierstr += "," + *iter;
	  }

	  // Output record
	  std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << rec->pos + 1 << "\t" << rec->d.allele[0] << "\t";
	  std::cout << rec->d.allele[1] << "\t" << exvar << "\t" << carrierstr << "\t" << gtstr << "\t";
	  std::cout << symb << "\t" << exon << "\t" << strand << "\t" << biotyp << "\t";
	  std::cout << cons << "\t" << clinsig << "\t";
	  std::cout << popmax << "\t" << hgvsc << "\t" << hgvsp << "\t";
	  std::cout << impact << "\t" << polyphen << "\t" << sift << "\t" << loftool << "\t" << mescan << "\t" << canonical << "\t";
	  std::cout << carrier.size() << "\t" << af << "\t" << missingRate << std::endl;
	}
      }
    }
  }

  // Clean-up
  if (gt != NULL) free(gt);
  if (csq != NULL) free(csq);
  if (htxt != NULL) free(htxt);
  
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}

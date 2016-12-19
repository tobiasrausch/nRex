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
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include "htslib/vcf.h"


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
	      for(;itC != itCEnd; ++itC, ++i) cmap.insert(std::make_pair(*itC, i));
	    }
	  }
	}
      }
    }
  }
  return true;
}


template<typename TStrParts, typename TMap>
bool
candidateVar(TStrParts const& vals, TMap const& cmap) {
  // Filter sets
  std::string tmp[] = {"transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", "regulatory_region_ablation"};
  //std::string tmp[] = {"TF_binding_site_variant"};
  std::set<std::string> consequence_filter(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
      
  // Filter Consequence
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
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <snps.vep.bcf> " << std::endl;
    return 1; 
  }

  // Load bcf file
  htsFile* ifile = bcf_open(argv[1], "r");
  if (ifile == NULL) {
    std::cerr << "Fail to load " << argv[1] << std::endl;
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

  std::cout << "chr\tpos\tref\talt\texisting_variation\tsingleton\tgt\tsymbol\tbiotype\tconsequence\tclin_sig\tpopmax\thgvsc\thgvsp\timpact\tcarrier\taf\tmissingrate" << std::endl;
    
  // Parse VCF records
  bcf1_t* rec = bcf_init();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    // Only bi-allelic
    if (rec->n_allele == 2) {
      bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
      int32_t ac[2];
      ac[0] = 0;
      ac[1] = 0;
      std::string rareCarrier = "NA";
      int32_t uncalled = 0;
      int32_t singlecarrier = 0;
      int32_t carrier = 0;
      for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
	if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	  int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	  ++ac[bcf_gt_allele(gt[i*2])];
	  ++ac[bcf_gt_allele(gt[i*2 + 1])];
	  if (gt_type != 0) {
	    if (singlecarrier) rareCarrier = "NA";
	    else {
	      rareCarrier = hdr->samples[i];
	      singlecarrier = gt_type;
	    }
	    ++carrier;
	  }
	}
      }
      if (carrier) {
	// Compute GT stats
	float af = (float) ac[1] / (float) (ac[0] + ac[1]);
	float missingRate = (float) uncalled / (float) bcf_hdr_nsamples(hdr);
	std::string gtstr = "NA";
	if (rareCarrier != "NA") {
	  if (singlecarrier == 1) gtstr = "het";
	  else if (singlecarrier == 2) gtstr = "hom";
	}

	// Parse VEP
	typedef std::vector<std::string> TStrParts;
	std::string vep;
	if (bcf_get_info_string(hdr, rec, "CSQ", &csq, &ncsq) > 0) vep = std::string(csq);
	TStrParts mgenes;
	boost::split(mgenes, vep, boost::is_any_of(std::string(",")));
	int32_t prevImpact = 0;
	for(TStrParts::const_iterator mgIt = mgenes.begin(); mgIt != mgenes.end(); ++mgIt) {
	  TStrParts vals;
	  boost::split(vals, *mgIt, boost::is_any_of(std::string("|")));
	  float popmax = 0;
	  std::string tmp[] = {"AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_AF_AFR", "ExAC_AF_AMR", "ExAC_AF_EAS", "ExAC_AF_NFE", "ExAC_AF_SAS"};
	  typedef std::set<std::string> TAfSet;
	  TAfSet afs(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
	  for(typename TAfSet::const_iterator afIt = afs.begin(); afIt != afs.end(); ++afIt) {
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
		afval = boost::lexical_cast<float>(freqstr);
		if (afval > popmax) popmax = afval;
		//std::cerr << *afIt << "\t" << afval << "\t" << afstr << "\t" << popmax << std::endl;
	      }
	    }
	  }
	  std::string clinsig("unknown");
	  if (vals[cmap.find("CLIN_SIG")->second].size()) clinsig = vals[cmap.find("CLIN_SIG")->second];
	  std::string hgvsp("NA");
	  if (vals[cmap.find("HGVSp")->second].size()) hgvsp = vals[cmap.find("HGVSp")->second];
	  std::string impact("NA");
	  if (vals[cmap.find("IMPACT")->second].size()) impact = vals[cmap.find("IMPACT")->second];
	  std::string exvar("NA");
	  if (vals[cmap.find("Existing_variation")->second].size()) exvar = vals[cmap.find("Existing_variation")->second];
	  int32_t impscore = 0;
	  if (impact == "LOW") impscore = 1;
	  else if (impact == "MODIFIER") impscore = 2;
	  else if (impact == "MODERATE") impscore = 3;
	  else if (impact == "HIGH") impscore = 4;	  
	  else {
	    std::cerr << "Unknown IMPACT!" << std::endl;
	    return - 1;
	  }
	  //for (typename TColumnMap::const_iterator cIt = cmap.begin(); cIt != cmap.end(); ++cIt) std::cout << cIt->first << ',' << vals[cIt->second] << std::endl;

	  // Candidate variant ?
	  if ((popmax < 0.01) && (candidateVar(vals, cmap)) && (impscore > prevImpact)) {	    
	    std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << rec->pos + 1 << "\t" << rec->d.allele[0] << "\t" << rec->d.allele[1] << "\t" << exvar << "\t" << rareCarrier << "\t" << gtstr << "\t" << vals[cmap.find("SYMBOL")->second] << "\t" << vals[cmap.find("BIOTYPE")->second] << "\t" << vals[cmap.find("Consequence")->second] << "\t" << clinsig << "\t" << popmax << "\t" << vals[cmap.find("HGVSc")->second] << "\t" << hgvsp << "\t" << impact << "\t" << carrier << "\t" << af << "\t" << missingRate << std::endl;
	    prevImpact = impscore;
	  }
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
library(GenomicFeatures)

exonTable = function(id) {
	  db=makeTxDbFromUCSC(genome=id, tablename="ccdsGene")
	  ex=reduce(exons(db), ignore.strand=T)
	  ex=keepStandardChromosomes(ex)
	  ex=ex[width(ex)>1,]
	  df=data.frame(chr=seqnames(seqinfo(ex)), start=1, end=seqlengths(seqinfo(ex)), type="wgs")
	  gz = gzfile(paste0(id, ".wgs.bed.gz"), "w")
	  write.table(df, gz, quote=F, row.names=F, col.names=F, sep="\t")
	  close(gz)
	  gz = gzfile(paste0(id, ".wes.bed.gz"), "w")
	  df=data.frame(chr=seqnames(ex), start=start(ex), end=end(ex), type="wes")
	  write.table(df, gz, quote=F, row.names=F, col.names=F, sep="\t")
	  close(gz)
}

exonTable("hg19")
exonTable("mm10")

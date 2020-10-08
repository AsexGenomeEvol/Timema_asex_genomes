library('AsexStats')

tab <- read.table('tables/mapped_reseq_coverages.tsv', header = T, stringsAsFactors=T)
for(sp in timemas$codes){
  sp_name = substr(sp, 3, 5);
  sp_tab = tab[substr(tab$sample, 1, 3) %in% sp_name,];
  #    sample      cov
  sp_tab = data.frame(id = sp_tab$sample, path = paste0('data/mapped_reseq_reads/', sp_tab$sample ,'_to_b3v08_mapped_within_scfs.bam'), depth = round(sp_tab$cov, 2), 'read length' = 120, check.names = F)
  write.table(sp_tab, file = paste0('data/genotyping/', sp, '_samples.txt'), row.names = F, quote = F, sep = '\t')
}
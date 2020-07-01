
library(BUSpaRse)

createTr2g <- function(fasta_file, kallisto_out_path) {
  tr2g <- transcript2gene(fasta_file = fasta_file,
                          kallisto_out_path = kallisto_out_path)
  write.table(tr2g, paste0(kallisto_out_path, "/transcripts_to_genes.txt"), quote=F, row.names = F, col.names=F, sep="\t")
}

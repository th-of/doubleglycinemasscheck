# library("CNEr")
library("microseq")
# writeFasta(contigs[6,], "contig6.fasta", width = 80)
library("Peptides")
library("dplyr")
library("xlsx")

contig <- readFasta("/home/thomas/Documents/Genomics/KAI8_plus_KAI13/Assembly/contigs_ORFs.fasta")
#contig <- readFasta("contigs_orfs.fasta")

#contig6 <- cbind(contig6, mw=mw(contig6$Sequence, monoisotopic = TRUE))
#contig6 <- cbind(contig6, Interesting=between(contig6$mw, 1625, 1629))
#contig6 <- cbind(contig6, Interesting2=between(contig6$mw, 1852, 1856))

hits <- c()

for (i in contig$Sequence){
  for (j in strsplit(i, "GG")[[1]]){
    if (between(mw(j), 1625, 1629) == TRUE | between(mw(j), 1852, 1856) == TRUE){
      a <- c(j, mw(j))
      hits <- rbind(hits, a)
    } 
  }
}
colnames(hits) <- c("Sequence (|GG|)", "Mass (Da)")
hits <- as.data.frame(hits)
hits <- cbind(hits, 'Charge (Lehninger)'=charge(hits$Sequence, pH = 7, pKscale = "Lehninger"))
hits <- cbind(hits, 'Hydrophobicity (KyteDoolittle)'=hydrophobicity(hits$Sequence, scale = "KyteDoolittle"))

loc <- c()

for (k in hits$Sequence){
  loc <- append(loc, contig$Header[which(grepl(k, contig$Sequence))])
}

hits <- cbind(hits, 'ORF location'=loc)

write.xlsx(hits, "potential_bacteriocins_plus.xlsx")


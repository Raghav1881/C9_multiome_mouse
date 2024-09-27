library(ChIPseeker)
library(EnsDb.Mmusculus.v79)

edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
oligo <- subset(integrated_multiome, idents = 'Oligo')
DefaultAssay(oligo) <- 'peaks'

# Subset oligo peaks where DARs are
oligo_range <- oligo[cfeature$query_region, ]
oligoPeakAnno <- annotatePeak(granges(oligo_range), tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(oligoPeakAnno)
plotDistToTSS(oligoPeakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
pathway1 <- enrichKEGG(as.data.frame(oligoPeakAnno)$geneId, organism = 'mmu')

motifs <- li


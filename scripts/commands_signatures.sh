cd $HOME/ebicancerworkshop201607/vizu
mkdir -p signatureResults

tree data/signature/
less data/signature/S01.mutect.somatic.vcf

R

library(SomaticSignatures)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(Cairo)

files <- list.files("data/signature",pattern=".vcf$",recursive=T,full.names=TRUE)
files

vranges <- lapply(files, function(v) readVcfAsVRanges(v,"hs37d5"))

vranges.cat <- do.call(c,vranges)
vranges.cat

print(table(sampleNames(vranges.cat)))

mc <- mutationContext(vranges.cat, BSgenome.Hsapiens.1000genomes.hs37d5)

mc

mm <- motifMatrix(mc, group = "sampleNames", normalize=TRUE)
dim(mm)

gof_nmf <- assessNumberSignatures(mm, 2:7, nReplicates = 5)

Cairo(file="signatureResults/plotNumberOfSignatures.pdf", type="pdf", units="in", width=9, height=8, dpi=72)
plotNumberSignatures(gof_nmf)
dev.off()

sigs_nmf = identifySignatures(mm, 7, nmfDecomposition)

library(pheatmap)
Cairo(file="signatureResults/plot7Signatures_heatmat.pdf", type="pdf", units="in", width=9, height=6, dpi=72)
pheatmap(samples(sigs_nmf),cluster_cols=F, clustering_distance_cols = "correlation")
dev.off()

Cairo(file="signatureResults/plot7Signatures.pdf", type="pdf", units="in", width=10, height=8, dpi=72)
plotSignatures(sigs_nmf,normalize=TRUE, percent=FALSE) + ggtitle("Somatic Signatures: NMF - Barchart") + scale_fill_brewer(palette = "Set2")
dev.off()

Cairo(file="signatureResults/PlotSampleContribution7Signatures.pdf", type="pdf", units="in", width=9, height=6, dpi=72)
plotSamples(sigs_nmf, normalize=TRUE) + scale_y_continuous(breaks=seq(0, 1, 0.2), expand = c(0,0))+ theme(axis.text.x = element_text(size=6))
dev.off()

library(deconstructSigs)

sigs.input=as.data.frame(t(mm))
colnames(sigs.input)=c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G",
 "C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C",
 "T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A",
 "C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
 "T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G",
 "A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C",
 "G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A",
 "A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
 "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G",
 "T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C",
 "C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A",
 "T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
 "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G",
 "G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

Cairo(file=paste("signatureResults/PlotSample",i,"deconstructAlexandrov_pie.pdf",sep="_"), type="pdf", units="in", width=9, height=6, dpi=72)
layout(matrix(1:9,nrow=3,byrow=T))
for (i in rownames(sigs.input)) {
	output.sigs = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = i)
	makePie(output.sigs)
}
dev.off()

q("yes")

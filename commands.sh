
export APP_ROOT=/home/training/Applications/
export PICARD_JAR=$APP_ROOT/picard-tools/picard.jar
export SNPEFF_HOME=$APP_ROOT/snpEff/
export GATK_JAR=$APP_ROOT/gatk/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$APP_ROOT/bvatools-1.6/bvatools-1.6-full.jar
export REF=/home/training/ebicancerworkshop201507/reference



cd /home/training/ebicancerworkshop201507/visualization/


tree


less data/breakdancer.somatic.tsv
less data/mutec.somatic.vcf
less data/scones.somatic.30k.tsv


R


library(circlize)


par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.05)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))


circos.initializeWithIdeogram(species = "hg19")


snv_tmp=read.table("data/mutec.somatic.vcf",comment.char="#")
snv=cbind(snv_tmp[,1:2],snv_tmp[,2]+1)
circos.genomicTrackPlotRegion(snv,stack=TRUE, panel.fun = function(region, value, ...) {
	circos.genomicPoints(region, value, cex = 0.05, pch = 9,col='orange' , ...)
})


cnv=read.table("data/scones.somatic.30k.tsv",header=T)
dup=cnv[cnv[,5]>2,]
del=cnv[cnv[,5]<2,]
circos.genomicTrackPlotRegion(dup, stack = TRUE,panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "red",bg.border = NA, cex=1 , ...)
})
circos.genomicTrackPlotRegion(del, stack = TRUE,panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "blue",bg.border = NA, cex=1 , ...)
})


sv=read.table("data/breakdancer.somatic.tsv",header=T)
typeE=c("DEL","INS","INV")
colE=c("blue","black","green")
for (i in 1:3) { 
        bed_list=sv[sv[,8]==typeE[i],c(3,4,5,7)]
        circos.genomicTrackPlotRegion(bed_list,stack=TRUE, panel.fun = function(region, value, ...) {
                circos.genomicPoints(region, value, cex = 0.5, pch = 16, col = colE[i], ...)
        })
}

bed1=data.frame(chr=as.character(as.vector(sv[sv[,8]=="TRA",3])),start=as.integer(as.vector(sv[sv[,8]=="TRA",4]-5)),end=as.integer(as.vector(sv[sv[,8]=="TRA",4]+5 )),value1=as.numeric(as.vector(rnorm(dim(sv[sv[,8]=="TRA",])[1]))))
bed2=data.frame(chr=as.character(as.vector(sv[sv[,8]=="TRA",6])),start=as.integer(as.vector(sv[sv[,8]=="TRA",5]-5)),end=as.integer(as.vector(sv[sv[,8]=="TRA",5]+5 )),value2=as.numeric(as.vector(rnorm(dim(sv[sv[,8]=="TRA",])[1]))))


if (dim(bed1)[1] > 0 & dim(bed2)[1] > 0) {
        for (i in 1:dim(bed1)[1]) {
                circos.link(bed1[i,1],bed1[i,2],bed2[i,1],bed2[i,2])
        }
}


title("Somatic calls (SNV - SV - CNV) of sample LR376")
legend(0.7,1.4,legend=c("SNV", "CNV-DUPLICATION","CNV-DELETION","SV-DELETION","SV-INSERTION","SV-INVERSION"),col=c("orange","red","blue","blue","black","green","red"),pch=c(16,15,15,16,16,16,16,16),cex=0.75,title="Tracks:",bty='n')
legend(0.6,0.95,legend="SV-TRANSLOCATION",col="black",lty=1,cex=0.75,lwd=1.2,bty='n')


q("yes")


mkdir -p contamination


grep SOMATIC ../SNV/pairedVariants/mutect.vcf \
 | awk 'BEGIN {OFS="\t"} NR > 1 {print $1 , $2 , $4 , $5}' \
 > contamination/mutect.snpPos.tsv


for i in normal tumor
do
  mkdir -p contamination/${i}
  java -Xmx2G -jar $BVATOOLS_JAR basefreq \
    --pos contamination/mutect.snpPos.tsv \
    --bam ../SNV/alignment/${i}/${i}.sorted.dup.recal.bam \
    --out contamination/${i}/${i}.somaticSnpPos \
    --perRG
done


less  contamination/normal/normal.somaticSnpPos.normal_62DPDAAXX_8.alleleFreq.csv


for i in contamination/*/*.somaticSnpPos*_?.alleleFreq.csv
do
  NAME=`echo $i | sed 's/.*somaticSnpPos.\(.*\).alleleFreq.csv/\1/g'`
  echo "--freq $NAME $i";done | tr '\n' ' '
done


java -Xmx2G -jar $BVATOOLS_JAR clustfreq \
--snppos contamination/mutect.snpPos.tsv \
--threads 3 \
--prefix sampleComparison \
--outputFreq \
--freq normal_62DPDAAXX_8 contamination/normal/normal.somaticSnpPos.normal_62DPDAAXX_8.alleleFreq.csv \
--freq normal_62DVGAAXX_1 contamination/normal/normal.somaticSnpPos.normal_62DVGAAXX_1.alleleFreq.csv  \
--freq normal_62MK3AAXX_5 contamination/normal/normal.somaticSnpPos.normal_62MK3AAXX_5.alleleFreq.csv \
--freq normal_A81DF6ABXX_1 contamination/normal/normal.somaticSnpPos.normal_A81DF6ABXX_1.alleleFreq.csv \
--freq normal_A81DF6ABXX_2 contamination/normal/normal.somaticSnpPos.normal_A81DF6ABXX_2.alleleFreq.csv \
--freq normal_BC04D4ACXX_2 contamination/normal/normal.somaticSnpPos.normal_BC04D4ACXX_2.alleleFreq.csv  \
--freq normal_BC04D4ACXX_3 contamination/normal/normal.somaticSnpPos.normal_BC04D4ACXX_3.alleleFreq.csv \
--freq normal_BD06UFACXX_4 contamination/normal/normal.somaticSnpPos.normal_BD06UFACXX_4.alleleFreq.csv  \
--freq normal_BD06UFACXX_5 contamination/normal/normal.somaticSnpPos.normal_BD06UFACXX_5.alleleFreq.csv  \
--freq tumor_62DU0AAXX_8 contamination/tumor/tumor.somaticSnpPos.tumor_62DU0AAXX_8.alleleFreq.csv  \
--freq tumor_62DU6AAXX_8 contamination/tumor/tumor.somaticSnpPos.tumor_62DU6AAXX_8.alleleFreq.csv  \
--freq tumor_62DUUAAXX_8 contamination/tumor/tumor.somaticSnpPos.tumor_62DUUAAXX_8.alleleFreq.csv  \
--freq tumor_62DUYAAXX_7 contamination/tumor/tumor.somaticSnpPos.tumor_62DUYAAXX_7.alleleFreq.csv  \
--freq tumor_62DVMAAXX_4 contamination/tumor/tumor.somaticSnpPos.tumor_62DVMAAXX_4.alleleFreq.csv  \
--freq tumor_62DVMAAXX_5 contamination/tumor/tumor.somaticSnpPos.tumor_62DVMAAXX_5.alleleFreq.csv  \
--freq tumor_62DVMAAXX_6 contamination/tumor/tumor.somaticSnpPos.tumor_62DVMAAXX_6.alleleFreq.csv  \
--freq tumor_62DVMAAXX_7 contamination/tumor/tumor.somaticSnpPos.tumor_62DVMAAXX_7.alleleFreq.csv  \
--freq tumor_62DVMAAXX_8 contamination/tumor/tumor.somaticSnpPos.tumor_62DVMAAXX_8.alleleFreq.csv  \
--freq tumor_62JREAAXX_3 contamination/tumor/tumor.somaticSnpPos.tumor_62JREAAXX_3.alleleFreq.csv  \
--freq tumor_62JREAAXX_4 contamination/tumor/tumor.somaticSnpPos.tumor_62JREAAXX_4.alleleFreq.csv  \
--freq tumor_62JREAAXX_5 contamination/tumor/tumor.somaticSnpPos.tumor_62JREAAXX_5.alleleFreq.csv  \
--freq tumor_62JREAAXX_6 contamination/tumor/tumor.somaticSnpPos.tumor_62JREAAXX_6.alleleFreq.csv  \
--freq tumor_62JREAAXX_7 contamination/tumor/tumor.somaticSnpPos.tumor_62JREAAXX_7.alleleFreq.csv  \
--freq tumor_62JREAAXX_8 contamination/tumor/tumor.somaticSnpPos.tumor_62JREAAXX_8.alleleFreq.csv  \
--freq tumor_AC0756ACXX_4 contamination/tumor/tumor.somaticSnpPos.tumor_AC0756ACXX_4.alleleFreq.csv  \
--freq tumor_AC0756ACXX_5 contamination/tumor/tumor.somaticSnpPos.tumor_AC0756ACXX_5.alleleFreq.csv  \
--freq tumor_AD08C1ACXX_1 contamination/tumor/tumor.somaticSnpPos.tumor_AD08C1ACXX_1.alleleFreq.csv  \
--freq tumor_BD08K8ACXX_1 contamination/tumor/tumor.somaticSnpPos.tumor_BD08K8ACXX_1.alleleFreq.csv


R


colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol = c("blue");
        if(grepl("normal", a$label)) {
          labCol = c("red");
        }
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}


dataMatrix <- read.csv("sampleComparison.dist.csv", row.names=1, header=TRUE)
hc <- hclust(as.dist(dataMatrix));
hcd = as.dendrogram(hc)


clusDendro = dendrapply(hcd, colLab)


data <- read.csv("sampleComparison.freq.csv", header=FALSE,row.names=1, colClasses=c("character", rep("numeric",17)))
colLanes <- rownames(data)
colLanes[grep("normal", colLanes, invert=TRUE)] <- "blue"
colLanes[grep("normal", colLanes)] <- "red"
pca <- prcomp(data)


pdf("sampleComparison.laneMix.pdf")
par(mar=c(3,3,1,12))
cols <- c("red","blue")
plot(clusDendro, main = "Lane distances", horiz=TRUE)
legend("top", legend = c("Normal","Tumor"), fill = cols, border = cols)

par(mar=c(1,1,1,1))
screeplot(pca, type="lines")
plot(pca$x[,1:2])
text(pca$x[,1:2], rownames(data), col=colLanes)
plot(pca$x[,2:3])
text(pca$x[,2:3], rownames(data), col=colLanes)
dev.off()


q("yes")



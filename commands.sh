
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


## initiualize plot
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


# Generate a part of the command
for i in contamination/*/*.somaticSnpPos*_?.alleleFreq.csv
do
  NAME=`echo $i | sed 's/.*somaticSnpPos.\(.*\).alleleFreq.csv/\1/g'`
  echo "--freq $NAME $i";done | tr '\n' ' '
done

# Copy this output and paste it at the end of the command like so
java -Xmx2G -jar $BVATOOLS_JAR clustfreq \
--snppos contamination/mutect.snpPos.tsv \
--threads 3 \
--prefix sampleComparison \
--outputFreq \
--freq normal_C0LWRACXX_1 alignment/normal/normal.somaticSnpPos.normal_C0LWRACXX_1.alleleFreq.csv \
--freq normal_C0LWRACXX_6 alignment/normal/normal.somaticSnpPos.normal_C0LWRACXX_6.alleleFreq.csv \
--freq normal_C0PTAACXX_6 alignment/normal/normal.somaticSnpPos.normal_C0PTAACXX_6.alleleFreq.csv \
--freq normal_C0PTAACXX_7 alignment/normal/normal.somaticSnpPos.normal_C0PTAACXX_7.alleleFreq.csv \
--freq normal_C0PTAACXX_8 alignment/normal/normal.somaticSnpPos.normal_C0PTAACXX_8.alleleFreq.csv \
--freq normal_C0R2BACXX_6 alignment/normal/normal.somaticSnpPos.normal_C0R2BACXX_6.alleleFreq.csv \
--freq normal_C0R2BACXX_7 alignment/normal/normal.somaticSnpPos.normal_C0R2BACXX_7.alleleFreq.csv \
--freq normal_C0R2BACXX_8 alignment/normal/normal.somaticSnpPos.normal_C0R2BACXX_8.alleleFreq.csv \
--freq normal_D0YR4ACXX_1 alignment/normal/normal.somaticSnpPos.normal_D0YR4ACXX_1.alleleFreq.csv \
--freq normal_D0YR4ACXX_2 alignment/normal/normal.somaticSnpPos.normal_D0YR4ACXX_2.alleleFreq.csv \
--freq tumor_BC0TV0ACXX_8 alignment/tumor/tumor.somaticSnpPos.tumor_BC0TV0ACXX_8.alleleFreq.csv \
--freq tumor_C0LVJACXX_6 alignment/tumor/tumor.somaticSnpPos.tumor_C0LVJACXX_6.alleleFreq.csv \
--freq tumor_C0PK4ACXX_7 alignment/tumor/tumor.somaticSnpPos.tumor_C0PK4ACXX_7.alleleFreq.csv \
--freq tumor_C0PK4ACXX_8 alignment/tumor/tumor.somaticSnpPos.tumor_C0PK4ACXX_8.alleleFreq.csv \
--freq tumor_C0R29ACXX_7 alignment/tumor/tumor.somaticSnpPos.tumor_C0R29ACXX_7.alleleFreq.csv \
--freq tumor_C0R29ACXX_8 alignment/tumor/tumor.somaticSnpPos.tumor_C0R29ACXX_8.alleleFreq.csv \
--freq tumor_C0TTBACXX_3 alignment/tumor/tumor.somaticSnpPos.tumor_C0TTBACXX_3.alleleFreq.csv \
--freq tumor_D114WACXX_8 alignment/tumor/tumor.somaticSnpPos.tumor_D114WACXX_8.alleleFreq.csv


fileName <- "sampleComparison"
normalName <- "normal"


distFile <- paste(fileName, ".dist.csv",sep="")

#distName.noext = sub("[.][^.]*$", "", distName, perl=TRUE)

dataMatrix <- read.csv(distFile, row.names=1, header=TRUE)
hc <- hclust(as.dist(dataMatrix));
hcd = as.dendrogram(hc)

colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol = c("blue");
        if(grepl(normalName, a$label)) {
          labCol = c("red");
        }
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}

# using dendrapply
clusDendro = dendrapply(hcd, colLab)
cols <- c("red","blue")


freqFile <- paste(fileName, ".freq.csv",sep="")
data <- read.csv(freqFile, header=FALSE,row.names=1, colClasses=c("character", rep("numeric",4))
colLanes <- rownames(data)
colLanes[grep(normalName, colLanes, invert=TRUE)] <- "blue"
colLanes[grep(normalName, colLanes)] <- "red"
pca <- prcomp(data)

# make plot
pdfFile <- paste(fileName,".laneMix.pdf", sep="")
pdfFile
pdf(pdfFile)
par(mar=c(3,3,1,12))
plot(clusDendro, main = "Lane distances", horiz=TRUE)
legend("top", legend = c("Normal","Tumor"), fill = cols, border = cols)

par(mar=c(1,1,1,1))
plot(pca$x[,1:2])
text(pca$x[,1:2], rownames(data), col=colLanes)
screeplot(pca, type="lines")
plot(pca$x[,2:3])
text(pca$x[,2:3], rownames(data), col=colLanes)
dev.off()

pngFile <- paste(fileName,".laneMix.png", sep="")
pngFile
png(pngFile, type="cairo", width=1920, height=1080)
par(mar=c(3,3,1,12))
par(mfrow=c(2,1))
plot(clusDendro, main = "Lane distances", horiz=TRUE)
legend("top", legend = c("Normal","Tumor"), fill = cols, border = cols)

plot(pca$x[,1:2])
text(pca$x[,1:2], rownames(data), col=colLanes)
dev.off()



# Aligned or not, we want them all
samtools view alignment/normal/normal.sorted.bam | awk '{if($10 ~ /TTAGGGTTAGGGTTAGGG/) {SUM++}} END {print "NbTeloReads",SUM}'
samtools view alignment/tumor/tumor.sorted.bam | awk '{if($10 ~ /TTAGGGTTAGGGTTAGGG/) {SUM++}} END {print "NbTeloReads",SUM}'


less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv



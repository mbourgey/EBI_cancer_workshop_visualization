
# Introduction to DNA-Seq processing for cancer data - Visualization
***By Mathieu Bourgey, Ph.D***

# Analysis of cancer data
This workshop will be a mix of different methods to look and explore your data.

We will be working on the same BAMs you generated from the SNV part.
Again, for practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.
This leads to some strange results in this part.

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

### Environment setup
We will need an updated bvatools for these exercises

```{.bash}

export APP_ROOT=/home/training/Applications/
export PICARD_JAR=$APP_ROOT/picard-tools/picard.jar
export SNPEFF_HOME=$APP_ROOT/snpEff/
export GATK_JAR=$APP_ROOT/gatk/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$APP_ROOT/bvatools-1.6/bvatools-1.6-full.jar
export REF=/home/training/ebicancerworkshop201507/reference

```

### Software requirements
These are all already installed, but here are the original links.

  * [Genome MuSiC](http://gmt.genome.wustl.edu/genome-shipit/genome-music/current/)
  * [BVATools](https://bitbucket.org/mugqic/bvatools/downloads)
  * [SAMTools](http://sourceforge.net/projects/samtools/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)
  * [Genome Analysis Toolkit](http://www.broadinstitute.org/gatk/)
  * [Picard](http://picard.sourceforge.net/)
  * [SnpEff](http://snpeff.sourceforge.net/)
  * [R]()

# Circular represnetation of your calls
Many tools are available to do this the most common know is circos. But circos is a really not user friendly. In this tutoriel we show you an easy alternative to build circular representation of genomic data.

First we nee to go in the folder to do the analysis

```{.bash}
cd /home/training/ebicancerworkshop201507/visualization/
```

Let see what is in this folder

```{.bash}
tree
```

 data/
  -- breakdancer.somatic.tsv
  -- mutec.somatic.vcf
  -- scones.somatic.30k.tsv
 src/
  -- commands.sh

Take a look of the data files.

These are data of the same paired sample that we worked on during the SNV pratical. But this time the data are not limitated to a short piece of the chromosome 9.


```{.bash}
less data/breakdancer.somatic.tsv
less data/mutec.somatic.vcf
less data/scones.somatic.30k.tsv
```


**What can you see fron this data ?**
[solution](solution/_data1.md)

**Why don't we use the vcf format for all type of call?**
[solution](solution/_data2.md)

**Did you notice something different from the SNV pratical ?**
[solution](solution/_data3.md)


The analysis will be done using the R program

```{.bash}
R
```

We will use the circlize package from the cran R project. This package is dedicated to generate circular plot and had the advantage to provide pre-build function for genomics data. One of the main advantage of this tools is the use of bed format as input data.


```{.bash}
library(circlize)
```

We need to set-up the generic graphical parameters 


```{.bash}
## initiualize plot
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.05)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
```

Let's draw hg19 reference ideograms

```{.bash}
circos.initializeWithIdeogram(species = "hg19")
```


We can now draw 1 track for somatic mutations

```{.bash}
snv_tmp=read.table("data/mutec.somatic.vcf",comment.char="#")
snv=cbind(snv_tmp[,1:2],snv_tmp[,2]+1)
circos.genomicTrackPlotRegion(snv,stack=TRUE, panel.fun = function(region, value, ...) {
	circos.genomicPoints(region, value, cex = 0.05, pch = 9,col='orange' , ...)
})
```


Let's draw the 2 tracks for cnvs. One track for duplication in red and one blue track for deletion.

```{.bash}
cnv=read.table("data/scones.somatic.30k.tsv",header=T)
dup=cnv[cnv[,5]>2,]
del=cnv[cnv[,5]<2,]
circos.genomicTrackPlotRegion(dup, stack = TRUE,panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "red",bg.border = NA, cex=1 , ...)
})
circos.genomicTrackPlotRegion(del, stack = TRUE,panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "blue",bg.border = NA, cex=1 , ...)
})
```

We can cleary see a massive deletion in the chromosome 3, which is a very common observation for kidney cancer (85% of tumor)  


To  finsh we just need to draw 3 tracks + positional links to represent SVs

```{.bash}
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
```

A good graph needs title and legends

```{.bash}
title("Somatic calls (SNV - SV - CNV) of sample LR376")
legend(0.7,1.4,legend=c("SNV", "CNV-DUPLICATION","CNV-DELETION","SV-DELETION","SV-INSERTION","SV-INVERSION"),col=c("orange","red","blue","blue","black","green","red"),pch=c(16,15,15,16,16,16,16,16),cex=0.75,title="Tracks:",bty='n')
legend(0.6,0.95,legend="SV-TRANSLOCATION",col="black",lty=1,cex=0.75,lwd=1.2,bty='n')
```

you should obtain a plot like this one
![circular_view](img/somatic_circular_plot.pdf)

Exercice:


**Generate the graph and save it into a pdf file** 
[solution](solution/_image1.md)


Finally exit R

```{.bash}
q("yes")
```


# Finding contamination
This is more to QC but it can be very helpful to find strange patterns in your samples.

Create output folder contamination

```{.bash}
mkdir -p contamination
```

Extract positions of somatic variants from the SNV pratical 

```{.bash}
grep SOMATIC ../SNV/pairedVariants/mutect.vcf \
 | awk 'BEGIN {OFS="\t"} NR > 1 {print $1 , $2 , $4 , $5}' \
 > contamination/mutect.snpPos.tsv
```

Now we have our positions, we need the read counts *per lane* for these positions.
BVATools does this

```{.bash}
for i in normal tumor
do
  mkdir -p contamination/${i}
  java -Xmx2G -jar $BVATOOLS_JAR basefreq \
    --pos contamination/mutect.snpPos.tsv \
    --bam ../SNV/alignment/${i}/${i}.sorted.dup.recal.bam \
    --out contamination/${i}/${i}.somaticSnpPos \
    --perRG
done
```

We can look at one of the files to see what basefreq extracted

```{.bash}
less  contamination/normal/normal.somaticSnpPos.normal_62DPDAAXX_8.alleleFreq.csv
```

Now we need to extract and format the data so we can create a PCA and some hierarchical clusters

```{.bash}
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
```

Now you should have 2 files
sampleComparison.freq.csv
sampleComparison.dist.csv

One contains vectors of snp frequences, the other contains the pairwise Euclidean distance
Now plot the result in R
```{.bash}
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

```

Look at the graphs.

You could do this directly in R but
1. The basefreq format is not simple to parse
2. When you have thousands of somatics, and/or hundreds of samples, R struggles to build de pairwise distance and the PCA. This is why we precompute it in java before.


# Telomeres
In this first step we will try to qualitatively see if the normal and tumor have different
telomere lengths.

One way to do this is find the telomere motif.

A good link to get various telomere repeats is the [Telomerase Database](http://telomerase.asu.edu/sequences_telomere.html)

First step, count the number of reads with these repeats.

```{.bash}
# Aligned or not, we want them all
samtools view alignment/normal/normal.sorted.bam | awk '{if($10 ~ /TTAGGGTTAGGGTTAGGG/) {SUM++}} END {print "NbTeloReads",SUM}'
samtools view alignment/tumor/tumor.sorted.bam | awk '{if($10 ~ /TTAGGGTTAGGGTTAGGG/) {SUM++}} END {print "NbTeloReads",SUM}'
```

Why did we put multiple copied of the repeat in the search? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_telo.ex1.md)  

Next look in the alignments summary file we generated yesterday and extract the number of aligned reads.

```{.bash}
less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv
```

Now the rest can be done in good old excel.

Open a sheet up, load up your numbers and compute the fold change between normal and tumor.

This case is rather boring, there practically is no change.


## Aknowledgments
This tutorial is an adaptation of the one created by Louis letourneau [here](https://github.com/lletourn/Workshops/tree/ebiCancerWorkshop201407doc/01-SNVCalling.md). I would like to thank and acknowledge Louis for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta. I also want to acknowledge Joel Fillon, Louis Letrouneau (again), Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
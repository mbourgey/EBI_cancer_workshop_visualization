
# Introduction to DNA-Seq processing for cancer data - Interpretation and Visualization
***By Mathieu Bourgey, Ph.D***

================================

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

====================================

## Key Learning Outcomes
After completing this practical the trainee should be able to:
 
 * Perform mutional signature analysis using R
 * Generate circos like graphics using R


## Resources You'll be Using

### Tools Used

 * [R](https://cran.r-project.org/)
 * [R package SomaticSignatures](https://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html)
 * [R package deconstructSigs](https://cran.rstudio.com/web/packages/deconstructSigs/index.html)
 * [R package circlize](https://cran.r-project.org/web/packages/circlize/index.html)
 
-----------------------------
# Somatic mutational signature analysis 

## Introduction
The most common genetic model for cancer development is the accumulation of DNA mutations over time, eventually leading to the disruption or dysregulation of enough key genes that lead cells to uncontrolled growth. Cells in our bodies accumulate DNA mutations over time due to normal aging processes, through exposure to carcinogens or through defects in the cell’s ability to repair mistakes. Recently researchers found a method to take all the single nucleotide mutations identified in tumour cells (somatic SNVs) and cluster them together by the type of the mutation and also what the neighbouring bases are. This is commonly referred to as somatic mutational signatures. 

Common mutational processes that are regularly identified in cancer sequencing are:
 - Age: the aging process. These are high in C/T transitions due to deamination of methyl-cytidine.
 - Smoking: marks exposure to inhaled carcinogens and has high numbers of C/A transversions.
 - UV: UV exposure. These are also high in C/T transitions at di-pyrimidine sites.
 - BRCA: Indicates that the homologous recombination repair pathway is defective.
 - APOBEC: Thought to be marking dysregulated APOBEC enzyme activity on single stranded DNA produced during the repair processing of other lesions such as double stand breaks.
 - MMR: Mismatch repair pathway not working properly. These are high in C/T mutations too.

 
In cohort cancer analysis it is common to try to generate subtypes to group your data based on a particular molecular phenotype. A reason for doing may include finding sets of patients that have a similar form of the disease and therefore all might benefit from a particular treatment. We can use the somatic mutational signatures analysis to group the data from a cohort of patients to inform which genomes are most similar based on the pattern of exposures or processes that have contributed to their genome changes. 

The mathematical framework developed by Alexandrov et al to cluster the somoatic mutations was implemented in MATLAB. We are going to use a version implemented in
R by Gehring et al, called SomaticSignatures package, that is very quick and flexible but currently only accepts point mutations not insertions or deletions (indels). In tests on our data we have found that the Somatic Signatures package in R returns very similar results to the full implementation of Alexandrov’s framework.

## Data Source
We will be working on a seven cancer sample. Some of them come from the CageKid project which is part of ICGC and is focused on renal cancer in many of it's forms and the other come from colon cancer. 


For practical reasons we precomputed the mutect somatic mutations vcf of each sample.

## Environment setup
Everything is already installed on your machine and the analysis will be run using the`R` analysis

```{.bash}
cd $HOME/ebicancerworkshop201607/vizu
mkdir -p signatureResults

```



Let see what the data look like


```{.bash}
tree data/signature/
```

> data/signature/   
> ├── S01.mutect.somatic.vcf   
> ├── S02.mutect.somatic.vcf   
> ├── S03.mutect.somatic.vcf   
> ├── S04.mutect.somatic.vcf   
> ├── S05.mutect.somatic.vcf   
> ├── S06.mutect.somatic.vcf   
> └── S07.mutect.somatic.vcf   


we could explore one vcf file

```{.bash}
less data/signature/S01.mutect.somatic.vcf

```

**What can we see from this vcf compared to the one generated in the first practical ?** [solution](solution/_vcf1.md)


Now that we now what we are working on, we can start and prepare for the analysis

Just start by typing R onto the command line.

```{.bash}
R

```

Load all the package libraries needed for this analysis by running the commands.

```{.R}
library(SomaticSignatures)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(Cairo)

```
## Analysis


### loading mutation
Read in the mutations from the 7 vcf files

```{.R}
files <- list.files("data/signature",pattern=".vcf$",recursive=T,full.names=TRUE)
files

```

> [1] "data/signature/S01.mutect.somatic.vcf"   
> [2] "data/signature/S02.mutect.somatic.vcf"   
> [3] "data/signature/S03.mutect.somatic.vcf"   
> [4] "data/signature/S04.mutect.somatic.vcf"   
> [5] "data/signature/S05.mutect.somatic.vcf"   
> [6] "data/signature/S06.mutect.somatic.vcf"   
> [7] "data/signature/S07.mutect.somatic.vcf"   

Next we read in all the genomic positions of variants in the VCF files using the vranges class.

```{.R}
vranges <- lapply(files, function(v) readVcfAsVRanges(v,"hs37d5"))

```

Don't pay attention to the warnings

Now we can join all the lists of variant positions into one big data set so that it can be processed
together and look at what is contained in the concatenated vranges data

```{.R}
vranges.cat <- do.call(c,vranges)
vranges.cat

```
The first line of output of the vranges.cat shows us that in total we have put over 4,000 mutations. 
For each mutation we record between others:
 - the chromosome position
 - the mutation base changes 
 - the depths
 - the sample of origin 
 

We can print out how many mutations we have read in for each of the cancer samples we
are using by using the command.

```{.R}
print(table(sampleNames(vranges.cat)))

```

> S01 S02 S03 S04 S05 S06 S07   
> 921 485 233 846 967 793 539   


**Could you predict which sample belongs to kidney or colon cancers ?** [solution](solution/_vcf3.md)

now we can use the reference and the position of the mutation to look up the bases on either side of the mutation i.e. the mutation context.


### adding mutation context
Run the mutationContext function of SomaticSignatures.

```{.R}
mc <- mutationContext(vranges.cat, BSgenome.Hsapiens.1000genomes.hs37d5)

```

It is always important to select the correct reference for your data.

**why ?** [solution](solution/_vcf2.md)


We can inspect what information we had added to the vranges.cat objec

```{.R}
mc

```

Notice that the mutation and its context have been added to the last two columns


There are a total of 96 possible single base mutations and context combinations. We can calculate this by, first, listing out the six possible types of single nucleotide mutations: 

 - C/A the reverse compliment (G/T) is also in this group
 - C/G includes (G/C)
 - C/T includes (G/C)
 - T/A includes (A/T)
 - T/C includes (A/G)
 - T/G includes (A/C)


Then we can listing every of  mutation type using the neighbouring bases, on either side of a mutation, also referred to as the mutation context. 

There are 16 possible combinations of mutation contexts. Here [.] stands for one of the mutations listed above.

 - A[.]A A[.]C A[.]G A[.]T
 - C[.]A C[.]C C[.]G C[.]T
 - G[.]A G[.]C G[.]G G[.]T
 - T[.]A T[.]C T[.]G T[.]T

Now if we substitute the [.]’s with each of the 6 different mutations you will find there are
96 possible types of combined mutations and contexts (6 x 16).

**What about a mutation that looks like G[A/C]A, where should this go ?** [solution](solution/_vcf4.md)

 
Now we have all the information that is needed for each sample we can make a matrix that contains counts of mutations in each of the 96 possible combinations of mutations and contexts counting up the totals separately for each sample 

```{.R}
mm <- motifMatrix(mc, group = "sampleNames", normalize=TRUE)
dim(mm)

```

> [1] 96  7

The output of the command show us that there are 96 rows (these are the context
values) and 7 columns which are the 7 samples.

## Running the NMF analysis

Using the matrix we have made we can now run the non-negative matrix factorisation (NMF) process that attempts to find the most stable, grouping solutions for all of the combinations of mutations and contexts. It does this by trying to find similar patterns, or profiles, amongst the samples to sort the data into firstly just 2 groups. This is repeated to get replicate values for each attempt and then separating the data by 3 groups, and then 4 and so on.


```{.R}
gof_nmf <- assessNumberSignatures(mm, 2:10, nReplicates = 5)

```

These parameter choices have been made to keep running time short for this practical !

Visualise the results from the NMF processing by making a pdf of the plot

```{.R}
Cairo(file="signatureResults/plotNumberOfSignatures.pdf", type="pdf", units="in", width=9, height=8, dpi=72)
plotNumberSignatures(gof_nmf)
dev.off()

```

Open up the PDF and examine the curve.


The top plot shows the decreasing residual sum of squares for each increasing number of signatures and the bottom plot the increasing explained variance as the number of potential signatures increases. 


Ideally the best solution will be the lowest number of signatures with a low RSS and a high explained variance


Look at the y-axis scale on the bottom panel. The explained variance is already very high and so close to finding the correct solution for the number of signatures even with just 2. The error bars around each point are fairly small considering we have a very small sample set. Deciding how many signatures are present can be tricky but here let’s go for 8 This is where the gradient of both curves have become flat.

Now we can run the NMF again but this time stipulating that you want to group the data into 3 different mutational signatures.

```{.R}
sigs_nmf = identifySignatures(mm, 8, nmfDecomposition)

```

## Making sens of sample signature

Let's try to cluster samples based on the signture decomposition. 

```{.R}
library(pheatmap)
Cairo(file="signatureResults/plot8Signatures_heatmat.pdf", type="pdf", units="in", width=9, height=6, dpi=72)
pheatmap(samples(sigs_nmf),cluster_cols=F, clustering_distance_cols = "correlation")
dev.off()

```
Open up the `plot8Signatures.pdf` that will have been made.

**Are the coresponding cluster fiting with what we predict based on the number of mutation ?** 

Now, we can visualise the shape of the profiles for these 8 signatures

```{.R}
Cairo(file="signatureResults/plot8Signatures.pdf", type="pdf", units="in", width=10, height=8, dpi=72)
plotSignatures(sigs_nmf,normalize=TRUE, percent=FALSE) + ggtitle("Somatic Signatures: NMF - Barchart") + scale_fill_brewer(palette = "Set2")
dev.off()

```

Open up the `plot8Signatures.pdf` that will have been made.


The 96 possible mutation/context combinations are plotted along the x axis arranged in blocks of 6 lots of 16 (see information above). The height of the bars indicates the frequency of those particular mutation and context combinations in each signature.

Now we can plot out the results for the individual samples in our dataset to show what
proportion of their mutations have been assigned to each of the signatures.

```{.R}
Cairo(file="signatureResults/PlotSampleContribution8Signatures.pdf", type="pdf", units="in", width=9, height=6, dpi=72)
plotSamples(sigs_nmf, normalize=TRUE) + scale_y_continuous(breaks=seq(0, 1, 0.2), expand = c(0,0))+ theme(axis.text.x = element_text(size=6))
dev.off()

```

Open the resulting `PlotSampleContribution8Signatures.pdf`. This shows the results for the mutation grouping for each sample. The samples are listed on the x-axis and the proportion of all mutations for that sample is shown on the y-axis. The colours of the bars indicate what proportion of the mutations for that sample were grouped into each of the signatures. The colour that makes up most of the bar for each sample is called its ”major signature”.


## Interpreting the signature results
In their paper __Alexandrov et al__ used this analysis to generate profiles from the data for more than 7000 tumour samples sequenced through both exome and whole genome approaches. They were able to group the data to reveal which genomes have been exposed to similar mutational processes contributing to the genome mutations. More information can be found on the signatures page of the COSMIC website.


![Alexandrov signatures](img/alexandrov_signatures.png)  


**Can you match up, by eye, the profile shapes against a selection of known mutational signatures supplied ?** [solution](solution/_signatures1.md)


Unfortunately the SomaticSignatures package does not provide any autmated way to deconstruct the signal based on Alexandrov known signatures. To do this task we will need to use another R package, `deconstructSigs`, which implement that.

```{.R}
library(deconstructSigs)

```

First we need to reformat the data to  fit the `deconstructSigs` input format


```{.R}
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

```


We can now plot for each sample the contribution of known mutations
`
``{.R}
Cairo(file=paste("signatureResults/PlotSample",i,"deconstructAlexandrov_pie.pdf",sep="_"), type="pdf", units="in", width=9, height=6, dpi=72)
layout(matrix(1:9,nrow=3,byrow=T))
for (i in rownames(sigs.input)) {
	output.sigs = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = i)
	makePie(output.sigs)
}
dev.off()

```

------------------------------
# Circular representation of somtaic calls
 
Many tools are available to do this the most common know is circos. But circos is a really not user friendly. In this tutoriel we show you an easy alternative to build circular representation of genomic data.

First we nee to go in the folder to do the analysis

```{.bash}
cd /home/training/ebicancerworkshop201507/vizu
```

Let see what is in this folder

```{.bash}
tree  data/vizu/
```

 data/vizu
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

Now we need to extract and format the data so we can create a PCA and some hierarchical clusters. to do this we will use a specific command of bvatools: clustfreq  


The command is really complicated when the number of readgroup is large. So we will first generate a part of the command on the screen 


```{.bash}
for i in contamination/*/*.somaticSnpPos*_?.alleleFreq.csv
do
  NAME=`echo $i | sed 's/.*somaticSnpPos.\(.*\).alleleFreq.csv/\1/g'`
  echo "--freq $NAME $i";done | tr '\n' ' '
done
```

We will then copy this output and paste it at the end of the command bvatools. You should obtain that:

```{.bash}
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
```

Now you should have 2 files
sampleComparison.freq.csv
sampleComparison.dist.csv

One contains vectors of snp frequences, the other contains the pairwise Euclidean distance
Now let's plot the result in R

```{.bash}
R
```

First we need to implement a small function to assign different color to the library in function of which sample they belong to

```{.bash}
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
```

Load data, perform clutering and generate dendogram

```{.bash}
dataMatrix <- read.csv("sampleComparison.dist.csv", row.names=1, header=TRUE)
hc <- hclust(as.dist(dataMatrix));
hcd = as.dendrogram(hc)
```

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

Assign color using dendrapply and the function we created previously

```{.bash}
clusDendro = dendrapply(hcd, colLab)
```

Generate PCA on the data

```{.bash}
data <- read.csv("sampleComparison.freq.csv", header=FALSE,row.names=1, colClasses=c("character", rep("numeric",17)))
colLanes <- rownames(data)
colLanes[grep("normal", colLanes, invert=TRUE)] <- "blue"
colLanes[grep("normal", colLanes)] <- "red"
pca <- prcomp(data)
```

Make plot: a pdf file containing the hierchical clustering (dendogram) + Screeplot + PCA (COMP1 vs. COMP2) + PCA (COMP2 vs. COMP3)

```{.bash}
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
```

Finally exit R
```{.bash}
q("yes")
```

Look at the graphs.

Based on the subset of data we have here a potential library issue could be present.  


But when looking at the enitre set of somatic mutations we can this is not true.
![somatic library clustering](img/somatic_clustering.png)  



You could do this directly in R but
1. The basefreq format is not simple to parse
2. When you have thousands of somatics, and/or hundreds of samples, R struggles to build de pairwise distance and the PCA. This is why we precompute it in java before.


# Other visualizations

Many other visualizations of cancer data are possible. we will not go further in this pratical. But here is non-exhaustive list of other interesting visualization of DNA-seq cancer data:  

1. Somatic mutation distribution by type ![somatic mutation](img/somatic_mutations.png)
2. Genomic context of somatic mutations  ![lego plot](img/somatic_lego_plot.png)
3. Representation of a possible transcriptional bias for somatic mutation ![transcriptional bias](img/somatic_transciptional_bias.png) 



## Aknowledgments
This tutorial is an adaptation of the one created by Louis letourneau [here](https://github.com/lletourn/Workshops/tree/ebiCancerWorkshop201407doc/01-SNVCalling.md). I would like to thank and acknowledge Louis for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta. I also want to acknowledge Joel Fillon, Louis Letrouneau (again), Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.

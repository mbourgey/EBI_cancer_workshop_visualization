
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
We will be working on a CageKid sample pair, patient C0098.
The CageKid project is part of ICGC and is focused on renal cancer in many of it's forms.
The raw data can be found on EGA and calls, RNA and DNA, can be found on the ICGC portal. 
For more details about [CageKid](http://www.cng.fr/cagekid/)

For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.

------------------------------
# Circular represenation of somtaic calls
 
## Introduction

This short workshop will show you how to visualize your data.


We will be working on 3 types of somatic calls: 

 * SNV calls from MuTect (vcf)
 * SV calls from DELLY (vcf)
 * CNV calls from SCoNEs (tsv)


This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License. This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.


## The environement

We will use a dataset derived from the analysis of whole genome sequencing paired normal/tumour samples

The call files are contained in the folder visualization:

 * mutect.somatic.vcf
 * delly.somatic.vcf
 * scones.somatic.tsv
 


# Circular representation of your calls
Many tools are available to do this the most common know is circos. But circos is a really not user friendly. In this tutoriel we show you an easy alternative to build circular representation of genomic data.

First we nee to go in the folder to do the analysis

```{.bash}
cd visualization
```

Let see what is in this folder

```{.bash}
ls
```

circos.R delly.somatic.vcf mutect.somatic.vcf scones.somatic.tsv



Take a look of the data files.

These are data of the are not restricted to a short piece of the chromosome.

SNVs have already been filtered 


```{.bash}
less delly.somatic.vcf
less mutect.somatic.vcf
less scones.somatic.tsv
```


**What can you see fron this data ?**
[solution](solutions/_data1.md)

**Why don't we use the vcf format for all type of call?**
[solution](solutions/_data2.md)


The analysis will be done using the R program

```{.bash}
R
```

We will use the circlize package from the cran R project. This package is dedicated to generate circular plot and had the advantage to provide pre-build function for genomics data. One of the main advantage of this tools is the use of bed format as input data.


```{.bash}
library(circlize)
```

Let's import the variants


```{.bash}
snp=read.table("mutect.somatic.vcf")
sv=read.table("somatic.sv.vcf")
cnv=read.table("data/scones.somatic.tsv",header=T)
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

Unfortunately circlize does not support hg38 yet. So we will need to reformat our data to fit the hg19 standards
As we work only on autosomes we won't need to lift-over and we could simply add **chr** at the begin of the chromosome names


We can now draw 1 track for somatic mutations

```{.bash}
snv_tmp=read.table("data/mutec.somatic.vcf",comment.char="#")
snv=cbind(paste("chr",as.character(snp[,1]),sep=""),snp[2],snp[,2]+1)
circos.genomicTrackPlotRegion(snv,stack=TRUE, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, cex = 0.05, pch = 9,col='orange' , ...)
})
```


Let's draw the 2 tracks for cnvs. One track for duplication in red and one blue track for deletion.

```{.bash}
dup=cnv[cnv[,5]>2,]
dup[,1]=paste("chr",as.character(dup[,1]),sep="")
del=cnv[cnv[,5]<2,]
del[,1]=paste("chr",as.character(del[,1]),sep="")
circos.genomicTrackPlotRegion(dup, stack = TRUE,panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "red",bg.border = NA, cex=1 , ...)
})
circos.genomicTrackPlotRegion(del, stack = TRUE,panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "blue",bg.border = NA, cex=1 , ...)
})
```

We can cleary see a massive deletion in the chromosome 3.


To finish we just need to draw 3 tracks + positional links to represent SVs

Unfortunately the vcf format has not been designed for SVs. SVs are defined by 2 breakpoints and the vcf format store the second one in the info field. So we will need to extract this information to draw these calls.


```{.bash}
chrEnd=NULL
posEnd=NULL
for (i in 1:dim(sv)[1]) {
    addInfo=strsplit(as.character(sv[i,8]),split=";")
    chrInf=strsplit(addInfo[[1]][3],split="=")
    chrEnd=c(chrEnd,chrInf[[1]][2])
    posInf=strsplit(addInfo[[1]][4],split="=")
    posEnd=c(posEnd,posInf[[1]][2])
}
svTable=data.frame(paste("chr",sv[,1],sep=""),as.numeric(sv[,2]),as.numeric(posEnd),paste("chr",chrEnd,sep=""),as.character(sv[,5]))
```

Now that we reformat the SV calls, let's draw them


```{.bash}
typeE=c("<DEL>","<INS>","<INV>")
colE=c("blue","black","green")
for (i in 1:3) { 
        bed_list=svTable[svTable[,5]==typeE[i],]
        circos.genomicTrackPlotRegion(bed_list,stack=TRUE, panel.fun = function(region, value, ...) {
                circos.genomicPoints(region, value, cex = 0.5, pch = 16, col = colE[i], ...)
        })
}

bed1=cbind(svTable[svTable[,5]=="<TRA>",1:2],svTable[svTable[,5]=="<TRA>",2]+5)
bed2=cbind(svTable[svTable[,5]=="<TRA>",c(4,3)],svTable[svTable[,5]=="<TRA>",3]+5)

for (i in 1:dim(bed1)[1]) {
        circos.link(bed1[i,1],bed1[i,2],bed2[i,1],bed2[i,2])
}
```

A good graph needs title and legends

```{.bash}
title("Somatic calls (SNV - SV - CNV)")
legend(0.7,1.4,legend=c("SNV", "CNV-DUPLICATION","CNV-DELETION","SV-DELETION","SV-INSERTION","SV-INVERSION"),col=c("orange","red","blue","blue","black","green","red"),pch=c(16,15,15,16,16,16,16,16),cex=0.75,title="Tracks:",bty='n')
legend(0.6,0.95,legend="SV-TRANSLOCATION",col="black",lty=1,cex=0.75,lwd=1.2,bty='n')
```

you should obtain a plot like this one
![circular_view](img/circos.png)

Exercice:


**Generate the graph and save it into a pdf file** 
[solution](solutions/_image1.md)


Finally exit R

```{.bash}
q("yes")
```

## Aknowledgments
This tutorial is an adaptation of the one created by Louis letourneau [here](https://github.com/lletourn/Workshops/tree/ebiCancerWorkshop201407doc/01-SNVCalling.md). I would like to thank and acknowledge Louis for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta. I also want to acknowledge Joel Fillon, Louis Letrouneau (again), Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.

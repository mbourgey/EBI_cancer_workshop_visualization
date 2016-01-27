# Introduction to cancer genomics data visualization
***By Mathieu Bourgey, Ph.D***

## Introduction 
A major goal of cancer genomics analysis is to identify genetic changes which lead to the tumor evolution.

In a more technical view we ends-up with different types of genomics changes
Point mutations
Copy number variations
Structural variants (insertion, inversion, deletion, duplication, translocation)

A common visualization of these data is to give an overview of change for the 23 chromosome of the genome in circular view.


![PNAS_circos](img/circos-cancer-cell.png)


from PNAS journal
   

## Circular representation of your calls
Many tools are available to do this the most common know is circos. But circos is a really not user friendly. In this tutoriel we show you an easy alternative to build circular representation of genomic data. to do that we will use the circlize R package

## General presentation of circlize

Author: Zuguang Gu

Gu Z et. al. (2014) circlize implements and enhances circular visualization in R. Bioinformatics.


Circular layout is very useful to represent complicated information, especially for genomic data. It
has advantages to visualize data with long axes or large amount of categories, described with different
measurements. It is also effective to visualize relations between elements.

### Principle of design
Since most of the figures are composed of simple graphics, such as points, lines, polygon (for filled
colors) et al, circlize implements low-level graphic functions for adding graphics in circular layout, so
that more higher level graphics can be easily comprised by low-level graphics. This principle ensures
the generality that types of high-level graphics are not restricted by the software but determined by
users.

Currently there are following graphic functions that can be used for plotting, they are similar to
the functions without “circos.” prefix from the traditional graphic engine :

• circos.points: add points in a cell, similar as points.
• circos.lines: add lines in a cell, similar as lines.
• circos.rect: add rectangle in a cell, similar as rect.
• circos.polygon: add polygon in a cell, similar as polygon.
• circos.text: add text in a cell, similar as text.
• circos.axis: add axis in a cell, functionally similar as axis but with more features.
• circos.link: this maybe the unique feature for circular layout to represent relationships between
elements.

For adding points, lines and text in cells through the whole track (among several sectors), the
following functions are available:
• circos.trackPoints: this can be replaced by circos.points through a for loop.
• circos.trackLines: this can be replaced by circos.lines through a for loop.

For Genomics data the author provides a set of predefine functions (circos.genomic...) 


## Genomics data representation
First we need to get the data

```{.bash}
git clone git@github.com:mbourgey/EBI_cancer_workshop_visualization.git
```

Let see what is in the data folder

```{.bash}
ls data
```

- breakdancer.somatic.tsv
- mutec.somatic.vcf
- scones.somatic.30k.tsv


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

Let's draw human genome hg19 reference ideograms

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

## A generic circular plot for non-genomic data

this is an example coming from the main vignette of the package

Following is an example to show the basic feature and usage of circlize package. First let’s generate
some random data. 
There needs a factor to represent categories, values on x-axis, and values on y-axis.

```{.bash}
set.seed(999)
n = 1000
a = data.frame(factor = sample(letters[1:8], n, replace = TRUE),
x = rnorm(n), y = runif(n))


### First initialize the layout 

we need to reset the graphic parameters and internal variables, so that it will not mess up
your next plot. 

In this step, circos.initialize allocates sectors in the circle according
to ranges of x-values in different categories. 

```{.bash}
circos.clear()
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
```

### Draw the first track. 

Before drawing any track we need to know that all tracks should
firstly be created by circos.trackPlotRegion, then those low-level functions can be applied. Since x-lims for cells in the track have already been defined
in the initialization step, here we only need to specify the y-lim for each cell, either by y or ylim
argument.

We also add axes in the first track, The axis for each cell is added by panel.fun argument.
circos.trackPlotRegion creates plotting region cell by cell and the panel.fun is actually executed
immediately after the creation of the plotting region for a certain cell. So panel.fun actually means
adding graphics in the “current cell”. After that, we add points through the whole track by circos.trackPoints.
Finally, add two texts in a certain cell (the cell is specified by sector.index and track.index argument).
When adding the second text, we do not specify track.index because the package knows we
are now in the first track.


```{.bash}
circos.clear()
circos.trackPlotRegion(factors = a$factor, y = a$y,panel.fun = function(x, y) {
	circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
```


### Draw the second track.

We use circos.trackHist to add histograms in the track. The
function also creates a new track because drawing histogram is really high level, so we do not need to
call circos.trackPlotRegion here. 


```{.bash}
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
```

### Draw the third track. 

Here some meta data for the current cell can be obtained by
get.cell.meta.data. This function needs sector.index and track.index arguments, and if they are
not specified, it means it is the current sector index and the current track index.

```{.bash}
circos.trackPlotRegion(factors = a$factor, x = a$x, y = a$y,
panel.fun = function(x, y) {
grey = c("#FFFFFF", "#CCCCCC", "#999999")
sector.index = get.cell.meta.data("sector.index")
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), sector.index)
circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)
})
```


### Draw the fourth track. 

Here you can choose different line types which is similar as
type argument in lines.


```{.bash}
circos.trackPlotRegion(factors = a$factor, y = a$y)
circos.trackLines(a$factor[1:100], a$x[1:100], a$y[1:100], type = "h")
```


###Draw links. 

Links can be from point to point, point to interval or interval to interval.


```{.bash}
circos.link("a", 0, "b", 0, h = 0.4)
circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
border = "blue", h = 0.2)
circos.link("e", 0, "g", c(-1,1), col = "green", lwd = 2, lty = 2)
```

### get information about your plot
You can get a summary of your circular layout by circos.info.

```{.bash}
circos.info()
circos.info(sector.index = "a", track.index = 2)
```

## Aknowledgments
Louis Letourneau
Toby Hockings
Guillaume Bourque
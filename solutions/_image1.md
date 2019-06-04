Just redo the same analysis but previously open the connection to the pdf file and close it when the analysis is done.


```{.bash}
pdf(file="somatic_circular_plot.pdf", height=8, width=8, compress=TRUE)

...


dev.off()
```

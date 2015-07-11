The 1000 Genomes project try to use/include SVs call in the [vcf format](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20%28Variant%20Call%20Format%29%20version%204.0/encoding-structural-variants). This a good idea to try to include everything altogether but, to my point of view, this not the best way to handle SV and CNV.  


Why ?  


Due to the nature of these calls, you can not easily integrate the postional information of the two breakpoints (that could be located faraway or in an other chormosome) using a single position format. 



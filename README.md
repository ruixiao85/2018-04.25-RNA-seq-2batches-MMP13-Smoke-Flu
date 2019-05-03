---
title: "RNA-seq MMP13 Smoke Flu"
author: "Rui Xiao"
date: "May 3, 2019"
output: 
   html_document:
      keep_md: true
      toc: true
      number_sections: true
      toc_float: true
      toc_depth: 2
      df_print: paged
      theme: united
      highlight: textmate
---



# Read the raw count-table and sample annotation files

```r
if (!requireNamespace("BiocManager",quietly=T)) install.packages("BiocManager")
library(BiocManager)

c1=read.table("counts_table_raw.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
c2=read.table("counts_table_raw set 2.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
c=cbind(c1,c2) # knowing that the row names are matching perfectly
rm(c1,c2) # free memory if no longer needed
dim(c)
```

```
## [1] 23420    24
```

```r
head(c) # count table overview
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["JK001"],"name":[1],"type":["int"],"align":["right"]},{"label":["JK004"],"name":[2],"type":["int"],"align":["right"]},{"label":["JK005"],"name":[3],"type":["int"],"align":["right"]},{"label":["JK011"],"name":[4],"type":["int"],"align":["right"]},{"label":["JK012"],"name":[5],"type":["int"],"align":["right"]},{"label":["JK014"],"name":[6],"type":["int"],"align":["right"]},{"label":["JK017"],"name":[7],"type":["int"],"align":["right"]},{"label":["JK018"],"name":[8],"type":["int"],"align":["right"]},{"label":["JK021"],"name":[9],"type":["int"],"align":["right"]},{"label":["JK022"],"name":[10],"type":["int"],"align":["right"]},{"label":["JK023"],"name":[11],"type":["int"],"align":["right"]},{"label":["JK028"],"name":[12],"type":["int"],"align":["right"]},{"label":["JK031"],"name":[13],"type":["int"],"align":["right"]},{"label":["JK036"],"name":[14],"type":["int"],"align":["right"]},{"label":["JK039"],"name":[15],"type":["int"],"align":["right"]},{"label":["JK040"],"name":[16],"type":["int"],"align":["right"]},{"label":["JK215"],"name":[17],"type":["int"],"align":["right"]},{"label":["JK217"],"name":[18],"type":["int"],"align":["right"]},{"label":["JK218"],"name":[19],"type":["int"],"align":["right"]},{"label":["JK219"],"name":[20],"type":["int"],"align":["right"]},{"label":["JK228"],"name":[21],"type":["int"],"align":["right"]},{"label":["JK229"],"name":[22],"type":["int"],"align":["right"]},{"label":["JK230"],"name":[23],"type":["int"],"align":["right"]},{"label":["JK231"],"name":[24],"type":["int"],"align":["right"]}],"data":[{"1":"2","2":"0","3":"0","4":"0","5":"5","6":"2","7":"1","8":"1","9":"0","10":"0","11":"0","12":"5","13":"0","14":"0","15":"0","16":"0","17":"15","18":"7","19":"5","20":"10","21":"9","22":"11","23":"5","24":"15","_rn_":"Xkr4"},{"1":"76","2":"5","3":"89","4":"81","5":"34","6":"141","7":"100","8":"43","9":"25","10":"248","11":"173","12":"18","13":"106","14":"131","15":"2","16":"4","17":"154","18":"305","19":"56","20":"122","21":"20","22":"31","23":"162","24":"99","_rn_":"Rp1"},{"1":"1119","2":"912","3":"1169","4":"417","5":"1401","6":"433","7":"554","8":"201","9":"483","10":"355","11":"113","12":"427","13":"440","14":"612","15":"327","16":"616","17":"2046","18":"2288","19":"1202","20":"1614","21":"675","22":"663","23":"1514","24":"982","_rn_":"Sox17"},{"1":"678","2":"539","3":"586","4":"465","5":"765","6":"395","7":"428","8":"783","9":"612","10":"493","11":"515","12":"1110","13":"489","14":"428","15":"1160","16":"998","17":"926","18":"1180","19":"522","20":"943","21":"1101","22":"968","23":"1000","24":"717","_rn_":"Mrpl15"},{"1":"1434","2":"969","3":"1378","4":"963","5":"1746","6":"856","7":"719","8":"1280","9":"1200","10":"1190","11":"1173","12":"1788","13":"1110","14":"931","15":"1270","16":"1145","17":"2791","18":"3623","19":"1494","20":"2783","21":"2574","22":"2407","23":"2974","24":"2216","_rn_":"Lypla1"},{"1":"877","2":"866","3":"785","4":"688","5":"938","6":"593","7":"685","8":"739","9":"692","10":"772","11":"817","12":"1196","13":"813","14":"654","15":"1051","16":"857","17":"1773","18":"2519","19":"1132","20":"1853","21":"1705","22":"1694","23":"1841","24":"1354","_rn_":"Tcea1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
a=read.csv("SampleAnnotation2Batches.csv",header=T,stringsAsFactors=F,row.names=1)
a$SingleFactor=paste(a$Batch,a$MMP13,a$SMOKE,a$FLU,sep="_") # merge factors into one for some further analysis
head(a) # sample annotation overview
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Batch"],"name":[1],"type":["int"],"align":["right"]},{"label":["MMP13"],"name":[2],"type":["chr"],"align":["left"]},{"label":["SMOKE"],"name":[3],"type":["chr"],"align":["left"]},{"label":["FLU"],"name":[4],"type":["chr"],"align":["left"]},{"label":["DOB"],"name":[5],"type":["chr"],"align":["left"]},{"label":["SAC_DATE"],"name":[6],"type":["chr"],"align":["left"]},{"label":["SingleFactor"],"name":[7],"type":["chr"],"align":["left"]}],"data":[{"1":"1","2":"KO","3":"RA","4":"PBS","5":"2/27/2017","6":"6/26/2017","7":"1_KO_RA_PBS","_rn_":"JK001"},{"1":"1","2":"KO","3":"RA","4":"PBS","5":"2/27/2017","6":"6/26/2017","7":"1_KO_RA_PBS","_rn_":"JK004"},{"1":"1","2":"KO","3":"RA","4":"PBS","5":"2/27/2017","6":"6/26/2017","7":"1_KO_RA_PBS","_rn_":"JK005"},{"1":"1","2":"WT","3":"SM","4":"PBS","5":"2/27/2017","6":"6/29/2017","7":"1_WT_SM_PBS","_rn_":"JK011"},{"1":"1","2":"KO","3":"RA","4":"FLU","5":"2/27/2017","6":"6/30/2017","7":"1_KO_RA_FLU","_rn_":"JK012"},{"1":"1","2":"WT","3":"SM","4":"PBS","5":"2/27/2017","6":"6/29/2017","7":"1_WT_SM_PBS","_rn_":"JK014"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

# Preprocess the input data

```r
if (!require("edgeR",quietly=T)) BiocManager::install("edgeR")
library(edgeR)
y=DGEList(counts=c,samples=a)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
head(y$samples[,2:3])
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["lib.size"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["norm.factors"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"20100075","2":"0.9927551","_rn_":"JK001"},{"1":"13444148","2":"0.9744933","_rn_":"JK004"},{"1":"19755033","2":"1.0433833","_rn_":"JK005"},{"1":"17655062","2":"0.9687881","_rn_":"JK011"},{"1":"21767921","2":"1.0802167","_rn_":"JK012"},{"1":"14472144","2":"1.0178222","_rn_":"JK014"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
logcpm <- cpm(y, prior.count=2, log=TRUE)
logfile="logcpm.csv"
if (!file.exists(logfile)) write.csv(logcpm,logfile)

# setwd("2018 results") # allow results to be written into the subfolder
```

# Principal Component Analysis

```r
if (!require("FactoMineR",quietly=T)) install.packages("FactoMineR")
library(FactoMineR)
al=cbind(a,t(logcpm))
colnames(al)[1:10] # annotaion goes from 1st to 7th column
```

```
##  [1] "Batch"        "MMP13"        "SMOKE"        "FLU"         
##  [5] "DOB"          "SAC_DATE"     "SingleFactor" "Rp1"         
##  [9] "Sox17"        "Mrpl15"
```

```r
# multiple factors
m.pca=PCA(al[,-c(5:7)],scale.unit=T,ncp=5,quali.sup=1:4,graph=F)
plotellipses(m.pca)
```

![](README_files/figure-html/pca-1.png)<!-- -->

```r
# single factor
s.pca=PCA(al[,-c(1:6)],scale.unit=T,ncp=5,quali.sup=1,graph=F)
# plot.PCA(s.pca,axes=c(1,2),choix="ind",habillage = 1)
plotellipses(s.pca)
```

![](README_files/figure-html/pca-2.png)<!-- -->

# Fit generalized linear model using edgeR

```r
# design=model.matrix(~MMP13+SMOKE+FLU,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU,data=a)
design=model.matrix(~MMP13*SMOKE*FLU+Batch,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU*Batch,data=a); design=design[,-16]
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

allList=function(cf){ # select complete unordered test result regarding each coefficient
   dt=as.data.frame(topTags(glmQLFTest(fit,coef=cf),n=999999,sort.by="none"))
   colnames(dt)=paste0(colnames(design)[cf],"_",colnames(dt))
   return(dt)
}
write.csv(cbind(allList(2),allList(3),allList(4),allList(5),allList(6),allList(7),allList(8)),"GLM_All.csv") # combine and write to csv file
```

# Output filtered gene counts and/or their names to console or data files

```r
# functions that select stats on each factor
filterCount=function(cf){
   cat(paste0(colnames(design)[cf]," Counts:\n"))
   return(nrow(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),abs(logFC)>log2(2)&FDR<0.05)))
}
filterTop=function(cf,top=30){
   cat(paste0(colnames(design)[cf]," top",top," :\n"))
   return(topTags(glmQLFTest(fit, coef=cf),n=top))
}
filterName=function(cf){
   cat(paste0(colnames(design)[cf]," Name:\n"))
   return(row.names(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),abs(logFC)>log2(2)&FDR<0.05)))
   
}
filterWrite=function(cf){
   cat(paste0(colnames(design)[cf]," written to file:\n"))
   write.csv(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),abs(logFC)>log2(2)&FDR<0.05),
             paste0("",gsub(":","+",colnames(design)[cf]),".csv"))
}

for (i in 2:ncol(design)){ # uncomment following lines to enable the outputs
   cat(filterCount(i),end="\n")
   # print(filterTop(i))
   # print(filterName(i))
   # filterWrite(i)
}
```

```
## MMP13WT Counts:
## 896 
## SMOKESM Counts:
## 68 
## FLUPBS Counts:
## 0 
## Batch Counts:
## 1674 
## MMP13WT:SMOKESM Counts:
## 1 
## MMP13WT:FLUPBS Counts:
## 3 
## SMOKESM:FLUPBS Counts:
## 699 
## MMP13WT:SMOKESM:FLUPBS Counts:
## 2
```

# Plot heatmap with top genes found significant for each factor

```r
if (!require("pheatmap",quietly=T)) install.packages("pheatmap")
library(pheatmap)

order=a[order(a$Batch,a$MMP13,a$SMOKE,a$FLU),] # ordered annotation
logcpm=logcpm[,row.names(order)] # ordered logcpm

plotPheatmap=function(cf,clus_col=T,pdf=T){
   d=as.data.frame(logcpm[row.names(filterTop(cf)),])
   outfile=gsub(":","+",paste("heatmap","coef",i,colnames(design)[i],filterCount(cf),sep="_"))
   if (pdf) pdf(paste0(outfile,".pdf"),width=9,height=6)
   pheatmap(d,scale="row",annotation_col=order[,c(1,2:4)],cluster_cols=clus_col,main=outfile)
   if (pdf) dev.off()
}
for (i in 2:ncol(design)){
   plotPheatmap(i,pdf=F)
}
```

```
## MMP13WT top30 :
## MMP13WT Counts:
```

![](README_files/figure-html/heatmap-1.png)<!-- -->

```
## SMOKESM top30 :
## SMOKESM Counts:
```

![](README_files/figure-html/heatmap-2.png)<!-- -->

```
## FLUPBS top30 :
## FLUPBS Counts:
```

![](README_files/figure-html/heatmap-3.png)<!-- -->

```
## Batch top30 :
## Batch Counts:
```

![](README_files/figure-html/heatmap-4.png)<!-- -->

```
## MMP13WT:SMOKESM top30 :
## MMP13WT:SMOKESM Counts:
```

![](README_files/figure-html/heatmap-5.png)<!-- -->

```
## MMP13WT:FLUPBS top30 :
## MMP13WT:FLUPBS Counts:
```

![](README_files/figure-html/heatmap-6.png)<!-- -->

```
## SMOKESM:FLUPBS top30 :
## SMOKESM:FLUPBS Counts:
```

![](README_files/figure-html/heatmap-7.png)<!-- -->

```
## MMP13WT:SMOKESM:FLUPBS top30 :
## MMP13WT:SMOKESM:FLUPBS Counts:
```

![](README_files/figure-html/heatmap-8.png)<!-- -->

# Plot significant gene counts of each major factor in a Venn diagram

```r
if (!require("systemPipeR",quietly=T)) BiocManager::install("systemPipeR")
library(systemPipeR)

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
           Smoke=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
           SmokeFlu=row.names(topTags(glmQLFTest(fit, coef=7),p.value=0.05,n=999999)))

vennPlot(overLapper(setlist,type="vennsets"))
```

![](README_files/figure-html/venn-1.png)<!-- -->

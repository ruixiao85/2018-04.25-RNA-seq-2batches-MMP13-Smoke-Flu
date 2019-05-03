---
title: "RNA-seq MMP13 Smoke Flu"
author: "Rui Xiao"
date: "May 2, 2019"
output: 
  html_document:
    keep_md: true
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

```
##        JK001 JK004 JK005 JK011 JK012 JK014 JK017 JK018 JK021 JK022 JK023
## Xkr4       2     0     0     0     5     2     1     1     0     0     0
## Rp1       76     5    89    81    34   141   100    43    25   248   173
## Sox17   1119   912  1169   417  1401   433   554   201   483   355   113
## Mrpl15   678   539   586   465   765   395   428   783   612   493   515
## Lypla1  1434   969  1378   963  1746   856   719  1280  1200  1190  1173
## Tcea1    877   866   785   688   938   593   685   739   692   772   817
##        JK028 JK031 JK036 JK039 JK040 JK215 JK217 JK218 JK219 JK228 JK229
## Xkr4       5     0     0     0     0    15     7     5    10     9    11
## Rp1       18   106   131     2     4   154   305    56   122    20    31
## Sox17    427   440   612   327   616  2046  2288  1202  1614   675   663
## Mrpl15  1110   489   428  1160   998   926  1180   522   943  1101   968
## Lypla1  1788  1110   931  1270  1145  2791  3623  1494  2783  2574  2407
## Tcea1   1196   813   654  1051   857  1773  2519  1132  1853  1705  1694
##        JK230 JK231
## Xkr4       5    15
## Rp1      162    99
## Sox17   1514   982
## Mrpl15  1000   717
## Lypla1  2974  2216
## Tcea1   1841  1354
```

```r
a=read.csv("SampleAnnotation2Batches.csv",header=T,stringsAsFactors=F,row.names=1)
a$SingleFactor=paste(a$Batch,a$MMP13,a$SMOKE,a$FLU,sep="_") # merge factors into one for some further analysis
head(a) # sample annotation overview
```

```
##       Batch MMP13 SMOKE FLU       DOB  SAC_DATE SingleFactor
## JK001     1    KO    RA PBS 2/27/2017 6/26/2017  1_KO_RA_PBS
## JK004     1    KO    RA PBS 2/27/2017 6/26/2017  1_KO_RA_PBS
## JK005     1    KO    RA PBS 2/27/2017 6/26/2017  1_KO_RA_PBS
## JK011     1    WT    SM PBS 2/27/2017 6/29/2017  1_WT_SM_PBS
## JK012     1    KO    RA FLU 2/27/2017 6/30/2017  1_KO_RA_FLU
## JK014     1    WT    SM PBS 2/27/2017 6/29/2017  1_WT_SM_PBS
```

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

```
##       lib.size norm.factors
## JK001 20100075    0.9927551
## JK004 13444148    0.9744933
## JK005 19755033    1.0433833
## JK011 17655062    0.9687881
## JK012 21767921    1.0802167
## JK014 14472144    1.0178222
```

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

# Functions that selectively output stats for each factor in the model

```r
# stats on each factor
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
```

## Output filtered gene counts and/or their names to console or data files

```r
for (i in 2:ncol(design)){
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

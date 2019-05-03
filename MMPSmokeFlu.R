if (!requireNamespace("BiocManager",quietly=T)) install.packages("BiocManager")

c1=read.table("counts_table_raw.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
c2=read.table("counts_table_raw set 2.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
c=cbind(c1,c2) # knowing that the row names are matching perfectly
rm(c1,c2) # free memory if no long needed
dim(c)
head(c) # count table overview

a=read.csv("SampleAnnotation2Batches.csv",header=T,stringsAsFactors=F,row.names=1)
a$SingleFactor=paste(a$Batch,a$MMP13,a$SMOKE,a$FLU,sep="_") # join factors
head(a) # sample annotation overview

if (!require("edgeR",quietly=T)) BiocManager::install("edgeR")
library(edgeR)
y=DGEList(counts=c,samples=a)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

# design=model.matrix(~MMP13+SMOKE+FLU,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU,data=a)
design=model.matrix(~MMP13*SMOKE*FLU+Batch,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU*Batch,data=a); design=design[,-16]

y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

logcpm <- cpm(y, prior.count=2, log=TRUE)
logfile="logcpm.csv"
if (!file.exists(logfile)) write.csv(logcpm,logfile)

setwd("2018 results") # allow results to write into the subfolder

if (!require("FactoMineR",quietly=T)) install.packages("FactoMineR")
library(FactoMineR)

al=cbind(a,t(logcpm))
colnames(al)[1:10] # annotaion goes 1 : 7 column

# multiple factors
m.pca=PCA(al[,-c(5:7)],scale.unit=T,ncp=5,quali.sup=1:4,graph=F)
plotellipses(m.pca)

s.pca=PCA(al[,-c(1:6)],scale.unit=T,ncp=5,quali.sup=1,graph=F)
# plot.PCA(s.pca,axes=c(1,2),choix="ind",habillage = 1)
plotellipses(s.pca)

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
for (i in 2:ncol(design)){
   cat(filterCount(i),end="\n")
   # print(filterTop(i))
   print(filterName(i))
   # filterWrite(i)
}


allList=function(cf){
   dt=as.data.frame(topTags(glmQLFTest(fit,coef=cf),n=999999,sort.by="none"))
   colnames(dt)=paste0(colnames(design)[cf],"_",colnames(dt))
   return(dt)
}
write.csv(cbind(allList(2),allList(3),allList(4),allList(5),allList(6),allList(7),allList(8)),"GLM_All.csv")

write.csv(cbind(allList(2),allList(3),allList(7)),"GLM_3_MMP13_Smoke_SmokeFlu.csv")


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


if (!require("systemPipeR",quietly=T)) BiocManager::install("systemPipeR")
library(systemPipeR)

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
           Smoke=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
           SmokeFlu=row.names(topTags(glmQLFTest(fit, coef=7),p.value=0.05,n=999999)))

pdf("venn.pdf")
vennPlot(overLapper(setlist,type="vennsets"))
dev.off()


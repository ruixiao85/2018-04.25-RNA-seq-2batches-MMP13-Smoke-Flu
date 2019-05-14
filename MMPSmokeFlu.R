
if (!requireNamespace("BiocManager",quietly=T)) install.packages("BiocManager")
library(BiocManager)

c1=read.table("counts_table_raw.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
c2=read.table("counts_table_raw set 2.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
c=cbind(c1,c2) # knowing that the row names are matching perfectly
rm(c1,c2) # free memory if no longer needed
dim(c)
head(c) # count table overview
a=read.csv("SampleAnnotation2Batches.csv",header=T,stringsAsFactors=F,row.names=1)
a$SingleFactor=paste(a$Batch,a$MMP13,a$SMOKE,a$FLU,sep="_") # merge factors into one for some further analysis
a$Batch=factor(a$Batch,levels=c(1,2))
a$MMP13=factor(a$MMP13,levels=c("WT","KO"))
a$FLU=factor(a$FLU,levels=c("PBS","FLU"))
a$SMOKE=factor(a$SMOKE,levels=c("RA","SM"))
head(a) # sample annotation overview

if (!require("edgeR",quietly=T)) BiocManager::install("edgeR")
library(edgeR)
y=DGEList(counts=c,samples=a,group=a$SingleFactor)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
head(y$samples[,2:3])

logcpm <- cpm(y, prior.count=2, log=TRUE)
logfile="logcpm.csv"
if (!file.exists(logfile)) write.csv(logcpm,logfile)

# setwd("2018 results") # allow results to be written into the subfolder

if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)


if (!require("FactoMineR",quietly=T)) install.packages("FactoMineR")
library(FactoMineR)
al=cbind(a,t(logcpm))
colnames(al)[1:10] # annotaion goes from 1st to 7th column
# multiple factors
m.pca=PCA(al[,-c(5:7)],scale.unit=T,ncp=5,quali.sup=1:4,graph=F)
plotellipses(m.pca)
# single factor
s.pca=PCA(al[,-c(1:6)],scale.unit=T,ncp=5,quali.sup=1,graph=F)
# plot.PCA(s.pca,axes=c(1,2),choix="ind",habillage = 1)
plotellipses(s.pca)


y1 <- estimateDisp(y,model.matrix(~SingleFactor,data=a))
pair=c("2_WT_RA_PBS","1_KO_RA_PBS") #  c("A","B") = B - A

tt=topTags(exactTest(y1, pair),n=9999999,sort.by="none") # select all with original order
# write.csv(tt,paste0("Pair_",paste(pair,collapse=" vs "),".csv")) # save to local file

attach(tt$table)
tt$table$sig=ifelse(PValue<0.05&abs(logFC)>1,logFC,0)
p = ggplot(tt$table,aes(logFC,-log10(PValue)))+
  geom_point(aes(col=sig))+
  scale_colour_gradient2(low="blue",mid="black",high="red")
p
detach(tt$table)


# design=model.matrix(~MMP13+SMOKE+FLU,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU,data=a)
 design=model.matrix(~MMP13*SMOKE*FLU+Batch,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU*Batch,data=a); design=design[,-c(13:16)]
y2 <- estimateDisp(y,design)
fit <- glmQLFit(y2, design)

allList=function(cf){ # select complete unordered test result regarding each coefficient
   dt=as.data.frame(topTags(glmQLFTest(fit,coef=cf),n=999999,sort.by="none"))
   colnames(dt)=paste0(colnames(design)[cf],"_",colnames(dt))
   return(dt)
}
# write.csv(cbind(allList(2),allList(3),allList(4),allList(5),allList(6),allList(7),allList(8)),"GLM_All.csv") # combine and write to csv file


# functions that select stats on each factor
ifSig=function(logFC,PValue,FDR){
   # return(abs(logFC)>=1&PValue<0.05)
   return(abs(logFC)>=1&FDR<0.05)
   # return(abs(logFC)>=0.5849&FDR<0.1)
}
filterCount=function(cf){
   cat(paste0(colnames(design)[cf]," Counts:\n"))
   return(nrow(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),ifSig(logFC,PValue,FDR))))
}
filterTop=function(cf,top=24){
   cat(paste0(colnames(design)[cf]," top",top," :\n"))
   return(topTags(glmQLFTest(fit, coef=cf),n=top))
}
filterName=function(cf){
   cat(paste0(colnames(design)[cf]," Name:\n"))
   return(row.names(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),ifSig(logFC,PValue,FDR))))
   
}
filterWrite=function(cf){
   cat(paste0(colnames(design)[cf]," written to file:\n"))
   write.csv(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),ifSig(logFC,PValue,FDR)),
             paste0("",gsub(":","+",colnames(design)[cf]),".csv"))
}

for (i in 2:ncol(design)){ # uncomment following lines to enable the outputs
   cat(filterCount(i),end="\n")
   # print(filterTop(i))
   # print(filterName(i))
   # filterWrite(i)
}


if (!require("pheatmap",quietly=T)) install.packages("pheatmap")
library(pheatmap)

order=a[order(a$Batch,a$MMP13,a$SMOKE,a$FLU),] # ordered annotation
logcpm=logcpm[,row.names(order)] # ordered logcpm

plotPheatmap=function(cf,clus_col=T,pdf=T){
   d=as.data.frame(logcpm[row.names(filterTop(cf)),])
   outfile=gsub(":","+",paste("heatmap","coef",i,colnames(design)[i],filterCount(cf),sep="_"))
   if (pdf) pdf(paste0(outfile,".pdf"),width=9,height=6)
   pheatmap(d,scale="row",annotation_col=order[,1:4],cluster_cols=clus_col,main=outfile)
   if (pdf) dev.off()
}
for (i in 2:ncol(design)){
   plotPheatmap(i,pdf=F)
}


if (!require("systemPipeR",quietly=T)) BiocManager::install("systemPipeR")
library(systemPipeR)

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
           Smoke=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
           Flu=row.names(topTags(glmQLFTest(fit, coef=4),p.value=0.05,n=999999)))
vennPlot(overLapper(setlist,type="vennsets"))

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
           Smoke=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
           Flu=row.names(topTags(glmQLFTest(fit, coef=4),p.value=0.05,n=999999)),
           Batch=row.names(topTags(glmQLFTest(fit, coef=5),p.value=0.05,n=999999)))
vennPlot(overLapper(setlist,type="vennsets"))

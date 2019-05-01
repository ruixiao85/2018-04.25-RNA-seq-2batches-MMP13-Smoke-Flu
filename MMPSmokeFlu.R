setwd("D:/Cel files/2018-04.25 RNA-seq Kyle 2batches MMP13 Smoke Flu")

p=read.table("counts_table_raw.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
p=cbind(p,read.table("counts_table_raw set 2.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1))
dim(p)
head(p)
library(FactoMineR)
PCA(t(p))
dir()
a=read.csv("SampleNote2Batches.csv",stringsAsFactors=F,row.names=1)

pa=cbind(t(p),a[,3:5])
PCA(pa,quali.sup=c(23421:23423))

pa=cbind(t(p),Group=a[,2])
PCA(pa,quali.sup=23421)

head(a)
a$Batch=factor(a$Batch)
colnames(a)=c("Batch","SingleFactor","MMP13","SMOKE","FLU","DOB","SAC.DATE") 
# design=model.matrix(~MMP13+SMOKE+FLU,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU,data=a)
design=model.matrix(~MMP13*SMOKE*FLU+Batch,data=a)
# design=model.matrix(~MMP13*SMOKE*FLU*Batch,data=a); design=design[,-16]

library(edgeR)
y=DGEList(counts=p,samples=a)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

#need fit and design
filterCount=function(cf){
   cat(paste0(colnames(design)[cf]," Counts:"))
   n=nrow(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),abs(logFC)>log2(2)&FDR<0.05))
   print(n)
   return(n)
}
filterTop=function(cf,top=30){
   cat(paste0(colnames(design)[cf]," top",top," :"))
   return(topTags(glmQLFTest(fit, coef=cf),n=top))
}
filterName=function(cf){
   print(paste0(colnames(design)[cf]," Name:"))
   name=row.names(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),abs(logFC)>log2(2)&FDR<0.05))
   print(name)
   return(name)
}
filterWrite=function(cf){
   print(paste0(colnames(design)[cf]," Write to file:"))
   write.csv(subset(as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999)),abs(logFC)>log2(2)&FDR<0.05),
             paste0("",colnames(design)[cf],".csv"))
}
for (i in 2:ncol(design)){
   filterCount(i)
   # filterName(i)
   # filterWrite(i)
}


allList=function(cf){
   dt=as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999,sort.by = "none"))
   colnames(dt)=paste0(colnames(design)[cf],"_",colnames(dt))
   return(dt)
}
write.csv(cbind(allList(2),allList(3),allList(4),allList(5),allList(6),allList(7),allList(8)),"GLM_All.csv")

write.csv(cbind(allList(2),allList(3),allList(7)),"GLM_3_MMP13_Smoke_SmokeFlu.csv")

logcpm <- cpm(y, prior.count=2, log=TRUE)
# write.csv(logcpm,"logcpm.csv")


order=a[order(a$MMP13,a$SMOKE,a$ FLU),] # ordered annotation
logcpm=logcpm[,row.names(order)] # ordered logcpm

library(pheatmap)
library(grid)
plotPheatmap=function(cf){
   setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
   d=as.data.frame(logcpm[row.names(filterTop(cf)),])
   outfile=gsub(":",".",paste0("heatmap_coef",cf,"_",colnames(design)[cf],".pdf"))
   pdf(outfile,width=9,height=6)
   setHook("grid.newpage", NULL, "replace")
   pheatmap(d,scale="row",annotation_col=order[,c(1,3:5)],cluster_cols=F)
   grid.text(as.numeric(filterCount(cf)), y=-0.07, gp=gpar(fontsize=16))
   grid.text(outfile, x=-0.07, rot=90, gp=gpar(fontsize=16))
   dev.off()
}
for (i in 2:ncol(design)){
   plotPheatmap(i)
}


source("https://bioconductor.org/biocLite.R"); biocLite("systemPipeR")
library(systemPipeR)

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
           Smoke=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
           SmokeFlu=row.names(topTags(glmQLFTest(fit, coef=7),p.value=0.05,n=999999)))

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
             Smoke=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
             Flu=row.names(topTags(glmQLFTest(fit, coef=4),p.value=0.05,n=999999)),
             MMP13Smoke=row.names(topTags(glmQLFTest(fit, coef=5),p.value=0.05,n=999999)),
             MMP13Flu=row.names(topTags(glmQLFTest(fit, coef=6),p.value=0.05,n=999999)),
             SmokeFlu=row.names(topTags(glmQLFTest(fit, coef=7),p.value=0.05,n=999999)),
             MMP13SmokeFlu=row.names(topTags(glmQLFTest(fit, coef=8),p.value=0.05,n=999999)))

vennset <- overLapper(setlist,type="vennsets")
pdf("venn.pdf")
vennPlot(vennset)
dev.off()


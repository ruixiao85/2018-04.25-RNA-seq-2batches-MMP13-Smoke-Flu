setwd("D:/Cel files/2018-04.25 RNA-seq Kyle 2batches MMP13 Smoke Flu")

p=read.table("counts_table_raw.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
p=cbind(p,read.table("counts_table_raw set 2.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1))
dim(p)
head(p)

##PCA
library(FactoMineR)
PCA(t(p))
dir()
a=read.csv("SampleNote2Batches.csv",stringsAsFactors=F,row.names=1)

pa=cbind(t(p),a[,3:5])
PCA(pa,quali.sup=c(23421:23423))

pa=cbind(t(p),Group=a[,2])
PCA(pa,quali.sup=23421)

head(a)
a$SingleFactor=factor(a$SingleFactor)
a$Batch=factor(a$Batch)
a$MMP13.KO.WT=factor(a$MMP13.KO.WT,levels=c("WT","KO"))
a$SMOKE..ROOM.AIR=factor(a$SMOKE..ROOM.AIR,levels=c("RA","SM"))
a$FLU=factor(a$FLU,levels=c("PBS","FLU"))
colnames(a)=c("Batch","SingleFactor","MMP13","SMOKE","FLU","DOB","SAC.DATE") 

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
allList=function(cf){
   dt=as.data.frame(topTags(glmQLFTest(fit, coef=cf),n=999999,sort.by = "none"))
   colnames(dt)=paste0(colnames(design)[cf],"_",colnames(dt))
   return(dt)
}

library(pheatmap)
library(grid)
plotPheatmap=function(cf){
   setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
   d=as.data.frame(logcpm[row.names(filterTop(cf)),])
   outfile=gsub(":",".",paste0("heatmap_coef",cf,"_",colnames(design)[cf],".pdf"))
   setHook("grid.newpage", NULL, "replace")
   pheatmap(d,scale="row",annotation_col=a[colnames(logcpm),c(1,3:5)],cluster_cols=F)
   grid.text(as.numeric(filterCount(cf)), y=-0.07, gp=gpar(fontsize=16))
   grid.text(outfile, x=-0.07, rot=90, gp=gpar(fontsize=16))
}

library(edgeR)

sel=a$SMOKE=="RA"; design=model.matrix(~MMP13*FLU+Batch,data=a[sel,]); note="WhenSmokeRA"
sel=a$FLU=="PBS"; design=model.matrix(~MMP13*SMOKE+Batch,data=a[sel,]); note="WhenNoFlu"
sel=!((a$SMOKE=="RA"&a$FLU=="FLU")|(a$SMOKE=="SM"&a$FLU=="PBS")); design=model.matrix(~MMP13*FLU+Batch,data=a[sel,]); note="WhenSmokeFluComb"

y=DGEList(counts=p,samples=a)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples$group=y$samples$SingleFactor

y$counts=y$counts[,sel]
y$samples=y$samples[sel,]

y <- calcNormFactors(y)
y$samples
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

for (i in 2:ncol(design)){
   filterCount(i)
   # filterName(i)
   # filterWrite(i)
}

# write.csv(cbind(allList(2),allList(3),allList(4),allList(5)),paste0("GLM_2factors_",note,".csv"))
# write.csv(cbind(allList(2),allList(3),allList(4)),paste0("GLM_2factors_",note,".csv"))
# list=list() #first time only
# list=NULL #clear list
list[[note]]=cbind(allList(2),allList(3),allList(4),allList(5))


# write.csv(cbind(list$WhenSmokeRA,list$WhenNoFlu,list$WhenSmokeFluComb),"GLM_2Factor_3combined.csv",quote=F)
all=read.csv("GLM_2Factor_3combined.csv")
library(xlsx)
options(java.parameters = "-Xmx12G")
write.xlsx2(cbind(list$WhenSmokeRA,list$WhenNoFlu,list$WhenSmokeFluComb), file = "GLM_2Factor_3combined.xlsx",
           sheetName = "Sheet1", append = FALSE)


logcpm <- cpm(y, prior.count=2, log=TRUE)
# write.csv(logcpm,"logcpm.csv")

logcpm=logcpm[,order(a[sel,]$MMP13,a[sel,]$SMOKE,a[sel,]$FLU)] # ordered logcpm

pdf(paste0("Heatmap_GLM_2Factor_",note,".pdf"),height=6,width=9)
for (i in 2:ncol(design)){
   plotPheatmap(i)
}
dev.off()


## pairwise Only with exactTest
y <- calcNormFactors(y)
y$samples
y <- estimateDisp(y)
write.csv(topTags(exactTest(y, pair=c("WT RA PBS","KO RA PBS")),n=999999),"KO RA PBS vs WT RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("WT RA PBS","WT SM FLU")),n=999999),"WT SM FLU vs WT RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("WT RA PBS","KO SM FLU")),n=999999),"KO SM FLU vs WT RA PBS.csv")

write.csv(topTags(exactTest(y, pair=c("WT SM FLU","KO SM FLU")),n=999999),"KO SM FLU vs WT SM FLU.csv")
write.csv(topTags(exactTest(y, pair=c("WT SM PBS","KO SM FLU")),n=999999),"KO SM FLU vs WT SM PBS.csv")
write.csv(topTags(exactTest(y, pair=c("WT RA FLU","KO SM FLU")),n=999999),"KO SM FLU vs WT RA FLU.csv")

write.csv(topTags(exactTest(y, pair=c("WT RA PBS","WT SM PBS")),n=999999),"WT SM PBS vs WT RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("WT RA PBS","KO SM PBS")),n=999999),"KO SM PBS vs WT RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("WT RA PBS","WT RA FLU")),n=999999),"WT RA FLU vs WT RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("WT RA PBS","KO RA FLU")),n=999999),"KO RA FLU vs WT RA PBS.csv")

write.csv(topTags(exactTest(y, pair=c("KO RA PBS","WT SM PBS")),n=999999),"WT SM PBS vs KO RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("KO RA PBS","KO SM PBS")),n=999999),"KO SM PBS vs KO RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("KO RA PBS","WT RA FLU")),n=999999),"WT RA FLU vs KO RA PBS.csv")
write.csv(topTags(exactTest(y, pair=c("KO RA PBS","KO RA FLU")),n=999999),"KO RA FLU vs KO RA PBS.csv")

redf=function(baseline,affected){# b vs a
   dt=topTags(exactTest(y, pair=c(baseline,affected)),n=999999,sort="none")[[1]]
   dt$FC=ifelse(dt$logFC>0,2^dt$logFC,-1*(2^(-dt$logFC)))
   for (c in 1:ncol(dt)){
      colnames(dt)[c]=paste0(affected," vs ",baseline," _ ",colnames(dt)[c])
   }
   return(dt)
}

list=redf("WT RA PBS","KO RA PBS"); # list=cbind(list,redf("WT RA PBS","KO RA PBS"))
list=cbind(list,redf("WT RA PBS","WT SM FLU"))
list=cbind(list,redf("WT RA PBS","KO SM FLU"))

list=cbind(list,redf("WT SM FLU","KO SM FLU"))
list=cbind(list,redf("WT SM PBS","KO SM FLU"))
list=cbind(list,redf("WT RA FLU","KO SM FLU"))

list=cbind(list,redf("WT RA PBS","WT SM PBS"))
list=cbind(list,redf("WT RA PBS","KO SM PBS"))
list=cbind(list,redf("WT RA PBS","WT RA FLU"))
list=cbind(list,redf("WT RA PBS","KO RA FLU"))

list=cbind(list,redf("KO RA PBS","WT SM PBS"))
list=cbind(list,redf("KO RA PBS","KO SM PBS"))
list=cbind(list,redf("KO RA PBS","WT RA FLU"))
list=cbind(list,redf("KO RA PBS","KO RA FLU"))
write.csv(list,"All14pairs.csv")

# head(redf("WT RA PBS","KO RA PBS"))
# head(redf("WT RA FLU","KO SM FLU"))

# source("https://bioconductor.org/biocLite.R"); biocLite("systemPipeR")
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

setlist=list(MMP13=row.names(topTags(glmQLFTest(fit, coef=2),p.value=0.05,n=999999)),
             SmokeFlu=row.names(topTags(glmQLFTest(fit, coef=3),p.value=0.05,n=999999)),
             MMP13XSmokeFlu=row.names(topTags(glmQLFTest(fit, coef=5),p.value=0.05,n=999999)))


vennset <- overLapper(setlist,type="vennsets")
pdf(paste0("venn_",note,".pdf"))
vennPlot(vennset)
dev.off()


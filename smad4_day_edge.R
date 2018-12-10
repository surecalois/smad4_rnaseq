library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)

day_index = 4 # 1: hESC, 2: day1, 3: day3, 4:day6
lfc = 3
pvalue = 1e-3
#this script is adopted from the giwi thing by using day as the group

day=c("clone1","clone2","es03")
rep=c("A","B")
libs = c();
days = c();
reps = c();
for(ii in day){
  for(jj in rep){
    libs = c(libs,paste(ii,jj,sep = "."))
    days = c(days,ii)
    reps = c(reps,jj)
  }
}
expriment = data.frame(libs = libs, days = days, reps = reps)
design = model.matrix(~0+days,data = expriment)
colnames(design)=day
rownames(design)=libs

#make pairwise compare
L = length(day)
d=paste(day[-1],day[-L],sep="-")
my.contrasts = makeContrasts(contrasts=d,levels=design)

fds = matrix(0,nrow = L,ncol=L)
for(ii in seq(1,L)){
  for(jj in seq(1,L)){
    fds[ii,jj] = paste(day[ii],day[jj],sep = "-")
  }}

ds = c();
for(ii in seq(1,L)){
  for(jj in seq(1,L)){
    if(ii > jj){
      ds = c(ds,paste(day[ii],day[jj],sep = "-"))
    }
  }
}
my.contrasts2 = makeContrasts(contrasts = ds, levels = design)


#=====================end of prep for design information================

raw = load_smad4counts()
colselec = seq(day_index,24,4)
counttable = raw[, colselec]
#counttable = raw[keep,]


e = DGEList(counts=counttable,group = days)
e = calcNormFactors(e)
e = estimateGLMCommonDisp(e, design)
e = estimateGLMTrendedDisp(e,design)
e = estimateGLMTagwiseDisp(e,design)
fit = glmFit(e, design)

##================analysis any paired change===================##
r2=dim(counttable)[1]
c2=dim(my.contrasts2)[2]
d_map2=matrix(0,nrow = r2,ncol=c2)
upGs2=vector(mode="list",length=0)
dnGs2=vector(mode="list",length=0)
lrt_gene = rownames(fit$count);

for(ii in seq(1,c2)){
  lrt = glmLRT(fit,contrast = my.contrasts2[,ii])
  dt = decideTestsDGE(lrt,lfc=lfc,p.value = pvalue)
  d_map2[,ii]=dt
  upGs2[ii]=list(lrt_gene[dt ==1])
  dnGs2[ii]=list(lrt_gene[dt ==-1])
  
  detags = rownames(lrt$table)[as.logical(dt)]
  plotSmear(lrt,de.tags = detags)
}

B.neg = apply(d_map2,2,function(x) sum(x==-1))
B.pos = apply(d_map2,2,function(x) sum(x==1))

MB = matrix(0,L,L,dimnames = list(day,day));
a = 1
b = 1
for(ii in seq(1,L)){
  for(jj in seq(1,L)){
    if(ii>jj) {
      MB[ii,jj] = B.pos[b]
      b = b+1
    }
    if(ii<jj) {
      MB[ii,jj] = B.neg[a]
      a = a + 1
    }}}
##================end of analysis any paired change===================##

ens2sym(intersect(unlist(upGs2[3]),unlist(upGs2[2])))
ens2sym(intersect(unlist(dnGs2[3]),unlist(dnGs2[2])))
ens2sym(unlist(dnGs2[1]))
ens2sym(unlist(upGs2[1]))


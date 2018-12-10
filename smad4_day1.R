source('~/Dropbox/working/some_r/smad4_20170705/smad4_functions.R')

cell=c("clone1","clone2","es03")
rep=c("A","B")
libs = c();
groups = c();
for(ii in cell){
  for(kk in rep){
      libs = c(libs,paste(ii,kk,sep = "."))
      groups = c(groups,paste(ii,sep = "."))      
  }
}
expriment = data.frame(libs = libs, cells = cell, reps = rep, groups = groups)

design = model.matrix(~0+groups,data = expriment)
groups = factor(groups)
colnames(design)=levels(groups)

#make some compare
my.contrasts <- makeContrasts(
  es03smad4 = es03 - (clone1+clone2)/2,
  es03clone1 = es03 - clone1,
  es03clone2 = es03 - clone2,
  levels=design)

my_cmps = colnames(my.contrasts)

#=====================end of prep for design information================

counttable = load_smad4counts()#load_smad4clusters()
cluster = counttable#[counttable$cluster == 2,]
subcount = select(cluster,contains("d0"))

e = DGEList(counts=subcount,group = groups)
e = calcNormFactors(e)
e = estimateGLMCommonDisp(e, design)
e = estimateGLMTrendedDisp(e,design)
e = estimateGLMTagwiseDisp(e,design)
fit = glmFit(e, design)


lrt <- glmLRT(fit, contrast=my.contrasts[,'es03smad4'])

lfc = 1.5 #1.5
pvalue = 1e-3 #1e-3
gg_smear3(lrt,"es03smad4",lfc = lfc, pvalue = pvalue)

ltable = lrt$table
ltable$gene = rownames(ltable)

#ltable$symb = ens2sym(ltable$gene)
#write.table(ltable,'day1diff.txt',quote = F, sep = '\t',row.names = F)

dgenes = diff_genes(fit,my.contrasts[,'es03smad4'],'es03smad4',lfc = 1,pvalue = pvalue)
gg_volcano(ltable,dgenes)
    






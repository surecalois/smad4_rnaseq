source('~/Dropbox/working/some_r/smad4_20170705/smad4_functions.R')

cell=c("clone1","clone2","es03")
day=c("day0","day1","day3","day6")
rep=c("A","B")
libs = c();
groups = c();
for(ii in cell){
   for(kk in rep){
      for(jj in day){
        libs = c(libs,paste(ii,jj,kk,sep = "."))
        groups = c(groups,paste(ii,jj,sep = "."))      
    }
  }
}
expriment = data.frame(libs = libs, cells = cell, days = day, reps = rep, groups = groups)

design = model.matrix(~0+groups,data = expriment)
groups = factor(groups)
colnames(design)=levels(groups)

#make some compare
my.contrasts <- makeContrasts(
  day0.es03smad4 = es03.day0 - (clone1.day0+clone2.day0)/2,
  day1.es03smad4 = es03.day1 - (clone1.day1+clone2.day1)/2,
  day3.es03smad4 = es03.day3 - (clone1.day3+clone2.day3)/2,
  day6.es03smad4 = es03.day6 - (clone1.day6+clone2.day6)/2,
  day0.es03clone1 = es03.day0 - clone1.day0,
  day1.es03clone1 = es03.day1 - clone1.day1,
  day3.es03clone1 = es03.day3 - clone1.day3,
  day6.es03clone1 = es03.day6 - clone1.day6,
  day0.es03clone2 = es03.day0 - clone2.day0,
  day1.es03clone2 = es03.day1 - clone2.day1,
  day3.es03clone2 = es03.day3 - clone2.day3,
  day6.es03clone2 = es03.day6 - clone2.day6,
  levels=design)

my_cmps = colnames(my.contrasts)

#=====================end of prep for design information================

counttable = load_smad4counts()

e = DGEList(counts=counttable,group = groups)
e = calcNormFactors(e)
e = estimateGLMCommonDisp(e, design)
e = estimateGLMTrendedDisp(e,design)
e = estimateGLMTagwiseDisp(e,design)
fit = glmFit(e, design)

lfc = 3
pvalue = 1e-5

#ES03 vs SMAD4 mutants
up_counts = c();
down_counts = c();
for(k in seq(1,4)){
  dgenes = diff_genes(fit,my.contrasts[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change = data.frame(day = day,ups = up_counts,downs = down_counts)
z = gather(gene.change,change,value,-day)
p = ggplot(data = z,aes(x = day,y = value,fill = change))+
    geom_bar(stat="identity",position="dodge",width = 0.6)+
    #geom_path(group = day)+
    theme_light()+
    xlab("")+
    ylab("number of changed genes")+
    ggtitle("ES03 vs. SMAD4 mutants")
print(p)


#make a combined table
for(k in seq(1,4)){
  lrt <- glmLRT(fit, contrast=my.contrasts[,my_cmps[k]])
  table = lrt$table
  dt = decideTestsDGE(lrt,lfc=lfc,p.value = pvalue)
  table$dt = dt
  table$cmp = my_cmps[k]
  if(k == 1){big_table  = table}
  else big_table = rbind(big_table,table)
}

# a = big_table[big_table$dt != 0,]
# a$ens = substr(rownames(a),1,15)
# a$gene = ens2sym(a$ens)
# write.table(a, file = 'smad4_diff_table.txt',quote = F, sep = '\t',row.names =  F)

pdata = big_table#[big_table$cmp != "day0.es03smad4",]
p = ggplot(data = pdata,aes(x = logCPM,y = logFC,size = -log(PValue))) +
  geom_point(data = pdata[pdata$dt == 0,],alpha = I(1/21)) + 
  geom_point(data = pdata[pdata$dt > 0,],color = "green",alpha = I(1/3)) +
  geom_point(data = pdata[pdata$dt < 0,],color = "red",alpha = I(1/3)) +
  facet_grid(~cmp)+ theme_bw()
print(p)


#day1 vs. day0
my.contrasts2 <- makeContrasts(
  smad4c1.day1day0 = clone1.day1 - clone1.day0,
  smad4c2.day1day0 = clone2.day1 - clone2.day0,
     es03.day1day0 =   es03.day1 -   es03.day0,
  levels=design)
my_cmps = colnames(my.contrasts2)
up_counts = c();
down_counts = c();
for(k in seq(1,3)){
  dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change.d1d0 = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "d1d0")

#day3 vs. day0
my.contrasts2 <- makeContrasts(
  smad4c1.day3day0 = clone1.day3 - clone1.day0,
  smad4c2.day3day0 = clone2.day3 - clone2.day0,
  es03.day3day0 = es03.day3 - es03.day0,
  levels=design)
my_cmps = colnames(my.contrasts2)
up_counts = c();
down_counts = c();
for(k in seq(1,3)){
  dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change.d3d0 = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "d3d0")


#day6 vs. day0
my.contrasts2 <- makeContrasts(
  smad4c1.day1day0 = clone1.day6 - clone1.day0,
  smad4c2.day1day0 = clone2.day6 - clone2.day0,
     es03.day1day0 =   es03.day6 -   es03.day0,
  levels=design)
my_cmps = colnames(my.contrasts2)
up_counts = c();
down_counts = c();
for(k in seq(1,3)){
  dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change.d6d0 = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "d6d0")

gene.change = rbind(gene.change.d1d0,gene.change.d3d0,gene.change.d6d0)

z = gather(gene.change,change,value,-cell,-group)
p = ggplot(data = z,aes(x = cell,y = value,fill = change))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~group)+
  theme_light()+
  xlab("")+
  ylab("number of changed genes")+
  ggtitle("ES03 vs. SMAD4 mutants")
print(p)

#==============daily differential========

#day1 vs. day0
my.contrasts2 <- makeContrasts(
  smad4c1.day1day0 = clone1.day1 - clone1.day0,
  smad4c2.day1day0 = clone2.day1 - clone2.day0,
  es03.day1day0 =   es03.day1 -   es03.day0,
  levels=design)
my_cmps = colnames(my.contrasts2)
up_counts = c();
down_counts = c();
for(k in seq(1,3)){
  dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change.d1d0 = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "d1d0")

#day3 vs. day1
my.contrasts2 <- makeContrasts(
  smad4c1.day3day1 = clone1.day3 - clone1.day1,
  smad4c2.day3day1 = clone2.day3 - clone2.day1,
  es03.day3day1 = es03.day3 - es03.day1,
  levels=design)
my_cmps = colnames(my.contrasts2)
up_counts = c();
down_counts = c();
for(k in seq(1,3)){
  dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change.d3d1 = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "d3d1")

#day6 vs. day3
my.contrasts2 <- makeContrasts(
  smad4c1.day1day3 = clone1.day6 - clone1.day3,
  smad4c2.day1day3 = clone2.day6 - clone2.day3,
  es03.day1day3 =   es03.day6 -   es03.day3,
  levels=design)
my_cmps = colnames(my.contrasts2)
up_counts = c();
down_counts = c();
for(k in seq(1,3)){
  dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
  up_counts = c(up_counts,length(dgenes$up))
  down_counts = c(down_counts,length(dgenes$down))
}

gene.change.d6d3 = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "d6d3")

gene.change = rbind(gene.change.d1d0,gene.change.d3d1,gene.change.d6d3)

p = ggplot(data = gene.change, aes(x = ups, y = downs,shape = cell, color = group)) + 
    geom_point(size = 5)
print(p)

z = gather(gene.change,change,value,-cell,-group)
p = ggplot(data = z,aes(x = cell,y = value,fill = change))+
  geom_bar(stat="identity",position="dodge",width = 0.6)+
  facet_wrap(~group)+
  theme_light()+
  xlab("")+
  ylab("number of changed genes")+
  ggtitle("ES03 vs. SMAD4 mutants")
print(p)

#trying to find some common changes. 
#I don't think this will work. you should use cell group and do some intersect
# my.contrasts2 <- makeContrasts(
#   day1day0 = (clone1.day1+clone2.day1+2*es03.day1)/4 -  (clone1.day0+clone2.day0+2*es03.day0)/4,
#   day3day1 = (clone1.day3+clone2.day3+2*es03.day3)/4 -  (clone1.day1+clone2.day1+2*es03.day1)/4,
#   day6day3 = (clone1.day6+clone2.day6+2*es03.day6)/4 -  (clone1.day3+clone2.day3+2*es03.day3)/4,
#   levels=design)
# my_cmps = colnames(my.contrasts2)
# up_counts = c();
# down_counts = c();
# for(k in seq(1,3)){
#   dgenes = diff_genes(fit,my.contrasts2[,my_cmps[k]],my_cmps[k],lfc = lfc,pvalue = pvalue);
#   up_counts = c(up_counts,length(dgenes$up))
#   down_counts = c(down_counts,length(dgenes$down))
# }
# 
# gene.change.common = data.frame(cell = cells,ups = up_counts,downs = down_counts,group = "common")
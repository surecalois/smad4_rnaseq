source('~/Dropbox/working/some_r/smad4_20170705/smad4_functions.R')
#===========begin process=================#
day_index = 4 # 1: hESC, 2: day1, 3: day3, 4:day6
colselec = seq(day_index,24,4)

rpkmtable = load_smad4();
#rpkmtable = load_smad4cpm();

libs = colnames(rpkmtable)
A = grep('A$',libs)
B = grep('B$',libs)

#pat_table=data.frame(t(apply(rpkmtable,1,function(x) x/max(x))))
#pat_table=data.frame(t(apply(rpkmtable,1,function(x) x/sum(x))))
#rpkmtable=pat_table

file = "/Users/jiexu/Desktop/templist.txt"
#file = "/Users/jiexu/genome/annotation/list/FGFs.txt"
#file = "/Users/jiexu/genome/annotation/list/tgfbeta.txt"
#file = "/Users/jiexu/genome/annotation/list/wnt.txt"
#file = "/Users/jiexu/genome/annotation/list/my_cardiacDiff.txt"
#file = "/Users/jiexu/Desktop/smad4rnaseq_list1.txt"

listfile <- read.table(file, header=F);
inlist = listfile$V1
#YOU CAN PUT YOUR LIST HERE
# inlist = c("ISL1","NKX2-5","TNNT2","MYL7",
#            "PAX6","RAX","SIX3","FOXG1")

genelist=sym2ens(inlist);
genelist=gencode_clean(genelist)
if(anyNA(genelist)) print(paste0(inlist[is.na(genelist)]," not found"))

v=rpkmtable[genelist,]
v=gnorm(v,A,B)

v$gene = factor(inlist,levels = inlist,ordered = T)
#smad4_listbar(v)
#smad4_heatmap_long(v)
smad4_heatmap_list(v)






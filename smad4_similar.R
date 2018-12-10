source('~/Dropbox/working/some_r/smad4_20170705/smad4_functions.R')
#===========begin process=================#
day_index = 1 # 1: hESC, 2: day1, 3: day3, 4:day6
colselec = seq(day_index,24,4)

#raw = load_smad4();
raw = load_smad4cpm();
rpkmtable = raw
#rpkmtable = select(raw,contains("d0"))


libs = colnames(rpkmtable)
A = grep('A$',libs)
B = grep('B$',libs)

pat_table = gnorm(rpkmtable,A,B)

genelist = find_similar("T",pat_table,n = 20)
genelist=gencode_clean(genelist)
symbs = ens2sym(genelist)
inlist = symbs

if(anyNA(genelist)) print(paste0(inlist[is.na(genelist)]," not found"))

v=rpkmtable[genelist,]
v=gnorm(v,A,B)
v$gene = factor(inlist,levels = inlist,ordered = T)
#smad4_listbar(v)
smad4_heatmap_list(v)
#smad4_heatmap_long(v)






source('~/Dropbox/working/some_r/smad4_20170705/smad4_functions.R')

#==============begin process==========#
#raw = load_smad4();
raw = load_smad4cpm();
raw = rearrange_smad4(raw)

keep = apply(raw,1,function(x) sum(x) > 3 && sd2(x) > 0)
rpkmtable = raw[keep,]
ps =  apply(rpkmtable,1,ps2)
rpkmtable = rpkmtable[ps > 0.9,]

libs = colnames(rpkmtable)
A = grep('A$',libs)
B = grep('B$',libs)
pat_table = gnorm(rpkmtable,A,B)

smad4_mdsplot(pat_table)
smad4_pcaplot(pat_table)

counttable = load_smad4counts()
smad4_MDSedgeR(counttable)


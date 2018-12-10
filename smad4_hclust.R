rearrange_smad4 = function(raw,A,B){
  result = raw
  libs = colnames(raw)
  A = grep('A$',libs)
  B = grep('B$',libs)
  L  = dim(raw)[2]
  result[,seq(1,L,2)] = raw[,A]
  result[,seq(2,L,2)] = raw[,B]
  rlibs = libs
  rlibs[seq(1,L,2)] = libs[A]
  rlibs[seq(2,L,2)] = libs[B]
  colnames(result) = rlibs
  return(result)
}

#==============begin process==========#
#raw = load_smad4();
raw = load_smad4cpm();
raw = rearrange_smad4(raw)

keep = apply(raw,1,function(x) sum(x) > 3 && sd2(x) > 0)
rpkmtable = raw[keep,]
ps =  apply(rpkmtable,1,ps2)
rpkmtable = rpkmtable[ps > 0.98,]

libs = colnames(rpkmtable)
A = grep('A$',libs)
B = grep('B$',libs)
pat_table = gnorm(rpkmtable,A,B)

vv = pat_table[complete.cases(pat_table),] #remove NAs rows


vvk = t(apply(vv,1,function(x) x/max(x)))

dc = as.dist(1-cor(t(vv)))#dist(vv)#
hc = hclust(dc,method = "average")

hcc = cutree(hc,k = 12)
svv = data.frame(vvk,cluster = hcc,gene = rownames(vvk))
#just group the same cluster together, not yet order between clusters.
svv = svv[order(svv$cluster),] 

#cluster centeral mean
hc.centers = data.frame(group_by(svv,cluster) %>% select(-gene) %>% summarise_each(funs(mean)))

#cluster_orderbysize
hc.size = as.data.frame(table(hcc))[,"Freq"]
order.size = order(hc.size,decreasing = T)

#cluster_orderbytime
hc.time = apply(hc.centers[,-1],1,which.max)
order.time = order(hc.time)

#showing the cluster centers
gene = factor(hc.centers$cluster,levels = order.size,ordered = T)
z = t(apply(hc.centers[,-1],1,function(x) x/sum(x)))
show_data = data.frame(gene = gene, z)
show_data = show_data[order.size,]
smad4_heatmap(show_data)
smad4_listbar(show_data)

#show the heatmap,I think you should add cluster filter
svv$cluster = factor(svv$cluster,levels = order.time,ordered = T)
svv$gc = seq(dim(svv)[1],1)

vmat = gather(svv,day,value,c(-gene,-gc,-cluster))
p = lgeom_heatmap + geom_tile(data = vmat,aes(x = day, y = gc,fill = value))
print(p)

#visualized by cluster by PCA
cvv = select(svv,-gene,-gc)
giwi_cluster.pca(cvv)






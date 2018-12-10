source('~/Dropbox/working/some_r/smad4_20170705/smad4_functions.R')

#==============begin process==========#
#raw = load_smad4();
raw = load_smad4cpm();
raw = rearrange_smad4(raw)


keep = apply(raw,1,function(x) sum(x) > 3 && sd2(x) > 0)
#rpkmtable = raw[keep,seq(17,24)]
rpkmtable = raw[keep,]
ps =  apply(rpkmtable,1,ps2)
rpkmtable = rpkmtable[ps > 0.9,]
#rpkmtable = select(rpkmtable,starts_with("smad4"))

libs = colnames(rpkmtable)
A = grep('A$',libs)
B = grep('B$',libs)
pat_table = gnorm(rpkmtable,A,B)

vv = pat_table[complete.cases(pat_table),] #remove NAs rows

nc = 8 #12
#set.seed(12345)
km = kmeans(vv,nc, iter.max = 500000)

kmc = km$centers
kmc = t(apply(km$centers,1,function(x) x/sum(x)))
gene = factor(seq(1,nc),ordered = T)
kmc.df = data.frame(gene = gene,kmc) #use gene as column name because giwi_heatmap
smad4_heatmap(kmc.df) #show each cluster pattern
smad4_listbar(kmc.df)

#xi$ix keep the level order, xi$x is the size
xi = sort(km$size,index.return = T,decreasing = T)
#just to plot out the pattern
mmat = gather(kmc.df,day,value,-gene)
names(mmat)=c("cluster","day","value")
#factorized the size make alpha look better
mmat$size = factor(km$size[mmat$cluster]) 
#keep the plot order as decreasing size by factor the cluster
mmat$cluster = factor(mmat$cluster,levels = xi$ix,ordered = T)
p = ggplot(data = mmat,aes(x = day,group = cluster, y = value))+
  geom_line(aes(color = cluster,fill = cluster,alpha = size), size = 1)+
  facet_wrap(~ cluster,scales = "fix")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)#free fix free_x free_y


#here to rearrange the data to show the clusted heatmap
#basic idea is to group the same cluster together
#between the cluster order by the peak time order
unwanted = xi$ix[seq(1,2)] #usually the first two is alway not so cool, you can specify c(6,12)


wmi = apply(kmc,1,smad4_sd)
#wmi = apply(kmc,1,which.min) #according to time order
wi = rev(order(wmi)) #this keeps the new cluster level orders
#wi = c(9,12,7,3,6,5,1,10,11,4,8,2)#using manual order

#vvk = t(apply(vv,1,function(x) x/max(x)))
vvk = gnorm(vv)
svv = data.frame(vvk,
                 cluster = factor(km$cluster,levels = wi,ordered = T),
                 gene = rownames(vvk))
svv = svv[order(svv$cluster),]
discard = svv$cluster %in% unwanted
#svv = svv[!discard,]
svv$gc = seq(dim(svv)[1],1)

vmat = gather(svv,day,value,c(-gene,-gc,-cluster))
p = lgeom_heatmap + geom_tile(data = vmat,aes(x = day, y = gc,fill = value))
print(p)

#make a sider bar
# p = ggplot()+
#   geom_bar(data = data.frame(cluster = vmat$cluster),aes(x = 1,fill = cluster)) +
#   scale_fill_brewer(palette = "Set3") + theme_void() + scale_y_reverse()
# print(p)
# 
# #visualized cluster by pca
# cvv = select(svv,-gene,-gc)
# smad4_cluster.pca(cvv)

svv$sym = ens2sym(svv$gene)
#write.table(svv,file = "~/Desktop/smad4_cluster_kmeans.txt",quote = F, sep = "\t",row.names = F)


#tyring to extern the cluster to the whole set
# libs = colnames(raw)
# A = grep('A$',libs)
# B = grep('B$',libs)
# rww = gnorm(raw,A,B)
# ss = match(svv$gene,rownames(rww),nomatch = 0)
# ww = rww[ss,]
# ww$gc =  seq(dim(ww)[1],1)
# wmat = gather(ww,day,value,-gc)
# p = lgeom_heatmap + geom_tile(data = wmat,aes(x = day, y = gc,fill = value))
# print(p)


# SporeSeq 100% OTUs
library(metagenomeSeq)

# 2045 have same prevalence in whole and resistant
# 4329 are more prevalent in spore than whole
# 3855 are more prevalent in whole than spore

# C butyricum c(denovo2831, denovo500, denovo2976, denovo2732, denovo5872)
# cbut <- c('denovo2831', 'denovo500', 'denovo2976', 'denovo2732', 'denovo5872') 

poi.cts <- function(otu.tab,nboot=100){
  require(poilog)
  mu <- c()
  sig <- c()
  gof <- c()
  npoi <- matrix(NA,dim(otu.tab)[1],dim(otu.tab)[2])
  
  for (i in 1:dim(otu.tab)[2]){
    tmp <- poilogMLE(otu.tab[,i],nboot=nboot,zTrunc=F,startVals=c(mu=-3,sig=3),method="BFGS",control=list(maxit=10000))
    mu[i] <- tmp$par[1]
    sig[i] <- tmp$par[2]
    gof[i] <- tmp$gof
    npoi[,i] <- (otu.tab[,i]/exp(mu[i]))^1/sig[i]
  }
  
  rownames(npoi) <- rownames(otu.tab)
  colnames(npoi) <- colnames(otu.tab)
  
  return(list(npoi,cbind(mu,sig,gof)))
  
}

resist_index <- function(otus1,otus2,spore.i,whole.i,rng=1:length(spore.i),normalize=T,poiss=T,dev=2,rem.inf=T,qct=0){
  resist.index <- matrix(NA,dim(otus1)[1],length(spore.i))
  k <- 1
  for (i in rng){
    sp.a <- otus1[,spore.i[k]]; 
    sp.b <- otus2[,whole.i[k]]; 
    if (poiss){
      vps <- sqrt(sp.a + sp.b)
      aps <- abs(sp.a-sp.b)
      ips <- which(aps < dev*vps)
    }
    else{
      ips <- NULL
    }
    if (normalize){
      sp.a <- sp.a/sum(sp.a,na.rm=T)
      sp.b <- sp.b/sum(sp.b,na.rm=T)
    }
    sp.i <- sp.a/sp.b
    if (qct > 0){
      sp.qb <- which(sp.b < quantile(sp.b[sp.b > 0],qct))
      sp.i[sp.qb] <- NA
    }
    # exclude estimates which are bad b/c of low read counts in both samples
    sp.i[ips] <- NA
    resist.index[,k] <- sp.i
    k = k + 1
  }  
  colnames(resist.index) <- colnames(otus1)[whole.i]
  rownames(resist.index) <- rownames(otus1)
  # remove 0 over 0
  # used to remove INFINITE VALUES (ie when spores/whole -> Inf), so that it is possible to calculate resistance index
  if (rem.inf){
  resist.index[which(!is.finite(resist.index))] <- NA}
  # resist.index[which(resist.index==0)] <- NA
  resist.index[which(is.nan(resist.index))] <- NA
  return(resist.index)
}

# 100% OTUs
otu.100 <- read.table("/Users/skearney/Documents/Alm Lab/Spores/SS100/SS100_results/RDP/SporeSeq100.otu_table.100.denovo.rdp_assigned",sep="\t",header=T)
tax.100 <- as.character(otu.100[,1])
meta <- read.table("/Users/skearney/Documents/Alm Lab/Spores/SporeSeq Meta/spore_metadata.txt",sep="\t",header=T)
smp.names <- as.character(meta[,1])

contam <- grep("Alteromona|Oceanospir|Cyano",tax.100)
tax.100 <- tax.100[-contam]
tax.lst <- sapply(tax.100,function(x) strsplit(x,split=";d__"))
tax.phy <- sapply(tax.100,function(x) strsplit(x,split=";c__"))
tax.cls <- sapply(tax.100,function(x) strsplit(x,split=";o__"))
tax.gen <- sapply(tax.100,function(x) strsplit(x,split=";"))
tax.IDs <- c()
tax.pID <- c()
tax.cID <- c()
tax.names <- c()
for (i in 1:length(tax.lst)){tax.IDs[i] <- tax.lst[[i]][2]; tax.names[i] <- tax.lst[[i]][1]}
for (i in 1:length(tax.phy)){tax.pID[i] <- tax.phy[[i]][1]}
for (i in 1:length(tax.cls)){tax.cID[i] <- tax.cls[[i]][1]}

un.phy <- unique(tax.pID)
un.names <- unique(tax.names)
un.cls <- unique(tax.cID)
names(tax.names) <- tax.IDs

#tax.conv <- function(otu.table,names,tax){
#  taxa.mat <- matrix(NA,length(names),dim(otu.table)[2])
#  for (i in 1:length(names)){
#    tmp <- which(tax == names[i])
#    if (length(tmp) > 1){
#      taxa.mat[i,] <- colSums(otu.table[tmp,])
#    }
#    else{
#      taxa.mat[i,] <- otu.table[tmp,]
#    }
#  }
#  rownames(taxa.mat) <- names
#  return(taxa.mat)
#}

#tx.mt <- tax.conv(otu.cts,un.cls,tax.cID)

#idw.mn <- apply(tx.mt[,ss.idw],1,mean); idw.mn <- idw.mn/sum(idw.mn)
#ids.mn <- apply(tx.mt[,ss.ids],1,mean); ids.mn <- ids.mn/sum(ids.mn)

# plot top 9 classes from whole community
#w.class <- names(sort(idw.mn,decreasing=T)[1:9])

# plot top 7 classes from resist fraction
#r.class <- names(sort(ids.mn,decreasing=T)[1:7])

#vis.p <- cbind(idw.mn,ids.mn)
# 
#barplot(vis.p[w.class,])
otu.cts <- apply(otu.100[,2:dim(otu.100)[2]],2,as.numeric)
otu.cts <- otu.cts[,smp.names]
otu.cts <- otu.cts[-contam,]
rownames(otu.cts) <- tax.IDs
otu.css <- cumNormMat(newMRexperiment(otu.cts))
# otu.pod <- poi.cts(otu.cts) 
otu.rel <- apply(otu.cts,2,function(x) x/sum(x))
otu.bn <- apply(otu.cts,1,function(x) ifelse(x > 0,1,0))


# plot phyla greater than 100 counts



# separate whole samples (idw) from resistant samples (ids)
ts.ida <- c(49:72)
ts.idw <- c(6,9:24) + 48
ts.ids <- c(6,9:24) + 72
ss.idw <- c(1:24) 
ss.ids <- c(25:48) 
 

# calculate resist index

# this function looks good 
# finds statistically significant OTUs differing between samples
otu.rss <- resist_index(otu.css,otu.css,ss.ids,ss.idw,normalize=F,dev=2,rem.inf=T)

# counts the total # of times an OTU is seen
otu.fq1s <- apply(otu.cts[,ss.idw],1,function(x) sum(x > 0)) 
otu.fq2s <- apply(otu.cts[,ss.ids],1,function(x) sum(x > 0)) 
# otu100.fq <- apply(cbind(otu100.fq1,otu100.fq2),1,max)

# some OTUs are present predominantly as SPORES

# this step excludes ambiguous OTUs present in fewer than 2 people 
otu.ms <- apply(otu.rss,1,function(x) ifelse(length(which(!is.na(x))) > 0,median(x,na.rm=T),NA))
otu.ns <- apply(otu.rss,1,function(x) ifelse(length(which(!is.na(x))) > 0,length(x[!is.na(x)]),NA))

# 
rs.s <- names(which(otu.ms > 1 & otu.ns > 1))
nr.s <- names(which(otu.ms < 1 & otu.ns > 1))

res.tax <- tax.100[match(rs.s,tax.IDs)]
res.nam <- tax.names[rs.s]

# write.table(res.nam,file="ResistantRDPnames.txt")
# correlate OTU resist vs OTU overall
nres <- names(which(otu.ns > 11))
spo.cos <- rep(NA,dim(otu.css)[1])
spo.cot <- rep(NA,dim(otu.css)[1])

for (i in 1:dim(otu.css)[1]){
  if (otu.fq1s[i] > 11){
    spo.cos[i] <- cor.test(otu.css[i,ss.idw],otu.css['denovo236',ss.idw],method="spearman",use="complete.obs")$p.value
  }
  if (otu.fq1t[i] > 11){
    spo.cot[i] <- cor.test(otu.css[i,ts.ida],otu.css['denovo236',ts.ida],method="spearman",use="complete.obs")$p.value
  }
} 

# note denovo12 == GERM_otus denovo38 == 100% Faecalibacterium prausnitzii ATCC 



# plot unassigned -- don't lose Inf when doing this calculation (just assign -Inf < median Inf > median)
# otu100.med <- apply(otu100.rs,2,function(x) median(log(x[which(is.finite(x) & x > 0)],10)))
# otu100.log <- log(otu100.rs,10);
# otu100.cnt <- t(apply(otu100.log,1,function(x) x-otu100.med))
# otu100.ms <- apply(otu100.cnt,1,function(x) ifelse(sum(!is.na(x)) > 1,round(sum(x >= 0,na.rm=T)/sum(!is.na(x))),NA))
# top20.tax <- names(sort(otu100.ms,decreasing=T))[1:40]

# the only NAs present should be those for which the ratio of spore to whole was 0/0 or NaN
# or those which had "low" read counts in both samples such that difference couldn't be calculated
# ones with 0 spore and n > 0 whole will be 0s 

# find OTUs that don't have genus level info
gen.r <- grep("g__;",un.names)
un.names[-gen.r] -> un.names

# combine all taxa w same name
taxa.mat <- matrix(NA,length(un.names),dim(otu.cts)[2])

for (i in 1:length(un.names)){
  tmp <- which(tax.names == un.names[i])
  if (length(tmp) > 1){
    taxa.mat[i,] <- colSums(otu.css[tmp,])
  }
  else{
   taxa.mat[i,] <- otu.css[tmp,]
  }
}

rownames(taxa.mat) <- un.names
taxa.rel <- apply(taxa.mat,2,function(x) x/sum(x))

# plotting commands to make histograms for each phyla
hA <- hist(otu100.cnt,plot=F)
hAB <- hA$breaks
par(mfcol=c(2,2))
par(mar=c(0.5,0.5,0.5,0.5))
par(mar=c(0,0,0,0))
for (i in 1:40){
  boxplot(otu100.cnt[top20.tax[i],],outline=F,ylim=c(-2,2))
  stripchart(otu100.cnt[top20.tax[i],],add=T,pch=20,vertical=T)
  abline(h=0,lty=2)
}

hist(otu100.cnt[grep("p__",tax.pID),],xlim=c(-4,4),breaks=hAB,ylim=c(0,6e3),col=1,border=F)
hist(otu100.cnt[grep("p__Firmicutes",tax.pID),],xlim=c(-4,4),breaks=hAB,ylim=c(0,6e3),col=1,border=F)
hist(otu100.cnt[grep("p__Bacteroidetes",tax.pID),],xlim=c(-4,4),breaks=hAB,ylim=c(0,6e3),col=1,border=F)
hist(otu100.cnt[grep("p__Actinobacteria",tax.pID),],xlim=c(-4,4),breaks=hAB,ylim=c(0,6e3),col=1,border=F)

tax100.rs <- resist_index(taxa.mat,taxa.mat,ss.ids,ss.idw,normalize=F)
tax100.med <- apply(tax100.rs,2,function(x) median(log(x[which(is.finite(x) & x > 0)],10)))
tax100.log <- log(tax100.rs,10)
tax100.cnt <- t(apply(tax100.log,1,function(x) x-tax100.med))
tax100.ms <- apply(tax100.cnt,1,function(x) if(sum(!is.na(x) & is.finite(x)) > 7) {median(x,na.rm=T)} else{ return(NA)})
tax100.vr <- apply(tax100.cnt,1,function(x) if(sum(!is.na(x) & is.finite(x)) > 7) {IQR(x,na.rm=T)} else{ return(NA)})
top20.tax <- names(sort(tax100.ms,decreasing=T))[1:20]

gsp.ix <- which(tax100.ms > 0)

par(mfcol=c(1,20))
par(mar=c(0,0,0,0))
for (i in 1:20){
  boxplot(tax100.rs[top20.tax[i],],outline=F,ylim=c(-2,2),border=rgb(0.5,0.5,0.5,0.5))
  stripchart(tax100.rs[top20.tax[i],],add=T,pch=20,vertical=T,method="jitter")
  abline(h=0,lty=2)
}


# plot commands for resistant/non-resistant cells
par(mfcol=c(2,2))
colors <- c(rgb(1,1,1,0.2), rgb(0.8,0.8,0.8,0.2),rgb(0.1,0.1,0.1,0.2),rgb(0.2,0.2,0.2,0.2))
nrs.hist <- hist(otu.fq1s[nr.s],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[1],border="gray",ylim=c(0,0.15),right=F)
res.hist <- hist(otu.fq1s[rs.s],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[3],add=T,border="gray",ylim=c(0,0.15),right=F)
par(new=T)
plot(res.hist$mids,cumsum(nrs.hist$counts)/sum(nrs.hist$counts),type='o',xlim=c(0,25),ylim=c(0,1))
points(res.hist$mids,cumsum(res.hist$counts)/sum(res.hist$counts),type='o',xlim=c(0,25),pch=20,ylim=c(0,1))

nrs.hist <- hist(otu100.fq1[nrs.enr],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[2],border="gray",ylim=c(0,0.15),right=F)
res.hist <- hist(otu100.fq1[res.enr],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[4],add=T,border="gray",ylim=c(0,0.15),right=F)
par(new=T)
plot(res.hist$mids,cumsum(nrs.hist$counts)/sum(nrs.hist$counts),type='o',xlim=c(0,25),ylim=c(0,1))
points(res.hist$mids,cumsum(res.hist$counts)/sum(res.hist$counts),type='o',xlim=c(0,25),pch=20,ylim=c(0,1))

nrs.hist <- hist(otu100.fq1[res.enr],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[4],border="gray",ylim=c(0,0.15),right=F)
res.hist <- hist(otu100.fq[res.enr],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[3],add=T,border="gray",ylim=c(0,0.15),right=F)
par(new=T)
plot(res.hist$mids,cumsum(nrs.hist$counts)/sum(nrs.hist$counts),type='o',xlim=c(0,25),ylim=c(0,1))
points(res.hist$mids,cumsum(res.hist$counts)/sum(res.hist$counts),type='o',xlim=c(0,25),pch=20,ylim=c(0,1))

nrs.hist <- hist(otu100.fq1[nrs.enr],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[2],border="gray",ylim=c(0,0.15),right=F)
res.hist <- hist(otu100.fq[nrs.enr],freq=F,xlim=c(0,25),breaks=seq(0,25,1),col=colors[1],add=T,border="gray",ylim=c(0,0.15),right=F)
par(new=T)
plot(res.hist$mids,cumsum(nrs.hist$counts)/sum(nrs.hist$counts),type='o',xlim=c(0,25),ylim=c(0,1))
points(res.hist$mids,cumsum(res.hist$counts)/sum(res.hist$counts),type='o',xlim=c(0,25),pch=20,ylim=c(0,1))


# partition into Actinobacteria and Clostridia
act <- names(tax.names)[grep("Actino",tax.names)]
nr.s[match(act,nr.s)] -> nrs.act; nrs.act[which(!is.na(nrs.act))] -> nrs.act
rs.s[match(act,rs.s)] -> res.act; res.act[which(!is.na(res.act))] -> res.act

clo <- names(tax.names)[grep("Clostrid",tax.names)]
nr.s[match(clo,nr.s)] -> nrs.clo; nrs.clo[which(!is.na(nrs.clo))] -> nrs.clo
rs.s[match(clo,rs.s)] -> res.clo; res.clo[which(!is.na(res.clo))] -> res.clo

# plot OTUs partitioned by abundance (as spore?)
res.med <- apply(otu.css[rs.s,ss.idw],1,function(x) median(log(x[x!=0],10)))
nrs.med <- apply(otu.css[nr.s,ss.idw],1,function(x) median(log(x[x!=0],10)))

res.pmd <- apply(otu.css[res.enr,ss.idw],2,function(x) median(x[x!=0]))
nrs.pmd <- apply(otu.css[nrs.enr,ss.idw],2,function(x) median(x[x!=0]))

par(mfcol=c(1,2))
par(mar=c(0.5,0.5,0.5,0.5))
boxplot(nrs.med~otu100.fq1[nrs.enr],ylim=c(0,5),col="gray")
abline(h=seq(-1,4,1),col="gray",lty=2)
abline(h=median(nrs.med,na.rm=T))
boxplot(res.med~otu100.fq1[res.enr],ylim=c(0,5))
abline(h=seq(-1,4,1),col="gray",lty=2)
abline(h=median(res.med,na.rm=T))


# # identify OTUs within a person that are enriched (rather than looking in bulk)
# res.dat <- apply(tax100.rs,2,function(x) which(x > 1))
# nrs.dat <- apply(tax100.rs,2,function(x) which(x <= 1))
# 
# res.mat <- taxa.rel; res.mat[,] <- NA
# nrs.mat <- taxa.rel; nrs.mat[,] <- NA
# p.rn.w <- c()
# # make two matrices for resistant and non-resistant taxa
# for (i in 1:24){
#   res.mat[res.dat[[i]],i] <- taxa.mat[res.dat[[i]],i]
#   nrs.mat[nrs.dat[[i]],i] <- taxa.mat[nrs.dat[[i]],i]
#   tmp1 <- res.mat[,i]
#   tmp2 <- nrs.mat[,i]
#   p.rn.w[i] <- wilcox.test(tmp1,tmp2)$p.value
# }
# 
# nrs.m <- apply(nrs.mat,2,function(x) median(log(x,10),na.rm=T))
# res.m <- apply(res.mat,2,function(x) median(log(x,10),na.rm=T))

# combine independent p-values from non-parametric tests of distributions
res.dat <- apply(otu.rss,2,function(x) which(x >= 1))
nrs.dat <- apply(otu.rss,2,function(x) which(x < 1))

res.mat <- otu.css; res.mat[,] <- NA
nrs.mat <- otu.css; nrs.mat[,] <- NA
# make two matrices for resistant and non-resistant taxa
abd.comp <- c()
med.subt <- c()
for (i in 1:24){
  res.mat[res.dat[[i]],i] <- otu.css[res.dat[[i]],i]
  nrs.mat[nrs.dat[[i]],i] <- otu.css[nrs.dat[[i]],i]
  abd.comp[i] <- wilcox.test(res.mat[,i],nrs.mat[,i])$p.value
  med.subt[i] <- median(res.mat[,i],na.rm=T)-median(nrs.mat[,i],na.rm=T)
}


var.res <- apply(otu.css[res.enr,ss.idw],1,function(x) IQR(log(x[x > 0])))
med.res <- apply(otu.css[res.enr,ss.idw],1,function(x) median(log(x[x > 0])))
var.nrs <- apply(otu.css[nrs.enr,ss.idw],1,function(x) IQR(log(x[x > 0])))
med.nrs <- apply(otu.css[nrs.enr,ss.idw],1,function(x) median(log(x[x > 0])))

# when using variance and mean, resistant cells are more variant than non-resistant
# when using IQR and median, resistant cells are LESS variant than non-resistant cells

# plot # of cells counted only within resistant fraction


# time series methods
otu.rst <- resist_index(otu.css,otu.css,ts.ids,ts.idw,normalize=F,dev=2,rem.inf=T)

# counts the total # of times an OTU is seen
otu.fq1t <- apply(otu.cts[,ts.idw],1,function(x) sum(x > 0)) 
otu.fq2t <- apply(otu.cts[,ts.ids],1,function(x) sum(x > 0)) 
# otu100.fq <- apply(cbind(otu100.fq1,otu100.fq2),1,max)

# some OTUs are present predominantly as SPORES

# this step excludes ambiguous OTUs present in fewer than 2 people 
otu.mt <- apply(otu.rst,1,function(x) ifelse(length(which(!is.na(x))) > 0,median(x,na.rm=T),NA))
otu.nt <- apply(otu.rst,1,function(x) ifelse(length(which(!is.na(x))) > 0,length(x[!is.na(x)]),NA))

# 
rs.t <- names(which(otu.mt >= 1 & otu.nt > 1))
nr.t <- names(which(otu.mt < 1 & otu.nt > 1))

# here is a place for this bullshit so i don't have to keep rewriting it
half.p <- names(which(otu.fq1t > 11))

rs.c <- match(rs.t,half.p); rs.c[!is.na(rs.c)] -> rs.c
nr.c <- match(nr.t,half.p); nr.c[!is.na(nr.c)] -> nr.c
col.vec <- rep("white",length(half.p))
col.vec[nr.c] <- "gray"
col.vec[rs.c] <- "black"
col.vec[600] <- "green"

cor.half <- cor(t(otu.css[half.p,ts.ida]),method="spearman",use="complete.obs")

hc <- hclust(dist(cor.half),method="ward.D")
hcd <- as.dendrogram(hc)

k <- cutree(hc,2)
k.1 <- names(which(k==1))
k.2 <- names(which(k==2))
k.r1 <- intersect(k.1,rs.t)
k.r2 <- intersect(k.2,rs.t)
k.n1 <- intersect(k.1,nr.t)
k.n2 <- intersect(k.2,nr.t)

# dendrogram plotting command
library(dendextend)
hcd %>% set("leaves_pch", 22) %>%  # node point type
  set("leaves_cex", 0.5) %>%  # node point size
  set("leaves_col", col.vec[order.dendrogram(hcd)]) %>% # node point color
  set("labels",NULL) %>%
  plot()

# calculate correlations between metabolites and OTUs 
k = 0
otu.tsm <- matrix(NA,69,length(half.p))
colnames(otu.tsm) <- half.p
for (i in 1:69){
  k <- k + 1
  for (j in half.p){
    tmp <- cor(otu.css[j,ts.ida],mole.n[i,25:48],method="spearman",use="complete.obs") 
    otu.tsm[k,j] <- tmp
  }
}


# plot persistence of res/not-res
# extended ts.idw
ts.alt <- 49:72
rs.pt <- rowSums(otu.bn[ts.alt,rs.t])/length(rs.t)
nr.pt <- rowSums(otu.bn[ts.alt,nr.t])/length(nr.t)
rs.ps <- rowSums(otu.bn[ts.alt,rs.s])/length(which(otu.fq1t[rs.s] > 0))
nr.ps <- rowSums(otu.bn[ts.alt,nr.s])/length(which(otu.fq1t[nr.s] > 0))

# identify OTUs not shared between time series and overall as spores
# eg nr.s[!nr.s %in% nr.t] -- finds OTUs that are defined as not resistant in cross-section (nr.s)
# that are not defined as not resistant in the time series (nr.t)  
nr.pns <- rowSums(otu.bn[ts.alt,nr.s[!(nr.s %in% nr.t)]])/length(nr.s[!(nr.s %in% nr.t)])
rs.pns <- rowSums(otu.bn[ts.alt,rs.s[!(rs.s %in% rs.t)]])/length(rs.s[!(rs.s %in% rs.t)])
nr.pnt <- rowSums(otu.bn[ts.alt,nr.t[!(nr.t %in% nr.s)]])/length(nr.t[!(nr.t %in% nr.s)])
rs.pnt <- rowSums(otu.bn[ts.alt,rs.t[!(rs.t %in% rs.s)]])/length(rs.t[!(rs.t %in% rs.s)])

# plot organisms that do form spores in the time series but not in the 
rs.pns <- rowSums(otu.bn[ts.alt,rs.s[!(rs.s %in% rs.t)]])/length(rs.s[!(rs.s %in% rs.t)])
nr.pns <- rowSums(otu.bn[ts.alt,nr.s[!(nr.s %in% rs.t)]])/length(nr.s[!(nr.s %in% rs.t)])

# plot acf
res.acf <- apply(otu.css[res.enr,ts.ids],1,function(x) acf(log(x),plot=F)$acf)
nrs.acf <- apply(otu.css[nrs.enr,ts.idw],1,function(x) acf(log(x),plot=F)$acf)

ts.spd <- apply(otu.pod[[1]][,ts.ids],1,diff) 
ts.whd <- apply(otu.pod[[1]][,ts.idw],1,diff)

lag.tsp <- c()
lag.tsw <- c()
for (i in res.enr){
  lag.tsp[i] <- cor.test(otu.css[i,ts.ids[1:16]],otu.css[i,ts.idw[2:17]],method="spearman",use="complete.obs")$p.value
}

for (i in res.enr){
  lag.tsw[i] <- cor.test(otu.css[i,ts.idw[1:16]],otu.css[i,ts.idw[2:17]],method="spearman",use="complete.obs")$p.value
}

# make function to identify # of new organisms from start of time series plus number of lost organisms (never again present after x time point)
tsl <- 17
new.ts <- list()
new.ts[[1]] <- which(otu.css[,ts.ids[1]] > 0)
new.ct <- c()
new.ct[1] <- length(new.ts[[1]])

new.tt <- list()
new.tt[[1]] <- which(otu.css[,ts.idw[1]] > 0)
new.cw <- c()
new.cw[1] <- length(new.tt[[1]])

ovl.tt <- list(); ovl.tt[[1]] <- intersect(names(new.tt[[1]]),names(new.ts[[1]]))
ovl.ct <- c(); ovl.ct[1] <- length(ovl.tt[[1]])

for (i in 2:tsl){
  tmp1 <- which(otu.css[,ts.idw[i]] > 0)
  tmp2 <- which(!(tmp1 %in% unlist(new.tt)))
  new.tt[[i]] <- tmp1[tmp2]
  new.cw[i] <- length(new.tt[[i]])
  
  tmp1 <- which(otu.css[,ts.ids[i]] > 0)
  tmp2 <- which(!(tmp1 %in% unlist(new.ts)))
  new.ts[[i]] <- tmp1[tmp2]
  new.ct[i] <- length(new.ts[[i]])
  
  ovl.tt[[i]] <- intersect(names(new.tt[[i]]),names(new.ts[[i]]))
  ovl.ct[i] <- length(ovl.tt[[i]])
}

dis.ts <- list()
dis.ct <- c()

dis.tt <- list()
dis.cw <- c()

all.tax <- c()

for (i in 1:17){
  all.tax <- c(all.tax,names(unlist(new.ts[[i]])))
  all.tax <- unique(all.tax)
  if (i != 17){
    tmp1 <- apply(otu.pod[[1]][all.tax,ts.ids[i:17]],1,sum)
  }
  else {
    tmp1 <- otu.pod[[1]][all.tax,ts.ids[17]]
  }
  tmp2 <- which(tmp1 == 0)
  tmp3 <- which(!(tmp2 %in% unlist(dis.ts)))
  dis.ts[[i]] <- tmp2[tmp3]
  dis.ct[i] <- length(dis.ts[[i]])
}


all.tax <- c()


for (i in 1:17){
  all.tax <- c(all.tax,names(unlist(new.tt[[i]])))
  all.tax <- unique(all.tax)
  if (i != 17){
    tmp1 <- apply(otu.pod[[1]][all.tax,ts.idw[i:17]],1,sum)
  }
  else {
    tmp1 <- otu.pod[[1]][all.tax,ts.idw[17]]
  }
  tmp2 <- which(tmp1 == 0)
  tmp3 <- which(!(tmp2 %in% unlist(dis.tt)))
  dis.tt[[i]] <- tmp2[tmp3]
  dis.cw[i] <- length(dis.tt[[i]])
}





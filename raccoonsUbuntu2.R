library(adegenet)
library(vcfR)
library(wesanderson)
library(poppr)
library(hierfstat)
library(spdep)
library(dartR)
library(pcadapt )
library(qvalue)

vcf <- read.vcfR("./Fijarczyk/SNPs/SubSampl.maxMiss0.8.minGQ95.mac1.exclOverHetGenes.vcf.recode.vcf")

glPoly <- vcfR2genlight(vcf)
glPoly <- gl.compliance.check(glPoly,verbose = 0)
glPoly@ploidy <- rep(2L,length(glPoly@ind.names))
populations <- read.csv("populations.xls", sep="\t", row.names = 1) #zasysanie tabelki z danymi
popList <- data.frame(row.names = glPoly@ind.names)
for (i in rownames(popList))popList[i,1] <- populations[i,1]
pop(glPoly)<-popList[,1]
data.frame(glPoly@ind.names,glPoly@pop)

popsInNames <- c()
for (i in 1:nInd(glPoly)){
  popsInNames[i] <- sub("PL",paste0(pop(glPoly[i]),"_"),indNames(glPoly[i]))
  popsInNames[i] <- sub("-MP","x",popsInNames[i]) }
indNames(glPoly) <- popsInNames

myFreq <- glMean(glPoly)
hist(myFreq, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")

freqs <- data.frame(row.names = locNames(glPoly))
for (pop in popNames(glPoly)){
  freqs <- cbind(freqs,round(gl.alf(glPoly[pop(glPoly)==pop])[,1],4))
}
colnames(freqs)<- popNames(glPoly)

write.table(freqs,file = "1stAlleleFreq.xls",sep = "\t",quote = FALSE,col.names = NA,row.names = TRUE)


### PCADAPT
gl2plink(glPoly,outfile = "lowNAsPlink.csv",outpath = getwd())

data_raw <- read.csv("lowNAsPlink.csv",header = TRUE)
data <- data_raw[,2:length(data_raw[1,])]
PCAdapt <- read.pcadapt(t(data),type = "pcadapt")
rownames(PCAdapt) <- data_raw$X
x <- pcadapt(PCAdapt,K=30)
plot(x,option="screeplot")
k = 5
x <- pcadapt(PCAdapt,K=k)
plot(x,option = "scores",pop = pop(glPoly), col = c("red","blue","magenta","green", "black", "white"))
outliesFDR <- which(p.adjust(x$pvalues,method = "fdr") < alpha);outliesFDR
PCadaptRes <- glPoly$loc.names[outliesFDR]
PCadaptRes <- unique(gsub(("_[0-9]{1,}"),"",PCadaptRes));PCadaptRes
cat(PCadaptRes,file="PCadapt_Loci.txt",sep = "\n")
plot(1:length(x$pvalues), x$pvalues,log = "y",main = paste(levels(pop(glPoly)),collapse = " "),col="grey",pch=20,cex=0.5)+
  points(outliesQval,x$pvalues[outliesQval],col="blue",pch=20,cex=0.7)+
  points(outliesBon,x$pvalues[outliesBon],col="red",pch=20,cex=1)


locHet <- gl.basic.stats(glPoly)
write.table(locHet$perloc, col.names = NA, file = "heterozygsity.loc.txt",sep = "\t")
write.table(locHet$perloc, col.names = NA, file = "heterozygsity.loc.txt",sep = "\t")
write.table(locHet$Ho, col.names = NA, file = "Ho.locInPops.txt",sep = "\t")
write.table(locHet$Hs, col.names = NA, file = "Hs.locInPops.txt",sep = "\t")
write.table(locHet$Fis, col.names = NA, file = "Fis.locInPops.txt",sep = "\t")
glpoly.PCA <- glPca(glPoly,n.cores = 8,parallel = TRUE,nf = 3)
scatter(glpoly.PCA)
#loadingplot(glpoly.PCA, fac = pop(glPoly))
gl.pcoa.plot.3d(glpoly.PCA,glPoly)

#dartR--------------------------------------------------------
PCoA <- gl.pcoa(glPoly,nfactors = 6,parallel = TRUE,n.cores = 10)
barplot(PCoA$eig/sum(PCoA$eig)*100)

grp1 <- find.clusters(glPoly)
140
4

dapc1 <- dapc(glPoly,grp1$grp)
100
2
scatter(dapc1)
prop.table(table(pop(glPoly), grp1$grp),1)

mat <- tab(glPoly,NA.method = "mean")
grpGI <- pop(glPoly)

xval <- xvalDapc(mat,grpGI,n.pca.max = 80,training.set = 0.7,result = "groupMean", 
                 center = TRUE,scale = FALSE,n.pca = seq(5,50,by = 5),n.rep = 14000,
                 xval.plot = TRUE,parallel = "multicore",ncpus = 14L)

xval$`Number of PCs Achieving Highest Mean Success`
xvals <-as.integer(xval$`Number of PCs Achieving Lowest MSE`)
dapc2 <- dapc(glPoly, n.pca = xvals, n.da = 2, pop = pop(glPoly))

custom.cols <- c("blue","grey","yellow","orange","red")

pdf(paste0("DAPC Scatter",".pdf"), width = 10, height=10)
scatter(dapc2,scree.da = F,pch=19,cex=1,col = custom.cols)
dev.off()


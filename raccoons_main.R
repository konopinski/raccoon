suppressPackageStartupMessages({
  #library(QuasR)
  #library(BSgenome)
  #library(Rsamtools)
  #library(rtracklayer)
  #library(GenomicFeatures)
  #library(Gviz)
  library(adegenet)
  library(vcfR)
  library(wesanderson)
  #library(radiator)
  library(poppr)
  library(hierfstat)
  library(spdep)
  library(dartR)
  #library(LDna)
  library(igraph)
  library(ape)
  library(pegas)
  library(pcadapt)
  library(qvalue)
  # source("Modified.dartR.functions.R")
  }
)

vcf <- read.vcfR("./Fijarczyk/SNPs/SubSampl.maxMiss0.8.minGQ95.mac1.exclOverHetGenes.vcf.recode.vcf")

lociSNP <- as.data.frame(unique(vcf@fix[,1]))
for (i in 1:length(lociSNP[,1])){
  lociSNP[i,2] <- length(which(lociSNP[i,1]==vcf@fix[,1]))
  }
write.table(lociSNP, file="loci.txt",quote=FALSE,sep="\t",eol="\r\n")
gl <- vcfR2genlight(vcf)
gl <- gl.compliance.check(gl)
gl@ploidy <- rep(2L,length(gl@ind.names))

# file with population names
data <- read.csv("allSamples_SNP_VarFiltr_flt_biall_miss15_pops.adegenet.csv", sep=";", row.names = 1) #zasysanie tabelki z danymi
levels(data[,1]) <- c(levels(data[,1]),"D2")
populations <- as.data.frame(data[,1],row.names = rownames(data))
D2list <- c("PL288-MP","PL289-MP","PL290-MP","PL291-MP","PL292-MP","PL293-MP","PL294-MP")
for (D2 in D2list){
  populations[D2,1] <- "D2"}
# Putting population names in order of samples within genelight object
popList <- data.frame(row.names = gl@ind.names)
for (i in rownames(popList)){
  popList[i,1] <- populations[i,1]
}
popList[,1] <- sub("D1","D12",popList[,1])
popList[,1] <- sub("D2","D12",popList[,1])
popList[,1] <- sub("D3","D34",popList[,1])
popList[,1] <- sub("D4","D34",popList[,1])
popList[,1] <- sub("PL1","PL12",popList[,1])
popList[,1] <- sub("PL2","PL12",popList[,1])
pop(gl)<-popList[,1]

gl <- gl[,order(gl$loc.names)]
# Ordering the individuals
popOrder <- c("FL", "CZ", "D12", "D34", "PL12")
popOrdered <- c()
for (i in popOrder) popOrdered <- append(popOrdered,which(pop(gl)==i))
glPoly <- gl[popOrdered]
popsInNames <- c()
for (i in 1:nInd(glPoly)){
  popsInNames[i] <- sub("PL",paste0(pop(glPoly[i]),"_"),indNames(glPoly[i]))
  popsInNames[i] <- sub("-MP","x",popsInNames[i])

  }
indNames(glPoly) <- popsInNames

myFreq <- c(glMean(glPoly), 1-glMean(glPoly))
hist(myFreq, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)

gl_diversities_in_pops <- list()
for (pop in levels(glPoly$pop)){
  gl_diversities_in_pops <- genetic_diff(vcf,pops = as.factor(pop))
  write.table(gl_diversities_in_pops,file = paste0("diversity",pop,".xls"),sep="\t",eol="\r\n")
}


boot.fis <- function(gl,nboot = 9999,CI = 0.05){
  boots <- sapply(popNames(gl), function(pop){
    inds <- which(pop(gl)==pop)
    sapply(1:nboot,function(i)sample(inds, length(inds), TRUE))
  })
  loci <- genind2loci(gl2gi(gl))
  fises <- lapply(boots,function(i) apply(i,2,function(j){
    hets <- H(loci[j,],observed = TRUE,variance = FALSE)
    Fis <- apply(hets,1,function(k)round(1-(k[2]/k[1]),round(log10(nboot))))
    return(Fis)
  }))
    return(lapply(fises,function(l)t(apply(l,1,quantile,prob = c(CI/2,1-CI/2),na.rm = TRUE))))
}

fises <- boot.fis(glPoly)
for(i in names(fises)) write.table(fises[i],file = paste0(i,"Fis_CIs.xls"),sep = "\t",quote=FALSE,col.names =  NA)

### PCADAPT
gi <- gl2gi(glPoly)
data <- gi$tab[,1:ncol(gi@tab)%%2==1]
data[is.na(data[,])] <- 9
PCAdapt <- read.pcadapt(data,type = "lfmm")
x <- pcadapt(PCAdapt,K=30)
plot(x,option="screeplot")
k = 5
x <- pcadapt(PCAdapt,K=k)
#plot(x,option = "qqplot")
alpha=0.05
outliesFDR <- which(p.adjust(x$pvalues,method = "fdr") < alpha);outliesFDR
lociSelectedPCA <- glPoly$loc.names[outliesFDR]
lociSelectedPCAUNIQ <- unique(gsub(("_[0-9]{1,}"),"",lociSelectedPCA));lociSelectedPCAUNIQ
FDRy <- data.frame(ix = outliesFDR, ygrek = -log(x$pvalues[outliesFDR],base = 10))
nazwy <- read.csv("genenames.txt", header = TRUE, sep = "\t",row.names = 1)[,-1]
gene <- gsub(("_[0-9]{1,}"),"",lociSelectedPCA)
site <- gsub(("Gene.[0-9]{1,}_"),"",lociSelectedPCA)
selected <- data.frame(gene,site)
for (i in 1:nrow(selected)){
  nazwa_genu  <- selected[i,1]
  if(nazwy[nazwa_genu,"Ferret"]!="") selected[i,1] <- nazwy[nazwa_genu,"Ferret"] else 
    if(nazwy[nazwa_genu,"Dog"]!="") selected[i,1] <- nazwy[nazwa_genu,"Dog"] else
    if(nazwy[nazwa_genu,"Human"]!="") selected[i,1] <- nazwy[nazwa_genu,"Human"]
    }

labele <- data.frame(label = selected[,1], site = selected[,2], x=FDRy[,1], y=FDRy[,2])
labele[9,1] <- "SAA1"

for (i in dim(labele)[1]:1) if (labele[i,"label"] %in% labele[1:i-1,"label"]) labele <- labele[-i,]
p <- plot(x,option = "manhattan")
pdf("Fig. 9B. PCAdapt.pdf", width = 10,height = 7)
p+
  geom_point(aes(x=ix, y=ygrek), data = FDRy, size =2, col = "red")+
  geom_text(aes(x=x, y=y,label = label), data = labele, 
            hjust = -0.15)+
  labs(y = expression(italic(log[10])~p-value),
       x = "SNPs") +
  theme_bw()
dev.off()

### tree in ggplot
tre <- ape::nj(dist(as.matrix(glPoly)))
pdf(paste0("Drzewo.","pdf"), width = 15, height=15) # boxplot
plot(tre, type = "unrooted", cex=0.7)
title("NJ tree of SNP data")
dev.off()

#### heterozygosity i fisy
popHet <- gl.report.heterozygosity(glPoly)
write.table(popHet, col.names = NA, file = "heterozygsity.pop.txt")
locHet <- gl.basic.stats(glPoly)
write.table(locHet$perloc, col.names = NA, file = "heterozygsity.loc.txt")
write.table(locHet$perloc, col.names = NA, file = "heterozygsity.loc.txt")
write.table(locHet$Ho, col.names = NA, file = "Ho.locInPops.txt")
write.table(locHet$Hs, col.names = NA, file = "Hs.locInPops.txt")
write.table(locHet$Fis, col.names = NA, file = "Fis.locInPops.txt")

#dartR--------------------------------------------------------
dapc1 <- dapc(glPoly)
100
2
scatter(dapc1)

mat <- tab(glPoly,NA.method = "mean")
grpGI <- pop(glPoly)
save.image()
xval <- xvalDapc(mat,grpGI,n.pca.max = 60,training.set = 0.8,result = "groupMean", 
                 center = TRUE,scale = FALSE,n.pca = seq(10,40,by = 2),n.rep = 14000,
                 xval.plot = TRUE,parallel = "multicore",ncpus = 14L)

XVpca <- xval$`Number of PCs Achieving Lowest MSE`
dapc2 <- dapc(glPoly, n.pca = as.numeric(XVpca), n.da = 2, pop = pop(glPoly))
pdf(paste0("DAPC Scatter",".pdf"), width = 15, height=15)
scatter(dapc2,scree.da = F)
dev.off()


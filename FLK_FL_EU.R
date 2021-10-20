setwd("/media/sf_Sharing/R/MiPy/FLK")
source('FLK.R')
library(vcfR)
library(dartR)
vcf <- read.vcfR("../Fijarczyk/SNPs/SubSampl.maxMiss0.8.minGQ95.mac1.exclOverHetGenes.vcf.recode.vcf")
gl <- vcfR2genlight(vcf)
gl <- gl.compliance.check(gl)
gl@ploidy <- rep(2L,length(gl@ind.names))
#Old file that contains population names
data <- read.csv("../allSamples_SNP_VarFiltr_flt_biall_miss15_pops.adegenet.csv", sep=";", row.names = 1) #zasysanie tabelki z danymi
populations <- as.data.frame(data[,1],row.names = rownames(data))
populations["PL289-MP",] <- "CE"
#Putting population names in order of samples within genelight object
popList <- data.frame(row.names = gl@ind.names)
for (i in rownames(popList))  popList[i,1] <- populations[i,1]
popList[,1] <- sub("D.*","CE",popList[,1])
popList[,1] <- sub("PL.*","CE",popList[,1]);pop(gl)<-popList[,1]

#ułożenie po kolei populacjami
popOrder <- c("FL", "CZ", "CE")
popOrdered <- unlist(sapply(popOrder,function(i) which(pop(gl)==i)));gl <- gl[popOrdered]
gl <- gl.filter.callrate(gl, t=0.95, plot = FALSE,recalc = TRUE,mono.rm = TRUE)

prefix <- "FL_EU_2gr"
outgroup <- "FL"
folder <- paste0("./",prefix,"/")
system2("mkdir", args = prefix)
system2("cp", args = c("./FLKnulll", paste0(folder,"FLKnull")))
fname = paste0(folder,"freq.dat")
cat(file = fname)
sapply(levels(pop(gl)),FUN = function(pop){
  cat(paste(pop,"\t"),file = fname, sep = "\t",append = TRUE)
  cat(glMean(gl[pop(gl)==pop]),file = fname,append = TRUE)
  cat("\n",file = fname,append = TRUE)})
freq <- read.table(fname,row.names = 1)
colnames(freq) <- locNames(gl)

DR_MS <- read.table("Reynolds_FL_EU.txt",sep = "\t", row.names = 1L)
colnames(DR_MS) <- rownames(DR_MS)
pops <- levels(pop(gl))
DR_MS <- DR_MS[pops,pops]

workDir <- getwd()
setwd(folder)
FijMS=Fij(freq,outgroup=outgroup,D=DR_MS)
testsMS=FLK(freq,FijMS) 
write.table(testsMS,file = paste0("FLK_",prefix,".xls"),col.names = NA,sep = "\t")
simuPOP <- reticulate::import("simuPOP")
cat("\n######  ",folder,"  #####\n")
system2("python3", args = c(paste0("FLKnull"), "1000000"))
setwd(workDir)
system2("mv", args = c(paste0(folder,"FLK_",prefix,".xls"), paste0("./FLK_",prefix,".xls")))
system2("mv", args = c(paste0(folder,"envelope.txt"), paste0("./envelope_",prefix,".txt")))
system2("rm", args = c("-r", paste0("./", prefix)))

null <- read.table(paste0("envelope_", prefix, ".txt"),head=T)
maxFLK = max(c(testsMS$F.LK,null$q0.995))
p05 <- testsMS[testsMS$F.LK.p.val<0.05,]
UpperLim <- splinefun(x=null$Ht, y=null$q0.995)
LowerLim <- splinefun(x=null$Ht, y=null$q0.005)
above <- UpperLim(testsMS$Ht)<testsMS$F.LK
below <- LowerLim(testsMS$Ht)>testsMS$F.LK
above[is.na(above)] <- FALSE
below[is.na(below)] <- FALSE
aboveLoci <- round(testsMS[above,],3)
belowLoci <- round(testsMS[below,],3)
Dewianty <- cbind(aboveLoci,deviation = rep("ABOVE", nrow(aboveLoci)))
Dewianty <- rbind(Dewianty,cbind(belowLoci,deviation = rep("BELOW", nrow(belowLoci))))
write.table(Dewianty[,c(1,2,6)],sep = "\t",file=paste0("devs_", prefix,".txt"),
            quote = FALSE,col.names = NA)

pdf(paste0(prefix, ".pdf"))
plot(null$Ht,null$q0.005,xlim=c(0,0.5),ylim=c(0.001,max(unlist(maxFLK,null))),
     col='gray',xlab='Heterozygosity',ylab='FLK statistic', 
     main = paste("MS",prefix))+
  lines(null$Ht,null$q0.995,type='l',col='gray')+
  lines(null$Ht,null$q0.025,type='l',col='gray')+
  lines(null$Ht,null$q0.975,type='l',col='gray')+
  lines(null$Ht,null$q0.5,type='l',col='gray')+
  lines(smooth.spline(null$Ht,null$q0.995))+
  lines(smooth.spline(null$Ht,null$q0.005))+
  lines(smooth.spline(null$Ht,null$q0.025),lty=2)+
  lines(smooth.spline(null$Ht,null$q0.975),lty=2)+
  lines(smooth.spline(null$Ht,null$q0.5),lty=3)+
  points(testsMS$Ht,testsMS$F.LK,pch=16)+
  points(p05$Ht,p05$F.LK,pch=16,col = "blue")+
  points(aboveLoci$Ht, aboveLoci$F.LK, col = "gold",pch = 16)+
  points(belowLoci$Ht, belowLoci$F.LK, col = "red",pch = 16)
dev.off()


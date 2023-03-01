library(PopGenome)
library(PopGenome)
library(ggplot2)
library(reshape2)
setwd("./pg")
source("./PopGenomeFunctions.R")
source("./set.synnonsyn.R")

popfile <- read.table("../SubSampl.SG.txt")
pops <- unique(popfile[,2])
poplist <- sapply(pops,function(pop){
  inds <- popfile[popfile[,2]==pop,1]
  c(inds,paste0(inds,".2"))
})

PGfile <- readData(path="vcf",gffpath = "gff", 
                   format = "VCF", include.unknown = TRUE)
PGfile <- set.synnonsyn(PGfile, ref.chr=FALSE,
                        ref.fasta.folder="./fasta/")

popgenom <- set.populations(PGfile,poplist)
popgenom <- neutrality.stats(popgenom)
popgenom <- detail.stats(popgenom)
popgenom <- F_ST.stats(popgenom)
tajD <-popgenom@Tajima.D
rownames(tajD) <- popgenom@region.names
colnames(tajD) <- names(poplist)
FuLiD <-popgenom@Fu.Li.D
rownames(FuLiD) <- popgenom@region.names
colnames(FuLiD) <- names(poplist)
FuLiF <-popgenom@Fu.Li.F
rownames(FuLiF) <- popgenom@region.names
colnames(FuLiF) <- names(poplist)
nseg <- popgenom@n.segregating.sites
rownames(nseg) <- popgenom@region.names
colnames(nseg) <- names(poplist)
write.table(tajD,file = "../PG_TajimaD.xls",sep = "\t", col.names = NA)
write.table(FuLiD,file = "../PG_FuLiD.xls",sep = "\t", col.names = NA)
write.table(FuLiF,file = "../PG_FuLiF.xls",sep = "\t", col.names = NA)
write.table(nseg,file = "../PG_NsegSites.xls",sep = "\t", col.names = NA)

popgenom.syn <- set.populations(PGfile,poplist)
popgenom.syn <- neutrality.stats(popgenom,subsites = "syn")
popgenom.syn <- detail.stats(popgenom.syn,subsites = "syn")
popgenom.syn <- F_ST.stats(popgenom.syn,subsites = "syn")
tajDsyn <-popgenom.syn@Tajima.D
rownames(tajDsyn) <- popgenom.syn@region.names
colnames(tajDsyn) <- names(poplist)
FuLiDsyn <-popgenom.syn@Fu.Li.D
rownames(FuLiDsyn) <- popgenom.syn@region.names
colnames(FuLiDsyn) <- names(poplist)
FuLiFsyn <-popgenom.syn@Fu.Li.F
rownames(FuLiFsyn) <- popgenom.syn@region.names
colnames(FuLiFsyn) <- names(poplist)
nsegSyn <- popgenom.syn@n.segregating.sites
rownames(nsegSyn) <- popgenom.syn@region.names
colnames(nsegSyn) <- names(poplist)
write.table(tajDsyn,file = "../PG_TajimaD_syn.xls",sep = "\t", col.names = NA)
write.table(FuLiDsyn,file = "../PG_FuLiD_syn.xls",sep = "\t", col.names = NA)
write.table(FuLiFsyn,file = "../PG_FuLiF_syn.xls",sep = "\t", col.names = NA)
write.table(nsegSyn,file = "../PG_NsegSites_syn.xls",sep = "\t", col.names = NA)

popgenom.nonsyn <- set.populations(PGfile,poplist)
popgenom.nonsyn <- neutrality.stats(popgenom.nonsyn,subsites = "nonsyn")
popgenom.nonsyn <- detail.stats(popgenom.nonsyn,subsites = "nonsyn")
popgenom.nonsyn <- F_ST.stats(popgenom.nonsyn,subsites = "nonsyn")
tajDns <-popgenom.nonsyn@Tajima.D
rownames(tajDns) <- popgenom.nonsyn@region.names
colnames(tajDns) <- names(poplist)
FuLiDns <-popgenom.nonsyn@Fu.Li.D
rownames(FuLiDns) <- popgenom.nonsyn@region.names
colnames(FuLiDns) <- names(poplist)
FuLiFns <-popgenom.nonsyn@Fu.Li.F
rownames(FuLiFns) <- popgenom.nonsyn@region.names
colnames(FuLiFns) <- names(poplist)
nsegNS <- popgenom.nonsyn@n.segregating.sites
rownames(nsegNS) <- popgenom.nonsyn@region.names
colnames(nsegNS) <- names(poplist)
write.table(tajDns,file = "../PG_TajimaD_nonsyn.xls",sep = "\t", col.names = NA)
write.table(FuLiDns,file = "../PG_FuLiD_nonsyn.xls",sep = "\t", col.names = NA)
write.table(FuLiFns,file = "../PG_FuLiF_nonsyn.xls",sep = "\t", col.names = NA)
write.table(nsegNS,file = "../PG_NsegSites_nonsyn.xls",sep = "\t", col.names = NA)

exclGenes <- c("Gene.307249", "Gene.365012", "Gene.381666", "Gene.217376", "Gene.151995", "Gene207091")
TDall <- tajD[!rownames(tajD)%in%exclGenes,]
TDall <- melt(TDall)
TDall$Var2 <- as.factor(TDall$Var2)
TDall <-TDall[!is.nan(TDall$value),]
TDall <-TDall[!is.na(TDall$value),]
names(TDall) <- c("gene","pop","Tajima")

TDsyn <- tajDsyn[!rownames(tajDsyn)%in%exclGenes,]
TDsyn <- melt(TDsyn)
TDsyn$Var2 <- as.factor(TDsyn$Var2)
TDsyn <-TDsyn[!is.nan(TDsyn$value),]
TDsyn <-TDsyn[!is.na(TDsyn$value),]
names(TDsyn) <- c("gene","pop","Tajima")

TDnSyn <- tajDns[!rownames(tajDns)%in%exclGenes,]
TDnSyn <- melt(TDnSyn)
TDnSyn$Var2 <- as.factor(TDnSyn$Var2)
TDnSyn <-TDnSyn[!is.nan(TDnSyn$value),]
TDnSyn <-TDnSyn[!is.na(TDnSyn$value),]
names(TDnSyn) <- c("gene","pop","Tajima")
TDtab <- rbind(cbind(TDall,sites = rep("all",nrow(TDall))),
      cbind(TDnSyn,sites = rep("NS",nrow(TDnSyn))),
      cbind(TDsyn,sites = rep("Syn",nrow(TDsyn))))
p <- ggplot(TDtab,aes(x = sites, y = Tajima, fill = pop)) 
pdf("tajimaDviol.pdf", width = 15,height = 8)
p+  geom_violin(trim = FALSE)
dev.off()
TDlist <- list(tajD,tajDns,tajDsyn)
TDnames <- c("All","NS","Syn")
for (i in 1:3){
pdf(paste0("tajimaDensityPlot",TDnames[i],".pdf"))
plot(density(TDlist[[i]][,1],na.rm = TRUE),lty = 1,lwd = 2, main = paste0("Tajima's D\n",TDnames[i]))
lines(density(TDlist[[i]][,2],na.rm = TRUE),lty = 2,lwd = 2)
lines(density(TDlist[[i]][,3],na.rm = TRUE),lty = 3,lwd = 2)
lines(density(TDlist[[i]][,4],na.rm = TRUE),lty = 4,lwd = 2)
lines(density(TDlist[[i]][,5],na.rm = TRUE),lty = 5,lwd = 2)
legend("topright",legend = c("FL","D12","CZ","D34","PL12"),lty = c(1,2,3,4,5),lwd = 2)
dev.off()}
getwd()
legend = c("FL","D12","CZ","D34","PL12")
for (i in 1:5){
pdf(paste0("tajimaDensityPlot_",legend[i],".pdf"))
plot(density(tajD[,i],na.rm = TRUE),lty = 1,lwd = 2, main = paste0("Tajima's D\n",legend[i]))
lines(density(tajDns[,i],na.rm = TRUE),lty = 2,lwd = 2)
lines(density(tajDsyn[,i],na.rm = TRUE),lty = 3,lwd = 2)
legend("topright",legend = c("all","ns","syn"),lty = c(1,2,3),lwd = 2)
dev.off()}

#  permutacje Pi
perm <- 10000
popgenom


MSout <- MS(popgenom,niter=9999,thetaID="Watterson",params=FALSE,detail=FALSE,
   neutrality=TRUE,linkage=FALSE,F_ST=FALSE,MSMS=FALSE,big.data=FALSE)
MSout@prob.less[1]


######################
PGfile2 <- readData(path="./vcf/",gffpath = "./gff/", 
                   format = "VCF", include.unknown = TRUE)
PGfile2 <- set.synnonsyn(PGfile2, ref.chr=FALSE,
                        ref.fasta.folder="./fasta/")
PGfile2 <- set.populations(PGfile2, poplist)
PGfile2@Coding.region
sum(sapply(PGfile2@region.data@synonymous, function (i)sum(i==1)))
PGfile2@region.data@biallelic.sites[[188]]
PGfile2@region.data@synonymous[[188]]

cat(file = "synPositions.txt")
cat(file = "nsynPositions.txt")
sapply(1:length(PGfile2@region.data@synonymous),function(i){
  syns <- PGfile2@region.data@synonymous[[i]]
  posits <- PGfile2@region.data@biallelic.sites[[i]]
  regname <- PGfile2@region.names[[i]]
  for (posit in 1:length(posits)){
    if(syns[posit]==1) cat(regname,"\t",posits[posit],"\n",file = "synPositions.txt",append = TRUE)
    if(syns[posit]==0) cat(regname,"\t",posits[posit],"\n",file = "nsynPositions.txt",append = TRUE)
}})


infoRDS("./pg/MS_Tajima")
popgenom2 <- set.populations(PGfile2,poplist)
popgenom2 <- neutrality.stats(popgenom2)
popgenom2 <- diversity.stats(popgenom2,pi = TRUE)
popgenom2 <- detail.stats(popgenom2)
MSoutPops <- MS(popgenom2,niter=99999,thetaID="Watterson",neutrality=TRUE)

mean(MS_getStats(MSoutPops,locus = i)[,1],na.rm = TRUE)
lessThan <- sapply(MSoutPops@prob.less,function(i){
  i[,1]
})
colnames(lessThan) <- names(poplist)
rownames(lessThan) <- c(popgenom2@region.names,"average","variance")
write.table(lessThan,file = "../ProbTDsmallerThanSim.xls",
            quote = FALSE, row.names = TRUE,sep = "\t",col.names = NA,dec = ",")
equalTo <- sapply(MSoutPops@prob.equal,function(i){
  i[,1]
})
colnames(equalTo) <- names(poplist)
rownames(equalTo) <- c(popgenom2@region.names,"average","variance")
write.table(equalTo,file = "../ProbTDequalToSim.xls",
            quote = FALSE, row.names = TRUE,sep = "\t",col.names = NA,dec = ",")
pops <- names(poplist)
loci <- popgenom2@region.names
for (pop in 1){
  for (locus in 1){
    tajD <- MSoutPops@locus[[locus]]@obs.val[pop,1]
    kwantyle <- MSoutPops@locus[[locus]]@quantiles[[pop]][c(1,2,4,7,10,12,13),1]
    if (!is.na(tajD)){
      jpeg(paste0("./Locus.",loci[locus],"_pop.",pops[pop],"_TajimaTheta.jpg"),
           height = 1024,width = 1024)
      hist(MSoutPops@locus[[locus]]@stats[[pop]][,1],breaks = 100,
           main = paste(pops[pop], loci[locus]),xlab = "Tajima D")
      abline(v = kwantyle,lwd = 2,col = "green")
      abline(v = tajD,lwd = 1.5,col = "blue")
      dev.off()
    }
   }
}
MSoutPops@locus
greater <- sapply(1:5,function(j)sapply(MSoutPops@locus,function(i){
  tajD <- i@obs.val[j,1]
  simuls <- i@stats[[j]][,1]
  if(!is.na(tajD)){
  sum(simuls>=tajD,na.rm = TRUE)/sum(!is.na(simuls))
  } else NA
}
))
colnames(greater) <- pops
write.table(greater,row.names = loci,col.names = NA,
            file = "../TajimaGreater.xls",sep = "\t",quote = FALSE,dec = ",")

##### Tajima theta 
MSoutPopsTaj <- PopGenome::MS(popgenom2,niter=999,thetaID="Tajima", 
                neutrality=TRUE)
lessThan <- sapply(MSoutPopsTaj@prob.less,function(i)i[,1])
colnames(lessThan) <- names(poplist)
rownames(lessThan) <- c(popgenom2@region.names,"average","variance")
write.table(lessThan,file = "../ProbTDsmallerThanSimThetaTajima.xls",
            quote = FALSE, row.names = TRUE,sep = "\t",col.names = NA,dec = ",")

equalTo <- sapply(MSoutPopsTaj@prob.equal,function(i)i[,1])
colnames(equalTo) <- names(poplist)
rownames(equalTo) <- c(popgenom2@region.names,"average","variance")
write.table(equalTo,file = "../ProbTDequallToSimThetaTajima.xls",
            quote = FALSE, row.names = TRUE,sep = "\t",col.names = NA,dec = ",")
greaterT <- sapply(1:5,function(j)sapply(MSoutPopsTaj@locus,function(i){
  tajD <- i@obs.val[j,1]
  simuls <- i@stats[[j]][,1]
  if(!is.na(tajD)){
    sum(simuls>=tajD,na.rm = TRUE)/sum(!is.na(simuls))
  } else NA
}
))
colnames(greaterT) <- pops
write.table(greaterT,row.names = loci,col.names = NA,
            file = "../TajimaGreaterTajimasTheta.xls",sep = "\t",quote = FALSE,dec = ",")

#popgenom2 <- set.outgroup(popgenom2,poplist2[[1]],diploid = TRUE)
popgenom2 <- diversity.stats(popgenom2)
popgenom2 <- detail.stats(popgenom2)
popgenom2 <- neutrality.stats(popgenom2,detail = TRUE,do.R2 = TRUE)
get.neutrality(popgenom2)[1]
popgenom2 <- F_ST.stats(popgenom2)
FST_All<- get.F_ST(popgenom2)[,"nucleotide.F_ST"]
FSTs <- cbind(FST_All,popgenom2@nuc.F_ST.vs.all)
colnames(FSTs)[2:6] <- names(poplist2)
rownames(FSTs) <- popgenom2@region.names
write.table(FSTs, file="../PG_FSTs_total_and_popVSall_All.xls",sep = "\t", col.names = NA)
popgenom2 <- F_ST.stats(popgenom2,subsites = "nonsyn")
pair.fst <- get.F_ST(popgenom2,pairwise = TRUE,mode = "nucleotide")
write.table(pair.fst[[1]], file="../PG_FSTs_Pairwise_Nonsyn.xls",sep = "\t", col.names = NA)
popgenom2 <- detail.stats(popgenom2,biallelic.structure=TRUE,mismatch.distribution=TRUE,
                          site.spectrum=TRUE,site.FST=TRUE)

###########
popgenom2 <- calc.fixed.shared(popgenom2)
all.shared <- popgenom2@n.shared.sites
all.fixed <- popgenom2@n.fixed.sites
popgenom2 <- calc.fixed.shared(popgenom2, subsites = "syn")
syn.shared <- popgenom2@n.shared.sites
syn.fixed <- popgenom2@n.fixed.sites
popgenom2 <- calc.fixed.shared(popgenom2, subsites = "nonsyn")
nsyn.shared <- popgenom2@n.shared.sites
nsyn.fixed <-popgenom2@n.fixed.sites
rownames(all.shared) <- popgenom2@region.names
rownames(all.fixed) <- popgenom2@region.names
rownames(syn.shared) <- popgenom2@region.names
rownames(syn.fixed) <- popgenom2@region.names
rownames(nsyn.shared) <- popgenom2@region.names
rownames(nsyn.fixed) <- popgenom2@region.names
write.table(cbind(all.shared,syn.shared,nsyn.shared),"../Shared.xls",sep = "\t",col.names = NA)
write.table(cbind(all.fixed,syn.fixed,nsyn.fixed),"../Fixed.xls",sep = "\t",col.names = NA)

popgenom2@region.stats@minor.allele.freqs

sapply(slide@region.stats@minor.allele.freqs, function(x){if(length(x)==0){return(0)}return(mean(x[3,], na.rm=TRUE))})

# PopGplot()
popgenom2 <- linkage.stats(popgenom2,do.ZnS = TRUE,do.WALL = TRUE,detail = TRUE)
get.linkage(popgenom2)[1]


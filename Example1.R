
##############################################################################
# create genotypes in hap, legend, indv format
# creating CVs
##############################################################################
NCHR <- 3
NSNP <- rep(1000,NCHR) #number SNPs per chromosome
NCV <- rep(100,NCHR) #number CVs per chromosome
NIND <- 1000
VAR.A <- 1
VAR.D <- .1
INCLUDE.CHRS <- 1:NCHR

# We create this map in order to use the real genomic map distance
# genotypes and cvs should be in range of genomic map
map.pos <- matrix(ncol=2,nrow=22)
map.pos[ 1 ,1] <-  738555
map.pos[ 2 ,1] <-  1
map.pos[ 3 ,1] <-  1
map.pos[ 4 ,1] <-  1
map.pos[ 5 ,1] <-  1
map.pos[ 6 ,1] <-  105878
map.pos[ 7 ,1] <-  567276
map.pos[ 8 ,1] <-  64984
map.pos[ 9 ,1] <-  88894
map.pos[ 10 ,1] <-  26070
map.pos[ 11 ,1] <-  102856
map.pos[ 12 ,1] <-  91619
map.pos[ 13 ,1] <-  19198564
map.pos[ 14 ,1] <-  20326742
map.pos[ 15 ,1] <-  22684095
map.pos[ 16 ,1] <-  1263
map.pos[ 17 ,1] <-  1
map.pos[ 18 ,1] <-  12535
map.pos[ 19 ,1] <-  160912
map.pos[ 20 ,1] <-  1
map.pos[ 21 ,1] <-  15107860
map.pos[ 22 ,1] <-  17096300

map.pos[ 1 ,2] <-  249238555
map.pos[ 2 ,2] <-  242900000
map.pos[ 3 ,2] <-  197900000
map.pos[ 4 ,2] <-  1.91e+08
map.pos[ 5 ,2] <-  180750000
map.pos[ 6 ,2] <-  170955878
map.pos[ 7 ,2] <-  159167276
map.pos[ 8 ,2] <-  146364984
map.pos[ 9 ,2] <-  141088894
map.pos[ 10 ,2] <-  135526070
map.pos[ 11 ,2] <-  135002856
map.pos[ 12 ,2] <-  133841619
map.pos[ 13 ,2] <-  115148564
map.pos[ 14 ,2] <-  105826742
map.pos[ 15 ,2] <-  102484095
map.pos[ 16 ,2] <-  90201263
map.pos[ 17 ,2] <-  81100000
map.pos[ 18 ,2] <-  78062535
map.pos[ 19 ,2] <-  59160912
map.pos[ 20 ,2] <-  6.3e+07
map.pos[ 21 ,2] <-  48157860
map.pos[ 22 ,2] <-  51246300






NCHR <- length(NSNP)
cv.info <- c()

for (ichr in 1:NCHR)
{
    print("---------------------------------------")
    nsnp_chr <- NSNP[ichr]
    p <- runif(nsnp_chr,min=.05,max=.95)
    hap <- matrix(0, nrow=NIND, ncol=2*nsnp_chr)
    for (isnp in 1:nsnp_chr)
    {
        hap[,2*isnp-1] <- rbinom(NIND,1,p[isnp])
        hap[,2*isnp] <- rbinom(NIND,1,p[isnp])
    }
    # creating ref.hap file
    write.table(hap, file=paste("ref.chr",ichr,".hap",sep=''), quote=FALSE,row.names=FALSE, col.names=FALSE)
    # creating ref.legend file
    legend.id <- paste("rs",1:nsnp_chr,sep='')
    legend.pos <- sort(sample(map.pos[ichr,1]:map.pos[ichr,2], nsnp_chr, replace = FALSE))
    legend.al0 <- rep("C", nsnp_chr)
    legend.al1 <- rep("T", nsnp_chr)
    write.table(cbind(legend.id,legend.pos,legend.al0,legend.al1), file=paste("ref.chr",ichr,".legend",sep=''), quote=FALSE,row.names=FALSE, col.names=FALSE)
    # creating ref.indv file
    write.table(1:NIND, file=paste("ref.chr",ichr,".indv",sep=''), quote=FALSE,row.names=FALSE, col.names=FALSE)
    ####
    # creating cv.hap file
    ncv_chr <- NCV[ichr]
    cv_chr <- sort(sample(1:nsnp_chr, ncv_chr, replace = FALSE))
    cv_hap_index <- sort(c(2*cv_chr-1,2*cv_chr))
    cvs <- hap[,cv_hap_index]
    write.table(cvs, file=paste("cv.chr",ichr,".hap",sep=''), quote=FALSE,row.names=FALSE, col.names=FALSE)
    cv.pos <- legend.pos[cv_chr]
    cv.maf <- apply(matrix(c(cvs),nrow=2*NIND),2,mean)
    cv.maf[which(cv.maf>.5)] <- 1-cv.maf[which(cv.maf>.5)]
    cv.a <- rnorm(ncv_chr,mean=0,sd=sqrt(VAR.A))
    cv.d <- rnorm(ncv_chr,mean=0,sd=sqrt(VAR.D))
    cv.info.chr <- cbind(ichr,cv.pos,cv.maf,cv.a,cv.d)
    cv.info <- rbind(cv.info,cv.info.chr)
}

colnames(cv.info) <- c("chr","pos","maf","a","d")
write.table(cv.info, file=paste("cv.info",sep=''), quote=FALSE,row.names=FALSE, col.names=TRUE)


########### --file_cvs
b <- paste("cv.chr",INCLUDE.CHRS,".hap",sep='')
b <- cbind(INCLUDE.CHRS,b)
write.table(b,file=paste("par.pop1.cv_hap_files.txt",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE)


########### --file_hap_name
a1 <- paste("ref.chr",INCLUDE.CHRS,".hap",sep="")
a2 <- paste("ref.chr",INCLUDE.CHRS,".legend",sep="")
a3 <- paste("ref.chr",INCLUDE.CHRS,".indv",sep="")
a <- cbind(INCLUDE.CHRS,a1,a2,a3)
write.table(a,file=paste("par.pop1.hap_sample_address.txt",sep=''),quote=FALSE,row.names=FALSE,col.names=c("chr","hap","legend","sample"))



# creating population info

##############################################################################
#DEFINE WILDCARDS
#Must run this section. These need to be changed as you see fit
##############################################################################

NUM.GENERATIONS <- 20    #number of generations
POP.SIZE <- rep(3000,NUM.GENERATIONS)  #population size over time
PHENO.MATE.COR <- rep(.5,NUM.GENERATIONS)   #phenotypic correlation between mates at each generation
OFFSPRING_DIST <- rep("p",NUM.GENERATIONS) # p or P=Poisson distribution, f or F=fixed distribution
OFFSPRING_DIST[NUM.GENERATIONS] <- "f" # last generation is fixed
SELECTION.FUNCTION <- rep("logit",NUM.GENERATIONS) #Selection function for each generation
SELECTION.FUNCTION.PAR1 <- rep(20,NUM.GENERATIONS) #Selection function for each generation
SELECTION.FUNCTION.PAR2 <- rep(0,NUM.GENERATIONS) #Selection function for each generation

#How fine grained should the recombination map be? This has an important effect on speed & memory. The lower this number (in kb) the more RAM and longer the program takes. Numbers < 10 are too fine-grained to make any difference. Typical choices are between 10-50.



##############################################################################
#END WILDCARDS
##############################################################################




###########popinfo.txt
write.table(cbind(POP.SIZE,PHENO.MATE.COR,OFFSPRING_DIST,SELECTION.FUNCTION,SELECTION.FUNCTION.PAR1,SELECTION.FUNCTION.PAR2), file=paste("par.pop1.info.txt",sep=''),quote=FALSE,row.names=FALSE, col.names=c("pop_size","mat_cor","offspring_dist","selection_func","selection_func_par1","selection_func_par2"))










################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
#Run this command (change "path" to appropriate directory)
# you also need to unzip Recom_Map.zip

# /path/GeneEvolve --file_gen_info par.pop1.info.txt --file_hap_name par.pop1.hap_sample_address.txt --file_recom_map Recom.Map.b37.50KbDiff --file_cv_info cv.info --file_cvs par.pop1.cv_hap_files.txt --va 1 --ve 1 --no_output --prefix out1












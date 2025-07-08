library(qtl)
setwd(dir = '~/Documents/impute/')

############# 1. Read in the data from file & check it ############
## File with the GBSA type I genotypes, i.e. GBSA scored markers with alleles completely fixed in founder lines
f2cross <- read.cross(format="csv",dir="/Users/orjanc/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/",file = "F2_n707Rqtl_ID.csv",na.strings = "-",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")

## File with the GBSA type II genotypes, i.e. GBSA scored markers with alleles fixed between founders within F2 family
## To update
#f2cross <- read.cross(format="csv",dir="/Users/orjanc/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/",file = "F2_827_typeII_Rqtl.csv",na.strings = "-",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
#f2cross <- read.cross(format="csv",dir="/Users/orjanc/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/",file = "F2_827_typeII_15_Rqtl.csv",na.strings = "-",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
#f2cross <- read.cross(format="csv",dir="/Users/orjanc/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/",file = "F2_827_typeII_50_Rqtl.csv",na.strings = "-",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
f2cross <- read.cross(format="csv",dir="~/Documents/impute/data/",file = "F2_827.within.fam.permu.180312.txt",na.strings = "NA",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
f2cross <- read.cross(format="csv",dir="~/Documents/impute/data/",file = "F2_827.within.fam.mat.bin50.cut.bugcorrected.180315.cutzw.harsh.csv",na.strings = "NA",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
#F2_827.within.fam.mat.bin50.cut.bugcorrected.180315.cutzw.csv
#f2cross <- read.cross(format="csv",dir="~/Documents/impute/data/",file = "F2_827.within.fam.mat.bin50.cut.bugcorrected.180315.cutzw.csv",na.strings = "NA",genotypes = c("1", "0", "-1", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
f2cross <- read.cross(format="csv",dir="~/Documents/impute/data/",file = "F2_827.within.fam.bin20_180326_allF2.csv",na.strings = "NA",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")

############# 2. Have a look at the imported data #############

## First an overall summary
summary(f2cross)

## Generate a plot to visualize missing genotypes across the markers 
plotMissing(f2cross)

## Visual check for missing genotypes per individual & marker
par(mfrow=c(1,2), las=1)
plot(ntyped(f2cross), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(f2cross, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

############# 3. Do a first cleaning of the data ##########################

######## 3.A. First drop the markers that with genotype scores in few individuals as they complicates map estimation #######
## Here use a threshold of 100 markers based on the visual inspection above: 100 individuals enough to 
## estimate recombination with "reasonable" accuracy

nt.bymar <- ntyped(f2cross, "mar")
todrop <- names(nt.bymar[nt.bymar < 100])
mapthis <- drop.markers(f2cross, todrop)

#### Double-check that the markers have been dropped
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

############# 3.B. Identify individuals that have GBSA scores for few markers ################
#### Here, this is an indication that these individuals had low marker coverage
#### The individual-ID's are given in the phenotype-column of the F2-cross infile

id_lt200mrk <- mapthis$pheno[which(ntyped(mapthis)<300),1] 

## Check the sequence-coverage for the selected individuals
#### First read in the sequence coverage file
coverage_F2 <- read.table(file = './data/180208_F2_coverage.txt', row.names = 1)
#### Print out the coverages
coverage_F2[row.names(coverage_F2) %in% id_lt200mrk,1]
#### As expected, the scoring of few markers was due to the low coverage for these individuals.

## The individuals with few genotyped markers are dropped from the linkage map construction
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>300))
#### Note that this does not indicate that the genotypes are of low quality - can thus be included in later mapping!

#### Check dataset again after the low-coverage individuals have been dropped
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

############ 3.C. Check whether any individuals are duplicates that needs to be removed #############

## Visually inspect the relationships between the individuals in the population
cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
rug(cg[lower.tri(cg)])

##  Identify the individuals that are highly related 
#### Here, a threshold of 0.7 is used for caution, but not necessary
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

## Inspect their genotype similarity
#### This is done by checking the pairwise genotype tables; nothing here looks alarming
g <- pull.geno(mapthis)
table(g[223,], g[238,])
table(g[391,], g[492,])
table(g[391,], g[502,])
table(g[492,], g[502,])
table(g[492,], g[509,])

table(g[390,], g[512,])
table(g[492,], g[512,])
table(g[509,], g[512,])
table(g[390,], g[523,])
table(g[492,], g[523,])
table(g[509,], g[523,])
table(g[512,], g[523,])

table(g[492,], g[539,])
table(g[512,], g[539,])
table(g[523,], g[539,])

table(g[492,], g[676,])
table(g[509,], g[676,])
table(g[512,], g[676,])

table(g[492,], g[690,])
table(g[509,], g[690,])
table(g[512,], g[690,])
table(g[523,], g[690,])
table(g[539,], g[690,])
table(g[676,], g[690,])

########### 3.D. Check whether any markers are duplicates ############
## No duplicates are expected given the way the genotypes are scored & the size of the population
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))

########### 3.E. Check for distorted segregation patterns ###########
## This is an important quality control for this type of markers
## This as it could indicate an error in GBSA genotype assignment despite most individuals having reads for them
gt <- geno.table(mapthis)
## A Bonferroni corrected threshold for the number of markers is used to identify and remove markers
gt[gt$P.value < 0.05/totmar(mapthis),]
todrop <- rownames(gt[gt$P.value < 0.05/totmar(mapthis),])
mapthis <- drop.markers(mapthis, todrop)
## Note1: More stringent filtering could, however, be used to obtain an even higher quality dataset for map construction
## Note2: These markers do not necessarily have to be removed for mapping as a shortage of heterozygotes is expected for markers/bins with few reads

########### 3.F. Evaluate the recombination frequencies for the markers ##########
## This is a way to identify individuals that have skewed genotypes
#### Could be an indication of DNA contamination from a founder
g <- pull.geno(mydata)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))

## We note that some individuals have more genome-wide "AA/AB/BB" genotypes than expected
#### Which are these individuals?
id_gtfreqAA <- mapthis$pheno[which(gfreq[1,]>0.6),1]
id_gtfreqAB <- c(mapthis$pheno[which(gfreq[2,]>0.75),1],mapthis$pheno[which(gfreq[2,]<0.35),1])
id_gtfreqBB <- mapthis$pheno[which(gfreq[3,]>0.6),1]
id_gtfreqALL <- sort(unique(c(id_gtfreqAA,id_gtfreqAB,id_gtfreqBB)))

#### Is it in any way related to coverage? 
coverage_F2[row.names(coverage_F2) %in% id_gtfreqAA,1]
coverage_F2[row.names(coverage_F2) %in% id_gtfreqAB,1]
coverage_F2[row.names(coverage_F2) %in% id_gtfreqBB,1]
#### No, coverage was not low, so check them futher for other problems!
#### Genotype skew is for every chromosome, suggests sample mixup with parent
table(mapthis$geno$`1`$data[which(gfreq[3,]>0.6),])
table(mapthis$geno$`2`$data[which(gfreq[3,]>0.6),])
table(mapthis$geno$`3`$data[which(gfreq[3,]>0.6),])

## Drop these individuals due to the skewed genotype frequencies
#### Individuals to keep
id_gtfreq <- which(gfreq[1,]<0.6 & gfreq[2,]<0.75 & gfreq[2,]>0.35 & gfreq[3,]<0.6)
#### Drop
mapthis <- subset(mapthis, ind=(id_gtfreq))

########## 4. Estimate a first rough draft of the linkage-map ############

####  4.A. Estimate the recombination-fractions between pairs of markers #####
mapthis <- est.rf(mapthis)
mapthis_tst <- markerlrt(mapthis)

## Plot the obtained recombination fractions against LOD-scores
#### Do not want any high recombination-fractions together with high LOD-scores
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#### Or weird linkages indicating that markers are wrongly ordered in the input data
#### such as between chromosomes or orders within chromsomes

## Evaluate by plotting pairwise recombination-fractions and LOD scores across the genome
dev.off()
plotRF(mapthis, alternate.chrid=TRUE)

## Looks pretty OK from the graph - possibly a couple of linkage-groups that could be explored further 
#### later to see if there might be some structural rearrangements on some segregating haplotypes
#### or problematic genotypes in some areas

######## 4.B. Get the complete genetic linkage map  #########
genmap <- est.map(mapthis, error.prob=0.005)
## As numerical summary
summaryMap(genmap)
## As a plot
plotMap(genmap)
mapthis <- replace.map(genmap, newmap)

########### 5. Improve quality of linkage map ###########
## The raw linkage map is too long, suggesting that there might be too high recombination-rates scored in some individuals
## Therefore try to identify these

########### 5.A. Find & remove individuals with excessive recombination
## These could be ones with problems in assigning GBSA scores or sample mixups
par(mfrow=c(1,3))
#col <- rep(c("red","green","black"),65)
col <- "black"
plot(countXO(mapthis), ylab="Number of crossovers",col=col)
hist(countXO(mapthis), ylab="Number of crossovers",col=col)
mean(countXO(mapthis), ylab="Number of crossovers",col=col)

## hich individuals are they?
id_ex_co <- mapthis$pheno[which(countXO(mapthis) > 80),1]
not_ok_id <- mapthis$pheno[which(countXO(mapthis) > 150),1]
ok_id <- mapthis$pheno[which(countXO(mapthis) < 80),1]
okish_id <- mapthis$pheno[which(countXO(mapthis) >=80 & countXO(mapthis) <= 150),1]
save(ok_id,not_ok_id,okish_id,file = "~/Documents/impute/results/ok_and_not_ok_id.RData")
xo_pop <- cbind(mapthis$pheno[which(countXO(mapthis) > 0),1],countXO(mapthis)[which(countXO(mapthis) > 0)])


## Evaluate their sequence overage
coverage_F2[row.names(coverage_F2) %in% id_ex_co,1]
#### Have good overall coverage, hence likely not coverage related

## Bi-modal distribution of crossovers identifies individuals with approximately double the
## amount of recominations than the majority. This suggests mixup of F2 samples!

## Remove these problematic individuals
#### Threshold of 120 recombinations picked arbitrarily from visual inspection.
mapthis2 <- subset(mapthis, ind=(countXO(mapthis) < 120))
rf <- pull.rf(mapthis2)
lod <- pull.rf(mapthis2, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#### Alternatively, remove individuals with index lt 295 as this part of pop seems strange
keep <- logical(length = length(countXO(mapthis)))
keep[296:length(keep)]<- T
mapthis3 <- subset(mapthis, ind=keep)

plot(countXO(mapthis3), ylab="Number of crossovers")
hist(countXO(mapthis3), ylab="Number of crossovers")
mean(countXO(mapthis3), ylab="Number of crossovers")

mapthis3 <- subset(mapthis3, ind=(countXO(mapthis3) < 75))

rf <- pull.rf(mapthis3)
lod <- pull.rf(mapthis3, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#########  6. Create final map with clean dataset ###############
## Build map
newmap2 <- est.map(mapthis2, error.prob=0.05)
summaryMap(newmap2)
dev.off()
plotMap(newmap2)

newmap3 <- est.map(mapthis3, error.prob=0.05)
summaryMap(newmap3)
dev.off()
plotMap(newmap3)

#########  7. Output final map in text and graphics ##############
mapthis <- replace.map(mapthis, newmap2)
summaryMap(mapthis)
dev.off()
plotMap(mapthis)

mapthis_error5 <- mapthis
newmap_error5 <- newmap

#### Additional stuff ########

####### A.1. Test infering linkage groups from the data ########
##   Here this basically splits chromosomes with large distances between markers
lg <- formLinkageGroups(mapthis, max.rf=0.3, min.lod=6) 
table(lg[,2])

###### A.2. Try reordering markers where linkages extend over long regions - possible rearrangements? #######

#### For example chromosome 4 
mapthis <- orderMarkers(mapthis, chr=4)
pull.map(mapthis, chr=4)

###### A.3. Estimate maximum-likelihood estimate of the overall genotyping error from data ############
## Estimate likelihoods
loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02) 
for(i in seq(along=err)) { 
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(mapthis, error.prob=err[i])
  loglik[i] <- do.call(sum,sapply(tempmap, attr, "loglik")) 
}
lod <- (loglik - max(loglik))/log(10)

## Plot likelihoods to identify the most likely overall GBSA genotyping error rate
plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))

###### A.4. Look for genotyping errors in individual data ############
## Useful for final quality fix of data
mapthis <- calc.errorlod(mapthis, error.prob=0.05)
print(toperr <- top.errorlod(mapthis, cutoff=4))

## Manually go through chromosomes
plotGeno(mapthis, chr=1, ind=toperr$id[toperr$chr==1], cutoff=4, include.xo=TRUE,cex = .5)
plotGeno(mapthis, chr=14, ind=toperr$id[toperr$chr==14], cutoff=4, include.xo=TRUE)
mapthis$pheno[31,1]

## Revisit segregation distortion across genome to view these for markers
gt <- geno.table(mapthis, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value"))) > plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")

## Identify marker that seem problematic at Bf-threshold = 2 x 10-4
gt[which(gt$neglog10P > -log10(0.05/269)),]

# SD appears due to i) large number of missing genotypes and ii) missing heterozygotes. 
# Possibly due to difficulty in being able to call heterozytotes?? None overlaps with ind of error ind genotypes
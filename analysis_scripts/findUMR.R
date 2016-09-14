
setwd("/Users/sambuckberry/Desktop/Dropbox/postDoc/projects/hawkinsProject/aggregatedData")

library(data.table)
library(GenomicRanges)
library(rtracklayer)

# As there is random permutation in the analysis, set seed for reporducability
set.seed(123)

## Read in the methylation data and format in genomic ranges
# Function to load the CG data into genomicRanges
readCG <- function(file){
        
        # import the data                
        dat <- fread(file, header=TRUE,
                     colClasses=c("character", "integer", "integer"))
        
        # Split the locus string

        dat <- dat[ ,c("chr", "start") := tstrsplit(locus, "_", fixed=TRUE)]
        
        # Reduce to the required data
        dat <- dat[ ,c(4, 5, 2, 3), with=FALSE]
        
        GRanges(seqnames=Rle(dat$chr),
                ranges=IRanges(start=as.vector(as.numeric(dat$start)),
                               end=as.vector(as.numeric(dat$start))),
                T=as.vector(as.numeric(dat$CT_reads)),
                M=as.vector(as.numeric(dat$C_reads)))
        
}

# combine data for datasets with replicates
readCGreplicate <- function(file1, file2, nReads=5){
        
        # import the data                
        dat <- fread(file1, header=TRUE,
                     colClasses=c("character", "integer", "integer"))
        
        dat2 <- fread(file2, header=TRUE,
                      colClasses=c("character", "integer", "integer"))
        
        # Add the data for common loci
        dat <- merge(dat, dat2, by="locus", all=TRUE)
        dat[is.na(dat)] <- 0
        dat$C_reads.x <- dat$C_reads.x + dat$C_reads.y
        dat$CT_reads.x <- dat$CT_reads.x + dat$CT_reads.y
        dat <- dat[ ,c(1:3), with=FALSE]
        colnames(dat) <- c("locus", "C_reads", "CT_reads")
        
        # Split the locus string
        dat <- dat[ ,c("chr", "start") := tstrsplit(locus, "_", fixed=TRUE)]
        
        # Reduce to the required data
        dat <- dat[ ,c(4, 5, 2, 3), with=FALSE]
        
        #Make the GRanges object
        gr <- GRanges(seqnames=Rle(dat$chr),
                ranges=IRanges(start=as.vector(as.numeric(dat$start)),
                               end=as.vector(as.numeric(dat$start))),
                T=as.vector(as.numeric(dat$CT_reads)),
                M=as.vector(as.numeric(dat$C_reads)))
        
        # Sort the GRanges object
        gr <- sortSeqlevels(gr)
        gr <- sort(gr)
        
        return(gr)
}

# Read CG data
es <- readCG(file="ziller_ES.txt")
ne <- readCG(file="ziller_NE.txt")
erg <- readCG(file="ziller_ERG.txt")
mrg <- readCG(file="ziller_MRG.txt")
nr <- readCGreplicate(file1="hawkins_NR_rep1.txt", file2="hawkins_NR_rep2.txt")
npc <- readCGreplicate(file1="hawkins_project.xie_r1a.aggregatedCG.txt",
                       file2="hawkins_project.xie_r1b.aggregatedCG.txt")

# Run MethylSeekR
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
#biocLite("BSgenome")
#biocLite("MethylSeekR")
library(MethylSeekR)

# Get the chromosome lenghts and add to GRanges object
library('BSgenome.Hsapiens.UCSC.hg19')
sLengths <- seqlengths(Hsapiens)


# Add seqlenghts to objects
addLength <- function(x){
        sLengths <- seqlengths(Hsapiens)
        sLengths <- sLengths[names(sLengths) %in% names(seqlengths(x))]
        seqlengths(x) <- sLengths
        return(x)
}

es <- addLength(es)
ne <- addLength(ne)
erg <- addLength(erg)
mrg <- addLength(mrg)
nr <- addLength(nr)
npc <- addLength(npc)

# Check for PMD's

datasets <- list(es, ne, erg, mrg, nr, npc)

for(i in datasets){
        plotAlphaDistributionOneChr(m=i, chr.sel="chr22", num.cores=1)
}

# Calculate FDRs
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "cpgIslandExt")
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))

# Use the hawkins methylome for FDR calculation
#for(i in datasets){
        stats <- calculateFDRs(m=nr, CGIs=CpGislands.gr, num.cores=1)
#}

FDR.cutoff <- 5
m.sel <- 0.5
n.sel <- as.integer(names(stats$FDRs[as.character(m.sel), ]
 [stats$FDRs[as.character(m.sel), ] <FDR.cutoff]) [1])
n.sel


# GEt UMRs and LMRs

getUMRsLMRs <- function(x){
        UMRLMRsegments.gr <- segmentUMRsLMRs(m=x, 
                                             meth.cutoff=m.sel,
                                             nCpG.cutoff=n.sel,
                                             num.cores=1,
                                             myGenomeSeq=Hsapiens,
                                             seqLengths=sLengths)
        return(UMRLMRsegments.gr)
}

# Get the UMR's and LMR's for each dataset
allLmr <- lapply(X=datasets, getUMRsLMRs)

# Get the UMR and LMR for NRs
umr_NRs <- getUMRsLMRs(nr)
nrResults <- data.frame(chr=seqnames(umr_NRs), start=start(umr_NRs),
                      end=end(umr_NRs), nCG=umr_NRs$nCG, totalReads=umr_NRs$T,
                      methReads=umr_NRs$M, meanMeth=umr_NRs$pmeth,
                      medianMeth=umr_NRs$median.meth, type=umr_NRs$type)
write.table(x=nrResults, file = "NeuralRosette_UMR_LMR_regions_all.txt",
            quote = FALSE, sep = "\t")


# # Reduce to a common set of UMR's and LMRs
reducedLmr <- Reduce(intersect, allLmr)
reducedLmr <- reduce(reducedLmr)

# # Reduce data to loci >= 5 reads across all samples
# es <- es[es$T >= 5]
# ne <- ne[ne$T >= 5]
# erg <- erg[erg$T >= 5]
# mrg <- mrg[mrg$T >= 5]
# nr <- nr[nr$T >= 5]
# npc <- npc[npc$T >= 5]

# Intersect mCG levels with UMR/LMRs to get mCG levels for these regions
es_level <- subsetByOverlaps(query=es, subject=reducedLmr)
ne_level <- subsetByOverlaps(query=ne, subject=reducedLmr)
erg_level <-  subsetByOverlaps(query=erg, subject=reducedLmr)
mrg_level <-  subsetByOverlaps(query=mrg, subject=reducedLmr)
nr_level <-  subsetByOverlaps(query=nr, subject=reducedLmr)
npc_level <-  subsetByOverlaps(query=npc, subject=reducedLmr)


# Put data in 'bed-like' format to calculate the methylation for each region
# Make the data into bed 'like' files
makeBedfromGRanges <- function(data){
        data.table(chr=as.character(seqnames(data)), start=as.numeric(start(ranges(data))),
                   stop=as.numeric(end(ranges(data)))+1,
                   C=data$M, CT=data$T)
}

es_bed <- makeBedfromGRanges(data=es_level)
ne_bed <- makeBedfromGRanges(data=ne_level)
erg_bed <- makeBedfromGRanges(data=erg_level)
mrg_bed <- makeBedfromGRanges(data=mrg_level)
nr_bed <- makeBedfromGRanges(data=nr_level)
npc_bed <- makeBedfromGRanges(data=npc_level)

# Create bed-like table of reduced lmr's
lmr_bed <- data.table(chr=as.character(seqnames(reducedLmr)),
                      start=as.numeric(start(ranges(reducedLmr))),
                      stop=as.numeric(end(ranges(reducedLmr))))

lmr_bed <- lmr_bed[order(lmr_bed$chr, lmr_bed$start), ]


# Bed tools function to calculate methylation for each lmr/umr

### Function to call bedtools
bedTools.2in<-function(functionstring="mapBed",bed1,bed2,opt.string=""){
        
        #create temp files
        a.file=tempfile()
        b.file=tempfile()
        out=tempfile()
        options(scipen =99) # not to use scientific notation when writing out
        
        #write bed formatted dataframes to tempfile
        write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
        write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
        
        # create the command string and call the command using system()
        command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
        cat(command,"\n")
        try(system(command))
        
        # Set the column classes for for downstream calculations
        res=read.table(out,header=F, na.strings=".",
                       colClasses=c("character", "numeric","numeric",
                                    "numeric","numeric"))
        
        # unlink the temp files
        unlink(a.file);unlink(b.file);unlink(out)
        
        # return the results
        return(res)
}

### Function for calculating region average
calcMeth <- function(data){
        
        # call bedtools
        dat <- bedTools.2in(bed1=lmr_bed, b=data[order(data$chr, data$start), ], 
                            opt.string="-c 4,5 -o sum,sum")
        
        # Caclulate percent methylation
        dat$pc <- dat$V4/dat$V5
        
        # output the mCG/CG data
        return(dat$pc)
}

es_mC <- calcMeth(data=es_bed)
ne_mC <- calcMeth(data=ne_bed)
erg_mC <- calcMeth(data=erg_bed)
mrg_mC <- calcMeth(data=mrg_bed)
nr_mC <- calcMeth(data=nr_bed)
npc_mC <- calcMeth(data=npc_bed)

all_mC <- data.frame(es=es_mC, ne=ne_mC, erg=erg_mC, mrg=mrg_mC,
                     nr=nr_mC, npc=npc_mC)

all_mC <- all_mC[complete.cases(all_mC), ]

hist(all_mC$ne)

plot(hclust(dist(t(all_mC))), main="Cluster dendrom of mC levels in UMR/LMR's")

d.correlation <- as.dist(1 - cor(all_mC,method=c("spearman")))
fit <- hclust(d.correlation, method="complete")
plot(fit)


# PCA plot -----------------------------------------------
#install.packages('ggfortify')
library(ggfortify)
autoplot(prcomp(t(all_mC)), label=TRUE, main="PCA plot of mC levels in UMR/LMR's")

pdf(file="lmr_plots.pdf")
plot(hclust(dist(t(all_mC))), main="Cluster dendrom of mC levels in UMR/LMR's")
autoplot(prcomp(t(all_mC)), label=TRUE, main="PCA plot of mC levels in UMR/LMR's")
dev.off()

library(data.table)
library(R.utils)
library(GenomicRanges)
library(rtracklayer)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

file1 <- args[1]
file2 <- args[2]

# combine data for datasets with replicates
readCGreplicate <- function(file1, file2, nReads=1){
        
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
        dat <- dat[ ,c("chr", "start") := tstrsplit(locus, ":", fixed=TRUE)]
        
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
merged <- readCGreplicate(file1,file2)



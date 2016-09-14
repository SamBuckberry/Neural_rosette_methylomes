library(DSS)
library(bsseq)
library(R.utils)
library(data.table)

# Do not output scienfic notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args<-commandArgs(TRUE)

sample1 <- args[1]
sample2 <- args[2]
outFile <- args[3]

# Function for reading and formatting aggregatedCG files
inputFile <- function(file, minCoverage=3){

		# Create the temp directory for uncompressed CGmap file
        tmpFile <- tempfile(tmpdir = tempdir(), fileext = ".aggregatedCG.tmp")

        #gunzip(CGmap, remove=FALSE, destname = tmpFile)
        command <- paste("zcat ", file, " > ", tmpFile, sep = "")
        system(command)
         
        # Read and reformat the data     
        dat <- fread(tmpFile)
        dat <- dat[dat$CT_reads >= minCoverage, ]
        dat <- dat[ ,c("chr", "pos") := tstrsplit(locus, ":", fixed=TRUE)]
        dat <- dat[ ,c(4,5,3,2), with=FALSE]
        dat <- dat[ , pos:=as.integer(pos)]
        setDF(dat)
        colnames(dat) <- c("chr", "pos", "N", "X")

        # Delete the temp file 
        file.remove(tmpFile)

        # Return the DSS formatted data
        return(dat)
}

# Load the data
s1 <- inputFile(sample1)
s2 <- inputFile(sample2)

# Create the BSobject
BSobj <- makeBSseqData(dat = list(s1, s2), sampleNames = c("sample1","sample2"))
rm(s1, s2)

# garbage collection
gc()

# Perform the DML test
dmlTest <- DMLtest(BSobj, group1=c("sample1"), group2=c("sample2"), smoothing=TRUE)
rm(BSobj)

# garbage collection
gc()

# Call DMRs
dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=0.01, minCG=3, dis.merge=100)

# Write output file
write.table(x=dmrs, file=outFile, quote=FALSE, row.names=FALSE)


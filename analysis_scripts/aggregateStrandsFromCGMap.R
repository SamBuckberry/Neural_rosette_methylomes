library(data.table)
library(R.utils)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

CGmap <- args[1]
outFile <- args[2]

aggregateStrandsCG <- function(CGmap){
        
        # Create the temp directory for uncompressed CGmap file
        tmpFile <- tempfile(tmpdir = tempdir(), fileext = ".CGmap.tmp")

        #gunzip(CGmap, remove=FALSE, destname = tmpFile)
        command <- paste("zcat ", CGmap, " > ", tmpFile, sep = "")
        system(command)
              
        # Read the file
        dat <- fread(input = tmpFile, sep = "\t", select = c(1,2,3,5,7,8),
                     col.names = c("chr", "base", "position", "context", 
                                   "C_reads", "CT_reads"))
        
        # Subset to CG context only
        dat <- dat[dat$context == "CG", ]
        
        # Set up the locus id
        dat$locus <- NA
        
        # Get a locus id relative to forward strand
        dat$locus <- ifelse(test = dat$base == "G", 
                            yes = paste(dat$chr, dat$position - 1, sep = ":"),
                            no = paste(dat$chr, dat$position, sep = ":"))
        
        # Drop the unused columns
        dat <- dat[ ,c("chr", "base", "position", "context") := NULL]
        
        # Sum the read counts for + and - strand
        combined <- dat[, lapply(.SD, sum), by=.(locus), .SDcols=c("C_reads", "CT_reads")]
        rm(dat)
        
        # Delete the temp file 
        file.remove(tmpFile)
        
        # return the aggregated data object
        return(combined)
}        

# Apply the function
dat <- aggregateStrandsCG(CGmap = CGmap)

# Write the output file
out <- gzfile(outFile)
write.table(x = dat, file = out, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)


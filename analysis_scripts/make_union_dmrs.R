library(data.table)
library(GenomicRanges)
library(magrittr)
library(IRanges)
library(stringr)
library(pheatmap)

# Load data ---------------------------------------------------------------

# read in the DMR files
dmr_fls <- list.files(path = "processed_data/",
                      pattern = ".dmr", full.names = TRUE)

cg_fls <- list.files(path = "processed_data/",
                     pattern = "aggregated.tsv$", full.names = TRUE)


# Get union of all DMRs ---------------------------------------------------

dmr_to_gr <- function(dmr){
        dat <- fread(dmr)
        gr <- GRanges(seqnames = dat$chr,
                      ranges = IRanges(start = dat$start, end = dat$end))
        return(gr)
}

# Get all DMRs and reduce to union set
all_dmr <- lapply(dmr_fls, dmr_to_gr) %>% unlist() %>%
        GRangesList() %>% unlist() %>% reduce()

# Calc mCG/CG for DMRs ----------------------------------------------------

# GC file to GRanges function
cg_to_gr <- function(cg_file){
        
        # read the aggregated CG data
        dat <- fread(cg_file, header=TRUE,
                     colClasses=c("character", "integer", "integer"))
        
        # Reformat for GRanges object
        dat <- dat[ ,c("chr", "start") := tstrsplit(locus, ":", fixed=TRUE)]
        
        # make the GRanges object
        gr <- GRanges(seqnames = dat$chr,
                      ranges = IRanges(start = as.numeric(dat$start),
                                       end = as.numeric(dat$start)),
                      C_reads=as.numeric(dat$C_reads),
                      CT_reads=as.numeric(dat$CT_reads))
        return(gr)
}

# Calc mCG/CG for union DMRs
get_mC_level <- function(cg_file, intervals){
        
        CGdat <- cg_to_gr(cg_file)
        
        hits <- findOverlaps(query = CGdat, subject = intervals) %>% as.data.frame()
        
        # List unique subject hits (the intervals with data)
        uniqHits <- unique(hits$subjectHits)
        
        # For each subject hit, get the mean mC level
        get_interval_cg_mean <- function(x){
                cgIndex <- hits$queryHits[hits$subjectHits == uniqHits[x]]
                #mean(CGdat[cgIndex]$score) # mean methylation
                mean <- sum(CGdat[cgIndex]$C_reads) / sum(CGdat[cgIndex]$CT_reads) # mCG/CG
                return(mean)
        }
        
        # Get the mean for each interval with a hit
        interval_means <- lapply(X = 1:length(uniqHits),
                                 FUN = get_interval_cg_mean)
        
        
        # Get the intervals with hits and add methylation values
        intervals_with_hits <- intervals[uniqHits]
        score(intervals_with_hits) <- unlist(interval_means)
        
        # Return GRanges object with mean methylation values
        return(intervals_with_hits)
}

all_mcg <- lapply(cg_fls, get_mC_level, intervals=all_dmr)

save(all_mcg, file = "all_mcg_granges.Rda")
load("all_mcg_granges.Rda")

get_dmr_levels <- function(x){
        dat <- all_mcg[x] %>% GRangesList() %>% unlist() %>% data.frame()
        loci <- str_c(dat$seqnames, dat$start, sep = ":") %>%
                str_c(dat$end, sep = "-")
        df <- data.frame(dmr=as.character(loci), mCG=as.numeric(dat$score))
        df <- df[order(df$dmr), ]
        return(df)
}

erg <- get_dmr_levels(1)
nr <- get_dmr_levels(2)
ne <- get_dmr_levels(3)
npc <- get_dmr_levels(4)

# Get the DMRs with coverage in all samples
all_cov <- Reduce(intersect, list(erg$dmr, nr$dmr, ne$dmr, npc$dmr))

erg <- erg[erg$dmr %in% all_cov, ]
nr <- nr[nr$dmr %in% all_cov, ]
ne <- ne[ne$dmr %in% all_cov, ]
npc <- npc[npc$dmr %in% all_cov, ]

all(as.character(erg$dmr) == as.character(nr$dmr))
all(as.character(nr$dmr) == as.character(npc$dmr))
all(as.character(erg$dmr) == as.character(npc$dmr))

all_dmr_levels <- data.frame(row.names = erg$dmr, erg=erg$mCG, nr=nr$mCG,
                             ne=ne$mCG, npc=npc$mCG)

write.csv(all_dmr_levels, file = "processed_data/union_dmr_meth_levels.csv",
          quote = FALSE, row.names = TRUE)

pheatmap::pheatmap(as.matrix(all_dmr_levels), kmeans_k = 1000, 
                   show_rownames = FALSE)




library(pheatmap)
library(ggplot2)
library(reshape2)

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

all_dmr_levels <- data.frame(row.names = erg$dmr, NR=nr$mCG, ERG=erg$mCG, 
                             NE=ne$mCG, NPC=npc$mCG)

write.csv(all_dmr_levels, file = "processed_data/union_dmr_meth_levels.csv",
          quote = FALSE, row.names = TRUE)

pdf("plots/dmr_heatmap.pdf")
pheatmap(as.matrix(all_dmr_levels), kmeans_k = 2000, 
         show_rownames = FALSE)
dev.off()

png("plots/dmr_heatmap.png", width = 7, height = 7, res = 300, units = 'in')
pheatmap(as.matrix(all_dmr_levels), kmeans_k = 2000, 
         show_rownames = FALSE)
dev.off()

# Violin plot mC data

int_melt <- melt(all_dmr_levels)

pdf(file = "plots/dmr_violin_plot.pdf", height = 4)
gp <- ggplot(int_melt, aes(y=value, x=variable)) + 
        geom_violin(scale = "area", alpha=0.3, fill='black') + 
        geom_boxplot(width=.1, notch = TRUE,
                     outlier.size = 0, outlier.shape = 0, outlier.stroke = NA) +
        ylab("DMR mCG/CG") +
        xlab("")
gp
dev.off()

png(file = "plots/dmr_violin_plot.png", height = 4, width = 7, units = 'in', res = 300)
gp <- ggplot(int_melt, aes(y=value, x=variable)) + 
        geom_violin(scale = "area", alpha=0.3, fill='black') + 
        geom_boxplot(width=.1, notch = TRUE,
                     outlier.size = 0, outlier.shape = 0, outlier.stroke = NA) +
        ylab("DMR mCG/CG") +
        xlab("")
gp
dev.off()
#!/bin/bash

Rscript callDMR_DSS.R hawkins_NR_merged_reps_CGmap_aggregated.tsv.gz NE.CGmap_aggregated.tsv.gz NR_NE.dmr &
Rscript callDMR_DSS.R hawkins_NR_merged_reps_CGmap_aggregated.tsv.gz xie_NPC_merged_reps_CGmap_aggregated.tsv.gz NR_NPC.dmr &
Rscript callDMR_DSS.R hawkins_NR_merged_reps_CGmap_aggregated.tsv.gz ERG.CGmap_aggregated.tsv.gz NR_ERG.dmr &
Rscript callDMR_DSS.R NE.CGmap_aggregated.tsv.gz xie_NPC_merged_reps_CGmap_aggregated.tsv.gz NE_NPC.dmr &
Rscript callDMR_DSS.R NE.CGmap_aggregated.tsv.gz ERG.CGmap_aggregated.tsv.gz NE_ERG.dmr &
Rscript callDMR_DSS.R xie_NPC_merged_reps_CGmap_aggregated.tsv.gz ERG.CGmap_aggregated.tsv.gz NPC_ERG.dmr

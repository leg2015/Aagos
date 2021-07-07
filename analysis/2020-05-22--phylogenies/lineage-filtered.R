library(ggplot2)  # (Wickham, 2016)
library(tidyr)    # (Wickham and Henry, 2020)
library(dplyr)    # (Wickham et al., 2020)
library(reshape2) # (Wickham, 2007)
library(cowplot)  # (Wilke, 2019)
library(patchwork)
library(data.table)

theme_set(theme_cowplot())

####################################################################################
# Load/transform data
####################################################################################

chg_rate_label <- function(mag, interval, drift) {
  if (drift) { return("drift") }
  else if (interval == 0) { return("0") }
  else { return(paste(mag, interval, sep="/")) }
}

lineages <- fread("./data-lineages/filtered_lineage.csv")
lineages$BIT_FLIP_PROB <- as.factor(lineages$BIT_FLIP_PROB)
lineages$DRIFT <- lineages$TOURNAMENT_SIZE==1
lineages$chg_rate_label <- factor(mapply(chg_rate_label, 
                                         lineages$CHANGE_MAGNITUDE,
                                         lineages$CHANGE_FREQUENCY,
                                         lineages$DRIFT),
                                  levels=c("drift", "0", "1/256", "1/128",
                                           "1/64", "1/32", "1/16", "1/8",
                                           "1/4", "1/2", "1/1", "2/1", 
                                           "4/1", "8/1", "16/1", "32/1",
                                           "64/1", "128/1", "256/1", 
                                           "512/1", "1024/1", "2048/1",
                                           "4096/1"))
# Some convenient down-filtered data
nk_lineages <- filter(lineages, GRADIENT_MODEL==0)
gradient_lineages <- filter(lineages, GRADIENT_MODEL==1)
lineages <- NULL

# down-filter a bit more...
gradient_lineages_BF003 <- filter(gradient_lineages, BIT_FLIP_PROB==0.003)
nk_lineages_BF003 <- filter(nk_lineages, BIT_FLIP_PROB==0.003)
gradient_lineages_BF03 <- filter(gradient_lineages, BIT_FLIP_PROB==0.03)
nk_lineages_BF03 <- filter(nk_lineages, BIT_FLIP_PROB==0.03)


####################################################################################
# Coding sites over time
####################################################################################

# =========================
# Gradient model, mutation rate 0.003
# =========================
g_bf003 <-
ggplot(gradient_lineages_BF003, 
       aes(x=generation, y=coding_sites, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("Gradient Fitness Model") +
  ggsave("imgs/lineage_coding-sites_gradient_BF-0.003_full.png", width=20, height=10)

# =========================
# Gradient model, mutation rate 0.03
# =========================

g_bf03 <-
  ggplot(gradient_lineages_BF03, 
         aes(x=generation, y=coding_sites, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("Gradient Fitness Model") +
  ggsave("imgs/lineage_coding-sites_gradient_BF-0.03_full.png", width=20, height=10)

# =========================
# NK fitness model, mutation rate 0.003
# =========================

nk_bf003 <-
  ggplot(nk_lineages_BF003, 
         aes(x=generation, y=coding_sites, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("NK Fitness Model") +
  ggsave("imgs/lineage_coding-sites_nk_BF-0.003_full.png", width=20, height=10)

# =========================
# NK fitness model, mutation rate 0.03
# =========================

nk_bf03 <-
  ggplot(nk_lineages_BF03, 
         aes(x=generation, y=coding_sites, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("NK Fitness Model") +
  ggsave("imgs/lineage_coding-sites_nk_BF-0.03_full.png", width=20, height=10)

####################################################################################
# Genome Length Over Time
####################################################################################

# =========================
# Gradient model, mutation rate 0.003
# =========================
g_bf003 <-
  ggplot(gradient_lineages_BF003, 
         aes(x=generation, y=genome_length, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("Gradient Fitness Model") +
  ggsave("imgs/lineage_genome-length_gradient_BF-0.003_full.png", width=20, height=10)


# =========================
# Gradient model, mutation rate 0.03
# =========================
g_bf03 <-
  ggplot(gradient_lineages_BF03, 
         aes(x=generation, y=genome_length, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("Gradient Fitness Model") +
  ggsave("imgs/lineage_genome-length_gradient_BF-0.03_full.png", width=20, height=10)

# =========================
# NK fitness model, mutation rate 0.003
# =========================
n_bf003 <-
  ggplot(nk_lineages_BF003, 
         aes(x=generation, y=genome_length, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("Gradient Fitness Model") +
  ggsave("imgs/lineage_genome-length_nk_BF-0.003_full.png", width=20, height=10)

# =========================
# NK fitness model, mutation rate 0.03
# =========================
n_bf03 <-
  ggplot(nk_lineages_BF03, 
         aes(x=generation, y=genome_length, fill=chg_rate_label, color=chg_rate_label)) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data="mean_cl_boot", 
               fun.args=list(conf.int=0.95), 
               geom="ribbon", 
               alpha=0.2,
               linetype=0) +
  ggtitle("Gradient Fitness Model") +
  ggsave("imgs/lineage_genome-length_nk_BF-0.03_full.png", width=20, height=10)

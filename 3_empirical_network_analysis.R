library(igraph)
library(tidyverse)
library(aricode)
library(rjson)
library(ggplot2)

set.seed(108)

filepaths <- fromJSON(file = "params/filepaths.json")
reanalyze <- F

if(reanalyze) {
  #load("network_data/empirical_network.Rdata")
  load(filepaths$network_path)
  KM_master_tbl <- read_csv(filepaths$KM_master)
  
  E(g)$shared_audience <- E(g)$shared_audience/45
  #walk-trap with self-loop
  g2 <- g
  g2[from=V(g2), to=V(g2)] <- 1
  
  
  KM_master_total <- KM_master_tbl %>%
    select(Month, Media, UV) %>%
    group_by(Media) %>%
    summarize(MeanUV = mean(UV))
  
  for(v in V(g2)$name) {
    E(g2)[v %--% v]$shared_audience <- KM_master_total %>% 
      filter(Media == v) %>% 
      pull(MeanUV)
  }
  
  WT1 <- cluster_walktrap(g, weights = E(g)$shared_audience)
  L1 <- cluster_louvain(g, weights = E(g)$shared_audience)
  FG1 <- cluster_fast_greedy(g, weights = E(g)$shared_audience)
  EB1 <- cluster_edge_betweenness(g, weights = E(g)$shared_audience)
  IM1 <- cluster_infomap(g, e.weights = E(g)$shared_audience)
  LP1 <- cluster_label_prop(g, weights = E(g)$shared_audience)
  LE1 <- cluster_leading_eigen(g, weights = E(g)$shared_audience, options = list(maxiter=1000000))
  SL1 <- cluster_spinglass(g, weights = E(g)$shared_audience)
  
  WT2 <- cluster_walktrap(g2, weights = E(g2)$shared_audience)
  L2 <- cluster_louvain(g2, weights = E(g2)$shared_audience)
  FG2 <- cluster_fast_greedy(g2, weights = E(g2)$shared_audience)
  EB2 <- cluster_edge_betweenness(g2, weights = E(g2)$shared_audience)
  IM2 <- cluster_infomap(g2, e.weights = E(g2)$shared_audience)
  LP2 <- cluster_label_prop(g2, weights = E(g2)$shared_audience)
  LE2 <- cluster_leading_eigen(g2, weights = E(g2)$shared_audience, options = list(maxiter=1000000))
  SL2 <- cluster_spinglass(g2, weights = E(g2)$shared_audience)
  
  save(WT1, WT2, L1, L2, FG1, FG2, EB1, EB2, IM1, IM2, LP1, LP2, LE1, LE2, SL1, SL2, file = "network_data/empirical_network/results2.Rdata")

}

load("network_data/empirical_network/results2.Rdata")
media_types <- read_csv(filepaths$media_types)

# function to calculate Normalized Mutual Information content of a community structure
get_NMI <- function(m_t, com) {
  c_tbl <- tibble(
    Media = com$names,
    Community = com$membership
  )
  
  confusion_tbl <- m_t %>%
    merge(c_tbl)
  
  R <- confusion_tbl %>%
    pull(State) %>%
    unique() %>%
    length()
  
  S <- confusion_tbl %>%
    pull(Community) %>%
    unique() %>%
    length()
  
  scaling_factor <- exp(-abs(R-S)/R)
  
  return(list(
    NMI(confusion_tbl$State, confusion_tbl$Community),
    scaling_factor)
  )
}

media_types <- media_types %>%
  select(Media, State)

# get NMI scores
get_NMI(media_types, WT1)
get_NMI(media_types, WT2)

get_NMI(media_types, L1)
get_NMI(media_types, L2)

get_NMI(media_types, FG1)
get_NMI(media_types, FG2)

get_NMI(media_types, EB1)
get_NMI(media_types, EB2)

get_NMI(media_types, IM1)
get_NMI(media_types, IM2)

get_NMI(media_types, LP1)
get_NMI(media_types, LP2)

get_NMI(media_types, LE1)
get_NMI(media_types, LE2)

get_NMI(media_types, SL1)
get_NMI(media_types, SL2)

# algos <- c("WT", "L", "FG", "EB", "IM", "LP", "LE", "SL")

algos <- c("EB", "FG", "IM", "L", "LE", "LP", "SL", "WT")
algo_types <- c(1, 2)

algo_NMI <- NULL
for(algo in algos) {
  for(a_type in algo_types) {
    algo_NMI <- algo_NMI %>%
      rbind(tibble(
        algorithm = algo,
        algorithm_type = a_type,
        NMI_value = get_NMI(media_types, get(paste0(algo, a_type)))[[1]],
        scaling_factor = get_NMI(media_types, get(paste0(algo, a_type)))[[2]]
        )
      )
  }
}

algo_NMI <- algo_NMI %>%
  rename(algo = algorithm,
         type = algorithm_type) %>%
  mutate(algo = as_factor(algo),
         type = as_factor(type)) %>%
  mutate(SNMI_value = NMI_value * scaling_factor)

# algos_labels <- c("WalkTrap", "Multilevel", "Fast Greedy", "Edge Betweeneness",
#                   "Infomap", "Label Propagation", "Leading Eigenvector", "Spinglass")

algos_labels <- c("Edge Betweenness", "Fast Greedy", "Infomap", "Multilevel",
                  "Leading Eigenvector", "Label Propagation", "Spinglass", "WalkTrap")

names(algos_labels) <- algos

algo_NMI <- algo_NMI %>%
  mutate(type = ifelse(type == 1, "Baseline", "Augmented")) %>%
  mutate(type = factor(type, levels = c("Baseline", "Augmented")))

e_plot_nmi <- ggplot(data = algo_NMI,
       aes(x = type,
           y = NMI_value,
           group = 1)) +
  geom_line(size = 1, linetype = "solid") +
  geom_point(size = 3, aes(color = as_factor(type))) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "red") +
  ylab("NMI") +
  xlab("Network type") +
  facet_wrap(~algo,
             nrow = 2,
             labeller = labeller(algo = algos_labels)) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(
    strip.background = element_rect(
      color="black", fill="black", size=1.5, linetype="solid"
    ),
    strip.text.x = element_text(
      size = 10, color = "white"
    ),
    # panel.grid.major.x = element_line(colour="gray90",size = rel(0.5))
    # scale_x_continuous(breaks = seq(from=1, to=2, by = 0.5)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title=element_text(size=14,face="bold")
  )

ggsave(e_plot_nmi, filename = "plots/empirical_network_NMI.eps", device=cairo_ps)

e_plot_snmi <- ggplot(data = algo_NMI,
                     aes(x = type,
                         y = SNMI_value,
                         group = 1)) +
  geom_line(size = 1, linetype = "solid") +
  geom_point(size = 3, aes(color = as_factor(type))) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "red") +
  ylab("SNMI") +
  xlab("Network type") +
  facet_wrap(~algo,
             nrow = 2,
             labeller = labeller(algo = algos_labels)) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(
    strip.background = element_rect(
      color="black", fill="black", size=1.5, linetype="solid"
    ),
    strip.text.x = element_text(
      size = 10, color = "white"
    ),
    # panel.grid.major.x = element_line(colour="gray90",size = rel(0.5))
    # scale_x_continuous(breaks = seq(from=1, to=2, by = 0.5)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title=element_text(size=14,face="bold")
  ) 

ggsave(e_plot_snmi, filename = "plots/empirical_network_SNMI.eps", device=cairo_ps)
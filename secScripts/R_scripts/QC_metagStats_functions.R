### Metag stats based functions

#' Plots barplot of total reads among samples or individuals
#' @param metag_stats a data frame of metagStats
#' @param column a string indicating if reads should be plotted for each sample 
#' (\code{SMPLID}) or for each individual (\code{AssmblGrps})
#' @param hide_labels logical indicating if sample/individual labels should be hidden
#' @param log10 logical indicating if axis of total reads should be transfromed with log10
plot_totalrds <- function(metag_stats, column = "SMPLID", hide_labels = FALSE, log10 = FALSE) {
  xlabel <- if(column == "SMPLID") {
    "Sample"
  } else if(column == "AssmblGrps") {
    "Individual"
  } else {
    column
  }
  p <- if(log10) {
    ggplot(metag_stats, aes(x = reorder(get(column), log10(totRds)), y = log10(totRds))) +
      geom_col() +
      coord_flip() +
      ylab("log10 (Total reads)")
  } else {
    ggplot(metag_stats, aes(x = reorder(get(column), totRds), y = totRds)) +
      geom_col() +
      coord_flip() +
      ylab("Total reads")
  }
  p2 <- if(hide_labels) {
    p + xlab(paste0(column)) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  } else {
    p + xlab(paste0(column)) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 6))
  }
  p2
}


#' Plots basic assembly and binning statistics as boxplots 
#' @param metag_stats a data frame of metagStats
plot_ass_bin_stats <- function(metag_stats) {
  metag_stats %>% replace(is.na(.), 0) %>%
    select(all_of(c(
      "HQ_bins_SB", "MQ_bins_SB", 
      "ContigN50", "ScaffN50", 
      "CircCtgs", "CircCtgG1M"
    ))) %>%
    rename( 
      CircCtg = CircCtgs, 
      'CircCtg g.t. 1Mbp' = CircCtgG1M 
    ) %>% 
    pivot_longer(1:ncol(.), names_to = "stat", values_to = "value") %>% 
    mutate(Type = case_when(stat == "HQ_bins_SB" ~ "Bins",
                            stat == "MQ_bins_SB" ~ "Bins",
                            stat == "CircCtg" ~ "Circular Contigs",
                            stat == "CircCtg.g.t..1Mbp" ~ "Circular Contigs",
                            stat == "ContigN50" ~ "N50",
                            stat == "ScaffN50" ~ "N50")) %>% 
    ggplot(aes(x = stat, y = value)) +
    geom_boxplot() +
    facet_wrap(~Type, scales = "free", ncol = 3) +
    xlab("") +
    theme_bw()
}


#' Plots basic gene statistics as boxplots 
#' @param metag_stats a data frame of metagStats
plot_gene_stats <- function(metag_stats) {
  metag_stats %>% replace(is.na(.), 0) %>%
    select(all_of(c(
      "AvgGeneLength", "AvgComplGeneLength"
    ))) %>%
    pivot_longer(1:ncol(.), names_to = "stat", values_to = "value") %>% 
    mutate(Type = case_when(stat == "AvgGeneLength" ~ "Gene length",
                            stat == "AvgComplGeneLength" ~ "Gene length"
    )) %>% 
    ggplot(aes(x = stat, y = value)) +
    geom_boxplot() +
    facet_wrap(~Type, scales = "free", ncol = 1) +
    xlab("") +
    theme_bw()
}

#' Plots read statistics as boxplots 
#' @param metag_stats a data frame of metagStats
plot_read_stats <- function(metag_stats) {
  metag_stats %>% replace(is.na(.), 0) %>%
    select(all_of(c(
      "totRds", "FilteredContaRds", 
      "Accepted1", "Accepted2", 
      "Rejected1", "Rejected2"
    ))) %>%
    mutate( 
      Accepted = Accepted1 + Accepted2, 
      Rejected = Rejected1 + Rejected2 
    ) %>% 
    select(-one_of(c(
      "Accepted1", "Rejected1", 
      "Accepted2", "Rejected2"
    ))) %>%
    pivot_longer(1:ncol(.), names_to = "stat", values_to = "value") %>% 
    mutate(Type = case_when(stat == "totRds" ~ "N. Reads",
                            stat == "Accepted" ~ "N. Reads",
                            stat == "Rejected" ~ "N. Reads",
                            stat == "FilteredContaRds" ~ "N. Reads",
                            
    )) %>% 
    ggplot(aes(x = stat, y = value)) +
    geom_boxplot() +
    facet_wrap(~Type, scales = "free", ncol = 1) +
    xlab("") +
    theme_bw()
}



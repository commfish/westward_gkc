
# notes ----
# power analsysis on cl ch data for aigkc crab
# tyler jackson

# load ----

library(segmented)
library(tidyverse)
library(FNGr)

# data ----

aigkc <- read_csv("./data/size_at_maturity/aigkc_chela_complete.csv")

# power analysis function ---- 

# args
## data - data from containing fields specific to chela data of GKC (i.e., cl, ch1)
## sample_size - vector of sample sizes per 5 mm size bin
## effect_size - size of effect (i.e. reduction in size at maturity)
## bin_size - size bin intervals

f_power_analysis <- function(data, sample_size, effect_size, bin_size) {
  
  
  
  ## compute reference sm50
  ref = segmented(lm(ch1~cl, data), seg.Z = ~cl, psi = 130)
  ref_bp = ref$psi[2]
  ref_se = ref$psi[3]
  
  
  f_sig_test <- function(sample_size, effect_size){
    
    data %>%
      mutate(cl_bin = floor(cl / bin_size) * bin_size) %>%
      group_by(cl_bin) %>%
      sample_n(size = sample_size, replace = T) %>%
      ungroup %>%
      # add effect
      mutate(cl_eff = cl + effect_size) -> eff_dat
    # compute fit
    segmented(obj = lm(ch1~cl_eff, data = eff_dat), seg.Z = ~cl_eff, psi = 130) -> mod
    
    ## extract sm50 compute lower bound of 95% ci
    bp_est <- mod$psi[2]
    bp_se <- mod$psi[3]
    
    ## statistical test
    test <- pnorm((bp_est - ref_bp) / bp_se, lower.tail = T) < 0.05

    if(length(test) == 0) {test <- NA}
    return(test)
  }
  
  ## do analysis
  expand_grid(sample_size, effect_size, iteration = 1:1000) %>%
    
    ### two step process due to limits on purrr::map
    mutate(sig = purrr::map2_dbl(sample_size, effect_size, f_sig_test)) %>%
    ## summarise results
    group_by(sample_size, effect_size) %>%
    summarise(power = sum(sig, na.rm = T) / sum(!is.na(sig))) %>%
    ungroup -> res_tab
  
  # plot
  res_tab %>%
    mutate(effect_size_text = factor(paste0(effect_size, " mm"), 
                                     levels = paste0(sort(unique(effect_size)), " mm"))) %>%
    ggplot()+
    geom_point(aes(x = sample_size, y = power))+
    geom_line(aes(x = sample_size, y = power, linetype = effect_size_text))+
    scale_y_continuous(breaks = seq(0, 1, 0.2))+
    labs(x = paste0("Sample size per ", bin_size," mm bin"), y = "Power", linetype = NULL)+
    theme(legend.justification = c(1,0),
          legend.position = c(1,0)) -> plot
  
  return(list(out = res_tab,
              plot = plot))
  
}

# power analysis (special collections) ----

## power analysis only on special collections (20 mm size bins)
aigkc %>%
  filter(source %in% c("adfg_special_collection")) %>%
  f_power_analysis(data = ., sample_size = c(5, 10, 20, 30, 40, 60, 80, 100), effect_size = c(-1, -3, -5, -10), bin_size = 20) %>%
  saveRDS("./output/size_at_maturity/power_analysis_adfg_special_collection_20mm.RDS")

### read results
pa_special_20 <- readRDS("./output/size_at_maturity/power_analysis_adfg_special_collection_20mm.RDS")

### save plot
pa_special_20$plot + 
  scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 60, 80, 100)) +
  theme_sleek() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0)) -> x
ggsave("./figures/size_at_maturity/special_collection_power_20mm_bins.png", plot = x, height = 4, width = 5)


## power analysis only on special collections (5 mm size bins)
aigkc %>%
  filter(source %in% c("adfg_special_collection")) %>%
  f_power_analysis(data = ., sample_size = c(5, 10, 20, 30, 40, 60, 80, 100), effect_size = c(-1, -3, -5, -10), bin_size = 5) %>%
  saveRDS("./output/size_at_maturity/power_analysis_adfg_special_collection_5mm.RDS")

### read results
pa_special_5 <- readRDS("./output/size_at_maturity/power_analysis_adfg_special_collection_5mm.RDS")

### save plot
pa_special_5$plot + 
  scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 60, 80, 100)) +
  theme_sleek() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0)) -> x
ggsave("./figures/size_at_maturity/special_collection_power_5mm_bins.png", plot = x, height = 4, width = 5)



# power analysis (all data) ----

## power analysis only on special collections (20 mm size bins)
aigkc %>%
  # remove two obvious outliers
  filter(!(cl < 100 & ch1 > 30)) %>%
  f_power_analysis(data = ., sample_size = c(5, 10, 20, 30, 40, 60, 80, 100), effect_size = c(-1, -3, -5, -10), bin_size = 20) %>%
  saveRDS("./output/size_at_maturity/power_analysis_all_data_20mm.RDS")

### read results
all_data_20 <- readRDS("./output/size_at_maturity/power_analysis_all_data_20mm.RDS")

### save plot
all_data_20$plot + 
  scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 60, 80, 100)) +
  theme_sleek() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0)) -> x
ggsave("./figures/size_at_maturity/all_data_power_20mm_bins.png", plot = x, height = 4, width = 5)


## power analysis only on special collections (5 mm size bins)
aigkc %>%
  # remove two obvious outliers
  filter(!(cl < 100 & ch1 > 30)) %>%
  f_power_analysis(data = ., sample_size = c(5, 10, 20, 30, 40, 60, 80, 100), effect_size = c(-1, -3, -5, -10), bin_size = 5) %>%
  saveRDS("./output/size_at_maturity/power_analysis_all_data_5mm.RDS")

### read results
all_data_5 <- readRDS("./output/size_at_maturity/power_analysis_all_data_5mm.RDS")

### save plot
all_data_5$plot + 
  scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 60, 80, 100)) +
  theme_sleek() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0)) -> x
ggsave("./figures/size_at_maturity/all_data_power_5mm_bins.png", plot = x, height = 4, width = 5)



# notes ----
# explore existing cl ch data for aigkc crab
# tyler jackson

# load ----

library(segmented)
library(tidyverse)
library(FNGr)
library(sf); library(sp); library(spdplyr)

# wrapper function to call segmented package that works well in pipes
f_seg <- function(data, psi){
  l_mod <- lm(ch1~cl, data)
  seg_mod <- segmented(l_mod, seg.Z = ~cl, psi = psi)
  return(seg_mod)
}

# data ----

aigkc <- read_csv("./data/size_at_maturity/aigkc_chela_complete.csv")

# non-objective way of removing potentially erroneous data ----

aigkc %>%
# step 1: fit segmented regression to all of data, minus the two very obvious outliers
  filter(!(cl < 100 & ch1 > 30)) %>%
  f_seg(data = ., psi = 120) -> all_dat_fit
## extract psi
psi <- all_dat_fit$psi[2]

# step 2: filter initial dataset for observations within 1.5 x the IQR of log(ch1/cl) each side of the all data psi
aigkc %>%
  mutate(size = ifelse(cl <= psi, "imm", "mat")) %>%
  group_by(size) %>% nest %>% 
  
  mutate(data = purrr::map(data, function(data) {
    
    data %>%
      mutate(ch1cl = ch1/cl) %>%
      pull(ch1cl) %>%
      quantile(., c(0.25, 0.75)) %>%
      as.numeric -> quant
    
    data %>%
      filter(ch1/cl > (quant[1] - 1.5 * (quant[2] - quant[1])),
             ch1/cl < (quant[2] + 1.5 * (quant[2] - quant[1]))) -> trim
    
    return(trim)
    
  })) %>%
  
  unnest(data) %>%
  ungroup %>%
  filter(cl < 190) %>%
  dplyr::select(-size) -> aigkc_trim


# examine longitudinal clines in psi ----


raster::getData("GADM", country = c("USA"), level = 1, path = "./data") %>%
  filter(NAME_1 %in% c("Alaska")) %>%
  fortify() %>%
  mutate(long = ifelse(long >= 0, long + -360, long)) -> land

aigkc_trim %>%
  filter(!is.na(lon), !is.na(lat), lat <= 55) %>%
  mutate(lon = ifelse(lon >= 0, lon + -360, lon), 
         lon_bin = floor(lon / 5) * 5) -> points


# map of data
ggplot()+
  geom_polygon(data = land, aes(x = long, y = lat, group = group), fill = "grey60", size = 0.4)+
  geom_point(data = points, aes(x = lon, y = lat), color = "lightblue", alpha = 0.5)+
  coord_map("azequalarea", orientation = c(52.5, -180, 0), xlim = c(-190, -168), ylim = c(51, 55), clip = "on")+
  geom_vline(xintercept = unique(points$lon_bin), linetype = 2)+
  scale_x_continuous(breaks = seq(-190, -170, by = 5),
                     labels = c(expression(170*degree*E), expression(175*degree*E), expression(180*degree), 
                                expression(175*degree*W), expression(170*degree*W)))+
  scale_y_continuous(breaks = seq(51, 55, by = 1),
                     labels = c(expression(51*degree*N), expression(52*degree*N), expression(53*degree*N), 
                                expression(54*degree*N), expression(55*degree*N)))+
  labs(x = "Longitude", y = "Latitude")+
  theme_sleek() -> x

ggsave("./figures/size_at_maturity/eda_sample_map_re_model.png", height = 3, width = 7, units = "in")


# prep data for lme model
aigkc_trim %>%
  filter(!is.na(lon), !is.na(lat), lat <= 55) %>%
  mutate(lon_bin = (floor(lon / 5) * 5),
         lon_bin = ifelse(lon_bin > 0, -360 + lon_bin, lon_bin)) %>%
  arrange(lon_bin, cl)-> lme_data
  
# fit segmented lme
l_mod <- lme(ch1 ~ cl, random = ~1|lon_bin, lme_data)
seg_mod_lme <- segmented.lme(l_mod, seg.Z = ~cl, psi = psi, random=list(lon_bin = pdDiag(~1 + cl + U + G0)))
summary(seg_mod_lme) 

## issues with pulling the fitted values from the model, so need to join with data after sorting
fits <- fitted.segmented.lme(seg_mod_lme, sort = F)

## breakpoints
tibble(lon_bin = as.numeric(names(seg_mod_lme$psi.i)),
       psi = seg_mod_lme$psi.i) -> re_psi

## plot of psi as a function of longitude bin
re_psi %>%
  ggplot()+
  geom_point(aes(x = lon_bin, y = psi))+
  geom_smooth(aes(x = lon_bin, y = psi), method = "lm", color = 1, se = F)+
  scale_x_continuous(labels = c(expression(170*degree*E), expression(175*degree*E), expression(180*degree), 
                                expression(175*degree*W), expression(170*degree*W)))+
  theme_sleek()+
  labs(x = "Longitudinal Bin", y = "Size at Maturity (mm)") -> x
ggsave("./figures/size_at_maturity/re_model_psi_longitude_eda.png", height = 3, width = 4, units = "in")

## plot random effect model fits to data
lme_data %>%
  mutate(fit = fits) %>%
  dplyr::select(lon_bin, cl, ch1, fit) %>%
  ggplot()+
  geom_point(aes(x = cl, y = ch1), color = "grey70", alpha = 0.5)+
  geom_line(aes(x = cl, y = fit), color = 1)+
  facet_wrap(~lon_bin, 
             labeller =labeller(lon_bin = ~paste0(c(170, 175, 180, 175, 170), "°", c(rep("E", 2), "", rep("W", 2)))))+
  labs(x = "Carapace Length (mm)", y = "Standard Carapace Height (mm)")+
  theme_sleek() -> x
ggsave("./figures/size_at_maturity/re_model_fit_longitude_eda.png", height = 5, width = 6, units = "in")



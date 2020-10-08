#Illustrate basic detectability argument using CRC data

library(tidyverse)
library(patchwork)

base_detect_code<-"CCGA2"
#load clean manuscript sensitivity
iso_sens_joined<-read_tsv(sprintf("generated_data/%s_iso_seer_manuscript.tsv",input_date_code))

#remove everything not staged for initial analysis
incidence_sens_source<-iso_sens_joined %>% 
  filter(Stage!="NotStaged") %>%
  mutate(iso_sens=sens)

#keep not staged for adding back
incidence_excluded_source<-iso_sens_joined %>%
  filter(Stage=="NotStaged")

#choose Colon/Rectum to illustrate performance
tmp_sens<-incidence_sens_source %>% filter(Cancer=="Colon/Rectum")

#generate evenly spaced probability individuals
npers<-15
individuals_df<-tibble(id=1:npers,prob=seq(0.5/npers,1-0.5/npers,length=npers))

#join with sensitivities to determine who is detectable at each stage
isolated_plot<-individuals_df %>% 
  mutate(dummy="a") %>% 
  full_join(tmp_sens %>% mutate(dummy="a")) %>% 
  mutate(detected=prob<=sens) %>%
  group_by(id) %>%
  mutate(found=sum(detected)) %>%
  ungroup() %>% 
  mutate(found=factor(found)) %>%
  #  filter(detected==TRUE) %>%
  ggplot(aes(x=Stage,y=prob,group=id))+
  geom_point(aes(color=detected),size=10)+
  geom_line(aes(y=sens),size=2,color="green")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(
    axis.text = element_text(size=14),
    legend.text = element_text(size=16),
    title=element_text(size=18))+
  labs(y="Sensitivity")+
  ggtitle("Cases:distinct individuals")


individual_plot<-individuals_df %>% 
  mutate(dummy="a") %>% 
  full_join(tmp_sens %>% mutate(dummy="a")) %>% 
  mutate(detectable=prob<sens) %>%
  ggplot(aes(x=Stage,y=prob,group=id))+
  geom_line(size=2,color="purple")+
  geom_point(aes(color=detectable),size=10)+
  geom_line(aes(y=sens),size=2,color="green")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(
    axis.text = element_text(size=14),
    legend.text = element_text(size=16),
    title=element_text(size=18))+
  labs(y="Sensitivity")+
  ggtitle("Model: increasing detectability")

individual_explain_grid<-isolated_plot+individual_plot

#explain using crc detection
ggsave(sprintf("figs/%s_supplemental_explain_detectable_model_using_%s.eps",date_code,"crc"),
       individual_explain_grid,
       width=12,
       height=8)
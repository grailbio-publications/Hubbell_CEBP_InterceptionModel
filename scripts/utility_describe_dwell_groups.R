#generate useful summaries discussing the dwell groups
#useful in making tables
#also used for tables in supplemental data

library(tidyverse)
#get my functions
source("R/load_standard_dwell_model.R")

dwell_standard_model<-load_seer_dwell_model()

#human readable short description
text_levels<-c("MIS","VSlow","Slow","Fast","AggFast")


text_cancer_groups<-dwell_standard_model %>% 
  group_by(dwell_group) %>% 
  summarize(cancer=paste(unique(Cancer),collapse=",")) %>%
  ungroup()

text_duration_groups<-dwell_standard_model %>% 
  group_by(dwell_group,scenario,Stage) %>% 
  summarize(mdwell=dwell[1]) %>%
  ungroup() %>%
  pivot_wider(names_from=Stage,values_from=mdwell) %>%
  mutate(scenario=text_levels[scenario+1]) 

text_cancer_groups %>% write_tsv("reports/text_cancer_groups.tsv")
text_duration_groups %>% write_tsv("reports/text_duration_groups.tsv")

compact_duration_groups<-text_duration_groups %>%
  mutate(text_duration=paste(I,II,III,IV,sep=",")) %>%
  select(dwell_group,scenario,text_duration) %>%
  pivot_wider(names_from=scenario,values_from=text_duration)

compact_duration_groups %>% write_tsv("reports/compact_duration_groups.tsv")

summary_duration_groups<-dwell_standard_model %>%
  group_by(scenario,Stage) %>% 
  summarize(dwell_range=paste(c(min(dwell),max(dwell)),collapse="-")) %>%
  ungroup() %>%
  filter(Stage=="I") %>% 
  mutate(aggressive=text_levels[scenario+1]) %>%
  bind_rows(tibble(scenario=0,Stage="I",dwell_range="NONE",aggressive="MIS")) %>%
  arrange(scenario)

summary_duration_groups %>% write_tsv("reports/summary_duration_groups.tsv")


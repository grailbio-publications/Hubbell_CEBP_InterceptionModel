#this generates a modified summary table after excluding cancer cases that are covered by USPSTF screening
#that is, standard of care is the 'best that can be done' for this population in those cancers
#including interval cases and those presenting at late stage
# 100% of Breast, Colorectal, Cervical
# 27% of Lung [*reference*]

library(tidyverse)
source("R/utility_table_format.R")

detect_code<-"CCGA2"

#We're going to just modify the existing interception results rather than re-run interception using a modified sensitivity file
#which would be the other way to do this
rich_survival<-read_tsv(sprintf("reports/%s_%s_text_data_set.tsv",source_date_code,detect_code),guess_max=10000)
scenario_descriptor<-read_tsv("reports/summary_duration_groups.tsv")

text_levels<-c("MIS","VSlow","Slow","Fast","AggFast")

#population ratios for those covered by conventional screening already
conventional_screen<-tibble(Cancer=c("Breast","Colon/Rectum","Cervix","Lung"),
                            SOC_Ratio=c(1,1,1,0.27))

#method: use the no-screening scenario as 'soc'
no_screening_survival<-rich_survival %>% 
  filter(opt=="NO") %>%
  left_join(conventional_screen) %>% 
  mutate(SOC_Ratio=replace_na(SOC_Ratio,0)) %>%
  mutate_at(c("caught",
              "original_survivors",
              "shifted_survivors",
              "original_deaths",
              "shifted_deaths"),list(~ . * .data$SOC_Ratio)) %>%
  select(-SOC_Ratio)

#all other scenarios, changed appropriately for population performance
full_screening_survival<- rich_survival %>% 
  filter(opt!="NO") %>%
  left_join(conventional_screen) %>% 
  mutate(SOC_Ratio=replace_na(1-SOC_Ratio,1)) %>%
  mutate_at(c("caught",
              "original_survivors",
              "shifted_survivors",
              "original_deaths",
              "shifted_deaths"),list(~ . * .data$SOC_Ratio)) %>%
  select(-SOC_Ratio)

#just give me the cases that matter
decor_survival<-full_screening_survival %>%
  select(opt,screen_interval,dw_scenario,scan,aggressive) %>% 
  unique()

#rebuild the file one option at a time
fixed_full_screening<-sapply(decor_survival$opt,function(z){
  full_screening_survival %>% 
    filter(opt==z) %>%
    bind_rows(no_screening_survival) %>%
    group_by(Cancer,clinical,prequel,mode_found) %>%
    summarize_at(c("caught","original_survivors","shifted_survivors","original_deaths","shifted_deaths"),sum) %>%
    ungroup() %>%
    mutate(opt=z)
},simplify=FALSE)  %>% 
  bind_rows() %>%
  left_join(decor_survival) 

final_full_screening<-bind_rows(fixed_full_screening,rich_survival %>% filter(opt=="NO"))

conv_screen_list<-make_pretty_table_from_survival(final_full_screening,scenario_descriptor)

write_tsv(conv_screen_list$incident,
          sprintf("reports/%s_%s_incident_noconv_performance.tsv",date_code,detect_code))
write_tsv(conv_screen_list$prevalent,
          sprintf("reports/%s_%s_prevalent_noconv_performance.tsv",date_code,detect_code))

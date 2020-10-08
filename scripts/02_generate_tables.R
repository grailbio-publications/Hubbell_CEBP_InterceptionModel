#build specific tables for paper
#first table in paper gives format to compare multiple scenarios
library(tidyverse)
#function to generate table entries
source("R/utility_table_format.R")

detect_code<-"CCGA2"

rich_survival<-read_tsv(sprintf("reports/%s_%s_text_data_set.tsv",source_date_code,detect_code),guess_max=10000)

text_levels<-c("MIS","VSlow","Slow","Fast","AggFast")


#use the summary description of scenarios taken from dwell times
scenario_descriptor<-read_tsv("reports/summary_duration_groups.tsv")


text_table_list<-make_pretty_table_from_survival(rich_survival,scenario_descriptor)

#this is table 1 in the paper
write_tsv(text_table_list$incident,
          sprintf("reports/%s_%s_incident_full_performance.tsv",date_code,detect_code))

#this is supplemental
write_tsv(text_table_list$prevalent,
          sprintf("reports/%s_%s_prevalent_full_performance.tsv",date_code,detect_code))

#this is useful for computers
write_tsv(text_table_list$all,
          sprintf("reports/%s_%s_all_performance.tsv",date_code,detect_code))


#individual cancer percentages
table_MIS_proportion_cancers<-rich_survival %>%
  group_by(scan,aggressive,mode_found,Cancer) %>%
  summarize(caught=sum(caught)) %>%
  ungroup() %>%
  pivot_wider(names_from=mode_found,
              values_from=caught,
              values_fill=list(caught=0.0)) %>%
  group_by(scan,aggressive) %>%
  mutate(total_cancer=cfdna+soc,
         proportion_cfdna=proportion_to_percent(cfdna/sum(cfdna)),
         proportion_soc=proportion_to_percent(soc/sum(soc)),
         proportion_total=proportion_to_percent(total_cancer/sum(total_cancer))) %>%
  ungroup() %>%
  mutate(aggressive=factor(aggressive,levels=text_levels)) %>%
  arrange(scan,aggressive,desc(cfdna),desc(soc)) %>%
  mutate_if(is.numeric,function(z){round(z,0)})

#also supplemental
write_tsv(table_MIS_proportion_cancers,sprintf("reports/%s_%s_cancer_percentages.tsv",date_code,detect_code))

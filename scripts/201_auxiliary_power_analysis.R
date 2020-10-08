#extra analysis of approximate scale of experiment needed to detect effects
#this is not a full fledged experimental design analysis
#this is an estimate of scale only
library(pwr)

#read in summary all performance table
#in computer processable form
main_text_table<-read_tsv(sprintf("reports/%s_%s_all_performance.tsv",date_code,detect_code))

#Question: given the relatively robust predictions of incidence round stage shift
#how many individuals need to be sampled to likely detect the effect
my_significance<-0.05
my_needed_power<-0.95
my_population_scale<-1e-5

#using 2-sided test out of an abundance of caution
#in case mortality or late stage cancer could be increased

compute_needed_sample<-function(original,final){
  mapply(function(ox,fx){
  pwr.2p.test(ES.h(ox*my_population_scale,
                   fx*my_population_scale),
              sig.level=my_significance,
              power=my_needed_power)$n
  },original,final)
}

#how many people needed to detect change in late cancer?
late_power_table<-main_text_table %>% 
  filter(scan=="incident") %>%
  select(scan,aggressive,total_late_cancers_original, final_late_cancers) %>%
  mutate(number_needed_detect_late=compute_needed_sample(total_late_cancers_original,final_late_cancers))

#how many people needed to detect change in mortality?
#this ignores important real-world experimental design constraints 
#such as observing mortality over the next 5+ years 
#plus lead time in original diagnosis
#we're computing an approximate scale
mortality_power_table<-main_text_table %>% 
  filter(scan=="incident") %>%
  select(scan,aggressive,total_ied, final_ied) %>%
  mutate(number_needed_detect_mortality=compute_needed_sample(total_ied,final_ied))

total_power_table<-late_power_table %>%
  left_join(mortality_power_table)

write_tsv(total_power_table,
          sprintf("reports/%s_%s_auxilliary_power.tsv",date_code,detect_code))

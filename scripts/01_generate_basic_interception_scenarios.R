#this is to generate the basic scenarios
#models whatever cancer incidence is supplied: currently uses all cancer
#Idealized screening, dwell scenarios 1-4, with annual screening, prevalent and incident results
#generate the big data frame that allows slicing/dicing the data for figure production in the paper
#write the big data frame to a file to ensure summary stats start from a fixed item.

library(tidyverse)
#get my functions
source("R/load_standard_dwell_model.R")
source("R/slip_rate_from_dwell.R")
source("R/run_intercept_model.R")

#date code, input date codes provided from coordination file

#standard parameter set
dwell_standard_model<-load_seer_dwell_model()

#read in standard performance numbers
base_detect_code<-"CCGA2"

#load and clean manuscript sensitivity
iso_sens_joined<-read_tsv(sprintf("generated_data/%s_iso_seer_manuscript.tsv",input_date_code))

#remove everything not staged for initial analysis
incidence_sens_source<-iso_sens_joined %>% 
  filter(Stage!="NotStaged") %>%
  mutate(iso_sens=sens)

#keep not staged for adding back
incidence_excluded_source<-iso_sens_joined %>%
  filter(Stage=="NotStaged")


##dwell time
#annual screening examined here for sensitivity
#generate exact slip rates
dwell_slip_rate<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=1,weibull_shape=1)

#generate prevalent slip rate by clever use of very large interval and multiplying expectation
long_interval<-100
dwell_prevalent_rate<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=long_interval,weibull_shape=1)

#no screening is happening, therefore nothing is ever intercepted
dwell_no_rate<-dwell_slip_rate %>% 
  mutate(slip=1.0,
         slip_clinical=1.0)

#accumulate 4 dwell scenarios
# prevalent and incident results
# plus perfect screening
# plus no screening
my_list<-vector("list",4*2+2)
k<-1
for (dw_scenario in 1:4){
  print(k)
  local_performance<-run_intercept_model(incidence_sens_source,
                                         dwell_slip_rate %>% 
                                           filter(scenario==dw_scenario))
  local_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
  local_performance<-bind_rows(local_performance,local_excluded)
  
  incident_performance<-local_performance
  local_performance<-local_performance %>%
    mutate(screen_interval=1,
           dw_scenario=dw_scenario,
           scan="incident")
  
  my_list[[k]]<-local_performance
  k<-k+1
  
  #generate prevalent round starting off screening
  #only going to use caught by cfdna and combine with incident
  #because expected rates = average over all years, we can reverse identity
  #to obtain first-year screen by multiplying
  #rather than doing a special integral for prevalent rounds
  prevalent_performance<-run_intercept_model(incidence_sens_source,
                                             dwell_prevalent_rate %>% 
                                               filter(scenario==dw_scenario))
  
  prevalent_performance<-prevalent_performance %>% 
    filter(found_clinical==1) %>%
    mutate(caught=caught*long_interval,
           original_survivors=original_survivors*long_interval,
           shifted_survivors=shifted_survivors*long_interval,
           original_deaths=original_deaths*long_interval,
           shifted_deaths=shifted_deaths*long_interval) %>%
    bind_rows(incident_performance %>% 
                filter(found_clinical==2))
  
  prevalent_performance<-prevalent_performance %>%
    mutate(screen_interval=1,
           dw_scenario=dw_scenario,
           scan="prevalent")
  my_list[[k]]<-prevalent_performance
  k<-k+1
}

#this is the MIS scenario where schedule sensitivity is perfect so slip rates are 0
optimal_performance<-run_intercept_model(incidence_sens_source,
                                         dwell_slip_rate %>% 
                                           filter(scenario==0))
optimal_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
optimal_performance<-bind_rows(optimal_performance,optimal_excluded)
optimal_performance<-optimal_performance %>% mutate(opt="0",dw_scenario=0,scan="incident")

#no screening so nothing found by cfdna operations
no_screening_performance<-run_intercept_model(incidence_sens_source,
                                              dwell_no_rate %>%
                                                filter(scenario==0))
no_screening_performance<-bind_rows(no_screening_performance, optimal_excluded) %>%
  mutate(opt="NO",dw_scenario=0,scan="no")

all_options_df<-bind_rows(my_list,.id="opt") %>%
  bind_rows(optimal_performance) %>%
  bind_rows(no_screening_performance)


#Now we have the full, detailed data frame
#add some helper text fields to clarify states represented by each line of the file 
text_options_df<-all_options_df %>%
  select(opt,Cancer,clinical,prequel,found_clinical,
         caught,s_survival,c_survival, 
         original_survivors,shifted_survivors,
         original_deaths,shifted_deaths,
         screen_interval,dw_scenario,scan) %>%
  mutate(mode_found=c("cfdna","soc")[found_clinical],
         aggressive=c("MIS","VSlow","Slow","Fast","AggFast")[dw_scenario+1])

#write to permanent archive of current data
write_tsv(text_options_df,path=sprintf("reports/%s_%s_text_data_set.tsv",date_code,base_detect_code)) 

#this reruns interception using a special detailed model that tracks every flow
library(tidyverse)

source("R/load_standard_dwell_model.R")
source("R/slip_rate_from_dwell.R")
source("R/run_intercept_model.R")
source("R/reconstruct_flow_in_detail.R")
source("R/plot_flow_diagrams.R")


#standard parameter set
dwell_standard_model<-load_seer_dwell_model()
#dwell_standard_model<-load_standard_dwell_model()
#read in standard performance numbers
base_detect_code<-"CCGA2"
iso_sens_joined<-read_tsv(sprintf("generated_data/%s_iso_seer_manuscript.tsv",input_date_code))

#remove everything not staged for initial analysis
incidence_sens_source<-iso_sens_joined %>% 
  filter(Stage!="NotStaged") %>%
  mutate(iso_sens=sens)

#keep not staged for adding back
incidence_excluded_source<-iso_sens_joined %>%
  filter(Stage=="NotStaged")

#incidence_sens_source

slip_rate_df<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=1)

global_figure_scale<-7.5
global_box_scale<-0.045
global_text_scale<-1.35
color_caught=c("intercept"="plum","clinical"="grey80")

#all done
#make one combined figure
{
  diagram_setup<-set_up_diagram_clean(color_caught=color_caught)
  
  setEPS()
  postscript(sprintf("figs/%s_figure_1_ALL.eps",date_code),
             width=2*global_figure_scale,height=2*global_figure_scale)
  par(mfrow=c(2,2))
  
  #1A
  my_locus<-blank_intercept_para(diagram_setup,"State-Transition Graph",local_box_size=global_box_scale)
  
  #simply fill in the boxes with the identity
  text(diagram_setup$map$x,diagram_setup$map$y,diagram_setup$map$name,cex=global_text_scale)
  
  #1B
  dw_scenario<-0  #start with pure scenario
  
  local_slip_rate_df<-slip_rate_df %>% 
    filter(scenario==dw_scenario)
  
  my_title<-"No Interception"
  no_intercept_flow<-intercept_with_flow(incidence_sens_source,local_slip_rate_df,
                                         intercept_start_at_stage=0)
  #start plotting
  my_locus<-blank_intercept_para(diagram_setup,my_title,local_box_size=global_box_scale)
  plot_object_flow_tuned(no_intercept_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=global_text_scale)
  
  
  ##1C
  dw_scenario<-0  #start with pure scenario
  
  local_slip_rate_df<-slip_rate_df %>% 
    filter(scenario==dw_scenario)
  
  my_title<-"Interception Model"
  mis_flow<-intercept_with_flow(incidence_sens_source,local_slip_rate_df,
                                intercept_start_at_stage=4)
  my_locus<-blank_intercept_para(diagram_setup,my_title,local_box_size=global_box_scale)
  plot_object_flow_tuned(mis_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=global_text_scale)
  
  ## 1D
  flow_up_to_stage<-4
  my_title<-"Interception Model: Fast"
  #now use a finite slip rate scenario
  dw_scenario<-3  #start with pure scenario
  local_slip_rate_df<-slip_rate_df %>% 
    filter(scenario==dw_scenario)
  
  fast_flow<-intercept_with_flow(incidence_sens_source,local_slip_rate_df,
                                 intercept_start_at_stage=4)
  #start plotting
  #start plotting
  my_locus<-blank_intercept_para(diagram_setup,my_title,local_box_size=global_box_scale)
  plot_object_flow_tuned(fast_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=global_text_scale)
  
  
  dev.off()
}

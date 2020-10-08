#read in dwell group model 
#hardcoded for paper
load_seer_dwell_model<-function(origin_host_dir="."){

  xStage<-c("I","II","III","IV")
  
  dwell_model_group_df<-read_tsv(sprintf("%s/data/20200728_dwell_time_groups.tsv",origin_host_dir))
  dwell_model_timing_df<-read_tsv(sprintf("%s/data/20200728_dwell_group_timing.tsv",origin_host_dir))
  
  dwell_model_all_df<-dwell_model_group_df %>% 
    rename(dwell_group=group) %>%
    full_join(dwell_model_timing_df %>% rename(dwell_group=group)) %>%
    mutate(number_stage=match(Stage,xStage))
  
  dwell_model_all_df
}

#this date matches the data date from manuscript
iso_sens_joined<-read_tsv("generated_data/20200113_iso_ccga2_train_manuscript.tsv")
#this date comes from the date of SEER reprocessing
seer_draw<-read_tsv(sprintf("generated_data/%s_total_seer_draw.tsv",source_seer_date_code))

iso_join_seer<-seer_draw %>%
  left_join(iso_sens_joined %>%
              mutate(SEER_Draw=Cancer)) %>% 
  mutate(Cancer=replace_na(Cancer,"NotFound"),
         c=replace_na(c,0),
         n=replace_na(n,0),
         sens=replace_na(sens,0.0),
         original_sens=replace_na(original_sens,0.0)) %>%
  select(SEER_Draw,Stage,IR,Survival,c,n,sens) %>%
  rename(Cancer=SEER_Draw)

#data all joined and written using current date code output
write_tsv(iso_join_seer,sprintf("generated_data/%s_iso_seer_manuscript.tsv",date_code))
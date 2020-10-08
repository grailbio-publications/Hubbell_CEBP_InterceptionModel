#join SEER18 data
library(tidyverse)
library(readxl)

#date_code<-"20200520"

SEER_spreadsheet_path<-"data/SEER_Draw_Interception_Manuscript_20200728.xlsx"

#should make sure to output the raw data tables from the spreadsheet for future use
stage_aair_data<-read_excel(SEER_spreadsheet_path,
                            sheet="Incidence",
                            range="A2:E126",
                            col_names=c("SEER_Draw","Stage","IR","Count","Population"))

stage_css_data<-read_excel(SEER_spreadsheet_path,
                           sheet="CSS",
                           range="A4:E2503",
                           col_names=c("SEER_Draw","Stage","TIME","N","CSS"),
                           guess_max=5000)

#filter out "ERROR" codes in SEER results and turn them into NA values
#and return CSS to double type
stage_css_data<-stage_css_data %>%
  mutate(CSS=case_when(grepl("ERROR",CSS)~NA_character_,
                       TRUE ~ CSS)) %>%
  type_convert()


stage_css_filtered_data<-stage_css_data %>%
  mutate(Survival=CSS) %>%
  filter(TIME=="60 mo")

stage_joint_filtered<-stage_aair_data %>% 
  left_join(stage_css_filtered_data) %>% 
  select(SEER_Draw,Stage,IR,Survival)

#have to deal with these cancers specially
hard_cancers<-c("Lymphoid Leukemia","Myeloid Neoplasm","Plasma Cell Neoplasm","[OTHER]")

#okay, deal with typical cases where stage exists, and unknown/missing can be imputed sensibly
limited_joint_filtered<-stage_joint_filtered %>%
  filter(!(SEER_Draw %in% hard_cancers)) 

#impute these
unknown_joint_filtered<-limited_joint_filtered %>%
  filter(Stage=="Unknown/missing") %>%
  group_by(SEER_Draw) %>%
  summarize(UR=sum(IR,na.rm=TRUE)) %>%
  ungroup()

tStage=c("I","II","III","IV","Unknown/missing")

imputed_joint_filtered<-limited_joint_filtered %>%
  filter(Stage %in% tStage[1:4]) %>%
  left_join(unknown_joint_filtered) %>%
  group_by(SEER_Draw) %>%
  mutate(URX=UR*IR/sum(IR,na.rm=TRUE)) %>%
  mutate(URX=replace_na(URX,0.0)) %>%
  ungroup() %>%
  mutate(IR=IR+URX) %>%
  select(-UR,-URX)

#unstaged and expected not to be staged
#need to up-impute "staged" to "notstaged" for lymphoid leukemia
#because we don't have by-stage sensitivities that are relevant to those entries in SEER
#rate is relatively small
unstaged_joint_filtered<-stage_joint_filtered %>%
  filter(SEER_Draw %in% hard_cancers[1:3]) %>%
  group_by(SEER_Draw) %>%
  summarize(
    IR=sum(IR,na.rm=TRUE),
    Survival=Survival[Stage=="Unknown/missing"]) %>%
  ungroup() %>%
  mutate(Stage="NotStaged") %>%
  select(SEER_Draw,Stage,IR,Survival)

#other - heterogenous group so cannot impute unstaged to staged
#also no sensible group-level sensitivity
#but do need incidence and survival data
other_joint_filtered<-stage_joint_filtered %>%
  filter(SEER_Draw=="[OTHER]") %>%
  mutate(Stage=case_when(Stage!="Unknown/missing" ~ Stage,
                         TRUE ~ "NotStaged"))

total_joint_filtered<-bind_rows(imputed_joint_filtered,unstaged_joint_filtered,other_joint_filtered)

#output imputed SEER data
write_tsv(total_joint_filtered,sprintf("generated_data/%s_total_seer_draw.tsv",date_code))



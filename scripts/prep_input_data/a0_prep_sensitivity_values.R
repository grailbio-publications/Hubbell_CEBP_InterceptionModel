#retrieve Liu et al information from data files
#preprocess to generate isotonic sensitivity by stage
library(tidyverse)
library(Iso)


#use cross-validated data for estimates of individual cancer sensitivity
#this has sufficient numbers to be used for individual cancer sensitivity
raw_ccga2_train_sens<-read_csv("data/20200113_external_train_sens.csv")


tStage=c("I","II","III","IV","No Stage")
#fill in missing entries
full_ccga2_train_sens<-raw_ccga2_train_sens %>%
  select(Cancer=cancer_type_tfl,Stage=cstage,c=detec_cancers,n=total_cancers) %>%
  complete(Cancer,Stage=tStage,fill=list(c=0,n=0)) 

#isotone regression for sensitivity by stage
isotone_fix<-function(sens,Stage,num){
  out_val<-sens
  ndx<-match(tStage,Stage) #numbers from 1:5, natural ordering
  good_ndx<-ndx[!is.na(ndx)]
  if (length(good_ndx)>1){
    #need stages I-IV only in order
    y<-sens[good_ndx]
    w<-num[good_ndx]
    val<-pava(y,w)
    out_val[good_ndx]<-val #put back
  } 
  out_val
}

#fix sensitivity with isotonic regression within each cancer type
ccga2_train_manuscript_iso_sens<-full_ccga2_train_sens %>%
  mutate(sensitivity=case_when(n>0 ~ c/n,
                               TRUE ~ 0.0)) %>%
  group_by(Cancer) %>%
  mutate(original_sens = sensitivity,
         sens=isotone_fix(sensitivity,Stage,n)) %>%
  ungroup() %>%
  group_by(Cancer) %>%
  mutate(flag=(sum(n[Stage=="No Stage"])>0))  %>%
  ungroup() %>%
  filter((!flag & Stage!="No Stage") | (flag & Stage=="No Stage")) %>%
  select(Cancer,Stage,c,n,original_sens,sens) %>%
  mutate(Stage=case_when(Stage=="No Stage" ~ "NotStaged",
                         TRUE ~ Stage))

#write to file
#keep original data draw date
write_tsv(ccga2_train_manuscript_iso_sens,"generated_data/20200113_iso_ccga2_train_manuscript.tsv")

#specificity also retrieved from manuscript


# figures and tables for supplemental descriptions of SEER and sensitivity
# parallel plots with numbers: IR/survival (SEER)
# parallel plots with numbers: Sensitivity (isotone)

library(tidyverse)

#read in the basic input data
#read in standard performance numbers
base_detect_code<-"CCGA2"

#load clean manuscript sensitivity
iso_sens_joined<-read_tsv(sprintf("generated_data/%s_iso_seer_manuscript.tsv",input_date_code))

#remove everything not staged for initial analysis
incidence_sens_source<-iso_sens_joined %>% 
  filter(Stage!="NotStaged") %>%
  mutate(iso_sens=sens)

#keep not staged for adding back
incidence_excluded_source<-iso_sens_joined %>%
  filter(Stage=="NotStaged")
#iso_sens_joined contains all data
#incidence_sens_source contains stageable data
#incidence_excluded_source contains non-stageable data

#plot stageable SEER cancers in one block for IR/Survival
#same for sensitivity
#avoids over-plotting of stageI-IV

#generate cancer order required for pretty plots
# generate by expected mortality within 5 years
c_order<-incidence_sens_source %>%
  group_by(Cancer) %>%
  summarize(sum=sum(IR),
            death=sum(IR*(1-Survival))) %>%
  ungroup() %>%
  mutate(Cancer_Num=sprintf("%s(%s/%s)",Cancer,round(sum,0),round(death,0))) #get alternate name with number

#levels
fixed_cancer_order<-c_order$Cancer_Num[order(c_order$death,decreasing=FALSE)]

seer_ir_stageable_plot<-incidence_sens_source %>% 
  left_join(c_order) %>%
  filter(!is.na(Cancer_Num)) %>%
  mutate(Cancer=factor(Cancer_Num,levels=rev(fixed_cancer_order))) %>%
  ggplot(aes(x=Stage,y=IR,group=Cancer))+
  geom_point()+
  geom_blank(aes(y=IR*1.1))+
  geom_line()+
  geom_label(aes(x=Stage,y=IR,label=round(IR,0)),size=8)+
  expand_limits(y=c(-5,40))+
  facet_wrap(~Cancer,ncol=4,scales="free_y")+
  labs(x="Stage",y="Crude Incidence")+
  theme_gray()+
  theme(axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        title=element_text(size=18),
        strip.text=element_text(size=14))+
  ggtitle(sprintf("Incidence (SEER)"))

ggsave(sprintf("figs/%s_supplemental_SEER_IR_stageable.eps",date_code),
       seer_ir_stageable_plot,
       width=12,
       height=12)

seer_survival_stageable_plot<-incidence_sens_source %>% 
  left_join(c_order) %>%
  filter(!is.na(Cancer_Num)) %>%
  mutate(Cancer=factor(Cancer_Num,levels=rev(fixed_cancer_order))) %>%
  ggplot(aes(x=Stage,y=Survival,group=Cancer))+
  geom_point()+
  geom_line()+
  geom_label(aes(x=Stage,y=Survival,label=round(Survival*100,0)),size=8)+
  coord_cartesian(ylim=c(-0.1,1.1))+
  scale_y_continuous(breaks=seq(0,1,by=0.2))+
  facet_wrap(~Cancer,ncol=4)+
  theme_gray()+
  theme(axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        title=element_text(size=18),
        strip.text=element_text(size=14))+
  ggtitle(sprintf("Survival(SEER)",detect_code))

ggsave(sprintf("figs/%s_supplemental_SEER_Survival_stageable.eps",date_code),
       seer_survival_stageable_plot,
       width=12,
       height=12)

seer_sensitivity_stageable_plot<-incidence_sens_source %>% 
  left_join(c_order) %>%
  filter(!is.na(Cancer_Num)) %>%
  mutate(Cancer=factor(Cancer_Num,levels=rev(fixed_cancer_order))) %>%
  mutate(raw_sens=c/n) %>%
  ggplot(aes(x=Stage,y=sens,group=Cancer))+
  geom_point()+
  geom_line()+
  geom_point(aes(y=raw_sens))+
  geom_line(aes(y=raw_sens),lty="dashed")+
  geom_label(aes(x=Stage,y=sens,label=round(sens*100,0)),size=8)+
  coord_cartesian(ylim=c(-0.1,1.1))+
  scale_y_continuous(breaks=seq(0,1,by=0.2))+
  facet_wrap(~Cancer,ncol=4)+
  theme_gray()+
  theme(axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        title=element_text(size=18),
        strip.text=element_text(size=14))+
  ggtitle(sprintf("Sensitivity by stage",detect_code))

ggsave(sprintf("figs/%s_supplemental_SEER_sensitivity_stageable.eps",date_code),
       seer_sensitivity_stageable_plot,
       width=12,
       height=12)

#no plotting of not-stageable for now
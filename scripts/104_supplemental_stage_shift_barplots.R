#generate individual stage-shift plots for supplemental

library(tidyverse)
library(patchwork)

detect_code<-"CCGA2"
#pull in the precomputed interception set
rich_survival<-read_tsv(sprintf("reports/%s_%s_text_data_set.tsv",source_date_code,detect_code), guess_max=10000)

text_levels<-c("MIS","VSlow","Slow","Fast","AggFast")
cStage<-c("NS","I","II","III","IV")

#now use all the scenarios, except, of course the no-screening case
for (scenario_code in c("MIS","Slow","Fast","AggFast")) {
  #plot for incidence rounds
  a_intercept<-rich_survival %>% filter(scan=="incident",aggressive==scenario_code)
  
  #stage-shift
  intercept_shifted<-a_intercept %>% 
    filter(mode_found=="cfdna") %>%
    group_by(Cancer,prequel) %>% 
    summarize(caught=sum(caught)) %>% 
    ungroup() %>%
    mutate(Stage=cStage[prequel+1])
  
  intercept_original<-a_intercept %>% 
    filter(mode_found=="cfdna") %>%
    group_by(Cancer,clinical) %>% 
    summarize(caught=sum(caught)) %>% 
    ungroup() %>%
    mutate(Stage=cStage[clinical+1])
  
  
  #generate horizontal barplot to match figure 3 in main text
  
  stage_shift_in_intercepted_h_plot<-intercept_shifted %>% 
    select(Cancer,Stage,cfdna=caught) %>%
    left_join(intercept_original %>% 
                select(Cancer,Stage,clinical=caught)) %>%
    gather(key="case",value="caught",cfdna,clinical) %>%
    mutate(case=case_when(case=="clinical" ~ "pre-intercept",
                          TRUE ~ "intercepted")) %>%
    group_by(case,Stage) %>%
    summarize(caught=sum(caught)) %>%
    ungroup() %>%
    mutate(Stage=factor(Stage,levels=rev(cStage))) %>%
    ggplot(aes(x=Stage,fill=case,y=caught))+
    geom_col(width=0.9,position="dodge",color="black")+
    geom_text(aes(label=round(caught,0)),position=position_dodge(width=0.9),hjust=-0.1,color="black",size=6)+
    scale_fill_manual(values=c("plum","grey80"))+
    scale_x_discrete(breaks=rev(cStage),drop=FALSE)+
    scale_y_continuous(position="right")+
    theme_bw()+
    theme(
      axis.text = element_text(size=14),
      legend.text = element_text(size=14),
      title=element_text(size=18),
      legend.position=c(0.8,0.9),
      legend.title=element_blank())+
    labs(y="Number of diagnoses per 100K")+
    coord_flip(ylim=c(0,250))+
    ggtitle(sprintf("%s: Intercepted with Stage Shift",scenario_code))
  
  stage_same_in_not_caught_h_plot<-a_intercept %>% 
    filter(mode_found=="soc") %>%
    group_by(clinical) %>%
    summarize(caught=sum(caught)) %>%
    ungroup() %>%
    mutate(case="clinical",Stage=cStage[clinical+1]) %>%
    mutate(Stage=factor(Stage,levels=rev(cStage))) %>%
    mutate(case="usual care") %>%
    ggplot(aes(x=Stage,fill=case,y=caught))+
    geom_col(width=0.9,position="dodge",color="black")+
    geom_text(aes(label=round(caught,0)),position=position_dodge(width=0.9),hjust=-0.1,color="black",size=6)+
    scale_fill_manual(values=c("grey50"))+
    scale_x_discrete(breaks=rev(cStage),drop=FALSE)+
    scale_y_continuous(position="right")+
    theme_bw()+
    theme(
      axis.text = element_text(size=14),
      legend.text = element_text(size=14),
      title=element_text(size=18),
      legend.position=c(0.8,0.9),
      legend.title=element_blank())+
    labs(y="Number of diagnoses per 100K")+
    coord_flip(ylim=c(0,325))+
    ggtitle(sprintf("%s: Remaining Not Intercepted",scenario_code))
  
  #use patchwork to join, label 1 & 2 to match figure 3 in the main text
  joint_h_plot<-stage_shift_in_intercepted_h_plot+stage_same_in_not_caught_h_plot+ plot_annotation(tag_levels=c('1'))  & 
    theme(plot.tag = element_text(size = 16))
  
  ggsave(sprintf("figs/%s_supplemental_stage_shift_%s_h_barplot.eps",date_code,scenario_code),
         joint_h_plot,
         width=12,
         height=8)
}

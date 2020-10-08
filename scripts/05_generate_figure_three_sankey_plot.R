#figure 3
#exact tracking of stage shift
#including flow per stage and resulting change in 5 year mortality
library(tidyverse)
library(ggalluvial)
library(patchwork)

detect_code<-"CCGA2"
source("R/utility_sankey_plot.R")
#read in basic data set
rich_survival<-read_tsv(sprintf("reports/%s_%s_text_data_set.tsv",source_date_code,detect_code), guess_max=10000)

text_levels<-c("MIS","VSlow","Slow","Fast","AggFast")
cStage<-c("NS","I","II","III","IV")

for (scenario_code in c("MIS")) {
  #for each scenario looped over
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
  
  #generate horizontal barplot
  
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
  
  #flip reverse
  joint_h_plot<-stage_shift_in_intercepted_h_plot+stage_same_in_not_caught_h_plot
  

  #reduce to stage found, stage at, mode found
  all_a_summary<-a_intercept %>%
    group_by(mode_found,clinical,prequel) %>%
    summarize(IR=sum(caught),
              Deaths=sum(shifted_deaths),
              Delta=sum(original_deaths),
              scan=scan[1],
              aggressive=aggressive[1]) %>%
    mutate(Delta=Delta-Deaths,
           Survived=IR-Deaths-Delta) %>%
    ungroup()
  
  #generate some automatic scaling for option 0
  ee_plot_range<-all_a_summary %>% 
    group_by(mode_found) %>%
    summarize(v=sum(IR)) %>%
    ungroup() %>%
    summarize(v=max(v)) %>%
    mutate(height_positive=round(v/20+0.5,0)*20,
           height_ticks=pmin(100,round(pmax(height_positive/4,10)/10)*10))
  
  height_positive<-ee_plot_range$height_positive[1]
  height_ticks<-ee_plot_range$height_ticks[1]
  
  found_code<-"cfdna"
  all_b_summary<-all_a_summary %>% filter(mode_found==found_code)
  
  a_plot<-fancy_alluvial_plot_pct_more(all_b_summary,my_title=sprintf("%s: Intercepted With Mortality Shift",scenario_code),
                                       height_positive,height_ticks)
  

  found_code<-"soc"
  all_b_summary<-all_a_summary %>% filter(mode_found==found_code)
  
  b_plot<-fancy_alluvial_plot_pct_more(all_b_summary,my_title=sprintf("%s: Remaining With Mortality",scenario_code),
                                       height_positive,height_ticks)
  
 
  total_figure3_top<-(stage_shift_in_intercepted_h_plot | stage_same_in_not_caught_h_plot )+plot_layout(tag_level='new')
  total_figure3_bottom<-( a_plot | b_plot )+plot_layout(tag_level='new')
  total_figure3_plot<-total_figure3_top / total_figure3_bottom + plot_annotation(tag_levels=c('A','1'))  & 
    theme(
      plot.tag = element_text(size = 16))
  
  ggsave(sprintf("figs/%s_figure_3_combined_%s.pdf",date_code,scenario_code),
         total_figure3_plot,
         width=16,
         height=16)
  
}

#pick my big color set
color_set_sankey<-c("NS"="lightblue",
                    "I"="green",
                    "II"="yellow",
                    "III"="pink",
                    "IV"="red",
                    "Alive"="purple",
                    "Saved"="white",
                    "Dead"="grey")

#translate from numbers to stage
cStage<-c("NS","I","II","III","IV")

fancy_alluvial_plot<-function(my_summary,
                              height_positive=500,
                              height_ticks=50){
  
  #convert from numerical stage to desired symbol
  #build the number found in each unique situation
  #i.e. stage II->stage I -> Alive -> XX number
  #add unique id so sankey knows how to arrange
  #lengthen so that each number is on one line with the status it had in the appropriate column
  ee_sankey<-my_summary %>% 
    select(clinical,prequel,Dead=Deaths,Alive=Survived,Saved=Delta) %>%
    mutate(clinical=cStage[clinical+1],prequel=cStage[prequel+1]) %>%
    pivot_longer(cols=c("Dead","Alive","Saved"),names_to="status",values_to="value") %>%
    mutate(uid=1:length(status)) %>%
    select(uid,StageInSEER=clinical,StageInScenario=prequel,FiveYear=status,value) %>%
    pivot_longer(cols=c("StageInSEER","StageInScenario","FiveYear"),names_to="event",values_to="Status") %>%
    mutate(event=factor(event,levels=c("StageInSEER","StageInScenario","FiveYear"))) %>%
    arrange(uid,event,Status)
  
  saved_sankey_plot<-ee_sankey %>%
    filter(value>0.01) %>%
    mutate(Status=factor(Status,levels=c(cStage,"Alive","Saved","Dead"))) %>%
    ggplot(aes(x=event,stratum=Status,alluvium=uid,y=value,fill=Status,label=Status))+
    scale_x_discrete(expand=c(.1,.1),labels=c("SEER","Intercept","Outcome"))+
    geom_flow(color="black")+
    geom_stratum(alpha=0.65)+
    geom_text(stat="stratum",size=5)+
    scale_fill_manual(values=color_set_sankey)+
    guides(fill=FALSE)+
    labs(x="",y="Individuals Diagnosed With Cancer")+theme_bw()+
    coord_cartesian(ylim=c(0,height_positive))+
    scale_y_continuous(breaks=seq(0,height_positive,height_ticks))+
    theme(strip.text=element_text(size=18),
          axis.text = element_text(size=12),
          axis.title=element_text(size=16))+
    ggtitle(sprintf("Found: %s",found_code))
  
  saved_sankey_plot
}

#include percentage in label
fancy_alluvial_plot_pct<-function(my_summary,
                                  height_positive=500,
                                  height_ticks=50){
  ee_sankey<-my_summary %>% 
    select(clinical,prequel,Dead=Deaths,Alive=Survived,Saved=Delta) %>%
    mutate(clinical=cStage[clinical+1],prequel=cStage[prequel+1]) %>%
    pivot_longer(cols=c("Dead","Alive","Saved"),names_to="status",values_to="value") %>%
    mutate(uid=1:length(status)) %>%
    select(uid,StageInSEER=clinical,StageInScenario=prequel,FiveYear=status,value) %>%
    pivot_longer(cols=c("StageInSEER","StageInScenario","FiveYear"),names_to="event",values_to="Status") %>%
    mutate(event=factor(event,levels=c("StageInSEER","StageInScenario","FiveYear"))) %>%
    arrange(uid,event,Status)
  
  #compute percentages to supply as label
  ee_pct<-ee_sankey %>% 
    group_by(event,Status) %>% 
    summarize(tv=sum(value)) %>%
    ungroup() %>% 
    group_by(event) %>% 
    mutate(pct=round(100*tv/sum(tv))) %>%
    ungroup() %>% 
    mutate(text_pct=sprintf("%s(%s)",Status,pct))
  
  #percentages added
  saved_sankey_plot<-ee_sankey %>%
    filter(value>0.01) %>% 
    left_join(ee_pct %>% select(Status,event,text_pct)) %>%
    mutate(Status=factor(Status,levels=c(cStage,"Alive","Saved","Dead"))) %>% 
    ggplot(aes(x=event,stratum=Status,alluvium=uid,y=value,fill=Status,label=text_pct))+
    scale_x_discrete(expand=c(.1,.1),labels=c("SEER","Intercept","Outcome"))+
    geom_flow(color="black")+
    geom_stratum(alpha=0.65)+
    geom_text(stat="stratum",size=5)+
    scale_fill_manual(values=color_set_sankey)+
    guides(fill=FALSE)+
    labs(x="",y="Individuals Diagnosed With Cancer")+theme_bw()+
    coord_cartesian(ylim=c(0,height_positive))+
    scale_y_continuous(breaks=seq(0,height_positive,height_ticks))+
    theme(strip.text=element_text(size=18),
          axis.text = element_text(size=12),
          axis.title=element_text(size=16))+
    ggtitle(sprintf("Found: %s",found_code))
  
  saved_sankey_plot
}

#include percentage in label
fancy_alluvial_plot_pct_more<-function(my_summary,
                                       my_title="",
                                       height_positive=500,
                                       height_ticks=50){
  ee_sankey<-my_summary %>% 
    select(clinical,prequel,Dead=Deaths,Alive=Survived,Saved=Delta) %>%
    mutate(clinical=cStage[clinical+1],prequel=cStage[prequel+1]) %>%
    pivot_longer(cols=c("Dead","Alive","Saved"),names_to="status",values_to="value") %>%
    mutate(uid=1:length(status)) %>%
    select(uid,StageInSEER=clinical,StageInScenario=prequel,FiveYear=status,value) %>%
    pivot_longer(cols=c("StageInSEER","StageInScenario","FiveYear"),names_to="event",values_to="Status") %>%
    mutate(event=factor(event,levels=c("StageInSEER","StageInScenario","FiveYear"))) %>%
    arrange(uid,event,Status)
  
  #compute percentages to supply as label
  ee_pct<-ee_sankey %>% 
    group_by(event,Status) %>% 
    summarize(tv=sum(value)) %>%
    ungroup() %>% 
    group_by(event) %>% 
    mutate(pct=round(100*tv/sum(tv))) %>%
    ungroup() %>% 
    mutate(text_pct=sprintf("%s(%s)",Status,pct))
  
  #percentages added
  saved_sankey_plot<-ee_sankey %>%
    filter(value>0.01) %>% 
    left_join(ee_pct %>% select(Status,event,text_pct)) %>%
    mutate(Status=factor(Status,levels=c(cStage,"Alive","Saved","Dead"))) %>% 
    ggplot(aes(x=event,stratum=Status,alluvium=uid,y=value,fill=Status,label=text_pct))+
    scale_x_discrete(expand=c(.1,.1),labels=c("SEER","Intercept","Outcome"))+
    geom_flow(color="black")+
    geom_stratum(alpha=0.65)+
    geom_text(stat="stratum",size=5)+
    scale_fill_manual(values=color_set_sankey)+
    guides(fill=FALSE)+
    labs(x="",y="Individuals Diagnosed With Cancer")+theme_bw()+
    coord_cartesian(ylim=c(0,height_positive))+
    scale_y_continuous(breaks=seq(0,height_positive,height_ticks))+
    theme(strip.text=element_text(size=18),
          axis.text = element_text(size=14),
          axis.title=element_text(size=16),
          title=element_text(size=18))+
    ggtitle(my_title)
  
  saved_sankey_plot
}
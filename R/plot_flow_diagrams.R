#make pretty flow pictures
library(diagram)


#picture: 4 diagonal boxes I-IV
# 4 adjacent boxes detectable I-IV
# 4 vertical boxes in a column: detected I-IV
# 4 horizontal boxes in a row: clinical presented

set_up_diagram<-function(){
  #build the big diagonal matrix first
  preclinical_names <- c("NC", "I", "II", "III", "IV")
  detectable_names<-preclinical_names[2:5]
  intercept_names<-detectable_names
  clinical_names<-preclinical_names[2:5]
  
  all_names<-c(preclinical_names,detectable_names,intercept_names,clinical_names)
  
  preclinical_pos<-(matrix(c(c(1,5),c(2,5),c(4,5),c(6,5),c(8,5)),ncol=2,byrow=TRUE)-0.5)/9
  preclinical_pos[1,1]<-preclinical_pos[1,1]
  preclinical_pos[,2]<-preclinical_pos[,2]+0.05
  detectable_pos<-preclinical_pos[2:5,]
  detectable_pos[,1]<-detectable_pos[,1]+0.05
  detectable_pos[,2]<-detectable_pos[,2]+0.1
  intercept_pos<-detectable_pos
  intercept_pos[,1]<-intercept_pos[,1]+0.05
  intercept_pos[,2]<-mean(intercept_pos[,2])+0.15 #top
  
  preclinical_pos[2:5,2]<-preclinical_pos[2:5,2]-0.1
  
  clinical_pos<-detectable_pos
  clinical_pos[,1]<-clinical_pos[,1]+0.05
  clinical_pos[,2]<-preclinical_pos[2:5,2]
  clinical_pos[,2]<-mean(clinical_pos[,2])-0.20 #bottom
  
  all_pos<-rbind(preclinical_pos,detectable_pos,intercept_pos,clinical_pos)
  
  #box type
  box_vec<-c(rep("circle",length(preclinical_names)),
             rep("diamond",length(detectable_names)),
             rep("rect",length(intercept_names)),
             rep("rect",length(clinical_names)))
  
  bigM<-matrix(nrow=length(all_names),ncol=length(all_names),data=0)
  
  #fill in my arrows
  #bigM[1,1]<-"pop"
  #pre-clinical
  bigM[2, 1] <- bigM[3, 2] <- bigM[4, 3] <- bigM[5, 4] <- ""
  #detectable->detectable
  bigM[7,6]<-bigM[8,7]<-bigM[9,8]<-"slip"
  #not-detectable to detectable
  bigM[6,1]<-bigM[7,2]<-bigM[8,3]<-bigM[9,4]<-""
  #detectable->intercepted
  bigM[10,6]<-bigM[11,7]<-bigM[12,8]<-bigM[13,9]<-"intercept"
  #pre-clinical->clinical
  bigM[14,2]<-bigM[15,3]<-bigM[16,4]<-bigM[17,5]<-""
  #detectable->clinical
  bigM[14,6]<-bigM[15,7]<-bigM[16,8]<-bigM[17,9]<-""
  
  box_col<-c(rep("white",length(preclinical_names)),
             rep("white",length(detectable_names)),
             rep("lightgreen",length(intercept_names)),
             rep("pink",length(clinical_names)))
  
  
  #ideally, want to reverse this so it happens at the top
  diagram_map<-tibble(x=all_pos[,1],y=all_pos[,2],shape=box_vec,box_col=box_col,name=all_names)
  diagram_map<-diagram_map %>% 
    mutate(state=case_when(box_col=="lightgreen" ~ "detect",
                           box_col=="pink" ~ "clinical",
                           shape=="diamond" ~ "predetect",
                           shape=="circle" ~ "preclinical",
                           TRUE ~ "")) %>%
    mutate(prequel=match(name,xStage)) %>%
    replace_na(list(prequel=0)) %>%
    mutate(state=case_when(prequel==0 ~ "non_cancer",
                           TRUE ~ state)) %>% mutate(box_id=1:length(x))
  
  list(map=diagram_map,
       matrix=bigM,
       all_pos=all_pos,
       intercept_pos=intercept_pos,
       clinical_pos=clinical_pos)
}

set_up_diagram_clean<-function(color_caught=c("intercept"="lightgreen","clinical"="pink")){
  #build the big diagonal matrix first
  basic_names<-c("NC", "I", "II", "III", "IV")
  preclinical_names <- basic_names
  predetect_names<-basic_names[2:5]
  detect_names<-basic_names[2:5]
  clinical_names<-basic_names[2:5]
  
  all_names<-c(preclinical_names,predetect_names,detect_names,clinical_names)
  
  #set out the grid
  basic_horizontal=(c(1,2,4,6,8)-0.5)/9  #NC,I,II,III,IV
  basic_vertical = c(0.20,0.40,0.50,0.60,0.8) #clinical,preclinical,NC,predetect,detect
  predetect_delta<-0.05
  outcome_delta<-0.1
  
  #untangle: 
  preclinical_pos<-matrix(0,5,2)
  preclinical_pos[,1]<-basic_horizontal
  preclinical_pos[,2]<-basic_vertical[2]
  preclinical_pos[1,2]<-basic_vertical[3] # special case in exact middle
  
  predetect_pos<-matrix(0,4,2)
  predetect_pos[,1]<-basic_horizontal[2:5]+predetect_delta
  predetect_pos[,2]<-basic_vertical[4]
  
  detect_pos<-matrix(0,4,2)
  detect_pos[,1]<-basic_horizontal[2:5]+outcome_delta
  detect_pos[,2]<-basic_vertical[5]
  
  clinical_pos<-matrix(0,4,2)
  clinical_pos[,1]<-basic_horizontal[2:5]+outcome_delta
  clinical_pos[,2]<-basic_vertical[1]
  
  all_pos<-rbind(preclinical_pos,predetect_pos,detect_pos,clinical_pos)
  
  #box type
  box_vec<-c(rep("circle",length(preclinical_names)),
             rep("diamond",length(predetect_names)),
             rep("rect",length(detect_names)),
             rep("rect",length(clinical_names)))
  
  bigM<-matrix(nrow=length(all_names),ncol=length(all_names),data=0)
  
  #fill in my arrows
  #bigM[1,1]<-"pop"
  #pre-clinical
  bigM[2, 1] <- bigM[3, 2] <- bigM[4, 3] <- bigM[5, 4] <- ""
  #detectable->detectable
  bigM[7,6]<-bigM[8,7]<-bigM[9,8]<-"slip"
  #not-detectable to detectable
  bigM[6,1]<-bigM[7,2]<-bigM[8,3]<-bigM[9,4]<-""
  #detectable->intercepted
  bigM[10,6]<-bigM[11,7]<-bigM[12,8]<-bigM[13,9]<-"intercept"
  #pre-clinical->clinical
  bigM[14,2]<-bigM[15,3]<-bigM[16,4]<-bigM[17,5]<-""
  #detectable->clinical
  bigM[14,6]<-bigM[15,7]<-bigM[16,8]<-bigM[17,9]<-""
  
  box_col<-c(rep("white",length(preclinical_names)),
             rep("white",length(predetect_names)),
             rep(color_caught["intercept"],length(detect_names)),
             rep(color_caught["clinical"],length(clinical_names)))
  
  
  #ideally, want to reverse this so it happens at the top
  diagram_map<-tibble(x=all_pos[,1],y=all_pos[,2],shape=box_vec,box_col=box_col,name=all_names)
  diagram_map<-diagram_map %>% 
    mutate(state=case_when(box_col==color_caught["intercept"] ~ "detect",
                           box_col==color_caught["clinical"] ~ "clinical",
                           shape=="diamond" ~ "predetect",
                           shape=="circle" ~ "preclinical",
                           TRUE ~ "")) %>%
    mutate(prequel=match(name,xStage)) %>%
    replace_na(list(prequel=0)) %>%
    mutate(state=case_when(prequel==0 ~ "non_cancer",
                           TRUE ~ state)) %>% mutate(box_id=1:length(x))
  
  list(map=diagram_map,
       matrix=bigM,
       all_pos=all_pos,
       intercept_pos=detect_pos,
       clinical_pos=clinical_pos)
}


#just a blank diagram with boxes and arrows
#to fill in as appropriate
blank_intercept<-function(diagram_setup,my_title="Intercept Model"){
  
  par(mar=c(0,0,3,0))
  #dev.off()
  bigM<-diagram_setup$matrix
  bigV<-as.vector(bigM)
  bigV[nchar(bigV)>1]<-""
  nullM<-matrix(bigV,17,17)
  
  #try with void so counts can be plotted in boxes
  #need to add numbers by arrows as well
  my_locus<-plotmat(nullM, 
                    pos =diagram_setup$all_pos, 
                    curve = 0, 
                    name = "", lwd = 1, box.size=0.04,
                    box.col=diagram_setup$map$box_col,shadow.size=0,
                    box.lwd = 2, cex.txt = 0.8, 
                    box.type = diagram_setup$map$shape, 
                    box.prop = 1.0,
                    self.cex=1.5,
                    main=my_title,cex.main=2)
  
  my_top<-apply(diagram_setup$intercept_pos,2,mean)
  my_top[2]<-my_top[2]+0.1
  text(my_top[1],my_top[2]+0.05,"Intercepted by test",cex=1.5)
  
  my_bottom<-apply(diagram_setup$clinical_pos,2,mean)
  my_bottom[2]<-my_bottom[2]-0.1
  text(my_bottom[1],my_bottom[2],"Found by standard of care",cex=1.5)
  
  #decorate stages
  diagram_top<-diagram_setup$map %>% 
    group_by(name) %>% 
    summarize(x=mean(x),y=max(y)) %>% 
    ungroup() %>%
    mutate(y=y+0.075) %>% 
    filter(grepl("I",name))
  
  #top
  text(diagram_top$x,diagram_top$y,diagram_top$name,cex=1.5)
  segments(diagram_top$x-0.08,diagram_top$y-0.02,diagram_top$x+0.08,diagram_top$y-0.02,lwd=2)
  
  box()
  my_locus
}

#just a blank diagram with boxes and arrows
#to fill in as appropriate
blank_intercept_para<-function(diagram_setup,my_title="Intercept Model",local_box_size=0.04){
  
  par(mar=c(0,0,3,0))
  #dev.off()
  bigM<-diagram_setup$matrix
  bigV<-as.vector(bigM)
  bigV[nchar(bigV)>1]<-""
  nullM<-matrix(bigV,17,17)
  
  #diamonds occupy less space, so boost relative size
  boost_diamond<-(diagram_setup$map$shape=="diamond")*0.2+1
  
  #try with void so counts can be plotted in boxes
  #need to add numbers by arrows as well
  my_locus<-plotmat(nullM, 
                    pos =diagram_setup$all_pos, 
                    curve = 0, 
                    name = "", lwd = 1, box.size=local_box_size*boost_diamond,
                    box.col=diagram_setup$map$box_col,shadow.size=0,
                    box.lwd = 2, cex.txt = 0.8, 
                    box.type = diagram_setup$map$shape, 
                    box.prop = 1.0,
                    self.cex=1.5,
                    main=my_title,cex.main=2)
  
  my_top<-apply(diagram_setup$intercept_pos,2,mean)
  my_top[2]<-my_top[2]+0.1
  text(my_top[1],my_top[2]+0.05,"Intercepted by test",cex=1.5)
  
  my_bottom<-apply(diagram_setup$clinical_pos,2,mean)
  my_bottom[2]<-my_bottom[2]-0.1
  text(my_bottom[1],my_bottom[2],"Found by standard of care",cex=1.5)
  
  #decorate stages
  diagram_top<-diagram_setup$map %>% 
    group_by(name) %>% 
    summarize(x=mean(x),y=max(y)) %>% 
    ungroup() %>%
    mutate(y=y+0.085) %>% 
    filter(grepl("I",name))
  
  #top
  text(diagram_top$x,diagram_top$y,diagram_top$name,cex=1.5)
  segments(diagram_top$x-0.08,diagram_top$y-0.02,diagram_top$x+0.08,diagram_top$y-0.02,lwd=2)
  
  box()
  my_locus
}


#decorate a plot with the flows obtained from a model
# utility: plot up to a given clinical stage to demonstrate build
plot_object_flow<-function(total_flow,diagram_setup,my_locus,flow_up_to_stage=4){
  
  #trick: filter by clinical stage, filter by particular cancer(s)
  #2:11 are columns remaining describing flow - must change if format changes
  object_flow<-total_flow %>% 
    filter(clinical<=flow_up_to_stage) %>%
    #  filter(Cancer=="lung") %>%
    group_by(prequel) %>% 
    summarize_at(vars(-Cancer,-clinical),sum) %>%
    gather(key="object",value="count",2:11)
  
  #map into boxes on the map
  state_flow<-object_flow %>% 
    filter(!grepl("arrow",object)) %>%
    separate(object,into=c("shape","state"),sep="_") %>%
    mutate(shape=case_when(shape=="box" ~ "rect",
                           TRUE ~ shape),
           state=case_when(shape=="diamond" & state=="detect" ~"predetect",
                           TRUE ~ state)) 
  
  coord_flow<-state_flow %>% left_join(diagram_setup$map %>% mutate(box_id=1:length(x)))
  
  text(coord_flow$x,coord_flow$y,round(coord_flow$count,0))
  
  arrow_flow<-object_flow %>% 
    filter(grepl("arrow",object)) %>%
    separate(object,into=c("shape","arrow_state"),sep="_") %>%
    select(-shape)
  
  # how is arrow flow related to state objects
  #arrow slip = from detectable stage prequel to next state
  #arrow evolve = from preclinical stage arrive at prequel
  #arrow arrive = from preclinical state arrive at prequel
  
  arrow_flow<-arrow_flow %>%
    mutate(box_state_end=case_when(arrow_state=="arrive" ~ "predetect",
                                   arrow_state=="slip" & prequel<4 ~ "predetect",
                                   arrow_state=="slip" & prequel==4 ~ "clinical",
                                   arrow_state=="evolve" ~ "preclinical",
                                   arrow_state=="miss" ~ "clinical",
                                   arrow_state=="clinical" ~ "clinical",
                                   arrow_state=="detect" ~ "detect",
                                   TRUE ~ ""),
           box_state_start=case_when(arrow_state=="arrive" & prequel==1 ~ "non_cancer",
                                     arrow_state=="arrive" ~ "preclinical",
                                     arrow_state=="slip" ~ "predetect",
                                     arrow_state=="evolve" & prequel==1 ~ "non_cancer",
                                     arrow_state=="evolve" ~ "preclinical",
                                     arrow_state=="miss" ~ "predetect",
                                     arrow_state=="clinical" ~ "preclinical",
                                     arrow_state=="detect" ~ "predetect",
                                     TRUE ~ ""),
           prequel_start=case_when(arrow_state=="arrive" ~ as.integer(prequel-1),
                                   arrow_state=="slip" ~ prequel,
                                   arrow_state=="evolve" ~ as.integer(prequel-1),
                                   TRUE ~ prequel),
           prequel_end=case_when(arrow_state=="arrive" ~ prequel,
                                 arrow_state=="slip" ~ as.integer(pmin(prequel+1,4)),
                                 arrow_state=="evolve" ~ prequel,
                                 TRUE ~ prequel))
  
  
  arrow_flow<-arrow_flow %>% 
    left_join(diagram_setup$map  %>% select(prequel_start=prequel,box_state_start=state,box_id_start=box_id)) %>%
    left_join(diagram_setup$map %>% select(prequel_end=prequel,box_state_end=state,box_id_end=box_id)) %>%
    filter(!(arrow_state=="slip" & box_state_end=="clinical")) #last slip arrow is actually the miss arrow
  
  
  tidy_arrow<-as_tibble(my_locus$arr)
  
  
  arrow_flow<-arrow_flow %>% 
    left_join(tidy_arrow %>% select(box_id_start=col,box_id_end=row,TextX,TextY,Angle)) %>%
    mutate(muddleY=box_state_start=="preclinical" & box_state_end=="clinical") %>%
    group_by(muddleY) %>%
    mutate(meanY=case_when(muddleY==TRUE ~ mean(TextY),
                           muddleY==FALSE ~ 0.0)) %>%
    ungroup() %>% 
    mutate(meanY=max(meanY)) %>%
    mutate(locX=case_when(Angle==0 ~ TextX,
                          box_state_start=="predetect" & box_state_end=="clinical" ~ TextX+0.03,
                          TRUE ~ TextX-0.04),
           locY=case_when(Angle==0 ~ TextY+0.01,
                          box_state_start=="predetect" & box_state_end=="clinical" ~ meanY+0.0,
                          TRUE ~ TextY+0.0)) %>%
    mutate(locY=case_when(Angle>0 & box_state_start=="non_cancer" ~ TextY+0.02,
                          Angle<0 & box_state_start=="non_cancer" ~ TextY-0.02,
                          TRUE ~ locY))
  
  text(arrow_flow$locX,arrow_flow$locY,round(arrow_flow$count))
  
}

#decorate a plot with the flows obtained from a model
# utility: plot up to a given clinical stage to demonstrate build
plot_object_flow_tuned<-function(total_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=1.0){
  
  #trick: filter by clinical stage, filter by particular cancer(s)
  #2:11 are columns remaining describing flow - must change if format changes
  object_flow<-total_flow %>% 
    filter(clinical<=flow_up_to_stage) %>%
    group_by(prequel) %>% 
    summarize_at(vars(-Cancer,-clinical),sum) %>%
    gather(key="object",value="count",2:11)
  
  #map into boxes on the map
  state_flow<-object_flow %>% 
    filter(!grepl("arrow",object)) %>%
    separate(object,into=c("shape","state"),sep="_") %>%
    mutate(shape=case_when(shape=="box" ~ "rect",
                           TRUE ~ shape),
           state=case_when(shape=="diamond" & state=="detect" ~"predetect",
                           TRUE ~ state)) 
  
  coord_flow<-state_flow %>% left_join(diagram_setup$map %>% mutate(box_id=1:length(x)))
  
  text(coord_flow$x,coord_flow$y,round(coord_flow$count,0),cex=cex_scale*1.2) 
  
  arrow_flow<-object_flow %>% 
    filter(grepl("arrow",object)) %>%
    separate(object,into=c("shape","arrow_state"),sep="_") %>%
    select(-shape)
  
  # how is arrow flow related to state objects
  #arrow slip = from detectable stage prequel to next state
  #arrow evolve = from preclinical stage arrive at prequel
  #arrow arrive = from preclinical state arrive at prequel
  
  arrow_flow<-arrow_flow %>%
    mutate(box_state_end=case_when(arrow_state=="arrive" ~ "predetect",
                                   arrow_state=="slip" & prequel<4 ~ "predetect",
                                   arrow_state=="slip" & prequel==4 ~ "clinical",
                                   arrow_state=="evolve" ~ "preclinical",
                                   arrow_state=="miss" ~ "clinical",
                                   arrow_state=="clinical" ~ "clinical",
                                   arrow_state=="detect" ~ "detect",
                                   TRUE ~ ""),
           box_state_start=case_when(arrow_state=="arrive" & prequel==1 ~ "non_cancer",
                                     arrow_state=="arrive" ~ "preclinical",
                                     arrow_state=="slip" ~ "predetect",
                                     arrow_state=="evolve" & prequel==1 ~ "non_cancer",
                                     arrow_state=="evolve" ~ "preclinical",
                                     arrow_state=="miss" ~ "predetect",
                                     arrow_state=="clinical" ~ "preclinical",
                                     arrow_state=="detect" ~ "predetect",
                                     TRUE ~ ""),
           prequel_start=case_when(arrow_state=="arrive" ~ as.integer(prequel-1),
                                   arrow_state=="slip" ~ prequel,
                                   arrow_state=="evolve" ~ as.integer(prequel-1),
                                   TRUE ~ prequel),
           prequel_end=case_when(arrow_state=="arrive" ~ prequel,
                                 arrow_state=="slip" ~ as.integer(pmin(prequel+1,4)),
                                 arrow_state=="evolve" ~ prequel,
                                 TRUE ~ prequel))
  
  
  arrow_flow<-arrow_flow %>% 
    left_join(diagram_setup$map  %>% select(prequel_start=prequel,box_state_start=state,box_id_start=box_id)) %>%
    left_join(diagram_setup$map %>% select(prequel_end=prequel,box_state_end=state,box_id_end=box_id)) %>%
    filter(!(arrow_state=="slip" & box_state_end=="clinical")) #last slip arrow is actually the miss arrow
  
  
  tidy_arrow<-as_tibble(my_locus$arr)
  
  
  arrow_flow<-arrow_flow %>% 
    left_join(tidy_arrow %>% select(box_id_start=col,box_id_end=row,TextX,TextY,Angle)) %>%
    mutate(muddleY=box_state_start=="preclinical" & box_state_end=="clinical") %>%
    group_by(muddleY) %>%
    mutate(meanY=case_when(muddleY==TRUE ~ mean(TextY),
                           muddleY==FALSE ~ 0.0)) %>%
    ungroup() %>% 
    mutate(meanY=max(meanY)) %>%
    mutate(locX=case_when(Angle==0 ~ TextX+0.01, #horizontal at middle of box
                          box_state_start=="preclinical" & box_state_end=="predetect" ~ TextX-0.035,
                          box_state_start=="predetect" & box_state_end=="clinical" ~ TextX+0.035,
                          box_state_start=="preclinical" & box_state_end=="clinical" ~ TextX-0.05,
                          box_state_start=="predetect" & box_state_end=="detect" ~ TextX-0.045,
                          box_state_start=="non_cancer" & box_state_end=="predetect" ~ TextX-0.05,
                          box_state_start=="non_cancer" & box_state_end=="preclinical" ~ TextX-0.05,
                          TRUE ~ TextX-0.06),
           locY=case_when(Angle==0 ~ TextY+0.02,
                          box_state_start=="preclinical" & box_state_end=="predetect" ~ TextY+0.03,
                          box_state_start=="predetect" & box_state_end=="clinical" ~ meanY, #special: match horizontal
                          box_state_start=="preclinical" & box_state_end=="clinical" ~ meanY,
                          box_state_start=="predetect" & box_state_end=="detect" ~ TextY+0.015,
                          box_state_start=="non_cancer" & box_state_end=="predetect" ~ TextY+0.035,
                          box_state_start=="non_cancer" & box_state_end=="preclinical" ~ TextY-0.035,
                          TRUE ~ TextY+0.0)) 
  
  #scale text up
  text(arrow_flow$locX,arrow_flow$locY,round(arrow_flow$count),cex=cex_scale)
  
}
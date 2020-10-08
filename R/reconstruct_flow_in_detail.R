#reconstruct the flow to make sure everything is fine
#explicitly record graph entries

#constructs the proportional flow for one path leading to clinical presentation normally but intercepted along the way
reconstruct_flow<-function(cumulative_sens, dwell_detect_rate){
  nn<-length(cumulative_sens)
  box_detect<-rep(0,nn)
  diamond_detect<-box_detect
  circle_preclinical<-box_detect
  arrow_slip<-box_detect
  arrow_arrive<-box_detect
  arrow_evolve<-box_detect
  arrow_clinical<-box_detect
  miss<-0
  remain<-1
  i<-1
  
  arrive<-cumulative_sens[i]  #diamond i new content from previous non-detected
  arrow_arrive[i]<-arrive #arrive at diamond i from previous non_detectable state
  live<-arrive+miss  #content of diamond 1
  diamond_detect[i]<-live
  box_detect[i]<-live*dwell_detect_rate[i] #flow to box i
  miss<-live*(1-dwell_detect_rate[i]) #flow to diamond i+1
  remain<-remain-arrive # proportion undetectable at box i
  arrow_evolve[i]<-remain
  circle_preclinical[i]<-remain
  arrow_slip[i]<-miss  #flow to diamond i+1 or slip to clinical
  
  if (nn>1){
    for (i in 2:nn)
    {
      #newly detectable cases
      arrive<-cumulative_sens[i]-cumulative_sens[i-1] 
      arrow_arrive[i]<-arrive
      live<-arrive+miss #currently detectable is newly detectable + missed at earlier stages
      diamond_detect[i]<-live
      box_detect[i]<-live*dwell_detect_rate[i] # flow to box i
      miss<-live*(1-dwell_detect_rate[i]) # flow to diamond i+1
      arrow_slip[i]<-miss
      remain = remain-arrive
      arrow_evolve[i]<-remain
      circle_preclinical[i]<-remain
    }
  }
  
  #status at 'clinical presentation' which is always the last item
  # arrow_evolve + arrow_slip
  box_clinical<-rep(0,nn)
  arrow_clinical<-box_clinical
  arrow_miss=box_clinical
  box_clinical[nn]<-remain+miss #everyone not detectable yet goes to clinical presentation, plus everyone "missed"
  arrow_clinical[nn]<-remain
  arrow_miss[nn]<-miss
  arrow_slip[nn]<-0 #this is the stage at clinical presentation, so no slipping can happen to a higher stage
  #detect[1-4] is detected at each stage
  #clinical detect= miss[4]
  tibble(prequel=1:nn,
             box_detect=box_detect,
             box_clinical=box_clinical,
             diamond_detect=diamond_detect,
             circle_preclinical=circle_preclinical,
             arrow_arrive=arrow_arrive,
             arrow_slip=arrow_slip,
             arrow_evolve=arrow_evolve,
             arrow_miss=arrow_miss,
             arrow_clinical=arrow_clinical, #two ways to get to clinical box
             arrow_detect=box_detect)
}

#get detailed information on flow
compute_effective_detection_flow<-function(incidence_sens_source,dwell_slip_df, 
                                           active_slip_clinical,
                                           intercept_start_at_stage){
  #just detection rate
  just_detection<-incidence_sens_source %>%
    mutate(number_stage=match(Stage,xStage),
           prequel=number_stage,
           detect=iso_sens) %>%
    select(Cancer,prequel,detect)
  
  #differences - marginal detection rate of remaining cases 
  #given that cases already detectable at earlier stage were removed or treated separately
  just_delta<- just_detection %>%
    group_by(Cancer) %>%
    arrange(prequel,.by_group=TRUE) %>%
    mutate(delta_detect=diff(c(0,detect))) %>%
    ungroup() %>%
    arrange(Cancer)
  
  just_slip_delta_detailed<-just_delta %>%
    left_join(dwell_slip_df %>% 
                select(Cancer,prequel=number_stage,slip,slip_clinical), by=c("Cancer","prequel")) %>%
    filter(!is.na(prequel)) %>%
    mutate(unroll=4) %>%
    uncount(unroll,.id="clinical") %>%
    filter(clinical>=prequel) %>%
    mutate(modified_slip=case_when(prequel<clinical ~ slip,
                                   prequel==clinical & active_slip_clinical ~ slip_clinical,
                                   prequel==clinical & !active_slip_clinical ~ slip,
                                   TRUE ~ 1.0)) %>%
    mutate(modified_slip=case_when(prequel<=intercept_start_at_stage ~ modified_slip,
                                   TRUE ~ 1.0)) %>%
    arrange(Cancer,clinical,prequel) %>%
    group_by(Cancer,clinical) %>%
    summarize(sens_slip=list(reconstruct_flow(detect,1-modified_slip))) %>% 
    unnest(cols=c(sens_slip))
  
  just_slip_delta_detailed
}


#run interception model
#but keep detailed flow information
#do not decorate with 5 year mortality from stage shift
#convenience option to turn off interception entirely to get a null model
intercept_with_flow<-function(incidence_sense_source,dwell_slip_df,
                              active_slip_clinical=TRUE,intercept_start_at_stage=0){
  #tricks: set dwell_detect to 0 to obtain the flow in the absence of running screening
  #aggregate individual flows by incidence to get final graph entries
  #set up all previous stages where cases could be intercepted given clinical detection
  incidence_set<-incidence_sens_source %>% 
    filter(Stage %in% xStage) %>%
    select(Cancer,Stage,IR) %>%
    mutate(number_stage=match(Stage,xStage),
           clinical=number_stage,
           unroll=number_stage) %>%
    uncount(unroll,.id="prequel")
  
  #compute effective detection by stage conditional on slip rate model
  just_slip_delta_flow<-compute_effective_detection_flow(incidence_sens_source,
                                                         dwell_slip_df, 
                                                         active_slip_clinical,
                                                         intercept_start_at_stage)
  
  #note 4:13 are flow columns and need to change if flow diagram changes
  #note start with slip which is effectively unrolled
  #apply incidence set - this is the reverse of the operation done in the intercept
  total_flow<-just_slip_delta_flow %>%
    gather(key="state",value="flow",4:13) %>%
    left_join(incidence_set %>% 
                select(Cancer,clinical,prequel,IR),
              by=c("Cancer","clinical","prequel")) %>%
    mutate(flow=flow*IR) %>%
    select(-IR) %>%
    spread(state,flow)
  
  total_flow
}

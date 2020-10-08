xStage<-c("I","II","III","IV")

#okay: this multiplies a final destination IR: all the people who wind up at some final destination
#so everything scales within itself
#assume stage IV destination unless the cases are stopped earlier
effective_sens<-function(cumulative_sens, dwell_detect_rate){
  detect<-rep(0,length(cumulative_sens))
  miss<-0
  i<-1
  arrive<-cumulative_sens[i]
  live<-arrive+miss
  detect[1]<-cumulative_sens[i]*dwell_detect_rate[i]
  miss<-cumulative_sens[i]*(1-dwell_detect_rate[i])
  if (length(cumulative_sens)>1){
    for (i in 2:length(cumulative_sens))
    {
      #newly detectable cases
      arrive<-cumulative_sens[i]-cumulative_sens[i-1]
      live<-arrive+miss #currently detectable is newly detectable + missed at earlier stages
      detect[i]<-live*dwell_detect_rate[i] #would detect all of them, but miss some of them due to timing
      miss<-live*(1-dwell_detect_rate[i])
    }
  }
  arrive<-1-cumulative_sens[i] #final miss
  live=arrive+miss
  #detect[1-4] is detected at each stage
  #clinical detect= miss[4]
  list(intercept=detect,clinical=live)
}

add_survival_to_stage_shift<-function(incidence_sens_source,incidence_intercepted){
  #add survival: original before shift, and shifted survival
  #using 5 year survival as "crude estimate of statistical cure"
  #unlikely to be affected by 1-2 year lead time significantly
  #extract survival by 'stage_at_detection'
  just_survival<-incidence_sens_source %>%
    mutate(prequel=match(Stage,xStage)) %>%
    select(Cancer,prequel,Survival) %>%
    filter(!is.na(prequel))
  
  intercept_survival<-incidence_intercepted %>%
    left_join(just_survival %>%
                select(Cancer,prequel,s_survival=Survival),by=c("Cancer","prequel")) %>% 
    left_join(just_survival %>% 
                select(Cancer,clinical=prequel,c_survival=Survival),by=c("Cancer","clinical")) 
  
  #compute absolute numbers rather than local rates
  intercept_survival<-intercept_survival %>%
    mutate(original_survivors=c_survival*caught,
           shifted_survivors=s_survival*caught,
           original_deaths=(1-c_survival)*caught,
           shifted_deaths=(1-s_survival)*caught)
  intercept_survival
}

compute_effective_detection_with_slip<-function(incidence_sens_source,dwell_slip_df, active_slip_clinical){
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
  
  #modify using slip rate
  #intercept using slip rate
  #slip to next
  
  #include modification of slip rate by clinical stage of detection
  #extra 'parameter'
  just_slip_delta_extra<-just_delta %>%
    left_join(dwell_slip_df %>% 
                #filter(scenario==dw_scenario) %>% 
                select(Cancer,prequel=number_stage,slip,slip_clinical),by=c("Cancer","prequel")) %>%
    filter(!is.na(prequel)) %>%
    mutate(unroll=4) %>%
    uncount(unroll,.id="clinical") %>%
    filter(clinical>=prequel) %>%
    mutate(modified_slip=case_when(prequel<clinical ~ slip,
                                   prequel==clinical & active_slip_clinical ~ slip_clinical,
                                   prequel==clinical & !active_slip_clinical ~ slip,
                                   TRUE ~ 1.0)) %>%
    arrange(Cancer,clinical,prequel) %>%
    group_by(Cancer,clinical) %>%
    mutate(sens_slip=effective_sens(detect,1-modified_slip)$intercept) %>%
    ungroup()
  
  just_slip_delta_extra
}




run_intercept_model<-function(incidence_sens_source, dwell_slip_df, active_slip_clinical=TRUE){
  
  #set up all previous stages where cases could be intercepted given clinical detection
  incidence_set<-incidence_sens_source %>% 
    filter(Stage %in% xStage) %>%
    select(Cancer,Stage,IR) %>%
    mutate(number_stage=match(Stage,xStage),
           clinical=number_stage,
           unroll=number_stage) %>%
    uncount(unroll,.id="prequel")
  
  #compute effective detection by stage conditional on slip rate model
  just_slip_delta_extra<-compute_effective_detection_with_slip(incidence_sens_source,
                                                               dwell_slip_df, 
                                                               active_slip_clinical)
  
  #updated: split effective slip rate in 2 for last stage
  #as the "time spent" should be halved
  #this involves a more elaborate model
  #note that "lives saved" is not affected, because those individuals are not stage shifted
  #this assumes that 'stage 4' is just automatically halved anyway
  incidence_intercepted<-incidence_set %>%
    left_join(just_slip_delta_extra,by=c("Cancer","clinical","prequel")) %>% 
    mutate(unroll=1+(number_stage==prequel)) %>%
    uncount(unroll,.id="found_clinical") %>%
    group_by(Cancer,clinical) %>%
    mutate(c_slip=cumsum(sens_slip),
           delta_detect=case_when(
             found_clinical==2 ~ 1-c_slip+sens_slip, #anyone not caught by new screening must be found clinically
             TRUE ~ sens_slip)) %>%
    mutate(caught=IR*delta_detect) %>%
    ungroup()
  
  intercept_survival<-add_survival_to_stage_shift(incidence_sens_source,incidence_intercepted)
  
  intercept_survival
}


run_excluded_model<-function(excluded_source){
  #fills out individuals not staged as they are not modeled
  excluded_survival<-excluded_source %>%
    mutate(number_stage=0,
           clinical=0,
           prequel=0,
           detect=0.0,
           delta_detect=0.0,
           slip=1.0,
           slip_clinical=1.0,
           modified_slip=1.0,
           sens_slip=0.0,
           found_clinical=2,
           c_slip=1.0,
           caught=IR,
           s_survival=Survival,
           c_survival=Survival) %>%
    mutate(original_survivors=c_survival*caught,
           shifted_survivors=s_survival*caught,
           original_deaths=(1-c_survival)*caught,
           shifted_deaths=(1-s_survival)*caught) %>%
    select(Cancer,Stage,IR,
           number_stage,clinical,prequel,
           detect,delta_detect,
           slip,slip_clinical,modified_slip,sens_slip,
           found_clinical,c_slip,
           caught,
           s_survival,c_survival,
           original_survivors,shifted_survivors,
           original_deaths,shifted_deaths)
  
  excluded_survival
 }

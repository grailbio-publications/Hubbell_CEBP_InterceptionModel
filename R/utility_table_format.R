#updated table: reformat for clarity

proportion_to_percent<-function(z){100*z}

pretty_string_plus_percent<-function(a_number,a_percent,number_sig=0,percent_sig=0){
  sprintf("%s(%s)",round(a_number,number_sig),round(a_percent,percent_sig))
}

#my_survival = rich survival, note specific dwell scenarios used here
#scenario_descriptor comes from dwell time, used for getting years of dwell time
make_pretty_table_from_survival<-function(my_survival,
                                          scenario_descriptor,
                                          population_scale=100000,
                                          fp_rate=0.007){
  
  #late stage is III+IV
  #number_stage 3+4
  #basic performance measures of mced under different circumstances
  
  fixed_descriptor<-scenario_descriptor %>% 
    mutate(scenario_type=sprintf("%s(%s)",aggressive,dwell_range)) %>%
    select(aggressive,scenario_type)
  
  #three parts of the table
  basic_performance<-my_survival %>%
    group_by(scan,aggressive) %>%
    summarize(total_cancers=sum(caught),
              found_usual_care=sum(caught[mode_found=="soc"]),
              found_intercepted=sum(caught[mode_found=="cfdna"])) %>%
    ungroup() %>%
    mutate(population=population_scale,
           fp_rate=fp_rate,
           total_fp = case_when(scan!="no" ~ (population-total_cancers)*fp_rate,
                                TRUE ~ 0.0), # no fp if no scans
           ppv=case_when(scan!="no" ~ found_intercepted/(total_fp+found_intercepted),
                         TRUE ~ NA_real_),
           scenario_zero=case_when(scan=="no" ~ "No MCED",
                                   TRUE ~ aggressive))
  
  text_basic_performance<-basic_performance %>%
    mutate(text_found_usual_care=pretty_string_plus_percent(found_usual_care,100*found_usual_care/total_cancers),
           text_found_intercepted=pretty_string_plus_percent(found_intercepted,100*found_intercepted/total_cancers),
           text_tp_ppv=pretty_string_plus_percent(found_intercepted,100*ppv)) %>%
    left_join(fixed_descriptor) %>%
    mutate(scenario_type=case_when(scan=="no" ~ "No MCED",
                                   TRUE ~ scenario_type)) %>%
    select(scenario_zero,
           scenario_type,
           scan,
           total_cancers,
           text_found_usual_care,
           text_found_intercepted,
           total_fp,
           text_tp_ppv) %>%
    mutate_if(is.numeric,function(z){round(z,0)})
  
  #unstaged data, not included in final table
  unstaged_data<-my_survival %>%
    group_by(scan,aggressive) %>%
    summarize(total_unstaged=sum(caught[clinical==0]),
              ied_unstaged=sum(original_deaths[clinical==0]))
  
  
  #stage shift
  stage_shift_performance<-my_survival %>%
    mutate(late_original=clinical %in% c(3,4),
           late_shifted=prequel %in% c(3,4)) %>%
    group_by(scan,aggressive) %>%
    summarize(total_incidence=sum(caught),
              total_late_cancers_original=sum(caught[late_original]),
              final_late_cancers=sum(caught[late_shifted]),
              #breakout
              total_found_usual_care=sum(caught[mode_found=="soc"]),
              total_found_intercepted=sum(caught[mode_found=="cfdna"]),
              total_late_usual_care=sum(caught[late_original & mode_found=="soc"]),
              total_late_pre_intercept=sum(caught[late_original & mode_found=="cfdna"]),
              total_late_post_intercept=sum(caught[late_shifted & mode_found=="cfdna"])) %>%
    ungroup() %>%
    mutate(scenario_zero=case_when(scan=="no" ~ "No MCED",
                                   TRUE ~ aggressive))
  
  text_stage_shift_performance<-stage_shift_performance %>%
    mutate(text_total_late_cancers=pretty_string_plus_percent(total_late_cancers_original,
                                                              100*total_late_cancers_original/total_incidence),
           text_reduction_late_cancers=pretty_string_plus_percent(final_late_cancers-total_late_cancers_original,
                                                                  100*(1-final_late_cancers/total_late_cancers_original)),
           text_final_late_cancers=pretty_string_plus_percent(final_late_cancers,
                                                              100*final_late_cancers/total_incidence),
           #breakout
           text_late_usual_care=pretty_string_plus_percent(total_late_usual_care,
                                                           100*total_late_usual_care/total_found_usual_care),
           text_late_pre_intercept=pretty_string_plus_percent(total_late_pre_intercept,
                                                              100*total_late_pre_intercept/total_found_intercepted),
           text_late_post_intercept=pretty_string_plus_percent(total_late_post_intercept,
                                                               100*total_late_post_intercept/total_found_intercepted),
           text_reduction_late_intercepted=pretty_string_plus_percent(total_late_post_intercept-total_late_pre_intercept,
                                                                      100*(1-total_late_post_intercept/total_late_pre_intercept))) %>%
    select(scenario_zero,
           scan,
           text_total_late_cancers,
           text_final_late_cancers,
           text_reduction_late_cancers,
           text_late_usual_care,
           text_late_pre_intercept,
           text_late_post_intercept,
           text_reduction_late_intercepted)
  
  #deaths
  ied_performance<-my_survival %>%
    group_by(scan,aggressive) %>%
    summarize(total_cancers=sum(caught),
              total_ied=sum(original_deaths),
              final_ied=sum(shifted_deaths),
              #breakout
              ied_found_usual_care=sum(original_deaths[mode_found=="soc"]),
              ied_found_pre_intercept=sum(original_deaths[mode_found=="cfdna"]),
              ied_found_post_intercept=sum(shifted_deaths[mode_found=="cfdna"])) %>%
    ungroup() %>%
    mutate(scenario_zero=case_when(scan=="no" ~ "No MCED",
                                   TRUE ~ aggressive))
  
  text_ied_performance<-ied_performance %>%
    mutate(text_total_ied=pretty_string_plus_percent(total_ied,100*total_ied/total_cancers),
           text_final_ied=pretty_string_plus_percent(final_ied,100*final_ied/total_cancers),
           text_reduction_ied=pretty_string_plus_percent(final_ied-total_ied,100*(1-final_ied/total_ied)),
           #breakout
           text_ied_usual_care=pretty_string_plus_percent(ied_found_usual_care,100*ied_found_usual_care/total_ied),
           text_ied_pre_intercept=pretty_string_plus_percent(ied_found_pre_intercept,100*ied_found_pre_intercept/total_ied),
           text_ied_post_intercept=pretty_string_plus_percent(ied_found_post_intercept,100*ied_found_post_intercept/ied_found_pre_intercept),
           text_reduction_ied_intercept=pretty_string_plus_percent(ied_found_post_intercept-ied_found_pre_intercept,
                                                                   100*(1-ied_found_post_intercept/ied_found_pre_intercept))) %>%
    select(scenario_zero,
           scan,
           text_total_ied,
           text_final_ied,
           text_reduction_ied,
           text_ied_usual_care,
           text_ied_pre_intercept,
           text_ied_post_intercept,
           text_reduction_ied_intercept)
  
  #generate easy-for-computer table of summary data matching the pretty table
  all_performance<-basic_performance %>%
    left_join(stage_shift_performance) %>%
    left_join(ied_performance %>% select(-total_cancers)) %>%
    left_join(unstaged_data)
  
  
  #join together for text table
  text_all_together_performance<-text_basic_performance %>% 
    left_join(text_stage_shift_performance) %>%
    left_join(text_ied_performance) 
  
  happy_names<-c("scenario_zero", "header",
                 "scenario_type","Scenario(years in stage I)",
                 "scan","Screening Type",
                 "total_cancers","Total Cancer Incidence (100K person years),n",
                 "text_found_usual_care","Found Usual Care, n(% total cancer)",
                 "text_found_intercepted","Intercepted, n(% total cancer)",
                 "total_fp","False Positives, n",
                 "text_tp_ppv","True Positives, n (PPV = % of positives)",
                 "text_total_late_cancers","Total Late (III+IV),n (% total cancers)",
                 "text_final_late_cancers","Final Late (III+IV),n (% total cancers)",
                 "text_reduction_late_cancers","Reduction Late, n( % total late)",
                 "text_late_usual_care","Stage III+IV Usual Care, n (% found usual care)",
                 "text_late_pre_intercept","Stage III+IV Pre-Intercept, n (% intercepted)",
                 "text_late_post_intercept","Stage III+IV Post-Intercept,n (% intercepted)",
                 "text_reduction_late_intercepted","Reduction(III+IV), n (% late intercepted)",
                 "text_total_ied","Total IED, n(% total cancer)",
                 "text_final_ied","Final IED, n(% total cancer)",
                 "text_reduction_ied","Reduction IED, n, (% total IED)",
                 "text_ied_usual_care","IED Usual Care, n (% total IED)",
                 "text_ied_pre_intercept","IED Pre-Intercept, n (% total IED)",
                 "text_ied_post_intercept","IED Post-Intercept, n (% IED pre-intercept)",
                 "text_reduction_ied_intercept","Reduction IED, n (% IED pre-intercept)")
  
  pretty_name_table<-
    tibble(outcome=happy_names[seq(1,length(happy_names),by=2)],
           long=happy_names[seq(2,length(happy_names),by=2)])
  
  
  display_text_table_performance_incident<-
    text_all_together_performance %>% 
    filter(scan=="incident") %>%
    arrange(scenario_zero) %>%
    t() %>% (function(z){
      colnames(z)<-z[1,]
      z[-1,]}) %>%
    as_tibble(rownames="outcome")  
  
  display_text_table_performance_prevalent<-
    text_all_together_performance %>% 
    filter(scan=="prevalent") %>%
    arrange(scenario_zero) %>%
    t() %>% (function(z){
      colnames(z)<-z[1,]
      z[-1,]}) %>%
    as_tibble(rownames="outcome")  
  
  display_text_table_performance_no<-
    text_all_together_performance %>% 
    filter(scan=="no") %>%
    arrange(scenario_zero) %>%
    t() %>% (function(z){
      colnames(z)<-z[1,]
      z[-1,,drop=FALSE]}) %>%
    as_tibble(rownames="outcome")  
  
  #select columns for paper, add spacer labels for table
  display_text_table_performance_incident_final<-display_text_table_performance_incident %>%
    bind_cols(display_text_table_performance_no %>% select('No MCED')) %>%
    left_join(pretty_name_table) %>%
    select(Measure=long,
           No='No MCED',
           AggFast,
           Fast,
           Slow,
           MIS) %>%
    mutate(row=1:length(Measure)) %>%
    bind_rows(tibble("Measure"="Basic Performance","No"="","AggFast"="","Fast"="","Slow"="","MIS"="","row"=2.5)) %>%
    bind_rows(tibble("Measure"="Stage Shift","No"="","AggFast"="","Fast"="","Slow"="","MIS"="","row"=7.5)) %>%
    bind_rows(tibble("Measure"="Individuals Expected to Die (IED) of Cancer in 5 years","No"="","AggFast"="","Fast"="","Slow"="","MIS"="","row"=14.5)) %>%
    bind_rows(tibble("Measure"="Breakout By Mode Found","No"="","AggFast"="","Fast"="","Slow"="","MIS"="","row"=10.5)) %>%
    bind_rows(tibble("Measure"="Breakout By Mode Found","No"="","AggFast"="","Fast"="","Slow"="","MIS"="","row"=17.5)) %>% 
    arrange(row) %>%
    select(-row)
  
  display_text_table_performance_prevalent_final<-display_text_table_performance_prevalent %>%
    bind_cols(display_text_table_performance_no %>% select('No MCED')) %>%
    left_join(pretty_name_table) %>%
    select(Measure=long,
           No='No MCED',
           AggFast,
           Fast,
           Slow) %>%
    mutate(row=1:length(Measure)) %>%
    bind_rows(tibble("Measure"="Basic Performance","No"="","AggFast"="","Fast"="","Slow"="","row"=2.5)) %>%
    bind_rows(tibble("Measure"="Stage Shift","No"="","AggFast"="","Fast"="","Slow"="","row"=7.5)) %>%
    bind_rows(tibble("Measure"="Individuals Expected to Die (IED) of Cancer in 5 years","No"="","AggFast"="","Fast"="","Slow"="","row"=14.5)) %>%
    bind_rows(tibble("Measure"="Breakout By Mode Found","No"="","AggFast"="","Fast"="","Slow"="","row"=10.5)) %>%
    bind_rows(tibble("Measure"="Breakout By Mode Found","No"="","AggFast"="","Fast"="","Slow"="","row"=17.5)) %>% 
    arrange(row) %>%
    select(-row)
  
  
  
  list(prevalent=display_text_table_performance_prevalent_final,
       incident=display_text_table_performance_incident_final,
       all=all_performance)
}

#slip rate computation

#exact: integrate over those missed during a screening interval
#see supplemental information for argument that this function is appropriate
integrate_slip_rate<-function(screen_interval, weibull_shape, dwell){
  # slip rate
  # yield of escape = integral cumululate distribution function F(t), 0<t<screen_interval
  # mean of weibull = scale*gamma(1+1/shape)
  dwell_scale=dwell/gamma(1+1/weibull_shape)
  tiny_delta<-365
  #trapezoidal integration
  days_low<-seq(0,screen_interval*tiny_delta-1,by=1)/tiny_delta
  days_hi<-days_low+1/tiny_delta
  F_by_day<-0.5*(pweibull(days_low,shape=weibull_shape,scale=dwell_scale)+pweibull(days_hi,shape=weibull_shape,scale=dwell_scale))
  escaped_yield<-sum(F_by_day)*(1/tiny_delta) #day width in years
  #convert to slip rate: how many missed
  total_yield<-screen_interval #total incidence is just duration
  slip_rate<-escaped_yield/total_yield
  slip_rate
}

#use integrated weibull cumulative distribution functions
#"exact" solution to slip rate
exact_slip_rate_from_dwell<-function(dwell_model_all_df,screen_interval=1,weibull_shape=1){
  
   #slip is "before clinical"
  #slip_clinical = "at stage of clinical detection"
  #assume expected is half-duration of stage of clinical detection
  #completeness in modeling
  dwell_slip_df<-dwell_model_all_df %>%
    mutate(slip = sapply(dwell,function(z){integrate_slip_rate(screen_interval,weibull_shape,z)}),
           slip_clinical=sapply(dwell*0.5,function(z){integrate_slip_rate(screen_interval,weibull_shape,z)}),
           screen_interval=screen_interval)
  
  dwell_ideal_df <- dwell_slip_df %>%
    filter(scenario==1) %>%
    mutate(dwell=10000,slip=0,slip_clinical=0,scenario=as.integer(0),screen_interval=NA)
  
  #add "scenario 0" perfect interception
  dwell_slip_df<-bind_rows(dwell_slip_df,dwell_ideal_df)
  dwell_slip_df
}

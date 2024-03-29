#------------------------------------------------------------------------------
# Purpose: Simulation study for manuscript,
# "Blinded sample size re-estimation for registry-based randomized controlled
# trials with incomplete event adjudication."
# Author: Ann Marie Weideman
#
# Dataset: see script simulate_data.R for further details
#------------------------------------------------------------------------------

#set working directory to current location
script_path<-setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Sample function will produce undesirable output if the vector of IDs
# is of length 1. The sample function is modified to be safe as below. See
# first paragraph in "Details" in ?sample for further info.
safer_sample<-function(x, size, replace = FALSE, prob = NULL)
{
  if (missing(size))
    size <- length(x)
  x[sample.int(length(x), size, replace, prob)]
}

#------------------------------------------------------------------------------
# Function for regression imputation for incomplete adjudications
#------------------------------------------------------------------------------

impute_fun<-function(data, t_interim){

  #-----------------------------------------------------
  # Subset data
  #-----------------------------------------------------

  # A full dataset that excludes treatment arm
  W<-data

  # Note: in the below, Y0=-1 is adjudicated if follow up is complete and Y0=-1
  # is unadjudicated if patient is still under observation

  # A subset of W that includes all adjudicated events and non-events with
  # complete follow-up
  # X = x, Y0 in {-1, 0, 1}, Y1 in {-1, 1, 1}, Z=1
  W_adj<-W[W$adjudicated==1,]

  # A subset of W that includes incomplete adjudications with events
  # X = x, Y0 in {0, 1}, Y1 = NA, Z = 0
  W_unadj_event<-W[W$adjudicated==0 & W$y0!=-1,]

  # A subset of W that includes complete adjudications with events
  # X = x, Y0 in {0, 1}, Y1 in {0, 1}, Z = 1
  W_adj_event<-W[W$adjudicated==1 & W$event==1,]

  # A subset of W that includes incomplete adjudications with non-events
  # X = x, Y0 = -1, Y1 = NA, Z = 0
  W_unadj_noevent<-W[W$adjudicated==0 & W$event==0,]

  #-----------------------------------------------------
  # Make predictions
  #-----------------------------------------------------

  #############################################
  # Method 1: logit(y1) ~ x1 + x2 + t_event
  #############################################

  # Step 1: Fit logistic regression to adjudicated patients with events
  binom1<-glm(y1_obs~x1+x2+t_event, family=binomial, data=W_adj_event)
  # Predict y1=1 in unadjudicated patients with events
  prob.binom1<-predict(binom1, newdata=W_unadj_event, type="response")
  pred.binom1<-rbinom(nrow(W_unadj_event),1,prob.binom1) #new

  # Step 2: Fit Cox proportional hazards model to combined dataset that
  # includes adjudicated patients with events AND predictions from step 1.
  W_temp1<-W_adj #Adjudicated outcomes
  W_temp1$y1_obs[W_temp1$y1_obs==-1]<-0 # Treat any non-events as 0 for this prediction
  W_temp2<-W_unadj_event # Unadjudicated outcomes with events
  W_temp2$y1_obs<-pred.binom1 # Replace y1_obs with predicted values from step 1
  W_temp3<-rbind(W_temp1,W_temp2) #Bind adjudicated with predictions for unadjudicated

  library(survival)
  # Cox regression
  cox1 <- coxph(Surv(t_eligible, y1_obs) ~ x1 + x2, data = W_temp3)

  # Predict relative risks for W_unadj_noevent
  rr1 <- predict(cox1, newdata = W_unadj_noevent, type = "risk")

  # Estimate the baseline hazard function using Breslow's estimator
  baseline_h1 <- basehaz(cox1, centered = FALSE)

  # Function to interpolate baseline hazard for each time point
  interpolate_hazard <- function(time_points, baseline_h) {
    sapply(time_points, function(time) {
      if(any(baseline_h$time >= time)) {
        idx <- which.max(baseline_h$time[baseline_h$time <= time])
        return(baseline_h$hazard[idx])
      } else {
        return(NA)
      }
    })
  }

  # Interpolate the baseline hazard for each time point in W_unadj_noevent
  interpolated_baseline_h1 <- interpolate_hazard(W_unadj_noevent$t_eligible, baseline_h1)

  # Calculate absolute hazard rates
  h1 <- interpolated_baseline_h1 * rr1

  # Calculate conditional probability of an event between time t_eligible and t_interim:
  # P(t_eligible < T < t_interim | T > t_eligible)
  # = P(T > t_eligible | T > t_eligible) - P( T > t_interim | T > t_eligible)
  # = 1 - S(t_interim | T > t_eligible)
  # = 1 - exp(-h1*t_interim)/exp(-h1*t_eligible)
  # = 1 - exp(-h1*(t_interim-t_eligible))
  prob.cox1<-1-exp(-h1*(t_interim-W_unadj_noevent$t_eligible))
  pred.cox1<-rbinom(nrow(W_unadj_noevent),1,prob.cox1)

  ############################################################
  # Method 2: logit(y1) ~ x1 + x2 + t_event, stratified by y0
  ############################################################

  # Step 1: Fit logistic regression to adjudicated patients with events
  # model 1: subset to Y0 = 0
  # model 2: subset to Y0 = 1
  binom2.y0eq0<-glm(y1_obs~x1+x2+t_event, family=binomial, data=W_adj_event[W_adj_event$y0==0,])
  binom2.y0eq1<-glm(y1_obs~x1+x2+t_event, family=binomial, data=W_adj_event[W_adj_event$y0==1,])

  # Predict y1=1 given y0=0 in unadjudicated patients with events
  if(nrow(W_unadj_event[W_unadj_event$y0==0,])!=0){
    prob.binom2.y0eq0<-predict(binom2.y0eq0,
                             newdata=W_unadj_event[W_unadj_event$y0==0,],
                             type="response")
    pred.binom2.y0eq0<-rbinom(nrow(W_unadj_event[W_unadj_event$y0==0,]),
                        1,prob.binom2.y0eq0)
  }else{pred.binom2.y0eq0<-c()}
  # Predict y1=1 given y0=1 in unadjudicated patients with events
  if(nrow(W_unadj_event[W_unadj_event$y0==1,])!=0){
    prob.binom2.y0eq1<-predict(binom2.y0eq1,
                             newdata=W_unadj_event[W_unadj_event$y0==1,],
                             type="response")
    pred.binom2.y0eq1<-rbinom(nrow(W_unadj_event[W_unadj_event$y0==1,]),
                        1,prob.binom2.y0eq1)
  }else{pred.binom2.y0eq1<-c()}

  # Step 2: Fit Cox proportional hazards model to combined dataset that
  # includes adjudicated patients with events AND predictions from step 1.
  W_temp4<-W_adj #Adjudicated outcomes
  W_temp4$y1_obs[W_temp4$y1_obs==-1]<-0 # Treat any non-events as 0 for this prediction
  W_temp5<-W_unadj_event # Unadjudicated outcomes with events
  # Replace y1_obs with predicted values from step 1
  W_temp5[W_temp5$y0==0,]$y1_obs<-pred.binom2.y0eq0
  W_temp5[W_temp5$y0==1,]$y1_obs<-pred.binom2.y0eq1
  W_temp6<-rbind(W_temp4,W_temp5) #Bind adjudicated with predictions from unadjudicated

  library(survival)
  # Cox regression
  cox2 <- coxph(Surv(t_eligible, y1_obs) ~ x1 + x2, data = W_temp6)

  # Estimate the baseline hazard function using Breslow's estimator
  baseline_h2 <- basehaz(cox2, centered = FALSE)

  # Predict relative risks for W_unadj_noevent
  rr2 <- predict(cox2, newdata = W_unadj_noevent, type = "risk")

  # Interpolate the baseline hazard for each time point in W_unadj_noevent
  interpolated_baseline_h2 <- interpolate_hazard(W_unadj_noevent$t_eligible, baseline_h2)

  # Calculate absolute hazard rates
  h2 <- interpolated_baseline_h2 * rr2

  # Calculate conditional probability of an event between time t_eligible and t_interim:
  # P(t_eligible < T < t_interim | T > t_eligible)
  # = P(T > t_eligible | T > t_eligible) - P( T > t_interim | T > t_eligible)
  # = 1 - S(t_interim | T > t_eligible)
  # = 1 - exp(-h2*t_interim)/exp(-h2*t_eligible)
  # = 1 - exp(-h2*(t_interim-t_eligible))
  prob.cox2<-1-exp(-h2*(t_interim-W_unadj_noevent$t_eligible))
  pred.cox2<-rbinom(nrow(W_unadj_noevent),1,prob.cox2)

  #-----------------------------------------------------------------------------
  # Estimate probability of an event in all subjects irrespective of
  # assigned treatment arm for use in blinded sample size re-estimation
  # procedure
  #-----------------------------------------------------------------------------

  ##############################################################################
  # Method 0: logit(y1) ~ x1 + x2 + t_event
  # Approach: analyze all data to include adjudicated outcomes and unadjudicated
  # outcomes from predictions
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  #
  # This does NOT include:
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  # By not including these participants, we are effectively assuming
  # censoring at random.
  ##############################################################################
  # Include all observed y1
  yfit0<-W$y1_obs
  # Replace unadjudicated y1 that are events with predictions from binomial model
  yfit0[W$adjudicated==0 & W$event==1]<-pred.binom1
  # Remove unadjudicated y1 that are non-events since no predictions made
  yfit0<-yfit0[-which(W$adjudicated==0 & W$event==0)]
  # Average number of classifying events (y_obs=1)
  # (ymethod2>-1)*ymethod2 treats -1's as 0's when computing mean
  p_method0<-mean((yfit0>-1)*yfit0)
  p_method0

  ##############################################################################
  # Method 1: logit(y1) ~ x1 + x2 + t_event
  # Approach: analyze all data to include adjudicated outcomes and unadjudicated
  # outcomes from predictions
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  ##############################################################################
  # Include all observed y1
  ymethod1<-W$y1_obs
  # Replace unadjudicated y1 that are events with predictions from binomial model
  ymethod1[W$adjudicated==0 & W$event==1]<-pred.binom1
  # Replace unadjudicated y1 that are non-events with predictions from Cox model
  ymethod1[W$adjudicated==0 & W$event==0]<-pred.cox1
  # Average number of classifying events (y_obs=1)
  # (ymethod1>-1)*ymethod1 treats -1's as 0's when computing mean
  p_method1<-mean((ymethod1>-1)*ymethod1)
  p_method1

  # Proportion of incorrect predictions for unadjudicated events
  p_event_unadj_method1<-ymethod1[W$event==1 & W$adjudicated==0]
  event_unadj_method1<-rbinom(length(p_event_unadj_method1),1,p_event_unadj_method1)
  pred_error_method1<-event_unadj_method1-W$y1_true[W$event==1 & W$adjudicated==0]

  # Total error rate (TER)
  p_TER_method1<-mean(abs(pred_error_method1))
  # False negative rate (FNR)
  p_FNR_method1<-mean((pred_error_method1==-1)*1)
  # False positive rate (FPR)
  p_FPR_method1<-mean((pred_error_method1==1)*1)

  ##############################################################################
  # Method 2: logit(y1) ~ x1 + x2 + t_event, stratified by y0
  # Approach: analyze all data to include adjudicated outcomes and unadjudicated
  # outcomes from predictions
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  ##############################################################################
  # Include all observed y1
  ymethod2<-W$y1_obs
  # Replace unadjudicated y1 that are events with predictions from binomial model
  ymethod2[W$y0==0 & W$adjudicated==0]<-pred.binom2.y0eq0 #for y0=0
  ymethod2[W$y0==1 & W$adjudicated==0]<-pred.binom2.y0eq1 #for y0=1
  # Replace unadjudicated y1 that are non-events with predictions from Cox model
  ymethod2[W$adjudicated==0 & W$event==0]<-pred.cox2
  # Average number of classifying events (y_obs=1)
  # (ymethod2>-1)*ymethod2 treats -1's as 0's when computing mean
  p_method2<-mean((ymethod2>-1)*ymethod2)
  p_method2

  # Proportion of incorrect predictions for unadjudicated events
  p_event_unadj_method2<-ymethod2[W$event==1 & W$adjudicated==0]
  event_unadj_method2<-rbinom(length(p_event_unadj_method2),1,p_event_unadj_method2)
  pred_error_method2<-event_unadj_method2-W$y1_true[W$event==1 & W$adjudicated==0]

  # Total error rate (TER)
  p_TER_method2<-mean(abs(pred_error_method2))
  # False negative rate (FNR)
  p_FNR_method2<-mean((pred_error_method2==-1)*1)
  # False positive rate (FPR)
  p_FPR_method2<-mean((pred_error_method2==1)*1)

  ###########################################################################
  # Complete data analysis
  # Approach: analyze observed adjudications, y1_obs only (no predictions)
  #
  # Note: Complete follow-up ==> t_eligible = t_interim, for these simulations
  # assumed to be 1 yr of complete follow-up and interim analysis at 1 yr.
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  #
  # This does NOT include:
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  # By not including these participants, we are effectively assuming
  # censoring at random.
  ###########################################################################
  p_complete<-mean((W_adj$y1_obs==1)*1)
  p_complete

  ###########################################################################
  # Carry forward analysis
  # Approach: analyze all y1 but use y0 where y1 is missing (y1_obs=NA).
  #
  # Note: Complete follow-up ==> t_eligible = t_interim, for these simulations
  # assumed to be 1 yr of complete follow-up and interim analysis at 1 yr.
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  # E.g., y1_obs = 1 (CV death), y1_obs = 0 (non-CV death), y1_obs = -1 (alive).
  #
  # This does NOT include:
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  # By not including these participants, we are effectively assuming
  ###########################################################################
  # Include those with complete follow-up only
  W_completefu<-W[-which(W$y0==-1 & W$t_eligible<1),]
  y_carryforward<-W_completefu$y1_obs
  # Replace y1_obs=NA with y0
  y_carryforward[W_completefu$adjudicated==0]<-W_completefu$y0[W_completefu$adjudicated==0]
  # Average number of classifying events (y_obs=1)
  p_carryforward<-mean((y_carryforward==1)*1)
  p_carryforward

  ############################################################################
  # Unadjudicated outcomes treated as not a classifying event
  # Approach: analyze all y1 but replace y1 with 0 (not a classifying event)
  # when missing (y1_obs=NA). This is an extreme approach that will provide a
  # lower bound on the mean of y1_obs.
  #
  # Note: Complete follow-up ==> t_eligible = t_interim, for these simulations
  # assumed to be 1 yr of complete follow-up and interim analysis at 1 yr.
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  # E.g., y1_obs = 1 (CV death), y1_obs = 0 (non-CV death), y1_obs = -1 (alive).
  #
  # This does NOT include:
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  # By not including these participants, we are effectively assuming
  # censoring at random.
  ###########################################################################
  # Include those with complete follow-up only
  y_not_classifying_event<-W_completefu$y1_obs
  # Replace y_obs=NA with 0
  y_not_classifying_event[W_completefu$adjudicated==0]<-0
  # Average number of classifying events (y_obs=1)
  p_not_classifying_event<-mean((y_not_classifying_event==1)*1)
  p_not_classifying_event

  ############################################################################
  # Unadjudicated outcomes treated as classifying events
  # Approach: analyze all y1 but replace y1 with 1 (classifying event)
  # when missing (y1_obs=NA). This is an extreme approach that will provide an
  # upper bound on the mean of y1_obs.
  #
  # Note: Complete follow-up ==> t_eligible = t_interim, for these simulations
  # assumed to be 1 yr of complete follow-up and interim analysis at 1 yr.
  #
  # This includes:
  # Participants with classifying events (y1_obs=1) and complete follow-up
  # Participants with non-classifying events (y1_obs=0) and complete follow-up
  # Participants with non-events (y1_obs=-1) and complete follow-up
  # E.g., y1_obs = 1 (CV death), y1_obs = 0 (non-CV death), y1_obs = -1 (alive)
  #
  # This does NOT include:
  # Participants with non-events (y1_obs=-1) that lack complete follow-up
  # By not including these participants, we are effectively assuming
  # censoring at random.
  ###########################################################################
  # Include those with complete follow-up only
  y_classifying_event<-W_completefu$y1_obs
  # Replace y_obs=NA with 1
  y_classifying_event[W_completefu$adjudicated==0]<-1
  # Average number of classifying events (y_obs=1)
  p_classifying_event<-mean((y_classifying_event==1)*1)
  p_classifying_event

  return(list(p_method0=p_method0,
              p_method1=p_method1,
              p_method2=p_method2,
              p_complete=p_complete,
              p_carryforward=p_carryforward,
              p_not_classifying_event= p_not_classifying_event,
              p_classifying_event=p_classifying_event,
              p_TER_method1=p_TER_method1,
              p_FNR_method1=p_FNR_method1,
              p_FPR_method1=p_FPR_method1,
              p_TER_method2=p_TER_method2,
              p_FNR_method2=p_FNR_method2,
              p_FPR_method2=p_FPR_method2))
}

#-------------------------------------------------------------------------------
# Functions for stratified bootstrap of predicted probabilities
# Note: the stratification across adjudications and event status ensures balance
# in the resampled datasets
#-------------------------------------------------------------------------------

# Serial version
strat_boot_serial<-function(data, B, seedling){

  set.seed(seedling)

  # Create necessary indices for stratification
  id_adj_event<-which(data$adjudicated==1 & data$event==1)
  id_adj_noevent<-which(data$adjudicated==1 & data$event==0)
  id_unadj_event<-which(data$adjudicated==0 & data$event==1)
  id_unadj_noevent<-which(data$adjudicated==0 & data$event==0)

  boot_data<-list()
  boot_p<-matrix(NA,nrow=B,ncol=13)

  for (i in seq_len(B)){
    samp1<-data[safer_sample(id_adj_event, length(id_adj_event),replace=TRUE), ]
    samp2<-data[safer_sample(id_adj_noevent, length(id_adj_noevent),replace=TRUE), ]
    samp3<-data[safer_sample(id_unadj_event, length(id_unadj_event),replace=TRUE), ]
    samp4<-data[safer_sample(id_unadj_noevent, length(id_unadj_noevent),replace=TRUE), ]
    boot_data <- rbind(samp1, samp2, samp3, samp4)
    boot_p[i,]<-as.numeric(impute_fun(boot_data[[i]],t_interim))
  }

  colnames(boot_p)<-c("p_method0", "p_method1", "p_method2", "p_complete",
                      "p_carryforward", "p_y1eq0",  "p_y1eq1",
                      "p_TER_method1", "p_FNR_method1","p_FPR_method1",
                      "p_TER_method2", "p_FNR_method2","p_FPR_method2")
  return(boot_p)
}

# Parallel version
strat_boot_parallel <- function(data, B, seedling, cl) {

  set.seed(seedling)

  # Create necessary indices for stratification
  id_adj_event<-which(data$adjudicated==1 & data$event==1)
  id_adj_noevent<-which(data$adjudicated==1 & data$event==0)
  id_unadj_event<-which(data$adjudicated==0 & data$event==1)
  id_unadj_noevent<-which(data$adjudicated==0 & data$event==0)

  clusterExport(cl, varlist = c("data", "id_adj_event", "id_adj_noevent",
                                "id_unadj_event", "id_unadj_noevent",
                                "impute_fun", "t_interim","safer_sample"),
                 envir = environment())

  # Define the function to be applied in parallel
  boot_func <- function(i) {
    samp1<-data[safer_sample(id_adj_event, length(id_adj_event),replace=TRUE), ]
    samp2<-data[safer_sample(id_adj_noevent, length(id_adj_noevent),replace=TRUE), ]
    samp3<-data[safer_sample(id_unadj_event, length(id_unadj_event),replace=TRUE), ]
    samp4<-data[safer_sample(id_unadj_noevent, length(id_unadj_noevent),replace=TRUE), ]
    boot_data <- rbind(samp1, samp2, samp3, samp4)
    as.numeric(impute_fun(boot_data, t_interim))
  }

  # Run the loop in parallel
  boot_p <- parLapply(cl, seq_len(B), boot_func)
  boot_p <- do.call(rbind, boot_p)  # Convert the list to a matrix

  colnames(boot_p) <- c("p_method0", "p_method1", "p_method2", "p_complete",
                        "p_carryforward", "p_y1eq0",  "p_y1eq1",
                        "p_TER_method1", "p_FNR_method1","p_FPR_method1",
                        "p_TER_method2", "p_FNR_method2","p_FPR_method2")
  return(boot_p)
}

#-------------------------------------------------------------------------
# Function for bias-corrected CI
#-------------------------------------------------------------------------

ci_boot <- function(test_data, boots, alpha) {

  boots_corrected<-2*unlist(impute_fun(data=test_data, t_interim)[1:7])-colMeans(boots)

  # Extract the lower and upper confidence interval bounds
  CI_LB <- boots_corrected+qnorm(alpha/2,lower.tail=T)*apply(boots, 2, sd)
  CI_UB <- boots_corrected+qnorm(alpha/2,lower.tail=F)*apply(boots, 2, sd)

  return(list(CI_LB = CI_LB, CI_UB = CI_UB))
}

#-------------------------------------------------------------------------
# Benchmarking
#-------------------------------------------------------------------------

source(paste0(script_path,"/simulate_data.R"))

n_vec<-c(250,500,1000) #sample size available for interim analysis
error_y0_vec<-c(0.05, 0.1, 0.3) # error in y0: 5%, 10%, 30%
adj_rate<-0.75 #adjudication rate
B<-1000 # number of bootstraps
D<-100 # number of datasets
R<-10 # number of replicates
t_interim<-1 #time (in years) from first randomization to interim analysis

Pred_error_mat<-data.frame("pred_error"=0,"Method"=0,"n"=0,"error_y0"=0)
CP_mat<-data.frame("CP"=0,"Method"=0,"n"=0,"error_y0"=0)
Bias_mat<-data.frame("bias"=0,"Method"=0,"n"=0,"error_y0"=0)
RMSE_mat<-data.frame("RMSE"=0,"Method"=0,"n"=0,"error_y0"=0)

# Indices for incrementing rows in dataframes
i<-j<-k<-0

names_benchmark<-c("Method 0", "Method 1", "Method 2", "Complete Case",
                   "Carry Forward", "Assign Y1=0", "Assign Y1=1", "Truth")
names_error<-c("Method 1 total", "Method 1 false negatives", "Method 1 false positives",
               "Method 2 total", "Method 2 false negatives", "Method 2 false positives")

# Loop over sample size at interim analysis
for(n in n_vec){
  # Loop over proportion of error in Y0
  for(error_y0 in error_y0_vec){

    # Reset random number generation
    seedling<-12345

    # Loop over replicates
    for(r in 1:R){

      # Initialize
      coverage<-RMSE_boot<-matrix(nrow=D,ncol=7)
      p_boot<-matrix(nrow=D,ncol=8)
      pred_error<-matrix(nrow=D,ncol=6)

      # Loop over bootstraps
      for(d in 1:D){

        print(d)

        # Simulate data
        test_data<-sim_fun(n=n,
                            error_y0=error_y0,
                            adj_rate=adj_rate,
                            t_interim=t_interim,
                            seedling=seedling)

        # Compute estimate of truth
        p_true<-mean((test_data$y1_true==1)*1)

        # Stratified bootstrap
        parallel<-T
        # Only initiate cluster once if running in parallel (saves time)
        if(seedling==12345 & parallel==T){
          library(parallel)

          # Number of cores
          no_cores <- detectCores() - 1

          # Create a cluster
          cl <- makeCluster(no_cores)
        }
        boots_out<-strat_boot_parallel(data=test_data, B=B, seedling=seedling,
                                   cl=cl)

        # Store prediction error for unadjudicated events for methods 1 and 2
        pred_error[d,]<-colMeans(boots_out[,8:13])

        # Store percentile bootstrap CIs
        CI_boot<-ci_boot(test_data=test_data,boots=boots_out[,1:7],alpha=0.05)

        # Store indicator for coverage probability
        coverage[d,]<-(p_true > CI_boot$CI_LB & p_true < CI_boot$CI_UB)*1

        # Store estimate of pooled event rate
        p_boot[d,]<-c(colMeans(boots_out[,1:7]), p_true=p_true)

        # Estimate root mean squared error (RMSE)
        RMSE_boot[d,]<-sqrt(colMeans((boots_out[,1:7]-p_true)^2))

        seedling<-seedling+1
      }

      # Percent bootstrap bias
      boot_bias<-abs(p_boot[,-8] - p_boot[,8])
      boot_bias_long<-reshape2::melt(t(boot_bias))$value
      Bias_mat[(i+1):(i+7*D),]<-cbind(boot_bias_long,
                                      rep(names_benchmark[-8],D),
                                      rep(paste0("n=",n),D*7),
                                      rep(paste0("error=",error_y0),D*7))

      # Coverage probability (proportion of bootstrap CIs that contain truth)
      CP_mat[(j+1):(j+7),]<-cbind(colMeans(coverage),
                                  names_benchmark[-8],
                                  paste0("n=",rep(n,7)),
                                  paste0("error=",rep(error_y0,7)))

      # Average root mean squared error (RMSE)
      RMSE_mat[(j+1):(j+7),]<-cbind(colMeans(RMSE_boot),
                                    names_benchmark[-8],
                                    paste0("n=",rep(n,7)),
                                    paste0("error=",rep(error_y0,7)))

      # Prediction error for unadjudicated events for Methods 1 and 2
      Pred_error_mat[(k+1):(k+6),]<-cbind(colMeans(pred_error),
                                          names_error,
                                          paste0("n=",rep(n,2)),
                                          paste0("error=",rep(error_y0,2)))

      i<-i+7*D
      j<-j+7
      k<-k+6
      print(paste0("n=",n,", error=",error_y0,", replicate=",r))
    }
  }
}

# Stop the cluster if running in parallel
if(parallel==T){stopCluster(cl)}

save(CP_mat, Bias_mat, RMSE_mat, Pred_error_mat,
      file = paste0(script_path,"/output/sims_blinded_75pctadj.R"))

#-------------------------------------------------------------------------
# Graphing
#-------------------------------------------------------------------------

library(simhelpers)
library(ggplot2)
library(dplyr)
library(knitr)
library(kableExtra)
library(scales)

options(scipen=0)

custom_theme <- theme_bw(base_size = 10) +  # Set the base font size to 10
  theme(
    axis.text = element_text(size = rel(1.2)),    # Adjust axis text size
    axis.title = element_text(size = rel(1.4)),   # Adjust axis title size
    axis.title.y = element_text(size = rel(1.2)), # Enlarge y-axis title size specifically
    plot.title = element_text(size = rel(1.6)),   # Adjust plot title size
    strip.text = element_text(size = rel(1.3)),   # Adjust facet strip text size
    legend.text = element_text(size = rel(1.2)),  # Adjust legend text size
    legend.title = element_text(size = rel(1.4))  # Adjust legend title size
    # Add other theme adjustments here if needed
  )

# Panel plot for bias
Bias_mat$error_y0 <- factor(Bias_mat$error_y0,
                            levels = c("error=0.05", "error=0.1","error=0.3"),
                            labels = c('\u03f5==0.05', '\u03f5==0.1', '\u03f5==0.3')) 

Bias_mat %>%
  ggplot(aes(x = Method, y = as.numeric(bias), fill = Method)) +
  geom_boxplot() +
  facet_grid(factor(error_y0) ~ factor(n, levels=c('n=250', 'n=500', 'n=1000'),
             labels=c("n[int]==250", "n[int]==500", "n[int]==1000")), 
             scales = "free", labeller = label_parsed) +
  labs(x = "", y = "Bias") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file = paste0(script_path,"/output/Bias_75pctadj.jpg"),
       plot = last_plot(), width = 7.5, height = 5.82, dpi = 600)

# Panel plot for coverage prob
CP_mat$error_y0 <- factor(CP_mat$error_y0,
                            levels = c("error=0.05", "error=0.1","error=0.3"),
                            labels = c('\u03f5==0.05', '\u03f5==0.1', '\u03f5==0.3')) 
CP_mat %>%
  ggplot(aes(x = Method, y = as.numeric(CP)*100, fill = Method)) +
  geom_boxplot() +
  facet_grid(error_y0 ~ factor(n, levels=c('n=250', 'n=500', 'n=1000'),
                               labels=c("n[int]==250", "n[int]==500", "n[int]==1000")), 
             scales = "free", labeller = label_parsed) +
  labs(x = "", y = "Coverage probability of 95% CI (%)") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_hline(yintercept=95,linetype=2)
ggsave(file = paste0(script_path,"/output/Coverage_Prob_75pctadj.jpg"),
       plot = last_plot(), width = 7.5, height = 5.82, dpi = 600)

# Panel plot for average RMSE
RMSE_mat$error_y0 <- factor(RMSE_mat$error_y0,
                          levels = c("error=0.05", "error=0.1","error=0.3"),
                          labels = c('\u03f5==0.05', '\u03f5==0.1', '\u03f5==0.3')) 
RMSE_mat %>%
  ggplot(aes(x = Method, y = as.numeric(RMSE), fill = Method)) +
  geom_boxplot() +
  facet_grid(error_y0 ~ factor(n, levels=c('n=250', 'n=500', 'n=1000'),
                               labels=c("n[int]==250", "n[int]==500", "n[int]==1000")), 
             scales = "free", labeller = label_parsed) +
  labs(x = "", y = "Average root mean squared error (RMSE)") +
  custom_theme +
  scale_y_continuous(trans='log10',
                     breaks=c(0.01,0.03,0.10,0.30),
                     limits=c(0.01, 0.35)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file = paste0(script_path,"/output/RMSE_75pctadj.jpg"),
       plot = last_plot(), width = 7.5, height = 5.82, dpi = 600)

# Modify Pred_error_mat for plotting
Pred_error_mat$Method=rep(c(rep("Method 1", 3), rep("Method 2", 3)), nrow(Pred_error_mat)/6)
Pred_error_mat$Error_Type = rep(c("Total","FP","FN"),nrow(Pred_error_mat)/3)
Pred_error_mat$error_y0 <- factor(Pred_error_mat$error_y0,
                            levels = c("error=0.05", "error=0.1","error=0.3"),
                            labels = c('\u03f5==0.05', '\u03f5==0.1', '\u03f5==0.3')) 

# Panel plot for prediction error in unadjudicated events for Methods 1 and 2
# Save the plot to a variable
p <- ggplot(data = Pred_error_mat, aes(x = Method, y = as.numeric(pred_error))) +
  geom_boxplot() +
  facet_grid(error_y0 ~ factor(n, levels = c('n=250', 'n=500', 'n=1000'),
                               labels = c("n[int]==250", "n[int]==500", "n[int]==1000")) * Error_Type,
             scales = "fixed", labeller = label_parsed) +
  labs(x = "", y = "Prediction error for unadjudicated events") +
  custom_theme +
  scale_y_continuous(trans = 'log10', breaks = c(0.03, 0.10, 0.30), limits = c(0.03, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Convert the plot to grob an identify the location of the strips
library(tidyverse)
library(gtable)
library(grid)
g <- ggplot_gtable(ggplot_build(p))
stript <- grep("strip", g$layout$name)

# Merge common headings
grid_cols <- sort(unique(g$layout[stript,]$l))
t_vals <- rep(sort(unique(g$layout[stript,]$t)), each = length(grid_cols)/3)
l_vals <- rep(grid_cols[seq_along(grid_cols) %% 3 == 1], length = length(t_vals))
r_vals <- rep(grid_cols[seq_along(grid_cols) %% 3 == 0], length = length(t_vals))
labs   <- c(expression(n[int]==250), 
            expression(n[int]==500), 
            expression(n[int]==1000))

for(i in seq_along(labs))
{
  filler <- rectGrob(y = 0.72, height = 0.57, gp = gpar(fill = "gray85", col = "black"))
  tg    <- textGrob(label = labs[i], y = 0.75, gp = gpar(cex = 1.2))
  g     <- gtable_add_grob(g, filler, t = t_vals[i], l = l_vals[i], r = r_vals[i],
                           name = paste0("filler", i))
  g     <- gtable_add_grob(g, tg, t = t_vals[i], l = l_vals[i], r = r_vals[i],
                           name = paste0("textlab", i))
}
grid.newpage()
grid.draw(g)
ggsave(file = paste0(script_path,"/output/prediction_error_75pctadj.jpg"),
       plot = g, width = 9, height = 6, dpi = 600)

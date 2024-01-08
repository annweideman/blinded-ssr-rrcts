#------------------------------------------------------------------------
# Purpose: Function to simulate test data to pass to simulation_study.R
# Author: Ann Marie Weideman
#
# Details for final output dataset, test_data:
#
# y0: initial classification
# (-1 if no event, 0 if non-classifying event, 1 if classifying event)
#
# y1_obs: observed adjudicated classification
# (-1 if no event, 0 if non-classifying event, 1 if classifying event)
#
# y1_true: true adjudicated classification
# (-1 if no event, 0 if non-classifying event, 1 if classifying event)
#
# error: the proportion of y0 that differ from y1_true
# [0,1]
#
# event: the event status
# 0 if no event, 1 if event
#
# t_event: time of occurrence of any event (classifying or non-classifying)
# [0, inf)
#
# t_risk: time at risk of occurrence of any event (classifying or non-classifying)
# [0, inf)
#
# t_eligible: equal to t_event if event occurs before the interim analysis
#             equal to t_risk if no event occurs before the interim analysis
# [0, inf)
#
# t_interim: time between first randomization and planned interim analysis
# [0, inf)
#
# adjudicated: indicator of adjudication status
# 0 if adjudicated, 1 if adjudicated
#
# x: a vector of continuous baseline or time-dependent covariates
# (-inf, inf)
#
# arm: treatment arm
# 1 if treatment, 0 if control
#------------------------------------------------------------------------------

# Function arguments:
# n = sample size at interim analysis, e.g. 500
# error_y0 = proportion of error in y0 (compared to adjudicated y1), e.g. 0.1
# adj_rate = proportion of events adjudicated, e.g., 0.5
# t_interim = time between first randomization and interim analysis
# seedling = seed for random number generation

sim_fun<-function(n, error_y0, adj_rate, t_interim, seedling){

  check_data<-0

  # See end of function for explanation of while loop
  while(check_data==0){

    set.seed(seedling)+check_data

    # Initialize dataframe
    df<-data.frame(y1=rep(NA,n))

    # Generate baseline or time-dependent conts predictor x (e.g., age)
    df$x<-rnorm(n,0,1)

    # Generate treatment arm
    df$arm <- 1
    df$arm[sample(n, n/2)] <- 0 # set 50% to 0

    h0 <- 2  # baseline hazard
    beta1 <- -1
    beta2 <- 1
    hazard_rates <- h0 * exp(beta1 * df$arm + beta2 * df$x)

    # Generate event times
    df$t_event<-rexp(n,rate=hazard_rates)

    # Generate y1
    df$y1<-rep(NA,n)

    # Generate adjudicated outcome y1 via logistic regression to create
    # an association between the outcome y1 and predictors arm, x, and t_event
    a<-c(-1,runif(1,-2, 2),runif(1,-2, 2))
    sim_logit1 <-  a[1]*df$arm  +
                   a[2]*df$x +
                   a[3]*df$t_event
    sim_prob1 <- 1/(1+exp(-sim_logit1))
    df$y1 <- rbinom(length(sim_prob1), size = 1, sim_prob1)

    # Check simulation
    #model1 <- glm(y1~-1+arm+x+t_event, family="binomial", data=df)
    #model1$coefficients

    # Generate follow-up times
    # 1. Generate 2/3 of participants with complete follow-up
    # 2. Generate the other 1/3 as unif(0.5,t_interim)
    # 3. Randomly shuffle
    n_complete<-ceiling(2/3*n)
    n_incomplete<-n-n_complete
    df$t_risk<-sample(c(rep(t_interim,n_complete),
                        runif(n_incomplete,0.5,t_interim)))
    # If t_risk > t_event, then t_eligible=t_event, else t_eligible=t_risk
    df$t_eligible<-ifelse(df$t_risk>df$t_event,df$t_event,df$t_risk)
    # Set all non-events with complete follow-up to -1 (alive)
    df$y1[which(df$t_eligible<df$t_event & df$t_eligible==t_interim)]<--1

    # Generate event indicator
    df$event<-ifelse(df$t_event<t_interim & df$t_event==df$t_eligible,1,0)
    id_event<-which(df$event==1)
    id_noevent<-which(df$event==0)

    # Generate initial classifications y0 by incorporating error
    df$y0<-df$y1
    num_ones2 <- ceiling(length(id_event) * error_y0)
    error <- sample(c(rep(1, num_ones2),
                      rep(0, (length(id_event) - num_ones2))))
    df$error<-rep(NA,n)
    df$error[id_noevent]<-NA
    df$error[id_event]<-error
    df$y0[id_event]<-abs(error-df$y0[id_event])
    # Set all non-events  to -1
    df$y0[id_noevent]<--1

    # Generate indicator for adjudication status
    # Note that all non-events with complete follow-up are treated as
    # adjudicated b/c we assume we know their status with 100% accuracy,
    # although they do not count towards the proportion adjudicated (adj_rate)
    df$adjudicated<-rep(0,nrow(df))
    df$adjudicated[df$y0!=-1] <- rbinom(length(which(df$y0!=-1)),1,adj_rate)
    df$adjudicated[df$y1==-1 & df$t_eligible==1]<-1
    df$adjudicated[df$y1==-1 & df$t_eligible<1]<-0

    df$y1_obs<-df$y1_true<-df$y1
    # Set all unadjudicated observed y1 to NA (y1_obs=NA)
    df[which(df$adjudicated==0),"y1_obs"]<-NA
    # For those alive participants with incomplete follow-up, if event doesn't occur
    # before t_interim, then categorize their true y1 as alive (y1_true=-1) at
    # t_interim
    df$y1_true[which(df$t_event>t_interim)]<--1

    # Select variables of interest
    test_data<-subset(df,select=c(y0, y1_obs, y1_true, error, event, t_event, t_risk,
                                  t_eligible, adjudicated, x, arm))

    # Check that unadjudicated and adjudicated subsets contain at least
    # one y_0 = -1, 0, and 1
    # If not, repeat function and update increment seed by 1
    if(length(unique(test_data$y0[test_data$adjudicated==0]))==3 &
       length(unique(test_data$y0[test_data$adjudicated==1]))==3){
      check_data<-1
    }
  }

  return(test_data)
}

#Sivaram Cheruvu
#University of Texas at Dallas
#sivaram.cheruvu@utdallas.edu
#sivaramcheruvu.com

# ####Creating log file####
my_log <- file("my_log.txt")
sink(my_log, append = TRUE)
sink(my_log, append = TRUE, type = "message")

setwd("~/My Drive/Germany Paper")

####Loading Packages####

library(tidyverse)
library(lubridate)
library(haven)
library(lfe)
library(rdd)
library(rdrobust)
library(rddensity)
library(stargazer)
library(broom)
library(tidymodels)
library(rdlocrand)
library(rdpower)

####Putting Together Data####
ZA4588_v1_0_0 <- read_dta("ZA4588_v1-0-0.dta") #Cumulative ALLBUS data through 2016

a <- ZA4588_v1_0_0 %>% select(yborn, mborn, pt02, pt03, pt08, pt12, feduc,
                              meduc, ep03, dn07,dg03,dg04, year, sex,
                              xt02, xt03, educ, wghtpew)
a <- a %>% rename("trust_FCC" = pt02, "trust_parliament" = pt03, "trust_judicial" = pt08,  
                  "trust_federal" = pt12, "father_school" = feduc,
                  "mother_school" = meduc, 
                  "econ_sit" = ep03, 
                  "born_germany" = dn07, "youth_east_west" = dg03, "born_east_west" = dg04,
                  "year_survey" = year, "month_survey" = xt02, "date_survey" = xt03,
                   "weights" = wghtpew)
a$east <- ifelse(a$youth_east_west < 0, a$born_east_west, a$youth_east_west)
a <- a %>% filter(trust_FCC > 0, trust_judicial > 0, trust_parliament >0, trust_federal > 0,
                  east > 0, mborn > 0, yborn > 0, 
                  born_germany == 1, 
                  educ > -1) %>% 
  select(-c("born_east_west")) #Removing survey respondents that did not answer all questions
a$east <- ifelse(a$east == 1| a$east == 2,1,0) #Creating variable for being born in east germany
rm(ZA4588_v1_0_0)

#Adding 2018 ALLBUS data
ZA5272_v1_0_0 <- read_dta("ZA5272_v1-0-0.dta")
b <- ZA5272_v1_0_0 %>% select(yborn, mborn, pt02, pt03, pt08, pt12, feduc,
                              meduc, ep03,
                              dn07,dg03, sex,
                              xt02, xt03, educ, wghtpew)
b <- b %>% rename("trust_FCC" = pt02, "trust_parliament" = pt03, "trust_judicial" = pt08,   
                  "trust_federal" = pt12, "father_school" = feduc,
                  "mother_school" = meduc,
                  "econ_sit" = ep03,  
                  "born_germany" = dn07, "youth_east_west" = dg03,
                  "month_survey" = xt02, "date_survey" = xt03, "weights" = wghtpew)
b$year_survey <- "2018" #Adding year survey variable
b <- b %>% filter(trust_FCC > 0, trust_judicial > 0, trust_parliament > 0, trust_federal > 0,
                  youth_east_west > 0, mborn > 0, yborn > 0, 
                  born_germany == 1, 
                  educ > -1)
b$east <- ifelse(b$youth_east_west == 1| b$youth_east_west == 2,1,0) #Creating variable for being born in east germany
rm(ZA5272_v1_0_0)

#Full Data
dat <- rbind(a,b) #combining all of the data
rm(a,b)
dat <- dat %>% filter(educ > 2, educ != 6, year_survey !=1994) #Including only those that finished POS and excluding 1994 survey

####Coding variables####
dat$month_treat <- ifelse(dat$mborn >5,1,0) #Creating variable for being born June 1 or later
dat$trust_FCC <- as.numeric((dat$trust_FCC -1)/6) #Scaled from 1 - 7, subtracting 1 and dividing by six so variable is rescaled from 0 to 1
dat$trust_judicial <- as.numeric((dat$trust_judicial -1)/6) #Scaled from 1 - 7, subtracting 1 and dividing by six so variable is rescaled from 0 to 1
dat$trust_federal <- as.numeric((dat$trust_federal-1)/6) #Scaled from 1 - 7, subtracting 1 and dividing by six so variable is rescaled from 0 to 1
dat$trust_parliament <- as.numeric((dat$trust_parliament -1)/6) #Scaled from 1 - 7, subtracting 1 and dividing by six so variable is rescaled from 0 to 1
dat$mborn <- as.numeric(dat$mborn) 
dat$cohort <- ifelse(dat$yborn < 1983 & dat$yborn > 1972,1,0) #creating variable for those born within treated cohorts
dat$econ_sit <- ifelse(dat$econ_sit < 0, NA, dat$econ_sit) #Removing non-answers from financial situation variable
dat$econ_sit <- abs(dat$econ_sit - 5)/4 #rescaling from 0 to 1, and making 1 best economic situation

####Function for Bootstrap####

#Triple Diff bootstrap
btraps_blocks_weights <- function(reg,reg_dat){
  set.seed(12345) #Setting seed for replicability
  boots <-  reg_dat %>% nest(data = -c(mborn,cohort,east)) %>% #nesting data into clusters for block bootstrapping
    bootstraps(times = 1027) #Creating 1027 bootstrap replications. 27 replications to be dropped due to violation of full rank assumption when using full data
  a <- tidy(felm(reg, data = reg_dat, weights = reg_dat$weights)) #Running original model and putting in tidy format
  a <- a %>% select(term, estimate) #Selecting term and estimate and dropping std.error and pvalue columns
  a <- a %>% group_by(term) %>% summarize(beta.estimate=mean(estimate)) #putting estimates in same order for later and renaming the estimate column to avoid confusion
  boot_models <- boots %>%
    mutate(model = map(splits, ~ as_tibble(.) %>% unnest() %>% #Unnesting block data inside each bootstrap rplication
                         felm(reg, dat = ., weights = .$weights)), #running regression on each block within each bootstrap replication and creating variable with regression model
           coef_info = map(model, tidy)) #creating new variable with regression in tidy format
  boots_unnested <- boot_models %>% unnest(coef_info) #unnesting coefficients across all replications
  std.error_a <- left_join(boots_unnested,a, by = "term") #joining with estimates from original model
  drop <-  std.error_a %>% filter(is.na(estimate)) %>% select(id) #finding which bootstrap replications violate full rank assumption
  drop <- drop$id #turning the id numbers of the boostrap replications into a vector
  std.error_a <- std.error_a %>% filter(!id %in% drop) %>%  #dropping bootstrap replications that violate full rank assumption )
    group_by(term) %>% 
    summarize(critical_values = (estimate-beta.estimate)/sd(estimate), #calculating critical value t distribution
              estimate = estimate, #keeping bootstrap estimate
              std.error = sd(estimate)) #bootstrap standard error
  critical_values <- std.error_a %>% group_by(term) %>% 
    group_modify({~quantile(.x$critical_values, probs = c(0.05,0.95)) %>% #calulating critical values at 0.05 and 0.95 for alpha = 0.1
        tibble::enframe(name = "prob", value = "quantile")}) %>%
    pivot_wider(names_from = prob, values_from = quantile) %>%
    rename(lower_alpha="5%", upper_alpha="95%") #reshaping data and renaming variables for upper and lower alpha critical values
  std.errors <- std.error_a %>% distinct(term,.keep_all = T) %>% select(term,std.error) #separating out the bootstrapped standard errors
  a <- left_join(a, std.errors, by = "term") #combining standard errors with beta from original model
  a <- left_join(a, critical_values, by = "term") #joining critical values
  a <- a %>% select(term,beta.estimate, std.error, lower_alpha, upper_alpha) %>% rename(estimate = beta.estimate) %>%
    mutate(.lower = estimate + lower_alpha*std.error, .upper = estimate + upper_alpha*std.error) #calulating 90% confidence intervals
  return(a)
}

#Diff-in-Diff bootstrap
btraps_blocks_weights_DiD <- function(reg,reg_dat){
  set.seed(12345) #Setting seed for replicability
  boots <-  reg_dat %>% nest(data = -c(mborn,east)) %>% #nesting data into clusters for block bootstrapping
    bootstraps(times = 1027) #Creating 1027 bootstrap replications. 27 replications to be dropped due to violation of full rank assumption when using full data
  a <- tidy(felm(reg, data = reg_dat, weights = reg_dat$weights)) #Running original model and putting in tidy format
  a <- a %>% select(term, estimate) #Selecting term and estimate and dropping std.error and pvalue columns
  a <- a %>% group_by(term) %>% summarize(beta.estimate=mean(estimate)) #putting estimates in same order for later and renaming the estimate column to avoid confusion
  boot_models <- boots %>%
    mutate(model = map(splits, ~ as_tibble(.) %>% unnest() %>% #Unnesting block data inside each bootstrap rplication
                         felm(reg, dat = ., weights = .$weights)), #running regression on each block within each bootstrap replication and creating variable with regression model
           coef_info = map(model, tidy)) #creating new variable with regression in tidy format
  boots_unnested <- boot_models %>% unnest(coef_info) #unnesting coefficients across all replications
  std.error_a <- left_join(boots_unnested,a, by = "term") #joining with estimates from original model
  drop <-  std.error_a %>% filter(is.na(estimate)) %>% select(id) #finding which bootstrap replications violate full rank assumption
  drop <- drop$id #turning the id numbers of the boostrap replications into a vector
  std.error_a <- std.error_a %>% filter(!id %in% drop) %>%  #dropping bootstrap replications that violate full rank assumption )
    group_by(term) %>% 
    summarize(critical_values = (estimate-beta.estimate)/sd(estimate), #calculating critical value t distribution of wald statistics
              estimate = estimate, #keeping bootstrap estimate
              std.error = sd(estimate)) #bootstrap standard error
  critical_values <- std.error_a %>% group_by(term) %>% 
    group_modify({~quantile(.x$critical_values, probs = c(0.05,0.95)) %>% #calulating critical values at 0.05 and 0.95 for alpha = 0.1
        tibble::enframe(name = "prob", value = "quantile")}) %>%
    pivot_wider(names_from = prob, values_from = quantile) %>%
    rename(lower_alpha="5%", upper_alpha="95%") #reshaping data and renaming variables for upper and lower alpha critical values
  std.errors <- std.error_a %>% distinct(term,.keep_all = T) %>% select(term,std.error) #separating out the bootstrapped standard errors
  a <- left_join(a, std.errors, by = "term") #combining standard errors with beta from original model
  a <- left_join(a, critical_values, by = "term") #joining critical values
  a <- a %>% select(term,beta.estimate, std.error, lower_alpha, upper_alpha) %>% rename(estimate = beta.estimate) %>%
    mutate(.lower = estimate + lower_alpha*std.error, .upper = estimate + upper_alpha*std.error) #calulating 90% confidence intervals
  return(a)
}

####Pooled DiDiD Models (Table A7)####

reg_dat <- dat %>% filter(yborn > 1949) %>% #Filtering those born after 1949
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables

#Model 1 Trust FCC

mod_1 <- felm(trust_FCC~cohort*month_treat*east |year_survey|0|0, data  = filter(dat,yborn >1949), weights  = filter(dat,yborn >1949)$weights)
reg <- trust_FCC~cohort*month_treat*east |year_survey|0|0
a <- btraps_blocks_weights(reg,reg_dat)
a <- a %>% mutate(model = "DiDiD", survey_year = "Pooled", subset = "School Environment", dv = "Trust FCC")
DiDiD_pooled <- a %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 2 Trust Parliament

mod_3 <- felm(trust_parliament~cohort*month_treat*east|year_survey|0|0, data  = filter(dat,yborn >1949), weights  = filter(dat,yborn >1949)$weights)
reg <- trust_parliament~cohort*month_treat*east|year_survey|0|0
c <- btraps_blocks_weights(reg,reg_dat)
c <- c %>% mutate(model = "DiDiD", survey_year = "Pooled", subset = "School Environment", dv = "Trust Parliament")
DiDiD_parliament_pooled <- c %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 3 Trust Federal

mod_4 <- felm(trust_federal~cohort*month_treat*east|year_survey|0|0, data  = filter(dat,yborn >1949), weights  = filter(dat,yborn >1949)$weights)
reg <- trust_federal~cohort*month_treat*east|year_survey|0|0
d <- btraps_blocks_weights(reg,reg_dat)
d <- d %>% mutate(model = "DiDiD", survey_year = "Pooled", subset = "School Environment", dv = "Trust Federal")
DiDiD_federal_pooled <- d %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 4 Financial Situation

mod_5 <- felm(econ_sit~cohort*month_treat*east|year_survey|0|0, data  = filter(dat,yborn >1949), weights  = filter(dat,yborn >1949)$weights)
reg <- econ_sit~cohort*month_treat*east|year_survey|0|0
e <- btraps_blocks_weights(reg,reg_dat)
e <- e %>% mutate(model = "DiDiD", survey_year = "Pooled", subset = "School Environment", dv = "Financial Situation")
DiDiD_econsit_pooled <- e %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a <- a %>% arrange(factor(term, levels = c("cohort", "month_treat", "east", "cohort:month_treat",
                                           "cohort:east", "month_treat:east", "cohort:month_treat:east")))
c <- c %>% arrange(factor(term, levels = c("cohort", "month_treat", "east", "cohort:month_treat",
                                           "cohort:east", "month_treat:east", "cohort:month_treat:east")))
d <- d %>% arrange(factor(term, levels = c("cohort", "month_treat", "east", "cohort:month_treat",
                                           "cohort:east", "month_treat:east", "cohort:month_treat:east")))
e <- e %>% arrange(factor(term, levels = c("cohort", "month_treat", "east", "cohort:month_treat",
                                           "cohort:east", "month_treat:east", "cohort:month_treat:east")))

####Pooled DiDiD Models Table A7
stargazer(mod_1,mod_3,mod_4,mod_5,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Cohort}", "\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{Cohort:After May}", "\\textsc{Cohort:East}",
                               "\\textsc{After May:East}","\\textsc{After May:East:Cohort}"),
          title = "Pooled Difference-in-Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a$.lower,a$.upper),cbind(c$.lower,c$.upper),cbind(d$.lower,d$.upper),cbind(e$.lower,e$.upper)),
          add.lines = list(c("Survey-Year Fixed-Effects","Yes", "Yes", "Yes", "Yes"),
                           c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A7.tex")

####Pooled DiD Models (Table A13)####

#FCC
mod_1 <- felm(trust_FCC~month_treat*east |year_survey|0|0, data  = filter(dat,cohort == 1), weights  = filter(dat,cohort == 1)$weights)
reg <- trust_FCC~month_treat*east |year_survey|0|0
a_DiD_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_FCC <- a_DiD_FCC %>% mutate(model = "DiD", survey_year = "Pooled", subset = "School Environment", dv = "Trust FCC")
DiD_pooled_FCC <- a_DiD_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Parliament
mod_2 <- felm(trust_parliament~month_treat*east |year_survey|0|0, data  = filter(dat,cohort == 1), weights  = filter(dat,cohort == 1)$weights)
reg <- trust_parliament~month_treat*east |year_survey|0|0
a_DiD_Parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_Parliament <- a_DiD_Parliament %>% mutate(model = "DiD", survey_year = "Pooled", subset = "School Environment", dv = "Trust Parliament")
DiD_pooled_Parliament <- a_DiD_Parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Federal
mod_3 <- felm(trust_federal~month_treat*east |year_survey|0|0, data  = filter(dat,cohort == 1), weights  = filter(dat,cohort == 1)$weights)
reg <- trust_federal~month_treat*east |year_survey|0|0
a_DiD_Federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_Federal <- a_DiD_Federal %>% mutate(model = "DiD", survey_year = "Pooled", subset = "School Environment", dv = "Trust Federal")
DiD_pooled_Federal <- a_DiD_Federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Situation

mod_4 <- felm(econ_sit~month_treat*east |year_survey|0|0, data  = filter(dat,cohort == 1), weights  = filter(dat,cohort == 1)$weights)
reg <- econ_sit~month_treat*east |year_survey|0|0
a_DiD_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_econsit <- a_DiD_econsit %>% mutate(model = "DiD", survey_year = "Pooled", subset = "School Environment", dv = "Financial Situation")
DiD_pooled_econsit <- a_DiD_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_FCC <- a_DiD_FCC %>% arrange(factor(term, levels = c("month_treat", "east", "month_treat:east")))
a_DiD_Parliament <- a_DiD_Parliament %>% arrange(factor(term, levels = c("month_treat", "east", "month_treat:east")))
a_DiD_Federal <- a_DiD_Federal %>% arrange(factor(term, levels = c("month_treat", "east", "month_treat:east")))
a_DiD_econsit <- a_DiD_econsit %>% arrange(factor(term, levels = c("month_treat", "east", "month_treat:east")))

####Pooled DiD Models Table A13
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Pooled Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_FCC$.lower,a_DiD_FCC$.upper),cbind(a_DiD_Parliament$.lower,a_DiD_Parliament$.upper),
                           cbind(a_DiD_Federal$.lower,a_DiD_Federal$.upper),cbind(a_DiD_econsit$.lower,a_DiD_econsit$.upper)),
          add.lines = list(c("Survey-Year Fixed-Effects","Yes", "Yes", "Yes", "Yes"),
                           c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A13.tex")

####Pooled RD models (Table A1)####
east <- dat %>% filter(east == 1, cohort == 1) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_pooled <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                                model = "RD", survey_year = "Pooled", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat,east == 1, cohort == 1)) #place holder for regression table

#RD Parliament Pooled
east <- dat %>% filter(east == 1, cohort == 1) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_pooled <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                               model = "RD", survey_year = "Pooled", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat,east == 1, cohort == 1)) #place holder for regression table


#RD Federal Pooled
east <- dat %>% filter(east == 1, cohort == 1) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_pooled <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                            model = "RD", survey_year = "Pooled", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat,east == 1, cohort == 1))

#RD Financial Situation Pooled
east <- dat %>% filter(east == 1, cohort == 1) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_pooled <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                            model = "RD", survey_year = "Pooled", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat,east == 1, cohort == 1))


pooled_models <- rbind(DiDiD_pooled,DiDiD_parliament_pooled,DiDiD_federal_pooled,DiDiD_econsit_pooled,DiD_pooled_FCC, DiD_pooled_Federal, 
                       DiD_pooled_Parliament, DiD_pooled_econsit, rd_environment_pooled, rd_parliament_pooled, rd_federal_pooled, rd_econsit_pooled) #Combining all models into one dataframe

pooled_models_full <- rbind(a,b,c,d,e,a_DiD_FCC,a_DiD_Parliament, a_DiD_Federal,a_DiD_econsit)

#Pooled RD Models Table A1
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Pooled Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_pooled$estimate[1],rd_parliament_pooled$estimate[1],rd_federal_pooled$estimate[1],rd_econsit_pooled$estimate[1]),
          ci.custom = list(cbind(rd_environment_pooled$lower,rd_environment_pooled$upper), 
                           cbind(rd_parliament_pooled$lower,rd_parliament_pooled$upper),
                           cbind(rd_federal_pooled$lower,rd_federal_pooled$upper),
                           cbind(rd_econsit_pooled$lower,rd_econsit_pooled$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A1.tex")

####Survey Year 2000 DiDiD Models (Table A8)####

reg_dat <- dat %>% filter(yborn > 1949, year_survey == 2000) %>% #Filtering those born after 1949
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables

#Model 1 Trust FCC survey year 2000

mod_1 <- felm(trust_FCC~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2000), weights  = filter(dat,yborn >1949, year_survey == 2000)$weights)
reg <- trust_FCC~cohort*month_treat*east |0|0|0
a_2000 <- btraps_blocks_weights(reg,reg_dat)
a_2000 <- a_2000 %>% mutate(model = "DiDiD", survey_year = "2000", subset = "School Environment", dv = "Trust FCC")
DiDiD_2000 <- a_2000 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 2 Trust Parliament survey year 2000

mod_3 <- felm(trust_parliament~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2000), weights  = filter(dat,yborn >1949, year_survey == 2000)$weights)
reg <- trust_parliament~cohort*month_treat*east |0|0|0
c_2000 <- btraps_blocks_weights(reg,reg_dat)
c_2000 <- c_2000 %>% mutate(model = "DiDiD", survey_year = "2000", subset = "School Environment", dv = "Trust Parliament")
DiDiD_parliament_2000 <- c_2000 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 3 Trust Federal survey year 2000

mod_4 <- felm(trust_federal~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2000), weights  = filter(dat,yborn >1949, year_survey == 2000)$weights)
reg <- trust_federal~cohort*month_treat*east |0|0|0
d_2000 <- btraps_blocks_weights(reg,reg_dat)
d_2000 <- d_2000 %>% mutate(model = "DiDiD", survey_year = "2000", subset = "School Environment", dv = "Trust Federal")
DiDiD_federal_2000 <- d_2000 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 4 Financial Situation survey year 2000
mod_5 <- felm(econ_sit~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2000), weights  = filter(dat,yborn >1949, year_survey == 2000)$weights)
reg <- econ_sit~cohort*month_treat*east |0|0|0
e_2000 <- btraps_blocks_weights(reg,reg_dat)
e_2000 <- e_2000 %>% mutate(model = "DiDiD", survey_year = "2000", subset = "School Environment", dv = "Financial Situation")
DiDiD_econsit_2000 <- e_2000 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_2000 <- a_2000 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
c_2000 <- c_2000 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
d_2000 <- d_2000 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
e_2000 <- e_2000 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))

####2000 DiDiD Models Table A8
stargazer(mod_1,mod_3,mod_4,mod_5,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Cohort}", "\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{Cohort:After May}", "\\textsc{Cohort:East}",
                               "\\textsc{After May:East}","\\textsc{After May:East:Cohort}"),
          title = "Survey Year 2000 Difference-in-Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_2000$.lower,a_2000$.upper),cbind(c_2000$.lower,c_2000$.upper),
                           cbind(d_2000$.lower,d_2000$.upper),cbind(e_2000$.lower,e_2000$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = "(Intercept)",
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A8.tex")

####2000 DiD Models (Table A14)####

#Trust FCC 2000
reg_dat <- dat %>% filter(cohort == 1, year_survey == 2000) %>% #Filtering those within relevant cohort and survey year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2000), weights  = filter(dat,cohort == 1, year_survey == 2000)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_2000_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2000_FCC <- a_DiD_2000_FCC %>% mutate(model = "DiD", survey_year = "2000", subset = "School Environment", dv = "Trust FCC")
DiD_2000_FCC <- a_DiD_2000_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 2000
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2000), weights  = filter(dat,cohort == 1, year_survey == 2000)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_2000_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2000_parliament <- a_DiD_2000_parliament %>% mutate(model = "DiD", survey_year = "2000", subset = "School Environment", dv = "Trust Parliament")
DiD_2000_parliament <- a_DiD_2000_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 2000
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2000), weights  = filter(dat,cohort == 1, year_survey == 2000)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_2000_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2000_federal <- a_DiD_2000_federal %>% mutate(model = "DiD", survey_year = "2000", subset = "School Environment", dv = "Trust Federal")
DiD_2000_federal <- a_DiD_2000_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 2000
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2000), weights  = filter(dat,cohort == 1, year_survey == 2000)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_2000_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2000_econsit <- a_DiD_2000_econsit %>% mutate(model = "DiD", survey_year = "2000", subset = "School Environment", dv = "Financial Situation")
DiD_2000_econsit <- a_DiD_2000_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_2000_FCC <- a_DiD_2000_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2000_parliament <- a_DiD_2000_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2000_federal <- a_DiD_2000_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2000_econsit <- a_DiD_2000_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####2000 DiD Models Table A14
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Survey Year 2000 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_2000_FCC$.lower,a_DiD_2000_FCC$.upper),cbind(a_DiD_2000_parliament$.lower,a_DiD_2000_parliament$.upper),
                           cbind(a_DiD_2000_federal$.lower,a_DiD_2000_federal$.upper),cbind(a_DiD_2000_econsit$.lower,a_DiD_2000_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A14.tex")

####2000 RD Models (Table A2)####

#RD FCC 2000
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2000) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_2000 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", survey_year = "2000", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2000))

#RD Parliament 2000
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2000) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_2000 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", survey_year = "2000", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2000))

#RD Federal 2000
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2000) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_2000 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", survey_year = "2000", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2000))

#RD Financial Situation 2000
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2000) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_2000 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", survey_year = "2000", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2000))


#RD 2000 models Table A2
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Survey Year 2000 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_2000$estimate[1],rd_parliament_2000$estimate[1],rd_federal_2000$estimate[1],rd_econsit_2000$estimate[1]),
          ci.custom = list(cbind(rd_environment_2000$lower,rd_environment_2000$upper), 
                           cbind(rd_parliament_2000$lower,rd_parliament_2000$upper),
                           cbind(rd_federal_2000$lower,rd_federal_2000$upper),
                           cbind(rd_econsit_2000$lower,rd_econsit_2000$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A2.tex")

models_2000 <- rbind(DiDiD_2000,DiDiD_parliament_2000,DiDiD_federal_2000, DiDiD_econsit_2000, DiD_2000_FCC, DiD_2000_federal, DiD_2000_parliament, DiD_2000_econsit,
                     rd_environment_2000,rd_parliament_2000,rd_federal_2000, rd_econsit_2000)

models_full_2000 <- rbind(a_2000,c_2000,d_2000,e_2000,a_DiD_2000_FCC,a_DiD_2000_parliament, a_DiD_2000_federal, a_DiD_2000_econsit)

####Survey Year 2002 DiDiD Models (Table A9)####

reg_dat <- dat %>% filter(yborn > 1949, year_survey == 2002) %>% #Filtering those born after 1949
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables

#Model 1 Trust FCC survey year 2002

mod_1 <- felm(trust_FCC~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2002), weights  = filter(dat,yborn >1949, year_survey == 2002)$weights)
reg <- trust_FCC~cohort*month_treat*east |0|0|0
a_2002 <- btraps_blocks_weights(reg,reg_dat)
a_2002 <- a_2002 %>% mutate(model = "DiDiD", survey_year = "2002", subset = "School Environment", dv = "Trust FCC")
DiDiD_2002 <- a_2002 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 2 Trust Parliament survey year 2002

mod_3 <- felm(trust_parliament~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2002), weights  = filter(dat,yborn >1949, year_survey == 2002)$weights)
reg <- trust_parliament~cohort*month_treat*east |0|0|0
c_2002 <- btraps_blocks_weights(reg,reg_dat)
c_2002 <- c_2002 %>% mutate(model = "DiDiD", survey_year = "2002", subset = "School Environment", dv = "Trust Parliament")
DiDiD_parliament_2002 <- c_2002 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 3 Trust Federal survey year 2002

mod_4 <- felm(trust_federal~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2002), weights  = filter(dat,yborn >1949, year_survey == 2002)$weights)
reg <- trust_federal~cohort*month_treat*east|0|0|0
d_2002 <- btraps_blocks_weights(reg,reg_dat)
d_2002 <- d_2002 %>% mutate(model = "DiDiD", survey_year = "2002", subset = "School Environment", dv = "Trust Federal")
DiDiD_federal_2002 <- d_2002 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 4 Financial Situation survey year 2002
mod_5 <- felm(econ_sit~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2002), weights  = filter(dat,yborn >1949, year_survey == 2002)$weights)
reg <- econ_sit~cohort*month_treat*east|0|0|0
e_2002 <- btraps_blocks_weights(reg,reg_dat)
e_2002 <- e_2002 %>% mutate(model = "DiDiD", survey_year = "2002", subset = "School Environment", dv = "Financial Situation")
DiDiD_econsit_2002 <- e_2002 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_2002 <- a_2002 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
c_2002 <- c_2002 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
d_2002 <- d_2002 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
e_2002 <- e_2002 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))

####2002 DiDiD Models Table A9
stargazer(mod_1,mod_3,mod_4,mod_5,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Cohort}", "\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{Cohort:After May}", "\\textsc{Cohort:East}",
                               "\\textsc{After May:East}","\\textsc{After May:East:Cohort}"),
          title = "Survey Year 2002 Difference-in-Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_2002$.lower,a_2002$.upper),cbind(c_2002$.lower,c_2002$.upper),
                           cbind(d_2002$.lower,d_2002$.upper),cbind(e_2002$.lower,e_2002$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = "(Intercept)",
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A9.tex")

####DiD Models 2002 (Table A15)####

#Trust FCC 2002
reg_dat <- dat %>% filter(cohort == 1, year_survey == 2002) %>% #Filtering those within relevant cohort and survey year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2002), weights  = filter(dat,cohort == 1, year_survey == 2002)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_2002_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2002_FCC <- a_DiD_2002_FCC %>% mutate(model = "DiD", survey_year = "2002", subset = "School Environment", dv = "Trust FCC")
DiD_2002_FCC <- a_DiD_2002_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 2002
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2002), weights  = filter(dat,cohort == 1, year_survey == 2002)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_2002_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2002_parliament <- a_DiD_2002_parliament %>% mutate(model = "DiD", survey_year = "2002", subset = "School Environment", dv = "Trust Parliament")
DiD_2002_parliament <- a_DiD_2002_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 2002
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2002), weights  = filter(dat,cohort == 1, year_survey == 2002)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_2002_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2002_federal <- a_DiD_2002_federal %>% mutate(model = "DiD", survey_year = "2002", subset = "School Environment", dv = "Trust Federal")
DiD_2002_federal <- a_DiD_2002_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 2002
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2002), weights  = filter(dat,cohort == 1, year_survey == 2002)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_2002_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2002_econsit <- a_DiD_2002_econsit %>% mutate(model = "DiD", survey_year = "2002", subset = "School Environment", dv = "Financial Situation")
DiD_2002_econsit <- a_DiD_2002_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_2002_FCC <- a_DiD_2002_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2002_parliament <- a_DiD_2002_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2002_federal <- a_DiD_2002_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2002_econsit <- a_DiD_2002_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####2002 DiD Models Table A15
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Survey Year 2002 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_2002_FCC$.lower,a_DiD_2002_FCC$.upper),cbind(a_DiD_2002_parliament$.lower,a_DiD_2002_parliament$.upper),
                           cbind(a_DiD_2002_federal$.lower,a_DiD_2002_federal$.upper),cbind(a_DiD_2002_econsit$.lower,a_DiD_2002_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A15.tex")

####2002 RD Models Table A3####

#RD FCC 2002
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2002) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_2002 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", survey_year = "2002", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2002)) #model placeholder for tables

#RD Parliament 2002
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2002) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_2002 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", survey_year = "2002", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2002)) #model placeholder for tables


#RD Federal 2002
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2002) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_2002 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", survey_year = "2002", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2002)) #model placeholder for tables


#RD Financial Situation 2002
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2002) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_2002 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", survey_year = "2002", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2002)) #model placeholder for tables

####2002 RD Models (Table A3)
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Survey Year 2002 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_2002$estimate[1],rd_parliament_2002$estimate[1],rd_federal_2002$estimate[1],rd_econsit_2002$estimate[1]),
          ci.custom = list(cbind(rd_environment_2002$lower,rd_environment_2002$upper), 
                           cbind(rd_parliament_2002$lower,rd_parliament_2002$upper),
                           cbind(rd_federal_2002$lower,rd_federal_2002$upper),
                           cbind(rd_econsit_2002$lower,rd_econsit_2002$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A3.tex")

models_2002 <- rbind(DiDiD_2002,DiDiD_parliament_2002,DiDiD_federal_2002, DiDiD_econsit_2002, DiD_2002_FCC, DiD_2002_federal, DiD_2002_parliament, DiD_2002_econsit,
                     rd_environment_2002,rd_parliament_2002,rd_federal_2002, rd_econsit_2002)

models_full_2002 <- rbind(a_2002,c_2002,d_2002,e_2002,a_DiD_2002_FCC,a_DiD_2002_parliament, a_DiD_2002_federal, a_DiD_2002_econsit)


####Survey Year 2008 DiDiD Models (Table A10)####

reg_dat <- dat %>% filter(yborn > 1949, year_survey == 2008) %>% #Filtering those born after 1949
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables

#Model 1 Trust FCC survey year 2008

mod_1 <- felm(trust_FCC~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2008), weights  = filter(dat,yborn >1949, year_survey == 2008)$weights)
reg <- trust_FCC~cohort*month_treat*east |0|0|0
a_2008 <- btraps_blocks_weights(reg,reg_dat)
a_2008 <- a_2008 %>% mutate(model = "DiDiD", survey_year = "2008", subset = "School Environment", dv = "Trust FCC")
DiDiD_2008 <- a_2008 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 2 Trust Parliament survey year 2008

mod_3 <- felm(trust_parliament~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2008), weights  = filter(dat,yborn >1949, year_survey == 2008)$weights)
reg <- trust_parliament~cohort*month_treat*east|0|0|0
c_2008 <- btraps_blocks_weights(reg,reg_dat)
c_2008 <- c_2008 %>% mutate(model = "DiDiD", survey_year = "2008", subset = "School Environment", dv = "Trust Parliament")
DiDiD_parliament_2008 <- c_2008 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 3 Trust Federal survey year 2008

mod_4 <- felm(trust_federal~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2008), weights  = filter(dat,yborn >1949, year_survey == 2008)$weights)
reg <- trust_federal~cohort*month_treat*east|0|0|0
d_2008 <- btraps_blocks_weights(reg,reg_dat)
d_2008 <- d_2008 %>% mutate(model = "DiDiD", survey_year = "2008", subset = "School Environment", dv = "Trust Federal")
DiDiD_federal_2008 <- d_2008 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 4 Financial Situation survey year 2008
mod_5 <- felm(econ_sit~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2008), weights  = filter(dat,yborn >1949, year_survey == 2008)$weights)
reg <- econ_sit~cohort*month_treat*east|0|0|0
e_2008 <- btraps_blocks_weights(reg,reg_dat)
e_2008 <- e_2008 %>% mutate(model = "DiDiD", survey_year = "2008", subset = "School Environment", dv = "Financial Situation")
DiDiD_econsit_2008 <- e_2008 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_2008 <- a_2008 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
c_2008 <- c_2008 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
d_2008 <- d_2008 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
e_2008 <- e_2008 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))

####2008 DiDiD Models Table (Table A10)
stargazer(mod_1,mod_3,mod_4,mod_5,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Cohort}", "\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{Cohort:After May}", "\\textsc{Cohort:East}",
                               "\\textsc{After May:East}","\\textsc{After May:East:Cohort}"),
          title = "Survey Year 2008 Difference-in-Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_2008$.lower,a_2008$.upper),cbind(c_2008$.lower,c_2008$.upper),
                           cbind(d_2008$.lower,d_2008$.upper),cbind(e_2008$.lower,e_2008$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = "(Intercept)",
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A10.tex")

####Survey Year 2008 DiD Models (Table A16)####

#Trust FCC 2008
reg_dat <- dat %>% filter(cohort == 1, year_survey == 2008) %>% #Filtering those within relevant cohort and survey year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2008), weights  = filter(dat,cohort == 1, year_survey == 2008)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_2008_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2008_FCC <- a_DiD_2008_FCC %>% mutate(model = "DiD", survey_year = "2008", subset = "School Environment", dv = "Trust FCC")
DiD_2008_FCC <- a_DiD_2008_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 2008
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2008), weights  = filter(dat,cohort == 1, year_survey == 2008)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_2008_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2008_parliament <- a_DiD_2008_parliament %>% mutate(model = "DiD", survey_year = "2008", subset = "School Environment", dv = "Trust Parliament")
DiD_2008_parliament <- a_DiD_2008_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 2008
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2008), weights  = filter(dat,cohort == 1, year_survey == 2008)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_2008_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2008_federal <- a_DiD_2008_federal %>% mutate(model = "DiD", survey_year = "2008", subset = "School Environment", dv = "Trust Federal")
DiD_2008_federal <- a_DiD_2008_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 2008
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2008), weights  = filter(dat,cohort == 1, year_survey == 2008)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_2008_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2008_econsit <- a_DiD_2008_econsit %>% mutate(model = "DiD", survey_year = "2008", subset = "School Environment", dv = "Financial Situation")
DiD_2008_econsit <- a_DiD_2008_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_2008_FCC <- a_DiD_2008_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2008_parliament <- a_DiD_2008_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2008_federal <- a_DiD_2008_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2008_econsit <- a_DiD_2008_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####2008 DiD Models Table A16
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Survey Year 2008 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_2008_FCC$.lower,a_DiD_2008_FCC$.upper),cbind(a_DiD_2008_parliament$.lower,a_DiD_2008_parliament$.upper),
                           cbind(a_DiD_2008_federal$.lower,a_DiD_2008_federal$.upper),cbind(a_DiD_2008_econsit$.lower,a_DiD_2008_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A16.tex")

####Survey Year 2008 RD Models (Table A4)####

#RD FCC 2008
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2008) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_2008 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", survey_year = "2008", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2008)) #placeholder for table

#RD Parliament 2008
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2008) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_2008 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", survey_year = "2008", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2008)) #placeholder for table


#RD Federal 2008
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2008) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_2008 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", survey_year = "2008", subset = "School Environment", dv = "Trust Federal")

mod_3 <- lm(trust_federal~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2008)) #placeholder for table


#RD Financial Situation 2008
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2008) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_2008 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", survey_year = "2008", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2008)) #placeholder for table


stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Survey Year 2008 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_2008$estimate[1],rd_parliament_2008$estimate[1],rd_federal_2008$estimate[1],rd_econsit_2008$estimate[1]),
          ci.custom = list(cbind(rd_environment_2008$lower,rd_environment_2008$upper), 
                           cbind(rd_parliament_2008$lower,rd_parliament_2008$upper),
                           cbind(rd_federal_2008$lower,rd_federal_2008$upper),
                           cbind(rd_econsit_2008$lower,rd_econsit_2008$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A4.tex")

models_2008 <- rbind(DiDiD_2008,DiDiD_parliament_2008,DiDiD_federal_2008, DiDiD_econsit_2008, DiD_2008_FCC, DiD_2008_federal, DiD_2008_parliament, DiD_2008_econsit,
                     rd_environment_2008,rd_parliament_2008,rd_federal_2008, rd_econsit_2008)

models_full_2008 <- rbind(a_2008,c_2008,d_2008,e_2008,a_DiD_2008_FCC,a_DiD_2008_parliament, a_DiD_2008_federal, a_DiD_2008_econsit)

####Survey Year 2012 DiDiD Models (Table A11)####

reg_dat <- dat %>% filter(yborn > 1949, year_survey == 2012) %>% #Filtering those born after 1949
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables

#Model 1 Trust FCC survey year 2012

mod_1 <- felm(trust_FCC~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2012), weights  = filter(dat,yborn >1949, year_survey == 2012)$weights)
reg <- trust_FCC~cohort*month_treat*east |0|0|0
a_2012 <- btraps_blocks_weights(reg,reg_dat)
a_2012 <- a_2012 %>% mutate(model = "DiDiD", survey_year = "2012", subset = "School Environment", dv = "Trust FCC")
DiDiD_2012 <- a_2012 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 2 Trust Parliament survey year 2012

mod_3 <- felm(trust_parliament~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2012), weights  = filter(dat,yborn >1949, year_survey == 2012)$weights)
reg <- trust_parliament~cohort*month_treat*east|0|0|0
c_2012 <- btraps_blocks_weights(reg,reg_dat)
c_2012 <- c_2012 %>% mutate(model = "DiDiD", survey_year = "2012", subset = "School Environment", dv = "Trust Parliament")
DiDiD_parliament_2012 <- c_2012 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 3 Trust Federal survey year 2012

mod_4 <- felm(trust_federal~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2012), weights  = filter(dat,yborn >1949, year_survey == 2012)$weights)
reg <- trust_federal~cohort*month_treat*east|0|0|0
d_2012 <- btraps_blocks_weights(reg,reg_dat)
d_2012 <- d_2012 %>% mutate(model = "DiDiD", survey_year = "2012", subset = "School Environment", dv = "Trust Federal")
DiDiD_federal_2012 <- d_2012 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 4 Financial Situation survey year 2012
mod_5 <- felm(econ_sit~cohort*month_treat*east|0|0|0, data  = filter(dat,yborn >1949, year_survey == 2012), weights  = filter(dat,yborn >1949, year_survey == 2012)$weights)
reg <- econ_sit~cohort*month_treat*east|0|0|0
e_2012 <- btraps_blocks_weights(reg,reg_dat)
e_2012 <- e_2012 %>% mutate(model = "DiDiD", survey_year = "2012", subset = "School Environment", dv = "Financial Situation")
DiDiD_econsit_2012 <- e_2012 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_2012 <- a_2012 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
c_2012 <- c_2012 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
d_2012 <- d_2012 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
e_2012 <- e_2012 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))

####2012 DiDiD Models Table A11
stargazer(mod_1,mod_3,mod_4,mod_5,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Cohort}", "\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{Cohort:After May}", "\\textsc{Cohort:East}",
                               "\\textsc{After May:East}","\\textsc{After May:East:Cohort}"),
          title = "Survey Year 2012 Difference-in-Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_2012$.lower,a_2012$.upper),cbind(c_2012$.lower,c_2012$.upper),
                           cbind(d_2012$.lower,d_2012$.upper),cbind(e_2012$.lower,e_2012$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = "(Intercept)",
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A11.tex")

####Survey Year 2012 DiD Models (Table A17)####

#Trust FCC 2012
reg_dat <- dat %>% filter(cohort == 1, year_survey == 2012) %>% #Filtering those within relevant cohort and survey year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2012), weights  = filter(dat,cohort == 1, year_survey == 2012)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_2012_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2012_FCC <- a_DiD_2012_FCC %>% mutate(model = "DiD", survey_year = "2012", subset = "School Environment", dv = "Trust FCC")
DiD_2012_FCC <- a_DiD_2012_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 2012
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2012), weights  = filter(dat,cohort == 1, year_survey == 2012)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_2012_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2012_parliament <- a_DiD_2012_parliament %>% mutate(model = "DiD", survey_year = "2012", subset = "School Environment", dv = "Trust Parliament")
DiD_2012_parliament <- a_DiD_2012_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 2012
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2012), weights  = filter(dat,cohort == 1, year_survey == 2012)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_2012_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2012_federal <- a_DiD_2012_federal %>% mutate(model = "DiD", survey_year = "2012", subset = "School Environment", dv = "Trust Federal")
DiD_2012_federal <- a_DiD_2012_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 2012
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2012), weights  = filter(dat,cohort == 1, year_survey == 2012)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_2012_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2012_econsit <- a_DiD_2012_econsit %>% mutate(model = "DiD", survey_year = "2012", subset = "School Environment", dv = "Financial Situation")
DiD_2012_econsit <- a_DiD_2012_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_2012_FCC <- a_DiD_2012_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2012_parliament <- a_DiD_2012_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2012_federal <- a_DiD_2012_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2012_econsit <- a_DiD_2012_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####2012 DiD Models Table A17
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Survey Year 2012 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_2012_FCC$.lower,a_DiD_2012_FCC$.upper),cbind(a_DiD_2012_parliament$.lower,a_DiD_2012_parliament$.upper),
                           cbind(a_DiD_2012_federal$.lower,a_DiD_2012_federal$.upper),cbind(a_DiD_2012_econsit$.lower,a_DiD_2012_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A17.tex")

####Survey Year 2012 RD Models (Table A5)####

#RD FCC 2012
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2012) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_2012 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", survey_year = "2012", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2012)) #placeholder for table

#RD Parliament 2012
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2012) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_2012 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", survey_year = "2012", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2012)) #placeholder for table

#RD Federal 2012
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2012) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_2012 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", survey_year = "2012", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2012)) #placeholder for table

#RD Financial Situation 2012
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2012) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_2012 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", survey_year = "2012", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat,east == 1, cohort == 1, year_survey == 2012)) #placeholder for table

#RD Models 2012 table (Table A5)
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Survey Year 2012 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_2012$estimate[1],rd_parliament_2012$estimate[1],rd_federal_2012$estimate[1],rd_econsit_2012$estimate[1]),
          ci.custom = list(cbind(rd_environment_2012$lower,rd_environment_2012$upper), 
                           cbind(rd_parliament_2012$lower,rd_parliament_2012$upper),
                           cbind(rd_federal_2012$lower,rd_federal_2012$upper),
                           cbind(rd_econsit_2012$lower,rd_econsit_2012$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A5.tex")

models_2012 <- rbind(DiDiD_2012,DiDiD_parliament_2012,DiDiD_federal_2012, DiDiD_econsit_2012, DiD_2012_FCC, DiD_2012_federal, DiD_2012_parliament, DiD_2012_econsit,
                     rd_environment_2012,rd_parliament_2012,rd_federal_2012, rd_econsit_2012)

models_full_2012 <- rbind(a_2012,c_2012,d_2012,e_2012,a_DiD_2012_FCC,a_DiD_2012_parliament, a_DiD_2012_federal, a_DiD_2012_econsit)

####Survey Year 2018 DiDiD Models (Table A12)####

reg_dat <- dat %>% filter(yborn > 1949, year_survey == 2018) %>% #Filtering those born after 1949
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables

#Model 1 Trust FCC survey year 2018

mod_1 <- felm(trust_FCC~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2018), weights  = filter(dat,yborn >1949, year_survey == 2018)$weights)
reg <- trust_FCC~cohort*month_treat*east |0|0|0
a_2018 <- btraps_blocks_weights(reg,reg_dat)
a_2018 <- a_2018 %>% mutate(model = "DiDiD", survey_year = "2018", subset = "School Environment", dv = "Trust FCC")
DiDiD_2018 <- a_2018 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)
``
#Model 2 Trust Parliament survey year 2018

mod_3 <- felm(trust_parliament~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2018), weights  = filter(dat,yborn >1949, year_survey == 2018)$weights)
reg <- trust_parliament~cohort*month_treat*east |0|0|0
c_2018 <- btraps_blocks_weights(reg,reg_dat)
c_2018 <- c_2018 %>% mutate(model = "DiDiD", survey_year = "2018", subset = "School Environment", dv = "Trust Parliament")
DiDiD_parliament_2018 <- c_2018 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 3 Trust Federal survey year 2018

mod_4 <- felm(trust_federal~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2018), weights  = filter(dat,yborn >1949, year_survey == 2018)$weights)
reg <- trust_federal~cohort*month_treat*east |0|0|0
d_2018 <- btraps_blocks_weights(reg,reg_dat)
d_2018 <- d_2018 %>% mutate(model = "DiDiD", survey_year = "2018", subset = "School Environment", dv = "Trust Federal")
DiDiD_federal_2018 <- d_2018 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Model 4 Financial Situation survey year 2018
mod_5 <- felm(econ_sit~cohort*month_treat*east |0|0|0, data  = filter(dat,yborn >1949, year_survey == 2018), weights  = filter(dat,yborn >1949, year_survey == 2018)$weights)
reg <- econ_sit~cohort*month_treat*east |0|0|0
e_2018 <- btraps_blocks_weights(reg,reg_dat)
e_2018 <- e_2018 %>% mutate(model = "DiDiD", survey_year = "2018", subset = "School Environment", dv = "Financial Situation")
DiDiD_econsit_2018 <- e_2018 %>% filter(term=="cohort:month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_2018 <- a_2018 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
c_2018 <- c_2018 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
d_2018 <- d_2018 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))
e_2018 <- e_2018 %>% arrange(factor(term, levels = c("(Intercept)","cohort", "month_treat", "east", "cohort:month_treat",
                                                     "cohort:east", "month_treat:east", "cohort:month_treat:east")))

####2018 DiDiD Models Table A12
stargazer(mod_1,mod_3,mod_4,mod_5,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Cohort}", "\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{Cohort:After May}", "\\textsc{Cohort:East}",
                               "\\textsc{After May:East}","\\textsc{After May:East:Cohort}"),
          title = "Survey Year 2018 Difference-in-Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_2018$.lower,a_2018$.upper),cbind(c_2018$.lower,c_2018$.upper),
                           cbind(d_2018$.lower,d_2018$.upper),cbind(e_2018$.lower,e_2018$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = "(Intercept)",
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A12.tex")

####Survey Year 2018 DiD Models (Table A18)####

#Trust FCC 2018
reg_dat <- dat %>% filter(cohort == 1, year_survey == 2018) %>% #Filtering those within relevant cohort and survey year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,year_survey,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2018), weights  = filter(dat,cohort == 1, year_survey == 2018)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_2018_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2018_FCC <- a_DiD_2018_FCC %>% mutate(model = "DiD", survey_year = "2018", subset = "School Environment", dv = "Trust FCC")
DiD_2018_FCC <- a_DiD_2018_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 2018
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2018), weights  = filter(dat,cohort == 1, year_survey == 2018)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_2018_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2018_parliament <- a_DiD_2018_parliament %>% mutate(model = "DiD", survey_year = "2018", subset = "School Environment", dv = "Trust Parliament")
DiD_2018_parliament <- a_DiD_2018_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 2018
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2018), weights  = filter(dat,cohort == 1, year_survey == 2018)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_2018_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2018_federal <- a_DiD_2018_federal %>% mutate(model = "DiD", survey_year = "2018", subset = "School Environment", dv = "Trust Federal")
DiD_2018_federal <- a_DiD_2018_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 2018
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, year_survey == 2018), weights  = filter(dat,cohort == 1, year_survey == 2018)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_2018_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_2018_econsit <- a_DiD_2018_econsit %>% mutate(model = "DiD", survey_year = "2018", subset = "School Environment", dv = "Financial Situation")
DiD_2018_econsit <- a_DiD_2018_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_2018_FCC <- a_DiD_2018_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2018_parliament <- a_DiD_2018_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2018_federal <- a_DiD_2018_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_2018_econsit <- a_DiD_2018_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####2018 DiD Models Table A18
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Survey Year 2018 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_2018_FCC$.lower,a_DiD_2018_FCC$.upper),cbind(a_DiD_2018_parliament$.lower,a_DiD_2018_parliament$.upper),
                           cbind(a_DiD_2018_federal$.lower,a_DiD_2018_federal$.upper),cbind(a_DiD_2018_econsit$.lower,a_DiD_2018_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A18.tex")

####Survey Year 2018 RD Models (Table A6)####

#RD FCC 2018
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2018) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_2018 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", survey_year = "2018", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2018)) #placeholder for table

#RD Parliament 2018
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2018) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_2018 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", survey_year = "2018", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2018))

#RD Federal 2018
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2018) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_2018 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", survey_year = "2018", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2018))


#RD Financial Situation 2018
east <- dat %>% filter(east == 1, cohort == 1, year_survey == 2018) #subsetting data to relevant cohort and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_2018 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", survey_year = "2018", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, year_survey == 2018))


#RD Models 2018 table A6
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Survey Year 2018 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_2018$estimate[1],rd_parliament_2018$estimate[1],rd_federal_2018$estimate[1],rd_econsit_2018$estimate[1]),
          ci.custom = list(cbind(rd_environment_2018$lower,rd_environment_2018$upper), 
                           cbind(rd_parliament_2018$lower,rd_parliament_2018$upper),
                           cbind(rd_federal_2018$lower,rd_federal_2018$upper),
                           cbind(rd_econsit_2018$lower,rd_econsit_2018$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A6.tex")

models_2018 <- rbind(DiDiD_2018,DiDiD_parliament_2018,DiDiD_federal_2018, DiDiD_econsit_2018, DiD_2018_FCC, DiD_2018_federal, DiD_2018_parliament, DiD_2018_econsit,
                     rd_environment_2018,rd_parliament_2018,rd_federal_2018, rd_econsit_2018)

models_full_2018 <- rbind(a_2018,c_2018,d_2018,e_2018,a_DiD_2018_FCC,a_DiD_2018_parliament, a_DiD_2018_federal, a_DiD_2018_econsit)


####Making big plot for survey years (Figure 3)####
survey_year_coefs_full <- rbind(pooled_models_full,models_full_2000,models_full_2002,models_full_2008,models_full_2012,models_full_2018) #putting together coefficients for all models
survey_year_coefs <- rbind(pooled_models,models_2000,models_2002,models_2008,models_2012,models_2018) #putting together coefficients in plot
#saveRDS(survey_year_coefs_full, "survey_year_coefs_full.rds")
#write.csv(survey_year_coefs_full, "survey_year_coefs_full.csv")
#saveRDS(survey_year_coefs, "survey_year_coefs.rds")
#write.csv(survey_year_coefs, "survey_year_coefs.csv") #Code to save data 

survey_year_coefs$survey_year <- factor(survey_year_coefs$survey_year, levels = c("2018", "2012", "2008", "2002", "2000", "Pooled")) #ordering years for plot
survey_year_coefs$dv <- factor(survey_year_coefs$dv, levels = c("Trust FCC", "Trust Federal", "Trust Parliament", "Financial Situation")) #Ordering facets for plot
survey_year_coefs$model <- factor(survey_year_coefs$model, levels = c("DiDiD", "DiD", "RD")) #Ordering Models for plot

####Making Figure 3
all_dv_plot <- ggplot(survey_year_coefs, aes(x = survey_year, y = estimate, color = model, shape = model)) +
  geom_point(size=2,
             position = position_dodge(width = 0.5)) +
  scale_shape_manual(values=c(15, 16, 17)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1,
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#CC79A7","#0072B2", "#D55E00")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~dv) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position="bottom", legend.text = element_text(size=10)) + 
  xlab("Survey Year") + ylab("Estimate") +
  ggtitle("Coefficient Estimates by Survey Year") 
ggsave("Figure_3.pdf", plot= all_dv_plot)


####Year born models####

####1973 DiD Models (Table A29)####

#Trust FCC 1973
reg_dat <- dat %>% filter(cohort == 1, yborn == 1973) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1973), weights  = filter(dat,cohort == 1, yborn == 1973)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1973_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1973_FCC <- a_DiD_1973_FCC %>% mutate(model = "DiD", year_born = "1973", subset = "School Environment", dv = "Trust FCC")
DiD_1973_FCC <- a_DiD_1973_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1973
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1973), weights  = filter(dat,cohort == 1, yborn == 1973)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1973_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1973_parliament <- a_DiD_1973_parliament %>% mutate(model = "DiD", year_born = "1973", subset = "School Environment", dv = "Trust Parliament")
DiD_1973_parliament <- a_DiD_1973_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1973
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1973), weights  = filter(dat,cohort == 1, yborn == 1973)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1973_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1973_federal <- a_DiD_1973_federal %>% mutate(model = "DiD", year_born = "1973", subset = "School Environment", dv = "Trust Federal")
DiD_1973_federal <- a_DiD_1973_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1973
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1973), weights  = filter(dat,cohort == 1, yborn == 1973)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1973_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1973_econsit <- a_DiD_1973_econsit %>% mutate(model = "DiD", year_born = "1973", subset = "School Environment", dv = "Financial Situation")
DiD_1973_econsit <- a_DiD_1973_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1973_FCC <- a_DiD_1973_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1973_parliament <- a_DiD_1973_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1973_federal <- a_DiD_1973_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1973_econsit <- a_DiD_1973_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1973 DiD Models Table A29
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1973 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1973_FCC$.lower,a_DiD_1973_FCC$.upper),cbind(a_DiD_1973_parliament$.lower,a_DiD_1973_parliament$.upper),
                           cbind(a_DiD_1973_federal$.lower,a_DiD_1973_federal$.upper),cbind(a_DiD_1973_econsit$.lower,a_DiD_1973_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A29.tex")

####1973 RD Models (Table A19)####

#RD FCC 1973
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1973) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1973 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1973", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1973))

#RD Parliament 1973
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1973) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1973 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1973", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1973))

#RD Federal 1973
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1973) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1973 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1973", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1973))

#RD Financial Situation 1973
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1973) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1973 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1973", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1973))


#RD 1973 models table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1973 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1973$estimate[1],rd_parliament_1973$estimate[1],rd_federal_1973$estimate[1],rd_econsit_1973$estimate[1]),
          ci.custom = list(cbind(rd_environment_1973$lower,rd_environment_1973$upper), 
                           cbind(rd_parliament_1973$lower,rd_parliament_1973$upper),
                           cbind(rd_federal_1973$lower,rd_federal_1973$upper),
                           cbind(rd_econsit_1973$lower,rd_econsit_1973$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A19.tex")

models_1973 <- rbind(DiD_1973_FCC, DiD_1973_federal, DiD_1973_parliament, DiD_1973_econsit,
                     rd_environment_1973,rd_parliament_1973,rd_federal_1973, rd_econsit_1973)

models_full_1973 <- rbind(a_DiD_1973_FCC,a_DiD_1973_parliament, a_DiD_1973_federal, a_DiD_1973_econsit)

####1974 DiD Models (Table A30)####

#Trust FCC 1974
reg_dat <- dat %>% filter(cohort == 1, yborn == 1974) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1974), weights  = filter(dat,cohort == 1, yborn == 1974)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1974_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1974_FCC <- a_DiD_1974_FCC %>% mutate(model = "DiD", year_born = "1974", subset = "School Environment", dv = "Trust FCC")
DiD_1974_FCC <- a_DiD_1974_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1974
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1974), weights  = filter(dat,cohort == 1, yborn == 1974)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1974_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1974_parliament <- a_DiD_1974_parliament %>% mutate(model = "DiD", year_born = "1974", subset = "School Environment", dv = "Trust Parliament")
DiD_1974_parliament <- a_DiD_1974_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1974
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1974), weights  = filter(dat,cohort == 1, yborn == 1974)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1974_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1974_federal <- a_DiD_1974_federal %>% mutate(model = "DiD", year_born = "1974", subset = "School Environment", dv = "Trust Federal")
DiD_1974_federal <- a_DiD_1974_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1974
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1974), weights  = filter(dat,cohort == 1, yborn == 1974)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1974_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1974_econsit <- a_DiD_1974_econsit %>% mutate(model = "DiD", year_born = "1974", subset = "School Environment", dv = "Financial Situation")
DiD_1974_econsit <- a_DiD_1974_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1974_FCC <- a_DiD_1974_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1974_parliament <- a_DiD_1974_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1974_federal <- a_DiD_1974_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1974_econsit <- a_DiD_1974_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1974 DiD Models Table A30
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1974 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1974_FCC$.lower,a_DiD_1974_FCC$.upper),cbind(a_DiD_1974_parliament$.lower,a_DiD_1974_parliament$.upper),
                           cbind(a_DiD_1974_federal$.lower,a_DiD_1974_federal$.upper),cbind(a_DiD_1974_econsit$.lower,a_DiD_1974_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A30.tex")

####1974 RD Models (Table A20)####

#RD FCC 1974
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1974) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1974 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1974", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1974))

#RD Parliament 1974
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1974) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1974 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1974", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1974))

#RD Federal 1974
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1974) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1974 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1974", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1974))

#RD Financial Situation 1974
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1974) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1974 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1974", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1974))


#RD 1974 models table A20
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1974 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1974$estimate[1],rd_parliament_1974$estimate[1],rd_federal_1974$estimate[1],rd_econsit_1974$estimate[1]),
          ci.custom = list(cbind(rd_environment_1974$lower,rd_environment_1974$upper), 
                           cbind(rd_parliament_1974$lower,rd_parliament_1974$upper),
                           cbind(rd_federal_1974$lower,rd_federal_1974$upper),
                           cbind(rd_econsit_1974$lower,rd_econsit_1974$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A20.tex")

models_1974 <- rbind(DiD_1974_FCC, DiD_1974_federal, DiD_1974_parliament, DiD_1974_econsit,
                     rd_environment_1974,rd_parliament_1974,rd_federal_1974, rd_econsit_1974)

models_full_1974 <- rbind(a_DiD_1974_FCC,a_DiD_1974_parliament, a_DiD_1974_federal, a_DiD_1974_econsit)

####1975 DiD Models (Table A31)####

#DiD Models

#Trust FCC 1975
reg_dat <- dat %>% filter(cohort == 1, yborn == 1975) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1975), weights  = filter(dat,cohort == 1, yborn == 1975)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1975_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1975_FCC <- a_DiD_1975_FCC %>% mutate(model = "DiD", year_born = "1975", subset = "School Environment", dv = "Trust FCC")
DiD_1975_FCC <- a_DiD_1975_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1975
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1975), weights  = filter(dat,cohort == 1, yborn == 1975)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1975_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1975_parliament <- a_DiD_1975_parliament %>% mutate(model = "DiD", year_born = "1975", subset = "School Environment", dv = "Trust Parliament")
DiD_1975_parliament <- a_DiD_1975_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1975
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1975), weights  = filter(dat,cohort == 1, yborn == 1975)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1975_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1975_federal <- a_DiD_1975_federal %>% mutate(model = "DiD", year_born = "1975", subset = "School Environment", dv = "Trust Federal")
DiD_1975_federal <- a_DiD_1975_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1975
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1975), weights  = filter(dat,cohort == 1, yborn == 1975)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1975_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1975_econsit <- a_DiD_1975_econsit %>% mutate(model = "DiD", year_born = "1975", subset = "School Environment", dv = "Financial Situation")
DiD_1975_econsit <- a_DiD_1975_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1975_FCC <- a_DiD_1975_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1975_parliament <- a_DiD_1975_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1975_federal <- a_DiD_1975_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1975_econsit <- a_DiD_1975_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1975 DiD Models Table A31
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1975 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1975_FCC$.lower,a_DiD_1975_FCC$.upper),cbind(a_DiD_1975_parliament$.lower,a_DiD_1975_parliament$.upper),
                           cbind(a_DiD_1975_federal$.lower,a_DiD_1975_federal$.upper),cbind(a_DiD_1975_econsit$.lower,a_DiD_1975_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A31.tex")

####1975 RD Models (Table A21)####

#RD FCC 1975
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1975) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1975 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1975", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1975))

#RD Parliament 1975
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1975) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1975 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1975", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1975))

#RD Federal 1975
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1975) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1975 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1975", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1975))

#RD Financial Situation 1975
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1975) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1975 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1975", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1975))


#RD 1975 models table A21
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1975 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1975$estimate[1],rd_parliament_1975$estimate[1],rd_federal_1975$estimate[1],rd_econsit_1975$estimate[1]),
          ci.custom = list(cbind(rd_environment_1975$lower,rd_environment_1975$upper), 
                           cbind(rd_parliament_1975$lower,rd_parliament_1975$upper),
                           cbind(rd_federal_1975$lower,rd_federal_1975$upper),
                           cbind(rd_econsit_1975$lower,rd_econsit_1975$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A21.tex")

models_1975 <- rbind(DiD_1975_FCC, DiD_1975_federal, DiD_1975_parliament, DiD_1975_econsit,
                     rd_environment_1975,rd_parliament_1975,rd_federal_1975, rd_econsit_1975)

models_full_1975 <- rbind(a_DiD_1975_FCC,a_DiD_1975_parliament, a_DiD_1975_federal, a_DiD_1975_econsit)

####1976 DiD Models (Table A32)####

#DiD Models

#Trust FCC 1976
reg_dat <- dat %>% filter(cohort == 1, yborn == 1976) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1976), weights  = filter(dat,cohort == 1, yborn == 1976)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1976_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1976_FCC <- a_DiD_1976_FCC %>% mutate(model = "DiD", year_born = "1976", subset = "School Environment", dv = "Trust FCC")
DiD_1976_FCC <- a_DiD_1976_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1976
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1976), weights  = filter(dat,cohort == 1, yborn == 1976)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1976_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1976_parliament <- a_DiD_1976_parliament %>% mutate(model = "DiD", year_born = "1976", subset = "School Environment", dv = "Trust Parliament")
DiD_1976_parliament <- a_DiD_1976_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1976
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1976), weights  = filter(dat,cohort == 1, yborn == 1976)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1976_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1976_federal <- a_DiD_1976_federal %>% mutate(model = "DiD", year_born = "1976", subset = "School Environment", dv = "Trust Federal")
DiD_1976_federal <- a_DiD_1976_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1976
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1976), weights  = filter(dat,cohort == 1, yborn == 1976)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1976_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1976_econsit <- a_DiD_1976_econsit %>% mutate(model = "DiD", year_born = "1976", subset = "School Environment", dv = "Financial Situation")
DiD_1976_econsit <- a_DiD_1976_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1976_FCC <- a_DiD_1976_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1976_parliament <- a_DiD_1976_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1976_federal <- a_DiD_1976_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1976_econsit <- a_DiD_1976_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1976 DiD Models Table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1976 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1976_FCC$.lower,a_DiD_1976_FCC$.upper),cbind(a_DiD_1976_parliament$.lower,a_DiD_1976_parliament$.upper),
                           cbind(a_DiD_1976_federal$.lower,a_DiD_1976_federal$.upper),cbind(a_DiD_1976_econsit$.lower,a_DiD_1976_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A32.tex")

####1976 RD Models (Table A22)####

#RD FCC 1976
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1976) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1976 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1976", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1976))

#RD Parliament 1976
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1976) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1976 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1976", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1976))

#RD Federal 1976
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1976) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1976 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1976", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1976))

#RD Financial Situation 1976
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1976) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1976 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1976", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1976))


#RD 1976 models table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1976 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1976$estimate[1],rd_parliament_1976$estimate[1],rd_federal_1976$estimate[1],rd_econsit_1976$estimate[1]),
          ci.custom = list(cbind(rd_environment_1976$lower,rd_environment_1976$upper), 
                           cbind(rd_parliament_1976$lower,rd_parliament_1976$upper),
                           cbind(rd_federal_1976$lower,rd_federal_1976$upper),
                           cbind(rd_econsit_1976$lower,rd_econsit_1976$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A22.tex")

models_1976 <- rbind(DiD_1976_FCC, DiD_1976_federal, DiD_1976_parliament, DiD_1976_econsit,
                     rd_environment_1976,rd_parliament_1976,rd_federal_1976, rd_econsit_1976)

models_full_1976 <- rbind(a_DiD_1976_FCC,a_DiD_1976_parliament, a_DiD_1976_federal, a_DiD_1976_econsit)

####1977 DiD Models (Table A33)####

#DiD Models

#Trust FCC 1977
reg_dat <- dat %>% filter(cohort == 1, yborn == 1977) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1977), weights  = filter(dat,cohort == 1, yborn == 1977)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1977_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1977_FCC <- a_DiD_1977_FCC %>% mutate(model = "DiD", year_born = "1977", subset = "School Environment", dv = "Trust FCC")
DiD_1977_FCC <- a_DiD_1977_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1977
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1977), weights  = filter(dat,cohort == 1, yborn == 1977)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1977_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1977_parliament <- a_DiD_1977_parliament %>% mutate(model = "DiD", year_born = "1977", subset = "School Environment", dv = "Trust Parliament")
DiD_1977_parliament <- a_DiD_1977_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1977
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1977), weights  = filter(dat,cohort == 1, yborn == 1977)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1977_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1977_federal <- a_DiD_1977_federal %>% mutate(model = "DiD", year_born = "1977", subset = "School Environment", dv = "Trust Federal")
DiD_1977_federal <- a_DiD_1977_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1977
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1977), weights  = filter(dat,cohort == 1, yborn == 1977)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1977_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1977_econsit <- a_DiD_1977_econsit %>% mutate(model = "DiD", year_born = "1977", subset = "School Environment", dv = "Financial Situation")
DiD_1977_econsit <- a_DiD_1977_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1977_FCC <- a_DiD_1977_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1977_parliament <- a_DiD_1977_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1977_federal <- a_DiD_1977_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1977_econsit <- a_DiD_1977_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1977 DiD Models Table A33
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1977 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1977_FCC$.lower,a_DiD_1977_FCC$.upper),cbind(a_DiD_1977_parliament$.lower,a_DiD_1977_parliament$.upper),
                           cbind(a_DiD_1977_federal$.lower,a_DiD_1977_federal$.upper),cbind(a_DiD_1977_econsit$.lower,a_DiD_1977_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A33.tex")

####1977 RD Models (Table A23)####

#RD FCC 1977
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1977) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1977 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1977", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1977))

#RD Parliament 1977
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1977) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1977 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1977", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1977))

#RD Federal 1977
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1977) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1977 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1977", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1977))

#RD Financial Situation 1977
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1977) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1977 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1977", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1977))


#RD 1977 models table A23
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1977 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1977$estimate[1],rd_parliament_1977$estimate[1],rd_federal_1977$estimate[1],rd_econsit_1977$estimate[1]),
          ci.custom = list(cbind(rd_environment_1977$lower,rd_environment_1977$upper), 
                           cbind(rd_parliament_1977$lower,rd_parliament_1977$upper),
                           cbind(rd_federal_1977$lower,rd_federal_1977$upper),
                           cbind(rd_econsit_1977$lower,rd_econsit_1977$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A23.tex")

models_1977 <- rbind(DiD_1977_FCC, DiD_1977_federal, DiD_1977_parliament, DiD_1977_econsit,
                     rd_environment_1977,rd_parliament_1977,rd_federal_1977, rd_econsit_1977)

models_full_1977 <- rbind(a_DiD_1977_FCC,a_DiD_1977_parliament, a_DiD_1977_federal, a_DiD_1977_econsit)

####1978 DiD Models (Table A34)####

#Trust FCC 1978
reg_dat <- dat %>% filter(cohort == 1, yborn == 1978) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1978), weights  = filter(dat,cohort == 1, yborn == 1978)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1978_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1978_FCC <- a_DiD_1978_FCC %>% mutate(model = "DiD", year_born = "1978", subset = "School Environment", dv = "Trust FCC")
DiD_1978_FCC <- a_DiD_1978_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1978
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1978), weights  = filter(dat,cohort == 1, yborn == 1978)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1978_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1978_parliament <- a_DiD_1978_parliament %>% mutate(model = "DiD", year_born = "1978", subset = "School Environment", dv = "Trust Parliament")
DiD_1978_parliament <- a_DiD_1978_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1978
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1978), weights  = filter(dat,cohort == 1, yborn == 1978)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1978_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1978_federal <- a_DiD_1978_federal %>% mutate(model = "DiD", year_born = "1978", subset = "School Environment", dv = "Trust Federal")
DiD_1978_federal <- a_DiD_1978_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1978
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1978), weights  = filter(dat,cohort == 1, yborn == 1978)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1978_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1978_econsit <- a_DiD_1978_econsit %>% mutate(model = "DiD", year_born = "1978", subset = "School Environment", dv = "Financial Situation")
DiD_1978_econsit <- a_DiD_1978_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1978_FCC <- a_DiD_1978_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1978_parliament <- a_DiD_1978_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1978_federal <- a_DiD_1978_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1978_econsit <- a_DiD_1978_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1978 DiD Models Table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1978 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1978_FCC$.lower,a_DiD_1978_FCC$.upper),cbind(a_DiD_1978_parliament$.lower,a_DiD_1978_parliament$.upper),
                           cbind(a_DiD_1978_federal$.lower,a_DiD_1978_federal$.upper),cbind(a_DiD_1978_econsit$.lower,a_DiD_1978_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A34.tex")

####1978 RD Models (Table A24)####

#RD FCC 1978
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1978) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1978 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1978", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1978))

#RD Parliament 1978
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1978) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1978 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1978", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1978))

#RD Federal 1978
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1978) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1978 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1978", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1978))

#RD Financial Situation 1978
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1978) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1978 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1978", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1978))


#RD 1978 models table A24
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1978 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1978$estimate[1],rd_parliament_1978$estimate[1],rd_federal_1978$estimate[1],rd_econsit_1978$estimate[1]),
          ci.custom = list(cbind(rd_environment_1978$lower,rd_environment_1978$upper), 
                           cbind(rd_parliament_1978$lower,rd_parliament_1978$upper),
                           cbind(rd_federal_1978$lower,rd_federal_1978$upper),
                           cbind(rd_econsit_1978$lower,rd_econsit_1978$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A24.tex")

models_1978 <- rbind(DiD_1978_FCC, DiD_1978_federal, DiD_1978_parliament, DiD_1978_econsit,
                     rd_environment_1978,rd_parliament_1978,rd_federal_1978, rd_econsit_1978)

models_full_1978 <- rbind(a_DiD_1978_FCC,a_DiD_1978_parliament, a_DiD_1978_federal, a_DiD_1978_econsit)

####1979 DiD Models (Table A35)####

#Trust FCC 1979
reg_dat <- dat %>% filter(cohort == 1, yborn == 1979) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1979), weights  = filter(dat,cohort == 1, yborn == 1979)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1979_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1979_FCC <- a_DiD_1979_FCC %>% mutate(model = "DiD", year_born = "1979", subset = "School Environment", dv = "Trust FCC")
DiD_1979_FCC <- a_DiD_1979_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1979
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1979), weights  = filter(dat,cohort == 1, yborn == 1979)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1979_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1979_parliament <- a_DiD_1979_parliament %>% mutate(model = "DiD", year_born = "1979", subset = "School Environment", dv = "Trust Parliament")
DiD_1979_parliament <- a_DiD_1979_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1979
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1979), weights  = filter(dat,cohort == 1, yborn == 1979)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1979_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1979_federal <- a_DiD_1979_federal %>% mutate(model = "DiD", year_born = "1979", subset = "School Environment", dv = "Trust Federal")
DiD_1979_federal <- a_DiD_1979_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1979
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1979), weights  = filter(dat,cohort == 1, yborn == 1979)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1979_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1979_econsit <- a_DiD_1979_econsit %>% mutate(model = "DiD", year_born = "1979", subset = "School Environment", dv = "Financial Situation")
DiD_1979_econsit <- a_DiD_1979_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1979_FCC <- a_DiD_1979_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1979_parliament <- a_DiD_1979_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1979_federal <- a_DiD_1979_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1979_econsit <- a_DiD_1979_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1979 DiD Models Table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1979 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1979_FCC$.lower,a_DiD_1979_FCC$.upper),cbind(a_DiD_1979_parliament$.lower,a_DiD_1979_parliament$.upper),
                           cbind(a_DiD_1979_federal$.lower,a_DiD_1979_federal$.upper),cbind(a_DiD_1979_econsit$.lower,a_DiD_1979_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A35.tex")

####1979 RD Models (Table A25)####

#RD FCC 1979
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1979) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1979 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1979", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1979))

#RD Parliament 1979
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1979) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1979 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1979", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1979))

#RD Federal 1979
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1979) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1979 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1979", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1979))

#RD Financial Situation 1979
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1979) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1979 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1979", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1979))


#RD 1979 models table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1979 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1979$estimate[1],rd_parliament_1979$estimate[1],rd_federal_1979$estimate[1],rd_econsit_1979$estimate[1]),
          ci.custom = list(cbind(rd_environment_1979$lower,rd_environment_1979$upper), 
                           cbind(rd_parliament_1979$lower,rd_parliament_1979$upper),
                           cbind(rd_federal_1979$lower,rd_federal_1979$upper),
                           cbind(rd_econsit_1979$lower,rd_econsit_1979$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A25.tex")

models_1979 <- rbind(DiD_1979_FCC, DiD_1979_federal, DiD_1979_parliament, DiD_1979_econsit,
                     rd_environment_1979,rd_parliament_1979,rd_federal_1979, rd_econsit_1979)

models_full_1979 <- rbind(a_DiD_1979_FCC,a_DiD_1979_parliament, a_DiD_1979_federal, a_DiD_1979_econsit)

####1980 DiD Models (Table A36)####

#Trust FCC 1980
reg_dat <- dat %>% filter(cohort == 1, yborn == 1980) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1980), weights  = filter(dat,cohort == 1, yborn == 1980)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1980_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1980_FCC <- a_DiD_1980_FCC %>% mutate(model = "DiD", year_born = "1980", subset = "School Environment", dv = "Trust FCC")
DiD_1980_FCC <- a_DiD_1980_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1980
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1980), weights  = filter(dat,cohort == 1, yborn == 1980)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1980_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1980_parliament <- a_DiD_1980_parliament %>% mutate(model = "DiD", year_born = "1980", subset = "School Environment", dv = "Trust Parliament")
DiD_1980_parliament <- a_DiD_1980_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1980
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1980), weights  = filter(dat,cohort == 1, yborn == 1980)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1980_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1980_federal <- a_DiD_1980_federal %>% mutate(model = "DiD", year_born = "1980", subset = "School Environment", dv = "Trust Federal")
DiD_1980_federal <- a_DiD_1980_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1980
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1980), weights  = filter(dat,cohort == 1, yborn == 1980)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1980_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1980_econsit <- a_DiD_1980_econsit %>% mutate(model = "DiD", year_born = "1980", subset = "School Environment", dv = "Financial Situation")
DiD_1980_econsit <- a_DiD_1980_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1980_FCC <- a_DiD_1980_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1980_parliament <- a_DiD_1980_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1980_federal <- a_DiD_1980_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1980_econsit <- a_DiD_1980_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1980 DiD Models Table A36
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1980 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1980_FCC$.lower,a_DiD_1980_FCC$.upper),cbind(a_DiD_1980_parliament$.lower,a_DiD_1980_parliament$.upper),
                           cbind(a_DiD_1980_federal$.lower,a_DiD_1980_federal$.upper),cbind(a_DiD_1980_econsit$.lower,a_DiD_1980_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A36.tex")

####1980 RD Models (Table A26)####

#RD FCC 1980
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1980) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1980 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1980", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1980))

#RD Parliament 1980
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1980) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1980 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1980", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1980))

#RD Federal 1980
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1980) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1980 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1980", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1980))

#RD Financial Situation 1980
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1980) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1980 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1980", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1980))


#RD 1980 models table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1980 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1980$estimate[1],rd_parliament_1980$estimate[1],rd_federal_1980$estimate[1],rd_econsit_1980$estimate[1]),
          ci.custom = list(cbind(rd_environment_1980$lower,rd_environment_1980$upper), 
                           cbind(rd_parliament_1980$lower,rd_parliament_1980$upper),
                           cbind(rd_federal_1980$lower,rd_federal_1980$upper),
                           cbind(rd_econsit_1980$lower,rd_econsit_1980$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A26.tex")

models_1980 <- rbind(DiD_1980_FCC, DiD_1980_federal, DiD_1980_parliament, DiD_1980_econsit,
                     rd_environment_1980,rd_parliament_1980,rd_federal_1980, rd_econsit_1980)

models_full_1980 <- rbind(a_DiD_1980_FCC,a_DiD_1980_parliament, a_DiD_1980_federal, a_DiD_1980_econsit)

####1981 DiD Models (Table A37)####

#Trust FCC 1981
reg_dat <- dat %>% filter(cohort == 1, yborn == 1981) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1981), weights  = filter(dat,cohort == 1, yborn == 1981)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1981_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1981_FCC <- a_DiD_1981_FCC %>% mutate(model = "DiD", year_born = "1981", subset = "School Environment", dv = "Trust FCC")
DiD_1981_FCC <- a_DiD_1981_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1981
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1981), weights  = filter(dat,cohort == 1, yborn == 1981)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1981_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1981_parliament <- a_DiD_1981_parliament %>% mutate(model = "DiD", year_born = "1981", subset = "School Environment", dv = "Trust Parliament")
DiD_1981_parliament <- a_DiD_1981_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1981
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1981), weights  = filter(dat,cohort == 1, yborn == 1981)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1981_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1981_federal <- a_DiD_1981_federal %>% mutate(model = "DiD", year_born = "1981", subset = "School Environment", dv = "Trust Federal")
DiD_1981_federal <- a_DiD_1981_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1981
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1981), weights  = filter(dat,cohort == 1, yborn == 1981)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1981_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1981_econsit <- a_DiD_1981_econsit %>% mutate(model = "DiD", year_born = "1981", subset = "School Environment", dv = "Financial Situation")
DiD_1981_econsit <- a_DiD_1981_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1981_FCC <- a_DiD_1981_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1981_parliament <- a_DiD_1981_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1981_federal <- a_DiD_1981_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1981_econsit <- a_DiD_1981_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1981 DiD Models Table A37
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1981 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1981_FCC$.lower,a_DiD_1981_FCC$.upper),cbind(a_DiD_1981_parliament$.lower,a_DiD_1981_parliament$.upper),
                           cbind(a_DiD_1981_federal$.lower,a_DiD_1981_federal$.upper),cbind(a_DiD_1981_econsit$.lower,a_DiD_1981_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A37.tex")

####1981 RD Models (Table A27)####

#RD FCC 1981
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1981) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1981 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1981", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1981))

#RD Parliament 1981
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1981) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1981 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1981", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1981))

#RD Federal 1981
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1981) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1981 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1981", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1981))

#RD Financial Situation 1981
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1981) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1981 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1981", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1981))


#RD 1981 models table
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1981 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1981$estimate[1],rd_parliament_1981$estimate[1],rd_federal_1981$estimate[1],rd_econsit_1981$estimate[1]),
          ci.custom = list(cbind(rd_environment_1981$lower,rd_environment_1981$upper), 
                           cbind(rd_parliament_1981$lower,rd_parliament_1981$upper),
                           cbind(rd_federal_1981$lower,rd_federal_1981$upper),
                           cbind(rd_econsit_1981$lower,rd_econsit_1981$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A27.tex")

models_1981 <- rbind(DiD_1981_FCC, DiD_1981_federal, DiD_1981_parliament, DiD_1981_econsit,
                     rd_environment_1981,rd_parliament_1981,rd_federal_1981, rd_econsit_1981)

models_full_1981 <- rbind(a_DiD_1981_FCC,a_DiD_1981_parliament, a_DiD_1981_federal, a_DiD_1981_econsit)

####1982 DiD Models (Table A38)####

#Trust FCC 1982
reg_dat <- dat %>% filter(cohort == 1, yborn == 1982) %>% #Filtering those within relevant cohort and birth year
  select(trust_FCC, trust_parliament, trust_federal, econ_sit, cohort,month_treat,east,yborn,mborn,weights) #Selecting relevant independent and dependant variables
mod_1 <- felm(trust_FCC~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1982), weights  = filter(dat,cohort == 1, yborn == 1982)$weights)
reg <- trust_FCC~month_treat*east |0|0|0
a_DiD_1982_FCC <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1982_FCC <- a_DiD_1982_FCC %>% mutate(model = "DiD", year_born = "1982", subset = "School Environment", dv = "Trust FCC")
DiD_1982_FCC <- a_DiD_1982_FCC %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust parliament 1982
mod_2 <- felm(trust_parliament~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1982), weights  = filter(dat,cohort == 1, yborn == 1982)$weights)
reg <- trust_parliament~month_treat*east |0|0|0
a_DiD_1982_parliament <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1982_parliament <- a_DiD_1982_parliament %>% mutate(model = "DiD", year_born = "1982", subset = "School Environment", dv = "Trust Parliament")
DiD_1982_parliament <- a_DiD_1982_parliament %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Trust federal 1982
mod_3 <- felm(trust_federal~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1982), weights  = filter(dat,cohort == 1, yborn == 1982)$weights)
reg <- trust_federal~month_treat*east |0|0|0
a_DiD_1982_federal <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1982_federal <- a_DiD_1982_federal %>% mutate(model = "DiD", year_born = "1982", subset = "School Environment", dv = "Trust Federal")
DiD_1982_federal <- a_DiD_1982_federal %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

#Financial Sitution 1982
mod_4 <- felm(econ_sit~month_treat*east |0|0|0, data  = filter(dat,cohort == 1, yborn == 1982), weights  = filter(dat,cohort == 1, yborn == 1982)$weights)
reg <- econ_sit~month_treat*east |0|0|0
a_DiD_1982_econsit <- btraps_blocks_weights_DiD(reg,reg_dat)
a_DiD_1982_econsit <- a_DiD_1982_econsit %>% mutate(model = "DiD", year_born = "1982", subset = "School Environment", dv = "Financial Situation")
DiD_1982_econsit <- a_DiD_1982_econsit %>% filter(term=="month_treat:east") %>%
  select(-c(upper_alpha,lower_alpha, term, std.error)) %>% rename(lower = .lower, upper = .upper)

####reordering terms for regression table
a_DiD_1982_FCC <- a_DiD_1982_FCC %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1982_parliament <- a_DiD_1982_parliament %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1982_federal <- a_DiD_1982_federal %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))
a_DiD_1982_econsit <- a_DiD_1982_econsit %>% arrange(factor(term, levels = c("(Intercept)","month_treat", "east", "month_treat:east")))

####1982 DiD Models Table A38
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{After May}", "\\textsc{East}",
                               "\\textsc{After May:East}"),
          title = "Birth Year 1982 Difference-in-Differences Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          ci.custom = list(cbind(a_DiD_1982_FCC$.lower,a_DiD_1982_FCC$.upper),cbind(a_DiD_1982_parliament$.lower,a_DiD_1982_parliament$.upper),
                           cbind(a_DiD_1982_federal$.lower,a_DiD_1982_federal$.upper),cbind(a_DiD_1982_econsit$.lower,a_DiD_1982_econsit$.upper)),
          add.lines = list(c("Survey Weights","Yes", "Yes", "Yes", "Yes")),
          keep.stat = c("n", "rsq"),
          #omit = c("Constant"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals calculated from 1000 block bootstrap replications are in parentheses}",
          out = "Table_A38.tex")

####1982 RD Models (Table A28)####

#RD FCC 1982
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1982) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_FCC #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_FCC <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_environment_1982 <- tibble(estimate = tmp_FCC$obs.stat, lower = tmp_FCC$interf.ci[1], upper = tmp_FCC$interf.ci[2],
                              model = "RD", year_born = "1982", subset = "School Environment", dv = "Trust FCC")
mod_1 <- lm(trust_FCC~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1982))

#RD Parliament 1982
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1982) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_parliament #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_parliament <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_parliament_1982 <- tibble(estimate = tmp_parliament$obs.stat, lower = tmp_parliament$interf.ci[1], upper = tmp_parliament$interf.ci[2],
                             model = "RD", year_born = "1982", subset = "School Environment", dv = "Trust Parliament")
mod_2 <- lm(trust_parliament~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1982))

#RD Federal 1982
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1982) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$trust_federal #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_federal <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_federal_1982 <- tibble(estimate = tmp_federal$obs.stat, lower = tmp_federal$interf.ci[1], upper = tmp_federal$interf.ci[2],
                          model = "RD", year_born = "1982", subset = "School Environment", dv = "Trust Federal")
mod_3 <- lm(trust_federal~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1982))

#RD Financial Situation 1982
east <- dat %>% filter(east == 1, cohort == 1, yborn == 1982) #subsetting data to relevant birth year and east germany
R <- as.numeric(east$mborn - 6) #centering mborn variable at 0
Y <- east$econ_sit #dependent variable
D <- as.numeric(R>=0) #treatment when cutoff > 0
tmp_econsit <- rdrandinf(Y,R,wl=-5,wr=7,interfci=.1, p =1) #randomization inference with 90% confidence intervals
rd_econsit_1982 <- tibble(estimate = tmp_econsit$obs.stat, lower = tmp_econsit$interf.ci[1], upper = tmp_econsit$interf.ci[2],
                          model = "RD", year_born = "1982", subset = "School Environment", dv = "Financial Situation")
mod_4 <- lm(econ_sit~month_treat, data = filter(dat, east == 1, cohort == 1, yborn == 1982))


#RD 1982 models table A28
stargazer(mod_1,mod_2,mod_3,mod_4,
          style = "default", no.space = TRUE,
          covariate.labels = c("\\textsc{Estimate}"),
          title = "Birth Year 1982 RD Models",
          dep.var.labels  = c("Trust FCC", "Trust Parliament", "Trust Federal", "Financial Situation"),
          coef = list(rd_environment_1982$estimate[1],rd_parliament_1982$estimate[1],rd_federal_1982$estimate[1],rd_econsit_1982$estimate[1]),
          ci.custom = list(cbind(rd_environment_1982$lower,rd_environment_1982$upper), 
                           cbind(rd_parliament_1982$lower,rd_parliament_1982$upper),
                           cbind(rd_federal_1982$lower,rd_federal_1982$upper),
                           cbind(rd_econsit_1982$lower,rd_econsit_1982$upper)),
          add.lines = list(c("N",tmp_FCC$sumstats[1,1]+tmp_FCC$sumstats[1,2], tmp_parliament$sumstats[1,1]+tmp_parliament$sumstats[1,2], 
                             tmp_federal$sumstats[1,1]+tmp_federal$sumstats[1,2], tmp_econsit$sumstats[1,1]+tmp_econsit$sumstats[1,2])),
          omit.table.layout = "s",
          omit = c("month_treat"),
          star.cutoffs = NA,
          notes = "\\parbox[t]{0.8\\linewidth}{90 percent confidence intervals under interference calculated from the rdlocrand package in R are in parentheses.}",
          out = "Table_A28.tex")

models_1982 <- rbind(DiD_1982_FCC, DiD_1982_federal, DiD_1982_parliament, DiD_1982_econsit,
                     rd_environment_1982,rd_parliament_1982,rd_federal_1982, rd_econsit_1982)

models_full_1982 <- rbind(a_DiD_1982_FCC,a_DiD_1982_parliament, a_DiD_1982_federal, a_DiD_1982_econsit)

####Making big plot for birth years (Figure 2)####
year_born_coefs_full <- rbind(models_full_1973, models_full_1974, models_full_1975, models_full_1976,
                              models_full_1977, models_full_1978, models_full_1979, models_full_1980,
                              models_full_1981, models_full_1982)
year_born_coefs <- rbind(models_1973, models_1974, models_1975, models_1976,
                         models_1977, models_1978, models_1979, models_1980,
                         models_1981, models_1982)
#saveRDS(year_born_coefs_full, "year_born_coefs_full.rds")
#write.csv(year_born_coefs_full, "year_born_coefs_full.rds")
#saveRDS(year_born_coefs, "year_born_coefs.rds")
#write.csv(year_born_coefs, "year_born_coefs.rds")

year_born_coefs$year_born <- factor(year_born_coefs$year_born, levels = c("1982", "1981", "1980", "1979", "1978", "1977",
                                                                          "1976", "1975", "1974","1973")) #ordering years for plot
year_born_coefs$dv <- factor(year_born_coefs$dv, levels = c("Trust FCC", "Trust Federal", "Trust Parliament", "Financial Situation")) #Ordering facets for plot
year_born_coefs$model <- factor(year_born_coefs$model, levels = c("DiD", "RD")) #Ordering Models for plot


all_yb_plot <- ggplot(year_born_coefs, aes(x = year_born, y = estimate, color = model, shape = model)) +
  geom_point(size=2,
             position = position_dodge(width = 0.5)) +
  scale_shape_manual(values=c(16, 17)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1,
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~dv) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position="bottom", legend.text = element_text(size=10)) + 
  xlab("Birth Year") + ylab("Estimate") +
  ggtitle("Coefficient Estimates by Birth Year") 
ggsave("Figure_2.pdf", plot= all_yb_plot)


####Balance Plot (Figure A1)####
balplot_dat <- dat %>% filter(cohort ==1,father_school > 0, mother_school > 0) #removing non-answers from father_school and mother_school
a <- tidy(t.test(filter(balplot_dat, cohort ==1, month_treat == 0)$sex, 
                 filter(balplot_dat, cohort ==1, month_treat == 1)$sex)) %>% 
  mutate(Variable = "Sex", Model = "School Environment")
b <- tidy(t.test(filter(balplot_dat, cohort ==1, month_treat == 0)$father_school, 
                 filter(balplot_dat, cohort ==1, month_treat == 1)$father_school)) %>% 
  mutate(Variable = "Father Schooling", Model = "School Environment")
c <- tidy(t.test(filter(balplot_dat, cohort ==1, month_treat == 0)$mother_school, 
                 filter(balplot_dat, cohort ==1, month_treat == 1)$mother_school)) %>% 
  mutate(Variable = "Mother Schooling", Model = "School Environment")
d <- tidy(t.test(filter(balplot_dat, cohort ==1, month_treat == 0)$yborn, 
                 filter(balplot_dat, cohort ==1, month_treat == 1)$yborn)) %>% 
  mutate(Variable = "Year Born", Model = "School Environment")
balplot_environment_did <- rbind(a,b,c,d)

p1 <- ggplot(balplot_environment_did, aes(Variable, estimate)) + geom_point() +
  ggtitle("Balance on Pre-Treatment Covariates (Born 1973 - 1982)") + ylab("Difference in Means") + xlab("Variable") +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width = .3) +
  coord_flip(ylim = c(-1,1))
ggsave("Figure_A1.pdf", plot= p1)


####Balance plot for all data (Figure A2)####
balplot2_dat <- dat %>% filter(father_school > 0, mother_school > 0) #removing non-answers from father_school and mother_school
a <- tidy(t.test(filter(balplot2_dat, yborn > 1949, month_treat == 0)$sex, 
                 filter(balplot2_dat, yborn > 1949, month_treat == 1)$sex)) %>% 
  mutate(Variable = "Sex")
b <- tidy(t.test(filter(balplot2_dat, yborn > 1949, month_treat == 0)$father_school, 
                 filter(balplot2_dat, yborn > 1949, month_treat == 1)$father_school)) %>% 
  mutate(Variable = "Father Schooling")
c <- tidy(t.test(filter(balplot2_dat, yborn > 1949, month_treat == 0)$mother_school, 
                 filter(balplot2_dat, yborn > 1949, month_treat == 1)$mother_school)) %>% 
  mutate(Variable = "Mother Schooling")
balplot_all_did <- rbind(a,b,c)

p2 <- ggplot(balplot_all_did, aes(Variable, estimate)) + geom_point() +
  ggtitle("Balance on Pre-Treatment Covariates (All Birth Years)") + ylab("Difference in Means") + xlab("Variable") +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width = .3) +
  coord_flip(ylim = c(-1,1))
ggsave("Figure_A2.pdf", plot= p2)


####McCrary Test (Figure A3)####

mccrary <- DCdensity(east %>% filter(cohort == 1) %>% pull(mborn), bw = 6, cutpoint = 6, ext.out=T)
mccrary <- as.data.frame(mccrary$data) 
mccrary <- mccrary %>% filter(cellval > 0)
mccrary$cellmp <- round(mccrary$cellmp)
mccrary$cutoff <- ifelse(mccrary$cellmp > 5, "yes", "no")

McCraryPlot <- ggplot(mccrary, aes(x = as.factor(cellmp), y = cellval, group = factor(cutoff))) + 
  geom_smooth(method = "lm") + geom_point() + theme_bw() + xlab("Birth Month") +
  ylab("Density Estimate") +  ggtitle("McCrary Test") +
  theme(legend.position="none", strip.text = element_text(size=25), 
        text = element_text(size=13), plot.title = element_text(hjust = 0.5), 
        legend.text=element_text(size=18), panel.border = element_blank(),
        axis.line = element_line(colour = "black")) + 
  geom_vline(xintercept = 6, linetype = "dashed")
ggsave("Figure_A3.pdf", plot= McCraryPlot)


####RD Power Test (Figure A4)####
east <- dat %>% filter(east == 1, cohort == 1)
R <- as.numeric(east$mborn - 6)
Y <- east$trust_federal
Z <- cbind(Y,R)
pdf(file = "Figure_A4.pdf",   # The directory you want to save the file in
    width = 9.15, # The width of the plot in inches
    height = 6.5)
aux = rdpower(data=Z,masspoints="off",stdvars="on", h =6, plot = T,
              graph.range = c(-.3,.3), p=1, alpha = 0.1) 
dev.off()

####Parallel Trends Plot (Figure 1)####
ptrends_dat <- dat %>% filter(yborn > 1949, yborn < 1972)
ptrends_dat$east <- ifelse(ptrends_dat$east == 1, "East", "West")
ptrends_dat_FCC <- ptrends_dat %>% 
  group_by(east, yborn) %>%
  group_modify(~tidy(felm(trust_FCC~month_treat*mborn|0|0|0, data = ., weights = .$weights))) %>%
  filter(term == "month_treat") %>% mutate(model = "Trust FCC")
ptrends_dat_federal <- ptrends_dat %>% 
  group_by(east, yborn) %>%
  group_modify(~tidy(felm(trust_federal~month_treat*mborn|0|0|0, data = ., weights = .$weights))) %>%
  filter(term == "month_treat") %>% mutate(model = "Trust Federal")
ptrends_dat_parliament <- ptrends_dat %>% 
  group_by(east, yborn) %>%
  group_modify(~tidy(felm(trust_parliament~month_treat*mborn|0|0|0, data = ., weights = .$weights))) %>%
  filter(term == "month_treat") %>% mutate(model = "Trust Parliament")
ptrends_dat_econsit <- ptrends_dat %>% 
  group_by(east, yborn) %>%
  group_modify(~tidy(felm(econ_sit~month_treat*mborn|0|0|0, data = ., weights = .$weights))) %>%
  filter(term == "month_treat") %>% mutate(model = "Financial Situation")
ptrends_dat <- rbind(ptrends_dat_FCC,ptrends_dat_federal,ptrends_dat_parliament,ptrends_dat_econsit) #putting data togeher for plot
ptrends_dat$model <- factor(ptrends_dat$model, levels = c("Trust FCC", "Trust Federal", "Trust Parliament", "Financial Situation")) #Ordering facets for plot
interval <- -qnorm((1-0.90)/2)
ptrends_plot <- ggplot(ptrends_dat, aes(x = factor(yborn), y = estimate, color = factor(east), shape = factor(east))) +
  geom_point(size=2,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = estimate - std.error*interval, ymax = estimate + std.error*interval), width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_line(aes(group = factor(east), linetype = factor(east)),
            position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() + facet_wrap(~model) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  guides(guide_legend(reverse = TRUE)) + xlab("Birth Year") + ylab("Difference at Cutoff") +
  ggtitle("Parallel Trends East and West Germany") 
ggsave("Figure_1.pdf", plot= ptrends_plot,
       width = 11, height = 8.5)

closeAllConnections() # Close connection to log file

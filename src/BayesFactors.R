# LIBRARY IMPORTS =======================================================================
library(lme4)
library(emmeans)
library(brms)
library(tidyverse)
library(future)
library(future.apply)
plan(multiprocess)

# BAYESIAN POINT NULL HYPOTHESIS TESTING ================================================
H_naught.test <- function(N){
  # Generate random data
  # mean=0 for null, sd=.5 from real data
  new_old.sims <- future_replicate(100,
                                   list(tibble(ChanceArcsin = rnorm(N*2, 0, .5),
                                               Participant = rep(1:N, each = 2),
                                               Condition = factor(rep(c("Label", "No Label"),
                                                                      times = N/2, each = 2)),
                                               ContrastType = factor(rep(c("Head", "Tail"),
                                                                         times = N)))))
  
  # Define STB analysis function, returning emmeans analysis
  stb.analysis <- function(df.list){
    e <- future_lapply(df.list,
                       function(df){
                         m <- lmer(ChanceArcsin ~ ContrastType*Condition +
                                     (1 | Participant),
                                   data = df)
                         t <- emmeans(m, ~ ContrastType | Condition,
                                      options = list(infer = c(T, T), null = 0,
                                                     level = .95)) %>%
                           as_tibble()
                         return(t)
                       }) %>%
      bind_rows()
    return(e)
  }
  
  # Define Bayesian analysis function, returning hyp. test bf
  bayesian.analysis <- function(df.list){
    p <- c(set_prior("normal(0,.5)", class = "Intercept"),
           set_prior("normal(0,.5)", class = "b"))
    bf <- future_lapply(df.list,
                        function(df){
                          bf <- brm(ChanceArcsin ~ ContrastType + Condition +
                                      ContrastType:Condition +
                                      (1 | Participant),
                                    data = df, prior = p, family = gaussian(),
                                    chains = 4, iter = 2000,
                                    save_all_pars = T) %>%
                            hypothesis(c("Intercept > 0",
                                         "Intercept + ContrastTypeTail > 0",
                                         "Intercept + ConditionNo Label > 0",
                                         paste("Intercept +",
                                               "ConditionNo Label +",
                                               "ContrastTypeTail +",
                                               "ContrastTypeTail:ConditionNo Label",
                                               "> 0")))
                          return(as_tibble(bf$hypothesis))
                        }) %>%
      bind_rows()
    return(bf)
  }
  
  # Get evidence summary from simulations
  t <- proc.time()
  new_old.sims.stb.evid <- stb.analysis(new_old.sims)
  stb.time <- proc.time() - t
  print("STB analysis over")
  t <- proc.time()
  new_old.sims.bayesian.evid <- bayesian.analysis(new_old.sims) %>%
    mutate(Condition = factor(ifelse(grepl("NoLabel", Hypothesis), "No Label", "Label")),
           ContrastType = factor(ifelse(grepl("Tail", Hypothesis), "Tail", "Head")))
  bayesian.time <- proc.time() - t
  
  return(list(BayesianEvidence = new_old.sims.bayesian.evid,
              SampleTheoryEvidence = new_old.sims.stb.evid,
              BayesianTime = bayesian.time,
              SampleTheoryTime = stb.time))
}

new_old.sims.results.48 <- H_naught.test(N = 48)
bayesian.95 <- quantile(new_old.sims.results.200$BayesianEvidence$Evid.Ratio, .95)
stb.95 <- quantile(new_old.sims.results.200$SampleTheoryEvidence$p.value, .05)
bayesian.above3 <- sum(new_old.sims.results.200$BayesianEvidence$Evid.Ratio > 3) / 400
stb.bellow_dot05 <- sum(new_old.sims.results.200$SampleTheoryEvidence$p.value<.05) / 400
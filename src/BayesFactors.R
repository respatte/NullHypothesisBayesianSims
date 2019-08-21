# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(emmeans)
library(brms)
library(tidyverse)
library(scales)
library(future)
library(future.apply)
plan(multiprocess)

# BAYESIAN POINT NULL HYPOTHESIS TESTING ===========================================================
# Define simulation function
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
    future_lapply(seq_along(df.list),
                  function(i){
                    lmer(ChanceArcsin ~ ContrastType*Condition +
                           (1 | Participant),
                         data = df.list[[i]]) %>%
                      emmeans(~ ContrastType | Condition,
                              options = list(infer = c(T, T), null = 0,
                                             level = .95)) %>%
                      as_tibble() %>%
                      mutate(Sim = i) %>%
                      return()
                  }) %>%
      bind_rows() %>%
      return()
  }
  
  # Define Bayesian analysis function, returning hyp. test bf
  bayesian.analysis <- function(df.list){
    p <- c(set_prior("normal(0,.5)", class = "Intercept"),
           set_prior("normal(0,.5)", class = "b"))
    future_lapply(seq_along(df.list),
                  function(i){
                    bf <- brm(ChanceArcsin ~ ContrastType + Condition +
                                ContrastType:Condition +
                                (1 | Participant),
                              data = df.list[[i]], prior = p, family = gaussian(),
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
                    bf$hypothesis %>% as_tibble() %>%
                      mutate(Sim = i) %>%
                      return()
                  }) %>%
      bind_rows() %>%
      return()
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

# Get simulations (takes about one hour on 2.3GHz a quadcore with 16G RAM + swap)
run_sims <- F
if(run_sims){
  new_old.sims.results.48 <- H_naught.test(N = 48)
  saveRDS(new_old.sims.results.48, file = "../results/NewOld_48.rds")
}else{
  new_old.sims.results.48 <- readRDS("../results/NewOld_48.rds")
}

# Get quantiles for STB and Bayesian stats
bayesian.95 <- quantile(new_old.sims.results.48$BayesianEvidence$Evid.Ratio, .95)
stb.95 <- quantile(new_old.sims.results.48$SampleTheoryEvidence$p.value, .05)
bayesian.above3 <- sum(new_old.sims.results.48$BayesianEvidence$Evid.Ratio > 3) / 400
stb.bellow_dot05 <- sum(new_old.sims.results.48$SampleTheoryEvidence$p.value<.05) / 400

# Combine data for plot and analysis
bayesian.evidence <- new_old.sims.results.48$BayesianEvidence %>%
  select(c(Sim, Evid.Ratio, Condition, ContrastType, Hypothesis)) %>%
  mutate(Hypothesis = case_when(grepl("ContrastTypeTail:ConditionNoLabel",
                                      Hypothesis) ~ "(Interaction)",
                                grepl("ContrastTypeTail", Hypothesis) ~ "ContrastType",
                                grepl("ConditionNoLabel", Hypothesis) ~ "Condition",
                                T ~ "(Intercept)"),
         Hypothesis = parse_factor(Hypothesis, levels = c("(Intercept)",
                                                          "ContrastType",
                                                          "Condition",
                                                          "(Interaction)")))
stb.evidence <- new_old.sims.results.48$SampleTheoryEvidence %>%
  select(c(Sim, p.value, Condition, ContrastType))
global.evidence <- full_join(bayesian.evidence, stb.evidence)

# Plot b-values over p-values
generate_plots <- F
if(generate_plots){
  b_to_p.plot <- global.evidence %>%
    ggplot(aes(x = p.value,
               y = Evid.Ratio,
               colour = Hypothesis)) +
    scale_y_continuous(trans = log10_trans()) +
    ylab("Bayes Factor") + xlab("p-value") + theme_bw() +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    geom_point(alpha = .5) +
    geom_vline(xintercept = .05, colour = "black",
               linetype = "62", alpha = .5) +
    scale_colour_brewer(palette = "Dark2",
                        labels = c("Intercept",
                                   "Main Effect 1",
                                   "Main Effect 2",
                                   "Interaction"))
  ggsave("../results/b-value_p-value.pdf",
         plot = b_to_p.plot,
         width = 5, height = 5,
         dpi = 600)
}
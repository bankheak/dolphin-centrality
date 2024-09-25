# Code

This file contains all the code for the Dolphins repository. The following describes the analysis steps in the `main_code.R` file. 

## Data analysis process
<img src="https://github.com/user-attachments/assets/da400ccb-1dc0-4ff4-b630-4b3535e9d98e" align="middle" width="500px"/>

## PART 1: Data Wrangling

I start by fixing and combining data from 1993-2014 and seperating 6-year periods between 1995-2012. I then only include individuals that have been seen at least 10 times in all three study periods. I then look at how many individuals performed each human-centric foraging behavior.

## PART 2: Calculate Local Metrics

I calculated social centrality (degree and strength of connections) of each individual within each study period. I then calculate the proportion of time they spent engaging in each human-centric behavior.

## PART 3: Run Model

I built a linear mixed model (LMM) using a Markov chain Monte Carlo (MCMC) sampler under a Bayesian statistical framework. I compared leave-one-out cross-validation information criteria (LOOIC) between a null model, an additive model, and an interactive model. I selected the model with the lowest significant change in the expected log-likelihood predictive density (ELPD) of the approximated LOOIC. The most parsimonious model was the interaction model, which described each covariates’ effect on the dyadic association index:

Process model:
$Social Centrality ~ u_i+β1_{BG}+β2_{FG}}+β3_{SD}+β4_{BG}*During HAB+β5_{BG}*After HAB+β4_{FG}*During HAB+β5_{FG}*After HAB+β4_{SD}*During HAB+β5_{SD}*After HAB$

Observation model:
$SRI_{i,j,p} ~ Normal[μ_{SRI_{i,j,p}},ϕ_{SRI_{i,j,p}}]$

I then create the figures of the effect sizes of these predictors.

## PART 6: Display Networks

I create the figure for the network containing clusters on and off the Sarasota map.



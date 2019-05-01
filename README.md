# Parallel-Experimentation

This repository contains R code that generates data for the analysis of parallel experimentation between competing advertisers. Below you will find a more detailed explanation of the code's goals and approach.

There are J competing advertisers. The advertiser of interest, whose average treatment effects (ATEs) we aim to recover, is advertiser j. The dependent variable of interest, Y, is a Bernoulli random variable, which in our application is an indicator for whether the user visited j's product's detail page (you can find the latest draft at https://arxiv.org/abs/1903.11198).

The probability a user visits j's page is a function of which advertisers are actually advertising (i.e. present) and which ones are not advertising (i.e. absent). There are 2^J such present/absent combinations, whose corresponding probabilities are given in the vector theta. 

In the experimental setup, each of these combinations corresponds to a total treatment assignment, which is a vector of J indicators of whether the user is in each advertiser's treatment or control groups. Accordingly, there are also 2^J treatment/control combinations, which are given in the matrix grid_tot. From j's perspective, however, there are only 2^(J-1) relevant combinations, which reflect whether the remaining J-1 advertisers are present or not. Hence, in the experimental design there are 2^(J-1) such treatment assignments, given in the matrix grid.

Given theta, we can compute the ATE for j at each of the 2^(J-1) treatment assignments. They are the objects we wish to estimate.

The code proceeds to run S simulation rounds. The goal of running several simulations is to potentially compare results from simply running linear regressions with the kernel-based estimator proposed by Li, Ouyang and Racine (2013) (you can find the paper at https://onlinelibrary.wiley.com/doi/abs/10.1002/jae.1261).

For each simulation s, the code randomly assigns each of n users to one of the 2^J total treatment assignments with equal probability, that is, with probability 0.5/2^J. For each user, the code randomly draws Y from a Bernoulli distribution with parameter corresponding to that of the total treatment assignment given to that user. 

The code then implements two estimators. First, it runs a regression of Y on a constant and a treatment indicator for advertiser j separately for each of the 2^(J-1) treatment assignments. Second, it implements Li, Ouyang and Racine (2013)'s method to smooth over the treatment indicators for the J-1 remaining advertisers. The implementation makes use of the np package, created and mantained by Professor Jeffrey S. Racine.

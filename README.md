# Parallel-Experimentation

This repository contains R code that generates data for the analysis of parallel experimentation between competing advertisers. Below you will find a more detailed explanation of the code's goals and approach.

There are J competing advertisers. The advertiser of interest, whose average treatment effects (ATEs) we aim to recover, is advertiser j. The dependent variable of interest, Y, is a Bernoulli random variable, which in our application is an indicator for whether the user visited j's product's detail page (you can find the latest draft at https://arxiv.org/abs/1903.11198).

The probability a user visits j's page is a function of which advertisers are actually advertising (i.e. present) and which ones are not advertising (i.e. absent). There are $2^J$.

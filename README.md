# emews

This repository contains the functions and scripts used for analysing seal diving data using hidden Markov models. The specific implementation has been developed for adult Weddell seals and is presented in detail in the "Sex-specific variation in the use of vertical habitat by a resident Antarctic top predator" by Theoni Photopoulou, Karine Heerah, Jennifer Pohle and Lars Boehme (2020). 

The functions used to fit the model are in the file "hmm7Rcpp.R" and the Viterbi algorithm for calculating the most likely state sequence is in file "viterbi_hmm7Rcpp.R". The numerical maximisation of the likelihood is implemented in C++ in the file "mllk_cov.cpp". 

There are two script for fitting the model to data. One is for female Weddell seals "diveHMM_Female_hourwk.R" and one for male Weddell seals "diveHMM_Male_hourwk.R".

The dataset used in the analysis is available at the following URL <http://doi.org/10.5281/zenodo.3820359>

Please contact me at <theoni.photopoulou@gmail.com> with any questions.
# CYGNO Sensitivity Study with BAT

This repository contains the code developed to study the **CYGNO experiment sensitivity** using the **BAT library**.  
BAT is a C++ software package designed to address statistical problems within the framework of Bayesian inference. Utilizing Bayes' Theorem, it implements Markov Chain Monte Carlo (MCMC) techniques, allowing access to the full posterior distribution, enabling comprehensive parameter estimation, limit setting, and uncertainty quantification.

## Prerequisites
You need to have ROOT installed to be able to install BAT. Follow the instructions here to install BAT: https://github.com/bat/bat  
The following Python libraries are required to plot:
- uproot
- array 
- numpy
- matplotlib.pyplot 
- mpl_toolkits.axes_grid1

## Setup
Once BAT is installed, set the number of events in the input files:
- background events: input/background.txt
- data events: input/data.txt
- [background mean](https://github.com/antorita/cygno_sensitivity/blob/1009ffac2a68482bb6118f5ae11d5f63b60b7680/oldmodel.cxx#L19) in oldmodel.cxx
- [background sigma](https://github.com/antorita/cygno_sensitivity/blob/1009ffac2a68482bb6118f5ae11d5f63b60b7680/oldmodel.cxx#L20) in oldmodel.cxx

In this study, assuming no Dark Matter interactions are observed, **20% of the runs were randomly selected to estimate the background**, while the **remaining 80% were used as the main data set**. The background is then rescaled to the data.
This approach ensures that the **background level is properly accounted for** while maintaining the majority of the data for evaluating the exclusion limit.

## Compilation
Compile the code to generate the executable:

```
make
```

## Running the code
After compilation, run the code using:

```
./runfit
```
The output will be a **ROOT file** containing the distribution of $\mu_b$ (background) and $\mu_s$ (signal), where $\mu_b$ and $\mu_s$ are the posterior distribution.

Note: You may need to create the results folder (where the output is saved) before running the code.

To plot the results and evalualte the **90% Confidence Interval (C.I.)**, run:
```
python3 plot.py
```


## Statistical approach
The statistical analysis of the data uses a **Bayesian approach**.
When the true value of a parameter is consistent with zero, our case for the expected number of DM interactions, the **90% confidence interval (C.I.)** is used to define the **upper limit**. The definition is:  
![Bayesian apper limit](https://latex.codecogs.com/svg.image?\mu(90\%C.I.)=\int_0^{\mu_{90\%}}P(\mu\mid\vec{x},H)\,d\mu=0.9)  
The signal is the number of the interactions observed in the data. [Signal prior](https://github.com/antorita/cygno_sensitivity/blob/1009ffac2a68482bb6118f5ae11d5f63b60b7680/oldmodel.cxx#L81) is chosen to be a uniform distribution because we do not have any hint of the value of the DM cross-section. While the probability of each target element recoiling due to background events measured in the selected samples would follow a Poissonian distribution. However it has been rescaled to match the signal exposure time. This implies that the background uncertainty should be correctly scaled to the ratio data to background exposure. Given that for a larger number of events a Poissonian distribution can be approximated with a Gaussian distribution, it is used the latter for the [background prior](https://github.com/antorita/cygno_sensitivity/blob/1009ffac2a68482bb6118f5ae11d5f63b60b7680/oldmodel.cxx#L78) to set the proper sigma values.  
The [Likelihood function](https://github.com/antorita/cygno_sensitivity/blob/1009ffac2a68482bb6118f5ae11d5f63b60b7680/oldmodel.cxx#L66) for evaluating the confidence interval (C.I.) in this rare event analysis with a background contribution is defined as the Poissonina distribution:
![Likelihood](https://latex.codecogs.com/svg.image?&space;L(\vec{x}\mid\mu_s,\mu_b,H)=\frac{(\mu_b&plus;\mu_s)^{N_{evt}}}{N_{evt}!}e^{-(\mu_b&plus;\mu_s)})
where, $N_{evt}$ is the observed total number of events, $\mu_s$ is the signal contribution from WIMP interactions, and $\mu_b$ is the expected background contribution.  

To evaluate the **Dark Matter exclusion limit**, see: https://github.com/antorita/limit_evaluation

## References
For more details, see **Chapter 7 of my PhD thesis** https://arxiv.org/abs/2510.01646  
BAT manual: https://bat.github.io/bat-docs/master/manual/BAT-manual.pdf

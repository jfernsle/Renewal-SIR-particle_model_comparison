# Renewal-SIR-particle model comparison

This repository contains MATLAB 2020a code that executes particle simulations of an idealized epidemic, which are compared to SIR and Renewal models.

## Description

An epidemic model that can accurately respond to changing conditions is necessary for informed decision-making, and presently the SIR and SEIR models are the standard methods for predicting pandemic behavior. These models derive from the classic paper by Kermack and McKendrick (1927), but embedded in that paper is a second, more general approach in which infectivity is specified as a function of time after infection. The second approach, called the “renewal model” (RM) here, is closely related to an even older model by Lotka (1907) for population growth. To determine which model is more accurate, we compare results from both with those produced by particle simulations of highly idealized epidemics. Although all three models behaved similarly, the (smoothed) results produced by the particle simulations differed significantly from the SIR solutions but matched the RM solutions almost perfectly. Even more troubling, the SIR model responded differently to changes in the most fundamental parameter in epidemics, the basic reproduction number R0. Similar poor behavior is expected for the SEIR model as well. Models which fail simple tests like this are, of course, highly unlikely to treat more complicated cases better. Consequently, given the notably improved performance of RM, and given that its parameters are directly measurable and can be adjusted to accommodate changing conditions if needed, the renewal model appears to be the better choice. However, the renewal model as originally formulated is incomplete, because it fails to give the populations of the infectors I(t) and recovered people R(t). Fortunately, those populations can be determined by solving an additional, integral-differential equation using the RM solution for the infection rate. 

## Getting Started

Here are some resources for epidemic modeling that are useful for understanding this epidemic modeling project. 
The Centers for Disease Control and Prevention use an ensemble of models to forecast COVID infections found here: https://www.cdc.gov/coronavirus/2019-ncov/science/forecasting/forecasts-cases.html While each model uses different assumptions to generate their forecasts, the dominant class of model is SIR/SEIR, which uses a set of ordinary differential equations to solve for the Susceptible, Exposed, Infected and Recovered populations. As of 11/10/2021, three models on the CDC site are adaptive to changing social conditions to predict their forecast and all three are SEIR-type models.
Another class of epidemic model, the Renewal model, is used much more rarely but includes a function that calculates future infections based on how probable an infected person is to pass on a secondary infection a certain number of days after their original infection. This function is otherwise known as the generation-interval distribution, where the generation interval is the time between when an infected person becomes infected and when they pass on infection to anothe rperson. This distribution can be estimated from contact tracing and the measurement of the serial interval distribution, such as H. Nishiura, N. M. Linton, A. R. Akhmetzhanov, Serial interval of novel coronavirus (COVID-19) infections. Int. J. Infect. Dis. 93, 284–286 (2020). The Imperial College of London uses a Renewal model for their covid forecasthttps://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/covid-19-planning-tools/
A third class of model uses individual 'agents' which pass on infections to other agents in the simulation.  These are often visualized as particles representing people, which move within a simulation space and pass on an infection when an infected particle comes within the radius of infection with a susceptible particle. This type of model is stochastic because the particle positions are randomized and often include a Monte Carlo routine to determine whether the infection is passed. These models are less practical than the SIR/SEIR and Renewal models since they require large computing overhead for many particles, but they are more directly comparable to a real epidemic and are helpful to visualize epidemic dynamics. See, for example, a simple video by 3Blue1Brown https://www.youtube.com/watch?v=gxAaO2rsdIs. 

### Dependencies

Code included in this repository requires MATLAB. It was coded on MATLAB 2020a, but most functions are generic and should work with older versions of MATLAB.

### Installing

There are three primary pieces of code included for this project.
1) PNAS_paper_particle_code_v1.m This code allows the user to specify conditions for a particle simulation of a highly idealized epidemic, where infectors move through a sea of stationary susceptible particles. 
2) PNAS_paper_RM_and_SIR_models_v1.m This code inputs early data from a particle simulation of an epidemic to find initial conditions and parameters to run the Renewal model and the SIR model. The solutions to the models are plotted on top of data from particle simulations for comparison purposes.
3) Model_statistics_analysis_v1.m This code reads a series of particle simulation data run for the same scenario to explore the uncertainty in epidemic outcomes and the predictions of the RM and SIR models. Example data is presented in folder \Linear_pi_data\ which used a linearly decreasing probability of infection for 10 days to generate 10 simulations. 'Contact tracing' of particle infections was used to measure the generation-interval distribution over the first 15 days of the simulation. The linearly-fit distributoin was used as gi(tau) for each simulation, which was further used to solve for R0, used in RM and SIR models from Eq. 1 in model paper. This is similar to what is done for real epidemiology.

### Executing program

Code should be run in the order of 1) PNAS_paper_particle_code_v1.m, 2) PNAS_paper_RM_and_SIR_models_v1.m.  Model statistics from 3) Model_statistics_analysis_v1.m should be run after a series of particle code data sets using code 1) have been run and stored as .csv files.
* Particle code 1) can run many scenarios with variables includeing particle number, density, velocity, etc. all which can lead to a predicted R0 found from the equation in Table 1. The probability of passing on an infection can be altered for each particle as a function of the time after infection, or distributoins of infectious time can be chosen. Similar scenarios to Fig. 2-5 can be generated by running the code sections labeled by Figure at the beginning of the script.
* The Renewal modela and SIR model can be run after a particle simulation is either stored in memory on MATLAB, or has been read from a .csv file. Depending on the particle's time-dependent infectiousness and/or distributoin of infectious time, you should choose between options A-D to select gi(tau) and RM(tau).  For the SIR model, you can choose between the parametric and nonparametric parameter choices that use an exponential fit to early data. It is also possible to use the generation-interval distribution from the particle code as your gi(tau) function and solve for R0 using this function and Eq.1. This method is similar to what is done in real epidemiology.


## Help

Example data sets were run using the infection rate i(t) recorded as the number of infections at a particular time step.  These must be converted to an infetion rate by dividing by the time step dt.  The current version of particle code 1) records a the actual infection rate (per day).


## Authors

Contributors names and contact info

Jonathan Fernsler, Richard Fernsler, Alexis Mora Solick, Peter Lenz, Patrick Sullivan, and Madison Singleton
Corresponding author: Jonathan Fernsler
Address: Physics Dept., Cal Poly, San Luis Obispo, CA 93407
Email: jfernsle@calpoly.edu
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 1.0
    * Initial version of scripts


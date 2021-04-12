# HNM
Heterogenous Neural Model (HNM) for large scale modelling of brain activity.

Code here is to generate the simulations and the results for the manuscript titled "Dynamical consequences of regional heterogeneity in the brain’s transcriptional landscape"

Gustavo Deco*, Morten L. Kringelbach*, Aurina Arnatkevičiūtė , Stuart Oldham , Kristina Sabaroedin , Nigel C. Rogasch , Kevin Aquino+ and Alex Fornito+

*,+ Equal contribution

This model is an extension of the Deco et al. 2014 model, where each region's excitability is varied according to known external measures: such as the T1T2 ratio (marker for Myeline content), the first principle component from a vector all genetic coexpressions, or the ratio of genetic co-expressions of excitatory vs inhibitory coxpresions (see the preprint: https://www.biorxiv.org/content/10.1101/2020.10.28.359943v1 for more information).


## Model execution
The model is run in three seperate instances: Model Optimization, Model applications, and Measurement of model applications.

### Model Optimization  

Model optimization is under the subfolder GrandModelOptimization, which contains the code that is optimized for SLURM based cluster systems. It runs each model and finds the optimal balance of the Scaling (Z) and the Bias (B) factors, (See equations 12 and 13 for where they fit in) through multiple simulations to optimzed Node-FC, Edge-FC and FCD. It is an exhaustive search and takes ~2 months spread over multiple jobs and cores on our High Performance cluster.

This will generate the figure below

<img width="500" alt="image" src="https://user-images.githubusercontent.com/6628199/114129188-6a489e00-9941-11eb-8b36-2902015516f4.png">

### Model applications

Once the model has been optimzed for 5 seperate variations, we have stored the optimized results so that to run Ignition for a particular model, such as the Balanced Excitation model (BEI):

``` matlab
addpath(genpath(pwd))
run_id='1'
model='BEI'
perturb_range='1'
jbal_flag=0
run_ignition_decay(run_id,model,perturb_range,jbal_flag);
```

This runs the ignition simulation for run_id 1, this is coded in a way to work on a cluster (which we performed calculations on MASSIVE https://www.massive.org.au) and using the Balanced Excitation model. The parameter perturb_range sets the range of pertbuation values, it is set to sample the pertubation vector for ignition sparsely so it may be calculated more finely through parallelization. The last parameter jbal_flag is set in order to have a pre-calculation of the inhibitory feedback J_i in the model), if it is 1 it calculates the J then saves it and exits, if it is set to 0 it runs through. 

To calculate the Autocorrelation functions due to noisy fluctuations in V1:
``` matlab
addpath(genpath(pwd))
run_id='1'
model='BEI'
jbal_flag=0
SEED='10'
run_decay_models_ACF(run_id,model,jbal_flag,'SEED');
```
As the last code, the additional component is the ability to simulate the ACF due to different seeds SEED=10 corresponds to the lateral occipital cortex in the Deskian-Killany parcellation. 

<img width="500" alt="image" src="https://user-images.githubusercontent.com/6628199/114130921-f14b4580-9944-11eb-85ff-1c32f66066de.png">

<img width="280" alt="image" src="https://user-images.githubusercontent.com/6628199/114130975-06c06f80-9945-11eb-8908-c15a0d871710.png">

### Measurement of model applications

To measure metrics from the model simulations, two helper functions aid this:
```matlab
ignition_decay
```
and
```matlab
time_decay_ACF
```
Which calculate the ignition metrics, and the autocorrelation functions in the manuscript. 

<img width="339" alt="image" src="https://user-images.githubusercontent.com/6628199/114131009-13dd5e80-9945-11eb-8112-f1c452e9118b.png"><img width="228" alt="image" src="https://user-images.githubusercontent.com/6628199/114131023-1b046c80-9945-11eb-9925-8665171f9fd1.png">


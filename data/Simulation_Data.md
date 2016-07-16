# Variables in simulation data


## Y variables

### Y1 to Y10

Predictable variable by X:V1-V6, with increasing positives, from 50 to 500

### R1 to R100

Random variables with increasing positives, approximately from 0.1 to 10%, not predictable by any X variables

### F0 to F4

Predictable correspondingly by confoundering variables X:C0-C1. Some of them can also be predicted by X:V1-V6, but removing the confounders will reduce their predictivity. 


## X variables

### V1 to V6

Predictive variables, partially correlated to each other, some of them have different scales and directions. 

### C0 to C3

Confounding variables of V1 to V6, with increasing confounding effect. C0 has none. 

### F1 to F3

Highly correlated to confounders C1 to C3, correspondingly. Can predict Y10 as a single predictor variable. However, removing their correlation to corresponding confounders will get rid of the predictivity. In other words, there will be false positive results if the confounder is not inlucded in the model.

### N1 and W1

N1 is a confounder highly predictive of Y10, but impair the predictivity of the predictor variable W1. Removing confounding effect of N1 will increase the predictivity of W1. 


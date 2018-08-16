This folder contains simulation scenarios with time-to-event endpoint.
It is subdivided into 2 folders:

- 0_predictive (simulations for the null case with no predictive biomares)
- 2_predictive (simulations with 2 predictive biomarkers)

The simulate.R files should be run in BATCH mode using the following line:

`nohup R CMD BATCH "--args 2000" --vanilla simulate.R`

The args argument correspond to the number of subjects for the simulation. If it is not specified, we set it to 100.

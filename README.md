# MCD-variogram-estimator
R-Code for the simulation study and the application of  the paper "Nonparametric directional variogram estimation in the presene of outlier blocks".

Test

- Functions: an ordner containing the following documents
      - Auxiliary-functions.R: auxiliary functions needed for the variogram estimations, e.g. to build the vectors for the 
                               MCD variogram estimators
      - function-for-simulation.R: function to simulate one iteration, i.e. to calculate the different variogram estimators
                                   for one grf
      - functions-data-simulation.R: function to simulate the grf (without outliers, with a outlier block or with isolated outliers)
      - MCD-variogram-estimator.R: functions to calculate the different variogram estimators

 - Data: an ordner containing the following ordners
      - Block: simulated grf with an outlier block for different block sizes and different outlier distributions (needed for section 4.4)
      - Consistency: simulated grf without outliers for different grid sizes (neede for section 4.2 and 4.3)
      - Correctionfactors: simulated grf without outliers for different grid sizes needed for the simulation of the correction factors
                           (needed for section 4.1)
      - Isolated: simulated grf with isolated outliers for different amounts of outliers and differen outlier distributions
                  (needed for section 4.5)

 - Block: an ordner with the following .RData files
      - res_block_...: output of the simulations for grfs with an outlier block for different block sizes and different outlier distributions
                       (for section 4.4)
      - estimation_bias: an array with the bias per lag for each estimator for the different scenarios of the simulations
                         (needed for graphs in section 4.4 and Table 3 in the appendix)
      - esimtaion_mean: an array with the average estimation per lag for each estimator for the different scenarios of the simulations
      - estimation_sqrtMSE: an array with the sqrtMSE per lag for each estimator for the different scenarios of the simulations
                            (needed for Table 3 in the appendix)

 - Consistency: an ordner with the following .RData files
      - res_cons_...: output of the simulations for grfs without outliers for different grid sizes
                      (for section 4.2 and 4.3)
      - estimation_bias: an array with the bias per lag for each estimator for the different grid sizes
                         (needed for graphs in section 4.2)
      - estimation_errors: an array with the estimation errors per lag and per iteration for each estimator for the different grid sizes
                           (needed for the boxplots in section 4.3)
      - esimtaion_mean: an array with the average estimation per lag for each estimator for the different grid sizes
      - estimation_sqrtMSE: an array with the sqrtMSE per lag for each estimator for the different grid sizes
                            (needed for graphs in section 4.2)

 - Correctionfactors: an ordner with the following .RData files
      - res_corr_...: output of the simulations for grfs without outliers for different grid sizes and different true variogram models 
                      for simulation of the correctionfactors (for section 4.1)
      - correctionfactors: an array with the simulated correctionfactors for the different estimators, the different grid sizes and 
                           the different true variogram models (needed for the Table in section 4.1 and Table 2 in the appendix)

 
 - Isolated: an ordner with the following .RData files
      - res_Iso_...: output of the simulations for grfs with isolated outliers for different block sizes and different outlier distributions
                     (for section 4.5)
      - estimation_bias: an array with the bias per lag for each estimator for the different scenarios of the simulations
                         (needed for graphs in section 4.5)
      - esimtaion_mean: an array with the average estimation per lag for each estimator for the different scenarios of the simulations
      - estimation_sqrtMSE: an array with the sqrtMSE per lag for each estimator for the different scenarios of the simulations
   

 - Modified: an ordner with the following .RData files
      - res_block_...: output of the simulations for grfs with an outlier block for different block sizes and different outlier distributions
                       including the modified mcd variogram estimators defined in section 4.3 (for section 4.4)
      - res_corr_...: output of the simulations for grfs without outliers for different grid sizes for simulation of the correctionfactors 
                      including the modified mcd variogram estimators defined in section 4.3 
      - res_cons_...: output of the simulations for grfs without outliers for different grid sizes including the modified mcd variogram estimators
                      defined in section 4.3 (for section 4.3)
      - correctionfactors: an array with the simulated correctionfactors for the different estimators including the modified mcd estimators,
                           and different grid sizes 
      - estimation_bias: an array with the bias per lag for each estimator for the different scenarios of the simulations for grfs without outliers
      - esimtaion_mean: an array with the average estimation per lag for each estimator for the different scenarios of the simulations for grfs without outliers
      - estimation_sqrtMSE: an array with the sqrtMSE per lag for each estimator for the different scenarios of the simulations for grfs without outliers  
      - estimation_errors: an array with the estimation errors per lag and per iteration for each estimator for the different grid sizes for grfs without outliers
                           (needed for the boxplots in section 4.3) 
      - estimation_bias_block: an array with the bias per lag for each estimator for the different scenarios of the simulations for grfs with an outlier block
                               (needed for Table 3 in the appendix)
      - esimtaion_mean_block: an array with the average estimation per lag for each estimator for the different scenarios of the simulations for grfs with an outlier block
      - estimation_sqrtMSE_block: an array with the sqrtMSE per lag for each estimator for the different scenarios of the simulations for grfs with an outlier block 
                                  (needed for Table 3 in the appendix)

- Graphs: an ordner containting the pdf files with the graphs of the paper

- L8_cropped: an ordner containing the data of the application in section 5

- Application-Satellite-Data.R: R-Code for section 5 (Application to the satellite data)

- figure-tables.R: R-Code to obtain the graphs and tables of the paper (expect them of section 5)

- Parametercombinations.xlsx: an excel with all parameter combinations of the simulations for scenarios of the paper 

- Parameters.R: R-Code, which defines the parameters for the simulations

- simulation-block-outlier.R: R-Code of the simulations for section 4.4 (Contamination with block outliers)

- simulation-correctionfactors.R: R-Code of the simulations of the correctionfactors (section 4.1)

- simulation-grf-konsistency.R: R-Code of the simulations for section 4.2 and 4.3 (Gaussian Data and Consistency)

- simulation-isolated-outlier.R: R-Code of the simulations for section 4.5 (Contamination with isolated outliers)

- simulation-modified-estimators.R: R-Code for the simulation of the modified mcd variogram estimators defined in section 4.3 (needed for section 4.3 und 4.4)

- simulation-of-data.R: R-Code for the simulation of the grf (with and without outliers)              

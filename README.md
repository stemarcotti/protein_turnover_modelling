# Protein Turnover Modelling

See for reference [Matsubayashi, Y, Sánchez-Sánchez, BJ, et al. Developmental Cell (2020)](https://www.cell.com/developmental-cell/fulltext/S1534-5807(20)30455-X). An additional protocol paper is currently work in progress, a link will be added here when available. Copyright - Stefania Marcotti (Stramer lab) 2020, Code tested on MATLAB v.2018b. For issues please contact us at https://www.stramerlab.com/contact.html

#### Anterograde modelling `[anterograde_model.m]`
This code takes as input the fitted mRNA and the logistic parameters of the protein of choice. It calculates the turnover rates for synthesis and degradation. It calls the function `[anterograde_funct.m]`.

Solve numerically for P(t) by using as input the experimental mRNA data M(t) as input. The synthesis and degradation rates Sp and Dp are calculated by non-linear regression between the calculated solution and the corresponding experimental data.

#### Retrograde modelling `[retrograde_model.m]`
This code takes as input the fitted mRNA and the logistic parameters of the protein of choice. It calculates the turnover rates for synthesis and degradation. It calls the function `[retrograde_funct.m]`.

Solve analytically for M(t) by using as input the experimental protein expression P(t). The synthesis and degradation rates Sp and Dp are calculated by non-linear regression between the calculated solution and the corresponding experimental data.

##### Instructions
* Click on the green button ‘Code’ and select ‘Download ZIP’ to save the file locally on your machine
* Unzip the folder
* Open MATLAB and navigate to the unzipped folder where the code is located. Further assistance on how to do so can be found [here](https://uk.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html)
* Input files in csv format containing the interpolated mRNA profile and the fitted logistic parameters for the protein of interest should be saved in an empty folder. The first file should include two columns, one for time and one for mRNA values, and the number of rows will depend on the chosen time steps and interval; the second file should include three rows (K, ti, r) and one column for each biological replicate.
* In MATLAB, open the chosen model and hit Run. The code will request you to locate the two input files and to choose the time interval and span of the experiments.
* Both models will produce four output files in csv format, containing the values for the production and degradation rates and their confidence intervals. Each row represents a biological replicate.

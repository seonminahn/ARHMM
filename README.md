# ARHMM
Supplements files for "Autoregressive Statistical Modeling of 
a Peru Margin Multi-Proxy Holocene Record Shows Correlation Not Cause, Flickering Regimes and Persistence"

To run the code,
* Locate data under the newData folder
* Download the fminsearchbnd and fminsearchcon function files from File Exchange: 
    https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon Software
* Run the run_hmm_hmmar1001.m script file to run the algorithm. 

The data and results used this study are under the newData and results_hmm_ar_2_state folder.
* newData\new4D_normalized_flipped_doubleTime.mat includes a variable sampleD. Data for Figure 2.
    * first column: time step
    * second - fifth columns: SST, C37, del15N, %N
* results_hmm_ar_2_state\new4D_normalized_flipped_doubleTime_results_dat57.mat include results
    * dataFilled. Data for Figure 4.  
        * first row: time step
        * second - fifth rows: SST, C37, del15N, %N
    * Estimated parameters. Values in Table 2
        * a (transition probability), means (c), thetaE1 & thetaE2 (theta), sigmaE1 & sigmaE2 (sigma) 
        * variable name in the mat file (names in the paper)
    * The variable st include 1000 sampled states for using back sampling

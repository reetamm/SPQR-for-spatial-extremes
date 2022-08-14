# SPQR for spatial extremes

## HCDN Folder
This folder contains details of all the stations in the network (as a .csv and .xlsx file) as well as <code>HCDN_annual_max.RData</code> which has data on the 702 locations for 1950 -- 2021.

## <code>download_annual_max.R</code>
Contains code for downloading the annual extremal streamflow data used in this study. Will generate all the variables needed for the <code>HCDN_annual_max.RData</code> file.

## <code>gen_data.R</code>
Contains utility functions for generating data from the process mixture model. Also contains functions for generating data for the GP and process mixture model simulation studies that are presented in the paper.

<code>gev_mle.R</code>
Generates the initial maximum likelihood estimates for the GEV parameters, the empirical variograms for each parameter's MLEs, and their corresponding plots that are presented in Section 5 of the paper.



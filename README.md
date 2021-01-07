# Campbell_Salmon_2021
 matlab codes to re-run results from Campbell and Goodwell, 2021
 
 Publication: Dynamic Drivers of Adult Chinook Salmon in the Columbia River Basin
 Nicholas Campbell and Allison Goodwell
 submitted to PNAS, January 2021
 
 DOI for this repository: tbd
 
 Description of contents:
 The excel files masterfile10yr_b.xlsx and masterfile35yr.xlsx are input files that contain formatted datasets from various sources as described in manuscript and supplemental information.  This are used in the main analysis files for the 10 year and 35 year studies, respectively:
 TIPNet_Salmon_10yr.m - run this file first to perform the 10 year analysis, of daily salmon counts, flow rates, water temperatures as network variables
 TIPNet_Salmon_35yr.m - run this file to perform 35 year analysis of annual variables
 TIPNet_Salmon_AnalyzeResults10yr.m  - run this file to produce plots and result matrices associated with 10 year study of daily data.  This file also outputs a text file that can be input into Circos software for circular network plots.
 TIPNet_Salmon_AnalyzeResults35yr.m - run this file to produce plots and results presented for 35 year study, in addition to correlation analysis
 
 function files:
 EntropyFun.m and EntropyFun35year.m - functions that are run in the main files to compute IT measures.  35year case only computes mutual information, not transfer entropy
 compute_pdfGUI.m - function to estimate pdf from data, used in EntropyFun codes
 compute_info_measures.m - function to compute IT measures based on pdf estimated with compute_pdfGUI
 
 output files:
 RESULTS_10yr.mat
 RESULTS_35yr.mat
 contain results of TIPNet_Salmon_10yr and TIPNet_Salmon_35yr programs.

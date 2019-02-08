# McL_Climate-Competition

This is the repository for chapter 2 of my dissertation investigating the interactive effects of precipitation variability and competition from invasive annual grasses on native annual forbs that vary in drought strategies. The project manuscript titled "Invasive species reduce the relative success of drought-avoiding plant species under a variable climate." Coauthors include Dr. Andrew Latimer and Dr. Susan Harrison.

The project seeks to answer how competition with a novel group alters the ability of native forbs to cope with a variable climate, namely drought or increased rainfall. 

## Analyses

The script used to generate the project's analysis is titled "Final-Paper-Analysis-112018" and lives in the "Scripts" folder. To run this script, the following post-processing data files are needed:

* __grass-cover.csv__: data on grass cover per subplot; not ultimately used
* __Marina-Treatment-30.csv__: experimental design
* __dem-data-16.csv__: census data from 2016
* __dem-data-17.csv__: census data from 2017
* __final-flo-seed.csv__: seed set and flowering data, both years
* __seed-carryover-plot.csv__: seed bag data, both years
* __final-traits-w.csv__: trait data
* __final-sla-13c.csv__: trait data, includes more species; not ultimately used

Intermediates (these files are saved and read in to save time)
* Boostrapped confidence intervals for figures
  + CI-1.Rdata
  + CI-2.Rdata
  + CI-3.Rdata
  + CI-4.Rdata
* Bootstrapped parameters for Lambda
  + BS-m.Rdata
  + BS-F.Rdata
  + BS-g.Rdata
* Final Lambda constructed from bootstrapped parameters
  + 20190206-lambda-sim.Rdata

## Data Processing
For data processing workflow from raw data files, see "data-processing.docx" in the Data folder

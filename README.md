# Resilience modelling
Code to generate the plots in the Supplemental Material in the paper "Framing resilience to manage complex environmental systems", 2024, Townend I, French J, Nicholls R, One Earth, Cell Publications.

## Licence
The code is provided as Open Source code (issued under a BSD 3-clause License).

## Requirements
ModelUI is written in Matlab(TM) and requires Matlab R2016a or later and has been tested up to Matlab R2024a. To run the code data needs to be downloaded from other other sources, as detailed below.

## Summary
The paper explores the concept of resilience within a system dynamics framework and links this to the evaluation and selection of adaptation pathways and transitions within the constraints of a ‘safe operating space’ and shows how such a resilience-based approach could be used operationally. The code in this repository was used to generate the figures in the Supplemental Material and Figure 7 of the main text.

![FigS8](https://github.com/user-attachments/assets/919b812e-cd36-4729-acf7-468706cb2350)

## Data sources needed to run code
O'Neill et al data: https://static-content.springer.com/esm/art%3A10.1038%2Fs41893-018-0021-4/MediaObjects/41893_2018_21_MOESM2_ESM.xlsx  
Fanning et al data: https://static-content.springer.com/esm/art%3A10.1038%2Fs41893-021-00799-z/MediaObjects/41893_2021_799_MOESM3_ESM.xlsx  
Population data: https://ourworldindata.org/grapher/population  

## Usage
_Excel spreadsheets_  
* 41893_2018_21_MOESM2_ESM.xlsx - data from O'Neill et al, 2018  
* 41893_2021_799_MOESM3_ESM.xlsx - data from Fanning et al, 2021  
* population.xlsx - data from  Gapminder - Population v7, 2022  

_Matlab *.mat files_  
.mat files are used in the main plotting function and the following functions are provided to generate these from the source excel files
* read_oneill_data.m - creates the ONeillData.mat file  
* read_fanning_data.m - creates the FanningData.mat file  
* read_pop_data.m - creates the PopData1900-2021.mat and PopData2011.mat files  

_Matlab *.m files_  
* resilience_modelling.m - function to analyse and plot O'Neill et al, 2018 data and Fanning et al, 2021 data  
* general_logistic.m - return a curve defined by the generalised logisitc equation  
* mcolor.m - select a default Matlab colour definition from table  
* target_marker.m - generate a circle and cross "target" symbol  

_Run analysis function from the command line_  
*>>* resilience_modelling;  
* when prompted select any mat file (it is just getting the path)  
* select plot type required and provide additional information as prompted.  

## See Also
Townend, I.H., French, J., Nicholls, R.J., (2024), Framing resilience to manage complex environmental systems, One Earth.
 
O'Neill, D.W., Fanning, A.L., Lamb, W.F., and Steinberger, J.K. (2018). A good life for all within planetary boundaries. Nature Sustainability 1, 88-95. https://doi.org/10.1038/s41893-018-0021-4.

Fanning, A.L., O’Neill, D.W., Hickel, J., and Roux, N. (2022). The social shortfall and ecological overshoot of nations. Nature Sustainability 5, 26-36. https://doi.org/10.1038/s41893-021-00799-z.

Gapminder (2022). Gapminder - Population v7 and other sources – with major processing by Our World in Data. https://ourworldindata.org/grapher/population.

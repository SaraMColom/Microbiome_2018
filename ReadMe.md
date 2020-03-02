# Read me

Data and code used to produce the results and figures of our microbiome project, where we studied the effect of plant-plant competition on the rhizosphere microboime of _Ipomoea purpurea_. We also explored how rhizosphere microboime variables correlate with various plant traits (i.e., root traits and plant fitness).

This repository contains the data and code used to produce the results figures for a paper currently in preparation by Sara Colom and Regina Baucom PhD (University of Michigan, Ann Arbor, USA. See the markdown file in the 'Pipeline' folder for work flow and detailed instructions. 

## R script
* **Genus_Level_MicrobiomeAnalysis.RMD** : R mardkown file with analysis for data aggregated at the Genus level.
* Files with 'Preliminary' in the prefix consist of preliminary data assesment, prior to aggregating data at Genus level.
* **RandomForest_Microbe_Fitness.RMD** preliminary analysis of using random forest to predict plant fitness with microbiome variables.
* **SimilarityTestsMicrobiome.R** preliminary assesment of analysis between microbiome similarity between competitors and focal plant fitness.

## DataSets
Cleaned up data used for preliminary and final data analysis.

**/Microbiome_Output**
* Files with the 'stability' _prefix_  consist of Mothur output on experimental samples.
    * Files with 'tax; suffix are taxonomy Files
        *Files with 'dist' suffix are distance measurements
        *Files with 'shared' suffix includes OTU counts
* MetaDataTest.csv: Meta Data on samples

### Mock_Output
Mothur output of 'mock' samples supplied by sequencing core, these were not integrated into our study.

## Figures
Preliminary and final figures.

## Pipeline
Detailed pipeline of analysis.

## Runing_Mothur
General instructions and batch script used to run mothur for this study.

* **stability.batch** (Batch script file)
* ReadME_mothur (markdown file with general instructions)

## CreatingSILVA_Database
Silva File creation data and instructions. Note it has subsfolders (e.g., Rscript, RawData, CleanData) for organizational purposes, similar to above. 



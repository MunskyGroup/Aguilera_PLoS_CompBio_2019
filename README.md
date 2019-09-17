Codes for "Computational design and interpretation of single-RNA translation experiments"
=======

Luis U. Aguilera, William Raymond, Zachary Fox, Michael May, Elliot Djokic, Tatsuya Morisaki, Timothy J. Stasevich, and Brian Munsky. <br/>

Department of Biochemistry and Molecular Biology, Institute of Genome Architecture and Function, Colorado State University, Fort Collins, CO, 80523, USA <br/>
Department of Chemical and Biological Engineering and School of Biomedical Engineering, Colorado State University, Fort Collins, CO, 80523, USA <br/>
World Research Hub Initiative, Institute of Innovative Research, Tokyo Institute of Technology, Yokohama, Kanagawa, 226-8503, Japan <br/>

For questions about the codes, please contact:  Luis.aguilera@colostate.edu and brian.munsky@colostate.edu <br/>

---

This repository contains the codes necessary to reproduce figures from the above manuscript. All codes are implemented in Matlab 2018. <br/>

## Code organization <br/>

The codes are organized on the following categories. <br/>

* **Computes the sequence properties.** <br/>

  masterFunction_Kymograph.m <br/>
  masterFunction_RibosomeOccupancy.m <br/>
  masterFunction_Collisions_RibosomeDensity.m <br/>

* **Computes the parameter estimation routines and parameter uncertainty.**

  masterFunction_ParameterEstimation.m <br/>
  masterFunction_ParameterUncertainty.m <br/>
  
 * **Computes the elongation methods (FRAP, FCS, and ROA) for the example genes (H2b, beta-actin, KDM5B).** <br/>
 
 Mean values <br/>
 
  masterFunction_Mean_AC.m <br/>
  masterFunction_Mean_FRAP.m <br/>
  masterFunction_Mean_ROA.m <br/>

Mean values with standard deviation <br/>

  masterFunction_SD_AC.m <br/>
  masterFunction_SD_FRAP.m <br/>
  masterFunction_SD_ROA.m <br/>
  
  * **Computes the elongation methods for a library of human genes. (FRAP, FCS, and ROA).** <br/>

  masterFunction_TableOfGenes.m <br/>
  masterFunction_GeneLibrary_REP_AC.m <br/>
  masterFunction_GeneLibrary_REP_FRAP.m <br/>
  masterFunction_GeneLibrary_REP_ROA.m <br/>

 * **Computes the codon optimization experiments.** <br/>

  masterFunction_CodonOptimization.m <br/>
  
  * **Computes the tRNA depletion experiments.** <br/>

  masterFunction_tRNA_Depletion.m <br/>

---

## Experimental data. <br/>

* **Autocorrelation data** <br/>
actBIntensityData.xls <br/>
h2bIntensityData.xls <br/>
kdmIntensityData.xls <br/>

* **Intensities data** <br/>
intensities_Bact.xlsx <br/>
intensities_H2B.xlsx <br/>
intensities_KDM5B.xlsx <br/>

---

## Gene Sequences. <br/>

H2B_withTags.txt <br/>
Bactin_withTags.txt <br/>
KDM5B_withTags.txt <br/>

The list of human genes used in this study is given in: <br/>

DB_Genes.xlsx <br/>

---  

## Code implementation.

The codes performing the optimization and uncertainty (masterFunction_ParameterEstimation.m and masterFunction_ParameterUncertainty.m), rely upon expensive computations that were performed on the W.M. Keck Compute Cluster. 

## Additional notes

The codes are compiled to run in Mac OS and Linux OS.  If there is a conflict with the user operative system. Try the following steps: 1) Open Matlab. 2) Open the directory. 3) double click on the file SSA_longNames.prj 4) Generate a MEX-file with Matlab Coder.
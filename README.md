# TinyBirds
Code repository for Orkney and Hedrick 2024

This repository contains commented codes to enable readers to reproduce the analyses conducted in Orkney and Hedrick, 2024. 

Required data to run these analyses were downloaded from the open-access repository of 
Navalon et al., 2022, and can be accessed at:
https://doi.org/10.1038/s41586-022-05372-y
https://osf.io/wjk3m/

The script 'Read_Navalon_12_04_2023.R' will be available in this repository and must be run to 
produce an initial '.RData' object. 
Thereafter, the initial script 'loader.R' must be run to produce interrogable data objects for subsequent analysis.


Fig_var_stru_01_12_2024.R
>> Figure 1

Salami_slicer_01_11_2024.R
>> Figure 2 and 3
(this script is more computationally intensive and it may take some time to run)

Dm_plots_01_11_2024.R
>> Figure 4

Description of analysis: 
The scripts here interrogate a dataset of skeletal proportions/ bone sizes across a broad sample of 228 extant bird species, 
which vary in body mass by almost a factor of 10,000. 

Figure 1 produces a visualisation of how body mass is distributed across avian phylogeny (subplot a)
and demonstrates that the skeletal proportions of larger birds tend to have an increased variability (which is consistent with expectations; subplot b). 
We show that some modular regions of the avian body plan, such as the leg, exhibit a proximo-distal gradient, whereby the most distal bones are
the most variable, which is consistent with expectations. The avian wing a conspicuous exception to this rule. 

Figure 2 is a visualisation that shows how average levels of evolutionary integration- the propensity for bone proportions to evolve either in concert
rather than independently of one another- change within the head, wing, trunk and leg, over avian body mass. 
We show that there is no substantial change in integration within the head or leg between small and large birds. 
However, the wing is substantially more integrated in large birds, and the trunk is substantially more integrated in small birds. 

Figure 3 is a more detailed visualisation of this change, arranged to show the broad sense (increase or decrease of integration with body mass)
between all 78 pairwise combinations of bone sizes within the avian skeleton that we consier. Modular body regions (head, wing, trunk and leg)
are identified on this plot to facilitate viewing. This figure provides the additional insight that integration between the wing and trunk 
appears to be higher in smaller birds. 

Figure 4 is a visualisation that shows how individual taxa differ from major axes of trait integration between the carpometacarpus and humerus (wing), 
cranium and mandible (head), scapula and sternum (trunk) and the carpometacarpus and sternum (a proxy for integration between the wing and trunk).
This allows the reader to inspect whether distinctive structure in macroevolutionary trends across birds shapes changes in the strength of evolutionary
integration between structures. A significance assessor (inset flame-plots) is produced by sub-sampling the dataset many times and 
conducting statistical tests for unequal scatter of 'distance from trait integration' as a function of body mass. This is contrasted against a null distribution
in a permuted dataset. The negative logarithm of p-values of the test outcomes is plotted for both the real and null distributions. Significant 
results will show high values on this axis, with minimal overlap with the null distribution. 

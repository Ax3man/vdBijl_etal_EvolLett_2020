# vdBijl_etal_EvolLett_2020

[![DOI](https://zenodo.org/badge/296743020.svg)](https://zenodo.org/badge/latestdoi/296743020)


This repository accompanies the paper entitled "Butterfly dichromatism primarily evolved via Darwin’s, not Wallace’s, model", by Van der Bijl et al.

If you have questions or comments, file an issue or mail me at wouter@zoology.ubc.ca.

It contains the following data files:
1. centroids.csv: This contains data on average colors males of females from European butterfly species, based on drawings from the Tolman field guide. It is used in all analyses. It contains the follow columns:
- Species: Species names as used in the Tolmand field guide.
- L_mean_females: average across the L axis of the Lab color space for the female drawing.
- a_mean_females: average across the a axis of the Lab color space for the female drawing.
- b_mean_females: average across the b axis of the Lab color space for the female drawing.
- L_mean_males: average across the L axis of the Lab color space for the female drawing.
- a_mean_males: average across the a axis of the Lab color space for the female drawing.
- b_mean_males: average across the b axis of the Lab color space for the female drawing.

2. centroids_photos.csv: This contains data on average colors males of females from a subset of European butterfly species, based on photos of museum specimens. It is used to validate the data from drawings. It contains the follow columns:
- Tolman: Species names following the Tolman field guide (matches centroids.csv).
- Sex: male or female.
- L: average across the L axis of the Lab color space.
- a: average across the a axis of the Lab color space.
- b: average across the b axis of the Lab color space.

3. tree.nex: Phylogeny used in this study in Nexus format. Uses Tolman species names.

4. species_translation_guide.csv: Contains both the modern names from the Wiemers et al. phylogeny, as well as the matched names from the Tolman field guide. Can be used to match species name to other datasets.

It contains the following scripts:
1. compare_evol_rates.R: This contains the analysis to compare the evolutionary rates between males and females, and to create figure S2.

2. compare_evol_rates_with_dichromatism.R: This contains the analysis of the difference in male and female evolutionary rates in relation to the level of dichromatism along the same branches, and to create figure 3. 

3. compare_evol_angles_with_dichromatism_angles.R: This contains the analysis to to analyze male and female evolutionary vectors, and compare them to the direction of dichromatism. This allows us to relate changes in dichromatism to male or female specific color evolution. It also creates figure 4.

4. color_verification.R: This script contains the comparison between the drawings and the photos.

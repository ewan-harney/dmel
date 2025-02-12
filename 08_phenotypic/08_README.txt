Within this 08_interchromate_phenotypic folder are 7 text files used in phenotypic analysis (see InterChromaTE_phenotypic_analysis_final_250203.R script).

- ctmax_data_akaa_manz.txt
- Phenotypic_F3_eggs-pupa-adults-viab.txt
- Phenotypic_F6_eggs-pupa-adults-viab.txt
- f3_all_age-to-eclosion.txt
- f3_all_age-to-pupation.txt
- f6_all_age-to-eclosion.txt
- f6_all_age-to-pupation.txt


The first file provides information about the critical thermal maximum (CTMax) for flies from Akaa and Manz subjected to a heat-ramping experiment. CTMax was the temperature at which the individual fly was no longer responsive to agitation, but the table also contains CTbax, which is the temperture at which the individual first showed signs of fainting due to heat stress.

'Phenotypic_F3_eggs-pupa-adults-viab.txt' provides information about the number of eggs, pupae and imagoes counted in every tube following F3 heat shock. Individual females were allowed to lay for four 48 hour periods, leading to four separate cohorts (as indicated by the 'day' column). EMiscount is a record of whether there is a potential miscount of eggs (if the number of pupae was greater than the number of eggs counted, the difference is recorded in this column). PMiscount is a record of whether there was a miscount of Pupae (if the number of adults was greater than the number of pupae, the difference is recorded here).

'Phenotypic_F6_eggs-pupa-adults-viab.txt' provides information about the number of eggs, pupae and imagoes counted in every tube for offspring of F6 females, with the same data layout as for 'Phenotypic_F3_eggs-pupa-adults-viab.txt'

the '_age-to-eclosion.txt' and '_age-to-pupation.txt' files provide the measurements of what day all individual offspring of F3 or F6 flies were seen to pupate or eclose, with the population density ('Density' = adult density, 'Density_pup' = pupal density) of the tube and cohort recorded for inclusion as random effects in the statistical analysis. The column 'Eclose' is a binary variable signifying successful eclosion (1) or not (0).
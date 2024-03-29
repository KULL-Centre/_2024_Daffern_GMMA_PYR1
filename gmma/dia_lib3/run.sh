#!/bin/bash

# [1] "Found 1 groups of fully co-occuring substitutions"
#         i resi taa obs signal note gmma gmma_reparam
# V158M 153  158   M  15     NA       use            1
# G159F 154  159   F  15     NA       use            1
# V160A 155  160   A  15     NA       use            1
# "T156A:V158M:G159F:V160A:M178R",81,11,-3.19498770812703
# "T156A:V158M:G159F:V160A:A179T",36,9,-2.30992347497169
# "T156A:V158M:G159F:V160A:A179V",49,9,-2.74432629911746
# "L140A:T156A:V158M:G159F:V160A:M178R",73,18,-2.38392405641546
# "L140A:T156A:V158M:G159F:V160A:A179T",64,11,-2.8598035165374
# "L140A:T156A:V158M:G159F:V160A:A179V",51,8,-2.95291292092888
# "E141S:T156A:V158M:G159F:V160A:M178R",23,6,-2.20000578289365
# "V145A:T156A:V158M:G159F:V160A:M178R",46,7,-2.97698705590774
# "V145A:T156A:V158M:G159F:V160A:A179T",24,0,-5.06625439400482
# "V145A:T156A:V158M:G159F:V160A:A179V",21,0,-4.8818298228674
# "T156A:V158M:G159F:V160A:L171F:A179V",7,13,0.384956717827506
# "K131R:V145A:T156A:V158M:G159F:V160A:M178R",32,0,-5.46679232358855
# "E132D:V145A:T156A:V158M:G159F:V160A:M178R",32,0,-5.46679232358855
# "E132D:V145A:T156A:V158M:G159F:V160A:A179T",27,0,-5.2297531262877
# "E132D:V145A:T156A:V158M:G159F:V160A:A179V",25,6,-2.31548300031359

Rscript gmma01_structure.r ../../enrichment/pyr1_dia_lib3.csv > gmma01_structure.out
Rscript ../gmma02_fit_individual.r > gmma02_fit_individual.out
Rscript ../gmma03_graph.r gmma_structured.rda > gmma03_graph.out
Rscript ../gmma04_global_estimation.r > gmma04_global_estimation.out
Rscript ../gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out


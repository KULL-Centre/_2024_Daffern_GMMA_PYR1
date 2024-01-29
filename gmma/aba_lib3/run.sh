#!/bin/bash

# [1] "Found 1 groups of fully co-occuring substitutions"
#        i resi taa obs signal note gmma gmma_reparam
# M158V 36  158   V   3     NA       use            1
# F159G 37  159   G   3     NA       use            1
# A160V 38  160   V   3     NA       use            1
# "F108Y:M158V:F159G:A160V:M178R",13,9,-1.88048962668782
# "F108Y:M158V:F159G:A160V:A179T",9,13,-0.909635972347336
# "F108Y:M158V:F159G:A160V:A179V",13,18,-0.954490208131596

Rscript gmma01_structure.r ../../enrichment/pyr1_aba_lib3.csv > gmma01_structure.out
Rscript ../gmma02_fit_individual.r > gmma02_fit_individual.out
Rscript ../gmma03_graph.r gmma_structured.rda > gmma03_graph.out
Rscript ../gmma04_global_estimation.r > gmma04_global_estimation.out
Rscript ../gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out


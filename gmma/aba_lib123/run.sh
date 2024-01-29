#!/bin/bash

## Based on graph analysis output:
# [1] "Found 1 groups of fully co-occuring substitutions"
#         i resi taa obs signal note gmma gmma_reparam
# M158V 254  158   V   3     NA       use            1
# F159G 255  159   G   3     NA       use            1
# A160V 256  160   V   3     NA       use            1

## Remove fully co-occurring substitutions from data since these are few
# "F108Y:M158V:F159G:A160V:M178R",0,0,22,-1.88048962668782
# "F108Y:M158V:F159G:A160V:A179T",0,0,22,-0.909635972347336
# "F108Y:M158V:F159G:A160V:A179V",0,0,31,-0.954490208131596

Rscript gmma01_structure.r ../../enrichment/pyr1_aba_lib123.csv > gmma01_structure.out
Rscript ../gmma02_fit_individual.r > gmma02_fit_individual.out
Rscript ../gmma03_graph.r gmma_structured.rda > gmma03_graph.out
Rscript ../gmma04_global_estimation.r > gmma04_global_estimation.out
Rscript ../gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out


#!/bin/bash

Rscript gmma01_structure.r ../../enrichment/pyr1_dia_Y81V_D80E_lib12.csv > gmma01_structure.out
Rscript ../gmma02_fit_individual.r > gmma02_fit_individual.out
Rscript ../gmma03_graph.r gmma_structured.rda > gmma03_graph.out
Rscript ../gmma04_global_estimation.r > gmma04_global_estimation.out
Rscript ../gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out


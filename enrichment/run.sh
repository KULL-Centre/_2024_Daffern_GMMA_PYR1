#!/bin/bash

Rscript collect_aba.r               > collect_aba.out 
Rscript collect_diazi.r             > collect_diazi.out
Rscript collect_diazi_Y81V.r        > collect_diazi_Y81V.out
Rscript collect_diazi_Y81V_D80E.r   > collect_diazi_Y81V_D80E.out

Rscript numbers.r > numbers.out

# rm *.png *.out pyr*.csv *.rda

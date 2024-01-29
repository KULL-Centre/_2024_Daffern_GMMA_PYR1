options(width=160) 

##
## Number of D80E and Y81V subst in Diazi libs 1 and 2
##
dia1 = read.csv("pyr1_dia_lib1.csv")
dia2 = read.csv("pyr1_dia_lib2.csv")

n1_Y81V = sum(grepl("Y81V", dia1$var))
n2_Y81V = sum(grepl("Y81V", dia2$var))
print(sprintf("Y81V is in %d of %d (%.1f%%) variants of library 1", n1_Y81V, nrow(dia1), n1_Y81V/nrow(dia1)*100))
print(sprintf("Y81V is in %d of %d (%.1f%%) variants of library 2", n2_Y81V, nrow(dia2), n2_Y81V/nrow(dia2)*100))

n1_D80E = sum(grepl("D80E", dia1$var))
n2_D80E = sum(grepl("D80E", dia2$var))
print(sprintf("D80E is in %d of %d (%.1f%%) variants of library 1", n1_D80E, nrow(dia1), n1_D80E/nrow(dia1)*100))
print(sprintf("D80E is in %d of %d (%.1f%%) variants of library 2", n2_D80E, nrow(dia2), n2_D80E/nrow(dia2)*100))

n1_D80E_Y81V = sum(grepl("D80E:Y81V", dia1$var))
n2_D80E_Y81V = sum(grepl("D80E:Y81V", dia2$var))
print(sprintf("D80E:Y81V is in %d of %d (%.1f%%) variants of library 1", n1_D80E_Y81V, nrow(dia1), n1_D80E_Y81V/nrow(dia1)*100))
print(sprintf("D80E:Y81V is in %d of %d (%.1f%%) variants of library 2", n2_D80E_Y81V, nrow(dia2), n2_D80E_Y81V/nrow(dia2)*100))





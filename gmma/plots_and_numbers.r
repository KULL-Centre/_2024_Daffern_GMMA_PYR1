options(width=160) 

# read and process csv of expected substitutions
es_raw = read.table("expected_substitutions.csv", sep="\t", comment.char="#")
target_mut_list = lapply(es_raw$V3, function(s){ strsplit(s,"/")[[1]] })
expected_subst = unname(unlist( mapply(paste0, paste0(es_raw$V2,es_raw$V1), target_mut_list) ))
es = data.frame(s=expected_subst,
                expected=unlist(rep(es_raw$V4, sapply(target_mut_list, length))),
                source=unlist(rep(es_raw$V5, sapply(target_mut_list, length))))
stopifnot(length(expected_subst) == length(unique(expected_subst)))

es$wt = substr(es$s, 1, 1)
es$pos = as.numeric(substr(es$s, 2, nchar(es$s)-1))
es$mut = substr(es$s, nchar(es$s), nchar(es$s))
es = es[order(es$pos,es$mut),]


##
## Compare ABA individual libs with combined lib
##

a1 = read.table("aba_lib1/subst.csv", header=T, sep=";")
a2 = read.table("aba_lib2/subst.csv", header=T, sep=";")
a3 = read.table("aba_lib3/subst.csv", header=T, sep=";")
ac = read.table("aba_lib123/subst.csv", header=T, sep=";")

# estimated substitutions
s123 = ac[which(! is.na(ac$rank)),"s"]
s1 = a1[which(! is.na(a1$rank)),"s"]
s2 = a2[which(! is.na(a2$rank)),"s"]
s3 = a3[which(! is.na(a3$rank)),"s"]

# are all subst estimated in individual libs also estimated in the combined?
all(s1 %in% s123)  # TRUE
all(s2 %in% s123)  # TRUE
all(s3 %in% s123)  # TRUE

# are all subst estimated in the combined lib also estimated in the individual
ac[which(ac$s %in% a1$s & (! ac$s %in% s1) & (! is.na(ac$rank))), "s"]   # "F71V"
ac[which(ac$s %in% a2$s & (! ac$s %in% s2) & (! is.na(ac$rank))), "s"]   # "M73Q" "V75E"
ac[which(ac$s %in% a3$s & (! ac$s %in% s3) & (! is.na(ac$rank))), "s"]   # "K121H"
# these are all from the overlap regions

s123_indi = union(s1,union(s2,s3))
print(sprintf("ABA %d (combined) and %d (individual) confident effects of %d substitutions considered",length(s123),length(s123_indi),nrow(ac)))
# any effects estimated in individual libs nut not in combined or visa verse
print(sprintf("Confident subst in combined but not individual: %s",setdiff(s123_indi, s123)))
print(sprintf("Confident subst in individual but not combined: %s",setdiff(s123, s123_indi)))

print("Effects in combined ABA analysis")
print(table(ac$eff))

# check that the same subst are considered in individual and combined lib
setdiff(union(a1$s,union(a2$s,a3$s)), ac$s)
setdiff(ac$s, union(a1$s,union(a2$s,a3$s)))

aba = data.frame(s=s123, dfc=ac[match(s123,ac$s),"ddG_glob"], df1=NA, df1_sd=NA, df2=NA, df2_sd=NA, df3=NA, df3_sd=NA)
aba[match(s1,aba$s),"df1"] = a1[match(s1,a1$s),"ddG_glob"]
aba[match(s1,aba$s),"df1_sd"] = a1[match(s1,a1$s),"stderr_meas"]
aba[match(s2,aba$s),"df2"] = a2[match(s2,a2$s),"ddG_glob"]
aba[match(s2,aba$s),"df2_sd"] = a2[match(s2,a2$s),"stderr_meas"]
aba[match(s3,aba$s),"df3"] = a3[match(s3,a3$s),"ddG_glob"]
aba[match(s3,aba$s),"df3_sd"] = a3[match(s3,a3$s),"stderr_meas"]

jpeg("aba_combined_vs_lidividual.jpg", width=9, height=9, units="cm", res=300, pointsize=9, quality=90)
par(mar=c(5,4,3,2)+.1)
plot(0,0, col=0, xlim=c(-0.8,3.5), ylim=c(-0.8,3.5),
     xlab=expression(Delta*italic(F)~"combined library"), ylab=expression(Delta*italic(F)~"individual library"))
abline(0,1)
points(aba$dfc, aba$df1, pch=20, col=2)
points(aba$dfc, aba$df2, pch=20, col=3)
points(aba$dfc, aba$df3, pch=20, col=4)
legend("topleft", c("Library 1","Library 2","Library 3"), pch=20, col=c(2,3,4))
dev.off()


##
## Compare Diazi individual libs with combined lib
##

d1 = read.table("dia_lib1_reref/subst.csv", header=T, sep=";")
d2 = read.table("dia_lib2_reref/subst.csv", header=T, sep=";")
d3 = read.table("dia_lib3/subst.csv", header=T, sep=";")
dc12 = read.table("dia_lib12_reref/subst.csv", header=T, sep=";")
dc13 = read.table("dia_lib13/subst.csv", header=T, sep=";")

# estimated substitutions
s12 = dc12[which(! is.na(dc12$rank)),"s"]
s1 = d1[which(! is.na(d1$rank)),"s"]
s2 = d2[which(! is.na(d2$rank)),"s"]

# are all subst estimated in individual libs also estimated in the combined?
all(s1 %in% s12)  # TRUE
all(s2 %in% s12)  # TRUE

# are all subst estimated in the combined lib also estimated in the individual
dc12[which(dc12$s %in% d1$s & (! dc12$s %in% s1) & (! is.na(dc12$rank))), "s"]   # "M87L"
dc12[which(dc12$s %in% d2$s & (! dc12$s %in% s2) & (! is.na(dc12$rank))), "s"]   # 
# these are all from the overlap regions

s12_indi = union(s1,s2)
print(sprintf("Diazi:D80E:Y81V %d (combined) and %d (individual) lib 1+2 confident effects of %d substitutions considered",length(s12),length(s12_indi),nrow(dc12)))
# any effects estimated in individual libs nut not in combined or visa verse
print(sprintf("Confident subst in combined but not individual: %s",setdiff(s12_indi, s12)))
print(sprintf("Confident subst in individual but not combined: %s",setdiff(s12, s12_indi)))

# check that the same subst are considered in individual and combined lib
setdiff(union(d1$s,d2$s), dc12$s)
setdiff(dc12$s, union(d1$s,d2$s))

print("Effects in combined Diazi tile 1+2 analysis")
print(table(dc12$eff))

dia = data.frame(s=s12, dfc=dc12[match(s12,dc12$s),"ddG_glob"], df1=NA, df1_sd=NA, df2=NA, df2_sd=NA, df3=NA, df3_sd=NA)
dia[match(s1,dia$s),"df1"] = d1[match(s1,d1$s),"ddG_glob"]
dia[match(s1,dia$s),"df1_sd"] = d1[match(s1,d1$s),"stderr_meas"]
dia[match(s2,dia$s),"df2"] = d2[match(s2,d2$s),"ddG_glob"]
dia[match(s2,dia$s),"df2_sd"] = d2[match(s2,d2$s),"stderr_meas"]

jpeg("dia_combined_vs_lidividual.jpg", width=9, height=9, units="cm", res=300, pointsize=9, quality=90)
par(mar=c(5,4,3,2)+.1)
plot(0,0, col=0, xlim=c(-0.8,3.5), ylim=c(-0.8,3.5),
     xlab=expression(Delta*italic(F)~"combined library"), ylab=expression(Delta*italic(F)~"individual library"))
abline(0,1)
points(dia$dfc, dia$df1, pch=20, col=2)
points(dia$dfc, dia$df2, pch=20, col=3)
legend("topleft", c("Library 1","Library 2"), pch=20, col=c(2,3))
dev.off()


# plot correlation between lib 3 subst in lib 3 and the combined lib 1+3
quartz(width=10, height=5)
par(mfcol=c(1,2))
lim = c(-1,1.5)
common_subst13 = intersect(d3[which(! is.na(d3$rank)),"s"], dc13[which(! is.na(dc13$rank)),"s"])
x = d3[match(common_subst13,d3$s),"ddG_glob"]
y = dc13[match(common_subst13,dc13$s),"ddG_glob"]
rp = cor(x, y, method="pearson", use="complete.obs")
plot(x, y, xlim=lim, ylim=lim, xlab="Lib 3", ylab="Lib 1+3 combined", col="white",
     main=sprintf("%d common subst. Pearson %.2f",length(common_subst13),rp))
abline(0,1)
fit = lm(y ~ x)
abline(coef(fit), lty=2)

i = which(common_subst13 %in% d3$s & common_subst13 %in% expected_subst)
points(x[i], y[i], pch=16, col="orange")
points(x, y)
legend("topleft", c("Library 3","Unexpected"), pch=c(16,1), col=c("orange","black"))


# since these fully correlate, use this to align lib 1+3 to the re-referenced lib 1+2
common_subst12_13 = intersect(dc12[which(! is.na(dc12$rank)),"s"], dc13[which(! is.na(dc13$rank)),"s"])
x = dc12[match(common_subst12_13,dc12$s),"ddG_glob"]
y = dc13[match(common_subst12_13,dc13$s),"ddG_glob"]
rp = cor(x, y, method="pearson", use="complete.obs")
plot(x, y, xlim=lim, ylim=lim, xlab="Lib 1+2 combined", ylab="Lib 1+3 combined", col="white",
     main=sprintf("%d common subst. Pearson %.2f",length(common_subst12_13),rp))
abline(0,1)
fit = lm(y ~ x)
abline(coef(fit), lty=2)

i = which(common_subst12_13 %in% d1$s & common_subst12_13 %in% expected_subst)
points(x[i], y[i], pch=16, col="orange")
i = which(common_subst12_13 %in% d3$s & common_subst12_13 %in% expected_subst)
points(x[i], y[i], pch=16, col="steelblue")
points(x, y)

legend("topleft", c("Library 1","Library 3","Unexpected"), pch=c(16,16,1), col=c("orange","steelblue","black"))

quartz.save("lib13_scale_correlations.png", type="png")


##
## Fig. 3: Scatter plot of everything combined 
##

confident_ac = ac[which(! is.na(ac$rank)),"s"]
confident_d3 = d3[which(! is.na(d3$rank)),"s"]
confident_dc12 = dc12[which(! is.na(dc12$rank)),"s"]
confident_dc13 = dc13[which(! is.na(dc13$rank)),"s"]
print(sprintf("All confident d3 are also confident in the combined diazi lib 1+3: %s",all(confident_d3 %in% confident_dc13)))
print(sprintf("Confident in diazi combined 1+3 not in combined diazi:D80E:Y81V lib 1+2 plus diazi lib 3: %s",paste(setdiff(confident_dc13, union(confident_dc12,confident_d3)),collapse=",")))

confident_dc = union(confident_dc12, confident_dc13)
confident_common = intersect(confident_ac, confident_dc)
adc = data.frame(s=confident_common)

ia = match(confident_common,ac$s)
adc$aba_ddG = ac[ia,"ddG_glob"]
adc$aba_err = ac[ia,"stderr_meas"]
adc$aba_eff = ac[ia,"eff"]

id12 = match(confident_common,dc12$s)
adc$dia12_ddG = dc12[id12,"ddG_glob"]
adc$dia12_err = dc12[id12,"stderr_meas"]
adc$dia12_eff = dc12[id12,"eff"]

id13 = match(confident_common,dc13$s)
adc$dia13_ddG = dc13[id13,"ddG_glob"]
adc$dia13_err = dc13[id13,"stderr_meas"]
adc$dia13_eff = dc13[id13,"eff"]

id3 = match(confident_common,d3$s)
adc$dia3_ddG = d3[id3,"ddG_glob"]
adc$dia3_err = d3[id3,"stderr_meas"]
adc$dia3_eff = d3[id3,"eff"]

# for combined diazi 1+2 and 3, use most stabilizing effect where both dc12 and d3 are confident
use_12_or_3 = apply(adc[,c("dia12_ddG","dia3_ddG")], MARGIN=1, function(v){ if(all(is.na(v))) NA else which.min(v)})
adc$dia_ddG = NA
adc$dia_err = NA
adc$dia_eff = NA
i = which(use_12_or_3==1)
adc[i,"dia_ddG"] = adc[i,"dia12_ddG"]
adc[i,"dia_err"] = adc[i,"dia12_err"]
adc[i,"dia_eff"] = adc[i,"dia12_eff"]
i = which(use_12_or_3==2)
adc[i,"dia_ddG"] = adc[i,"dia3_ddG"]
adc[i,"dia_err"] = adc[i,"dia3_err"]
adc[i,"dia_eff"] = adc[i,"dia3_eff"]

# stats on combined analysis
adc$expect = es[match(adc$s,es$s),"expected"]
print(sprintf("Common confident analysis has %d unexpected and %d of %d (%.1f%%) expected subst",
              sum(! adc$s %in% es$s), nrow(adc), nrow(es), nrow(adc)/nrow(es)*100))

aba_conf = ac[which(! is.na(ac$rank)),"s"]
print(sprintf("Combined ABA lib 1+2+3 analysis has %d unexpected and %d of %d (%.1f%%) expected subst",
              sum(! aba_conf %in% es$s), sum(aba_conf %in% es$s), nrow(es), sum(aba_conf %in% es$s)/nrow(es)*100 ))

dia_conf = union(dc12[which(! is.na(dc12$rank)),"s"],  d3[which(! is.na(d3$rank)),"s"])
print(sprintf("Combined Diazi lib 1+2+3 analysis has %d unexpected and %d of %d (%.1f%%) expected subst",
              sum(! dia_conf %in% es$s), sum(dia_conf %in% es$s), nrow(es), sum(dia_conf %in% es$s)/nrow(es)*100 ))

# all stabilizing
stab_subst_indi = unique(c(ac[which(ac$eff=="stab"),"s"], dc12[which(dc12$eff=="stab"),"s"], d3[which(d3$eff=="stab"),"s"]))
all(stab_subst_indi %in% es$s)
stab_pos_indi = unique(c(ac[which(ac$eff=="stab"),"resi"], dc12[which(dc12$eff=="stab"),"resi"], d3[which(d3$eff=="stab"),"resi"]))
print(sprintf("From individual analyses, %d stabilizing subst at %d positions",length(stab_subst_indi),length(stab_pos_indi)))

stab_subst_any = adc[which(adc$aba_eff=="stab" | adc$dia_eff=="stab"),"s"] 
stab_pos_any = unique(substr(stab_subst_any, 2, nchar(stab_subst_any)-1))
print(sprintf("For combined analyses of any sensors, there are %d stabilizing subst at %d positions",length(stab_subst_any),length(stab_pos_any)))
print(paste(stab_pos_any,collapse="+"))

stab_subst_both = adc[which(adc$aba_eff=="stab" & adc$dia_eff=="stab"),"s"] 
stab_pos_both = unique(substr(stab_subst_both, 2, nchar(stab_subst_both)-1))
print(sprintf("For combined analyses of both sensors, there are %d stabilizing subst at %d positions: %s",length(stab_subst_both),length(stab_pos_both),paste(stab_subst_both,collapse=",")))
print(paste(stab_pos_any,collapse="+"))


# Single plot of combined analysis
# quartz(width=6, height=6)
jpeg("fig3_gmma_combined1_350.jpg", width=8, height=8, units="cm", res=350, pointsize=10)
par(mar=c(5,4,1,1)+.1)
x = adc$aba_ddG
y = adc$dia_ddG
size = .3/apply(adc[,c("aba_err","dia_err")], MARGIN=1, function(v){ sqrt(sum(v**2, na.rm=T)) })
lim = c(-0.6,0.07)
plot(x, y, col=0, xlim=lim, ylim=lim, xlab="ABA GMMA effects", ylab="Diazi GMMA effects")
abline(h=0,v=0, lty=2)

# arrows(x0=adc[i,"aba_ddG"], y0=adc[i,"dia_ddG"]-adc[i,"dia_err"], y1=adc[i,"dia_ddG"]+adc[i,"dia_err"], code=3, angle=90, length=.01)
# arrows(x0=adc[i,"aba_ddG"]-adc[i,"aba_err"], x1=adc[i,"aba_ddG"]+adc[i,"aba_err"], y0=adc[i,"dia_ddG"], code=3, angle=90, length=.01)

# white box behind label
lab = "Q24H"
i = which(adc$s == lab)
text(x[i], y[i]+.00, paste(rep("\U2588", nchar(lab)), collapse=""), cex=0.8, col="white", offset=0.2, pos=2)
lab = "D26S"
i = which(adc$s == lab)
text(x[i], y[i]+.00, paste(rep("\U2588", nchar(lab)), collapse=""), cex=0.7, col="white", offset=0.6, pos=1)

i = which(adc$s %in% c("S11E","Q24H","D26N"))                                        # small points
text(x[i], y[i]+.00, adc[i,"s"], cex=0.7, offset=0.25, pos=c(2,2,2))
i = which(adc$s %in% c("N15D","F20Y","D26G","S29Q","E43D","R104H","R134K","D146N"))  # medium points
text(x[i], y[i]+.00, adc[i,"s"], cex=0.7, offset=0.6, pos=c(2,1,3,2,2,2,1,2))
i = which(adc$s == "E4G")
text(x[i], y[i]+.00, adc[i,"s"], cex=0.7, offset=0.7, pos=2)
i = which(adc$s %in% c("E102D","T118R"))                                             # big points
text(x[i], y[i]+.00, adc[i,"s"], cex=0.7, offset=1.0, pos=c(1,2))

i = which(adc$s %in% c("R104H","E43D","T118R","S29Q","R134K"))
points(x[i], y[i], cex=size[i]**2, pch=16, col="orange")

points(x, y, cex=size**2)

i = which(adc$s == "H34L")
text(x[i], y[i]+.00, paste(rep("\U2588", 4), collapse=""), cex=0.65, col="white", offset=0.7, pos=3)
text(x[i], y[i]+.00, adc[i,"s"], cex=0.7, offset=0.65, pos=3)

# quartz.save("gmma_combined1.png", type="png")
dev.off()


# Two plots
quartz(width=10, height=5)
par(mfcol=c(1,2))
lim = c(-0.7,0.1)

size12 = .3/apply(adc[,c("aba_err","dia12_err")], MARGIN=1, function(v){ sqrt(sum(v**2, na.rm=T)) })
plot(adc$dia12_ddG, adc$aba_ddG, xlim=lim, ylim=lim, cex=size12**2, lwd=2,
     xlab="Diazi:D80E:Y81V lib 1+2 GMMA effects", ylab="ABA lib 1+2+3 GMMA effects", main=sprintf("%d common subst.",sum(! is.na(adc$dia12_ddG) & (! is.na(adc$aba_ddG)))))
abline(h=0,v=0, lty=2)
i = which(adc$s %in% c("R104H","E43D","T118R","S29Q","R134K"))
points(adc[i,"dia12_ddG"], adc[i,"aba_ddG"], cex=size12[i]**2, pch=16, col="orange")
text(adc$dia12_ddG, adc$aba_ddG, adc$s, pos=2, offset=0.7, cex=.8)

size13 = .3/apply(adc[,c("aba_err","dia13_err")], MARGIN=1, function(v){ sqrt(sum(v**2, na.rm=T)) })
plot(adc$dia13_ddG, adc$aba_ddG, xlim=lim, ylim=lim, cex=size13**2, lwd=2,
     xlab="Diazi lib 1+3 GMMA effects", ylab="ABA lib 1+2+3 GMMA effects", main=sprintf("%d common subst.",sum(! is.na(adc$dia13_ddG) & (! is.na(adc$aba_ddG)))))
abline(h=0,v=0, lty=2)
i = which(adc$s %in% c("R104H","E43D","T118R","S29Q","R134K"))
points(adc[i,"dia13_ddG"], adc[i,"aba_ddG"], cex=size13[i]**2, pch=16, col="orange")
text(adc$dia13_ddG, adc$aba_ddG, adc$s, pos=2, offset=0.7, cex=.8)

quartz.save("gmma_combined2.png", type="png")


# Three plots, one per tile
common_conf = function(df1,df2) {
    df1_conf = df1[which(! is.na(df1$rank)),"s"]
    df2_conf = df2[which(! is.na(df2$rank)),"s"]
    common_conf = intersect(df1_conf, df2_conf)
    print(sprintf("Confident in gmma1 but not in common 1: %s", paste(setdiff(df1_conf, common_conf),collapse=",")))
    print(sprintf("Confident in gmma2 but not in common 1: %s", paste(setdiff(df2_conf, common_conf),collapse=",")))

    df_ret = data.frame(s=common_conf)
    i1 = match(common_conf,df1$s)
    df_ret$ddG1 = df1[i1,"ddG_glob"]
    df_ret$err1 = df1[i1,"stderr_meas"]
    i2 = match(common_conf,df2$s)
    df_ret$ddG2 = df2[i2,"ddG_glob"]
    df_ret$err2 = df2[i2,"stderr_meas"]
    return(df_ret)
}

quartz(width=12, height=4.2)
par(mfcol=c(1,3))
lim = c(-0.7,0.1)

# tile 1
t1 = common_conf(a1,d1)

size1 = .4/apply(t1[,c("err1","err2")], MARGIN=1, function(v){ sqrt(sum(v**2, na.rm=T)) })
plot(0,0,col=0, xlim=lim, ylim=lim, xlab="ABA lib 1 GMMA effects", ylab="Diazi:D80E:Y81V lib 1 GMMA effects", main=sprintf("%d common subst.",nrow(t1)))
abline(h=0,v=0, lty=2)
i = which(t1$s %in% c("R104H","E43D","T118R","S29Q","R134K"))
points(t1[i,"ddG1"], t1[i,"ddG2"], cex=size1[i]**2, pch=16, col="orange")
points(t1$ddG1, t1$ddG2, cex=size1**2)
# text(t1$ddG1, t1$ddG2, t1$s, pos=2, offset=0.7, cex=.8)
i = which(t1$s %in% c("S11A","S11E","F20Y","Q24H","D26G","D26N","D26S","S29Q"))
text(t1[i,"ddG1"], t1[i,"ddG2"], t1[i,"s"], pos=c(4,2,2,2,2,2,3,2), offset=0.7, cex=.8)
i = which(t1$s %in% c("E4G","N15D","N15S","H34L","P42A","E43D"))
text(t1[i,"ddG1"], t1[i,"ddG2"], t1[i,"s"], pos=2, offset=1.1, cex=.8)

# tile 2
t2 = common_conf(a2,d2)

size2 = .4/apply(t2[,c("err1","err2")], MARGIN=1, function(v){ sqrt(sum(v**2, na.rm=T)) })
plot(0,0, col=0, xlim=lim, ylim=lim, cex=size2**2, xlab="ABA lib 2 GMMA effects", ylab="Diazi:D80E:Y81V lib 2 GMMA effects", main=sprintf("%d common subst.",nrow(t2)))
abline(h=0,v=0, lty=2)
i = which(t2$s %in% c("R104H","E43D","T118R","S29Q","R134K"))
points(t2[i,"ddG1"], t2[i,"ddG2"], cex=size2[i]**2, pch=16, col="orange")
points(t2$ddG1, t2$ddG2, cex=size2**2)
# text(t2$ddG1, t2$ddG2, t2$s, pos=2, offset=0.7, cex=.8)
i = which(t2$s %in% c("E102D","R104H"))
text(t2[i,"ddG1"], t2[i,"ddG2"], t2[i,"s"], pos=2, offset=0.8, cex=.8)
i = which(t2$s %in% c("E68K","Q69D","Q69E"))
text(t2[i,"ddG1"], t2[i,"ddG2"], t2[i,"s"], pos=2, offset=1.4, cex=.8)
i = which(t2$s == "T118R")
text(t2[i,"ddG1"], t2[i,"ddG2"], t2[i,"s"], pos=2, offset=2.3, cex=.8)

# tile 3
t3 = common_conf(a3,d3)

size3 = .4/apply(t3[,c("err1","err2")], MARGIN=1, function(v){ sqrt(sum(v**2, na.rm=T)) })
plot(0,0, col=0, xlim=lim, ylim=lim, xlab="ABA lib 3 GMMA effects", ylab="Diazi lib 3 GMMA effects", main=sprintf("%d common subst.",nrow(t3)))
abline(h=0,v=0, lty=3)
i = which(t3$s %in% c("R104H","E43D","T118R","S29Q","R134K"))
points(t3[i,"ddG1"], t3[i,"ddG2"], cex=size3[i]**2, pch=16, col="orange")
points(t3$ddG1, t3$ddG2, cex=size3**2)
# text(t3$ddG1, t3$ddG2, t3$s, pos=2, offset=0.7, cex=.8)
i = which(t3$s %in% c("R134K"))
text(t3[i,"ddG1"], t3[i,"ddG2"], t3[i,"s"], pos=2, offset=0.4, cex=.8)
i = which(t3$s %in% c("E102D","R104H","K131R","E132D","N133G","D146N"))
text(t3[i,"ddG1"], t3[i,"ddG2"], t3[i,"s"], pos=c(4,2,2,2,1,2), offset=0.8, cex=.8)
i = which(t3$s == "T118R")
text(t3[i,"ddG1"], t3[i,"ddG2"], t3[i,"s"], pos=2, offset=1.3, cex=.8)

quartz.save("gmma_combined3.png", type="png")


pcol = c("s","obs","active","inactive","init_ddG","ddG_glob","stderr_meas","stderr_subfit_est","rank","eff")


# Expected substitutions
es_pos_stab = unique(es[which(es$expected=="Potential stabilizing"),"pos"])
es_pos_destab = unique(es[which(es$expected=="Predicted destabilizing"),"pos"])
es_pos_rev = unique(es[which(es$expected=="DIAZI reversion"),"pos"])
print(sprintf("Positions with both expected stab. and destab. subst.: %s",paste0(intersect(es_pos_stab,es_pos_destab),collapse=",")))
print(sprintf("Positions with both expected stab. and diazi reversions subst.: %s",paste0(intersect(es_pos_stab,es_pos_rev),collapse=",")))
print(sprintf("Positions with both expected destab. and diazi reversions subst.: %s",paste0(intersect(es_pos_destab,es_pos_rev),collapse=",")))

# issue that D80E is in adc but not estimated in dia although stabilizing - put to dia_eff=stab?
i_stab = which(adc$expect=="Potential stabilizing" & ! is.na(adc$aba_eff) & ! is.na(adc$dia_eff))
i_destab = which(adc$expect=="Predicted destabilizing" & ! is.na(adc$aba_eff) & ! is.na(adc$dia_eff))

# all points of combined analysis
quartz(width=6, height=6)
lim=c(-.6,3.2)
stab_col = "green"
destab_col = "red"
plot(0,0,col=0, xlim=lim, ylim=lim, xlab="ABA GMMA effects", ylab="Diazi GMMA effects")
abline(0,1)
points(0, 0, pch=20, col=1)

fit = lm(dia_ddG~aba_ddG, adc)
abline(coef(fit), lty=2)
print(sprintf("Fit with intersept: dia = %.4f + %.4f x aba", coef(fit)[1], coef(fit)[2]))

# fit through 0,0
fit0 = lm(dia_ddG~aba_ddG-1, adc)
aba_scale = 1/coef(fit0)
# abline(0, coef(fit0), lty=2)
print(sprintf("Fit without intersept: dia = %.4f x aba (aba scale  %.4f)", coef(fit0)[1], aba_scale))

points(adc[i_stab,"aba_ddG"], adc[i_stab,"dia_ddG"], pch=16, col=stab_col)
points(adc[i_destab,"aba_ddG"], adc[i_destab,"dia_ddG"], pch=16, col=destab_col)
points(adc$aba_ddG, adc$dia_ddG, pch=1, lwd=1, col="black")
legend("topleft", c("Potential stabilizing","Predicted destabilizing","Diazi = ABA",sprintf("Diazi = %.2f + %.2f x ABA",coef(fit)[1], coef(fit)[2])),
       pch=c(16,16,NA,NA), lty=c(NA,NA,1,2), col=c(stab_col,destab_col,"black","black"))
quartz.save("gmma_combined1_full.png", type="png")


# Pross substitutions
print(sprintf("GMMA effect of %d confident and expected stablizing subst",length(i_stab)))
print(table( adc[i_stab,c("aba_eff","dia_eff")] ))

pross_subst = es[which(es$source=="PROSS"),"s"]
# pross_subst_est = pross_subst[which(pross_subst %in% adc$s)]
i_pross = which(adc$s %in% pross_subst  & (! is.na(adc$aba_eff)) & (! is.na(adc$dia_eff)))
print(sprintf("GMMA effects of %d confident and expected out of %d (%.2f%%) PROSS substitutions. Missing: %s", length(i_pross), length(pross_subst),
              length(i_pross)/length(pross_subst)*100, paste0(setdiff(pross_subst,adc[i_pross,"s"]),collapse=",")))
print(table( adc[i_pross,c("aba_eff","dia_eff")] ))

i_other_stab = which(adc$expect=="Potential stabilizing" & ! adc$s %in% pross_subst & ! is.na(adc$aba_eff) & ! is.na(adc$dia_eff))
print(sprintf("GMMA effect of %d confident and expected stablizing non-PROSS subst",length(i_other_stab)))
print(table( adc[i_other_stab,c("aba_eff","dia_eff")] ))

n_est = sum((! is.na(adc$aba_ddG)) & (! is.na(adc$dia_ddG)))
quartz(width=8, height=8)
par(bg="white")
plot(adc$aba_ddG, adc$dia_ddG, col=0, xlab="ABA GMMA effects", ylab="Diazi GMMA effects", main=sprintf("%d expected of %d",n_est,nrow(es)))

# arrows(x0=adc[,"aba_ddG"], y0=adc[,"dia_ddG"]-adc[,"dia_err"], y1=adc[,"dia_ddG"]+adc[,"dia_err"], code=3, angle=90, length=.01)
# arrows(x0=adc[,"aba_ddG"]-adc[,"aba_err"], x1=adc[,"aba_ddG"]+adc[,"aba_err"], y0=adc[,"dia_ddG"], code=3, angle=90, length=.01)

points(adc[i_destab,"aba_ddG"], adc[i_destab,"dia_ddG"], pch=16, col=2)
points(adc[i_pross,"aba_ddG"], adc[i_pross,"dia_ddG"], pch=16, col=3)
points(adc[i_other_stab,"aba_ddG"], adc[i_other_stab,"dia_ddG"], pch=16, col=4)
points(adc$aba_ddG, adc$dia_ddG, col=1, pch=1)
legend("topleft", c("Destabilizing","PROSS","Other stabilizing"), pch=16, col=c(2,3,4))
quartz.save("gmma_combined_pross.png", type="png")



jpeg("gmma_combined_full.jpg", width=10, height=10, units="cm", res=300, quality=90, pointsize=8)
par(bg="white", mar=c(5,4,2,2)+.1)
plot(0,0,col=0, xlim=lim, ylim=lim, xlab="ABA GMMA effects", ylab="Diazi GMMA effects")
abline(0,1)
abline(coef(fit), lty=2)

destab_col = 2
stabp_col = 3
stabo_col = 4
points(adc[i_pross,"aba_ddG"], adc[i_pross,"dia_ddG"], pch=16, col=stabp_col)
points(adc[i_other_stab,"aba_ddG"], adc[i_other_stab,"dia_ddG"], pch=16, col=stabo_col)
points(adc[i_destab,"aba_ddG"], adc[i_destab,"dia_ddG"], pch=16, col=destab_col)
points(adc$aba_ddG, adc$dia_ddG, pch=1, lwd=1, col="black")
legend("topleft", c("Potential stabilizing, PROSS","Potential stabilizing, other","Predicted destabilizing","Diazi = ABA",sprintf("Diazi = %.2f + %.2f x ABA",coef(fit)[1], coef(fit)[2])),
       pch=c(16,16,16,NA,NA), lty=c(NA,NA,NA,1,2), col=c(stabp_col,stabo_col,destab_col,"black","black"))
dev.off()



# plot fit
# load variant data
var_ac = read.csv2("aba_lib123/mutant.csv")
var_d12 = read.csv2("dia_lib12_reref/mutant.csv")
var_d3 = read.csv2("dia_lib3/mutant.csv")

# find main tile library for each variant
var_ac$lib = apply(var_ac[,c("reads1","reads2","reads3")], MARGIN=1, which.max)
var_d12$lib = apply(var_d12[,c("reads1","reads2")], MARGIN=1, which.max)

# plotting function
plot_var = function(sig, res, col=1) {
    sig = as.numeric(sig)
    res = as.numeric(res)
    pred = sig - res
    plot(pred, sig, col=col, pch=20, cex=.25, xlim=c(-7,0),
         xlab="Enrichment predicted", ylab="Enrichment screen")
    abline(0,1, lty=2)
    rp = cor(pred, sig, method="pearson", use="complete.obs")
    print(sprintf("Pearson %.4f",rp))
}


jpeg("predictions.jpg", width=18, height=6, units="cm", res=300, pointsize=8)
par(mfrow=c(1,3), bg="white", mar=c(5,4,2,2)+.1)
plot_var(var_ac$signal, var_ac$residual, col=var_ac$lib)
legend("topleft", c("ABA tile 1 variant","ABA tile 2 variant","ABA tile 3 variant"), pch=c(20,20,20), col=c(1,2,3))
mtext("A", cex=1.4, font=2, line=0, at=-8)

plot_var(var_d12$signal, var_d12$residual, col=var_d12$lib)
legend("topleft", c("Diazi tile 1 variant","Diazi tile 2 variant"), pch=c(20,20), col=c(1,2))
mtext("B", cex=1.4, font=2, line=0, at=-8)

plot_var(var_d3$signal, var_d3$residual, col=3		)
legend("topleft", c("Diazi tile 3 variant"), pch=20, col=3)
mtext("C", cex=1.4, font=2, line=0, at=-8)
dev.off()





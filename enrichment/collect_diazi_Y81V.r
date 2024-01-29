options(width=160)

re_reference = function(vars,substs) {
    # check that the new WT is present
    i_wt = which(vars == paste0(substs,collapse=":"))
    stopifnot(length(i_wt) == 1)
    
    # make a vector of reversed substitutions
    nc = nchar(substs)
    substs_resi = as.numeric(substr(substs,2,nc-1))
    substs_rev = paste0(substr(substs,nc,nc), substs_resi, substr(substs,1,1))

    # change reference one substitutions at a time
    var_splitlist = strsplit(vars, ":")
    for (is in seq_along(substs)) {
	has_subst = grepl(substs[is], vars)
	
	# since var_splitlist is updated, this needs to be updated
	resi_splitlist = lapply(var_splitlist, function(vl){ as.numeric(substr(vl,2,nchar(vl)-1)) })
	has_resi = sapply(resi_splitlist, function(ril,resi){ resi %in% ril }, resi=substs_resi[is])

	# remove wt subst where it is present
	i_hs = which(has_subst)
	if (length(i_hs) > 0) {
	    var_splitlist[i_hs] = lapply(var_splitlist[i_hs], function(vl,sn){ vl[vl != sn] }, sn=substs[is])
	}

	# where variants have another substitution at position change the WT identity
	i_hr = which(has_resi & (! has_subst))
	change_subst_wt = function(vl,sn){
	    resi = as.numeric(substr(sn,2,nchar(sn)-1))
	    resis = as.numeric(substr(vl,2,nchar(vl)-1))
	    i = which(resis==resi)
	    stopifnot( length(i)==1 )
	    vl[i] = paste0( substr(sn,nchar(sn),nchar(sn)), substr(vl[i],2,nchar(vl[i])) )	    
	    return( vl )
	}
	if (length(i_hr) > 0) {
            var_splitlist[i_hr] = lapply(var_splitlist[i_hr], change_subst_wt, sn=substs[is])
	}

	# add reverse substitution where new WT substitution is not present nor is the position substituted
	i_nhs = which((! has_subst) & (! has_resi))
	add_subst = function(vl,sn) {
	    resi = as.numeric(substr(sn,2,nchar(sn)-1))
	    resis = as.numeric(substr(vl,2,nchar(vl)-1))
	    return( c(vl[which(resis < resi)], sn, vl[which(resis > resi)]) )
	}
	if (length(i_nhs) > 0) {
            var_splitlist[i_nhs] = lapply(var_splitlist[i_nhs], add_subst, sn=substs_rev[is])
        }
    }

    # return
    new_vars = sapply(var_splitlist, paste0, collapse=":")
    new_vars[i_wt] = "WT"
    return(new_vars)
}

mutate = function(ref_vec, substs) {
    nc = nchar(substs)
    wts = substr(substs, 1, 1)
    resis = as.numeric(substr(substs, 2, nc-1))
    muts = substr(substs, nc, nc)
    stopifnot(ref_vec[resis] == wts)
    stopifnot(all(wts != muts))
    ref_vec[resis] = muts
    return(ref_vec)
}


#############################################################################################
###  PYR1
#############################################################################################
# PYR1 aba (pJS637):    H60P              N90S
# PYR1 diazi (pJS752):  H60P, V81Y, L87M, N90S, F108Y, M158V, F159G, A160V
#            4    10        20        30        40        50  54    60        70        80      88 90    97 100       110       120       130 134   140       150       160       170     179180      190
#            |     |         |         |         |         |   |     |         |         |       | |      |  |         |         |         |   |     |         |         |         |        ||         |        
# PYR1:   MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKHFIKSCSVEQNFEMRVGCTRDVIVISGLPANTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMARNSGDGSGSQVT
#                                                                    *                    #     #  *                 #                                                 ###
wt_dia = "MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKPFIKSCSVEQNFEMRVGCTRDYIVISGMPASTSTERLDILDDERRVTGYSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRVGVDTVVKLNLQKLATVAEAMARNSGDGSGSQVT"
wt_aba = "MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKPFIKSCSVEQNFEMRVGCTRDVIVISGLPASTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMARNSGDGSGSQVT"
# WT lib 1:  ELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKPFIKSCSVEQNFEMRVGCTRDVIVISGLP                                             |
# WT lib 2:                                                    KPQTYKPFIKSCSVEQNFEMRVGCTRDVIVISGLPASTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENR
# WT lib 3:                                                                                               DILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMA

wt = wt_dia

# Read list of expected substitutions
lib_csv = read.csv("expected_substitutions.csv")
#   position wt     mut expected                                               note
# 1       29  S     Q/T     stab Original list of stabilizing mutations, dark green
# 2       15  N D/E/S/Q     stab Original list of stabilizing mutations, dark green
# 3       26  D G/N/S/K     stab Original list of stabilizing mutations, dark green
# 4      104  R       H     stab Original list of stabilizing mutations, dark green

# make data frame of expected substitutions
mut_split_list = strsplit(lib_csv$mut, "/")
rep_times = sapply(mut_split_list,length)

expected_subst = data.frame(wt = rep(lib_csv$wt, times=rep_times),
                            resi = as.numeric(rep(lib_csv$position, times=rep_times)),
			    mut = unlist(mut_split_list))
expected_subst$subst = paste0(expected_subst$wt, expected_subst$resi, expected_subst$mut)			    
expected_subst$expected = rep(lib_csv$expected,times=rep_times)
expected_subst$lib1 = expected_subst$resi >= 4  & expected_subst$resi <= 82  # Nic writes 82 but some sequencing had up to 88?
expected_subst$lib2 = expected_subst$resi >= 54 & expected_subst$resi <= 134
expected_subst$lib3 = expected_subst$resi >= 97 & expected_subst$resi <= 179
expected_subst$note = note=rep(lib_csv$note,times=rep_times)

# remove M158R from expected because of mis-match with diazi background?!
expected_subst = expected_subst[-which(expected_subst$subst=="M158R"),]

# check that WT is consistent
wt_vec = strsplit(wt,"")[[1]]
stopifnot(all( wt_vec[expected_subst$resi] == expected_subst$wt)) 

# dump info
print(sprintf("Expected substitutions: %d on %d positions. Lib 1: %d on %d. Lib 2: %d on %d. Lib 3: %d on %d",
              nrow(expected_subst), length(unique(expected_subst$resi)),
	      sum(expected_subst$lib1), length(unique(expected_subst[which(expected_subst$lib1),"resi"])),
	      sum(expected_subst$lib2), length(unique(expected_subst[which(expected_subst$lib2),"resi"])),
	      sum(expected_subst$lib3), length(unique(expected_subst[which(expected_subst$lib3),"resi"]))))

i12 = which(expected_subst$lib1 & expected_subst$lib2)
print(sprintf("Expects %d subst in overlap between lib 1 and 2: %s", length(i12), paste(expected_subst[i12,"subst"],collapse=", ")))

i23 = which(expected_subst$lib2 & expected_subst$lib3)
print(sprintf("Expects %d subst in overlap between lib 2 and 3: %s", length(i23), paste(expected_subst[i23,"subst"],collapse=", ")))


#############################################################################################
###  Read counts
#############################################################################################
# function to parse counts files
read.counts = function(filename) {
    d = read.csv(filename)
    stopifnot( length(table(d$sel_total)) == 1 )
    stopifnot( all(d$sel_total == d$ref_total) )
    stopifnot( all(d$sel_counts == d$ref_counts) )    
    return( d[,c("variant","sel_counts")] )
}

# read data
pyr1_dia_lib1_ref = read.counts("../sequencing_counts/PYR1_Diazi_Lib1_ref_counts.csv")
pyr1_dia_lib1_bnd = read.counts("../sequencing_counts/PYR1_Diazi_Lib1_bound_counts.csv")
pyr1_dia_lib2_ref = read.counts("../sequencing_counts/PYR1_Diazi_Lib2_ref_counts.csv")
pyr1_dia_lib2_bnd = read.counts("../sequencing_counts/PYR1_Diazi_Lib2_bound_counts.csv")
# pyr1_dia_lib3_ref = read.counts("../sequencing_counts/PYR1_Diazi_Lib3_ref_counts.csv")
# pyr1_dia_lib3_bnd = read.counts("../sequencing_counts/PYR1_Diazi_Lib3_bound_counts.csv")

# reformat variants as a colon-separated sequence of substitutions
var_orig = unique(c(pyr1_dia_lib1_ref$variant, pyr1_dia_lib1_bnd$variant,
                    pyr1_dia_lib2_ref$variant, pyr1_dia_lib2_bnd$variant))
#		    pyr1_dia_lib3_ref$variant, pyr1_dia_lib3_bnd$variant))
var = gsub("(", "", var_orig, fixed=T)
var = gsub(")", "", var, fixed=T)
var = sapply(strsplit(var, ","), paste0, collapse=":")
var = gsub(" ","",var)

# make new reference
new_ref = c("Y81V")
var = re_reference(var, new_ref)
wt_vec = mutate(wt_vec, new_ref)
wt = paste0(wt_vec, collapse="")

subst_list = strsplit(var, ":")
d = data.frame(var=var, var_orig=var_orig, n_sub=sapply(subst_list, length))
d[which(d$var=="WT"),"n_sub"] = 0

# Order variants, remember that first element is WT with resi=NA and taa="T"
resi_list = lapply(subst_list, function(v){  nv=as.numeric(substr(v,2,nchar(v)-1)); stopifnot(all(diff(nv) > 0)); nv })
taa_list = lapply(subst_list, function(v){ nc=nchar(v); substr(v,nc,nc) })

sorting_df = data.frame(resi1 = sapply(resi_list, "[", 1),
                        resi2 = sapply(resi_list, function(v){if (length(v)>1) i=v[2] else i=0; i}),
                        resi3 = sapply(resi_list, function(v){if (length(v)>2) i=v[3] else i=0; i}),
                        resi4 = sapply(resi_list, function(v){if (length(v)>3) i=v[4] else i=0; i}),
                        resi5 = sapply(resi_list, function(v){if (length(v)>4) i=v[5] else i=0; i}),
                        resi6 = sapply(resi_list, function(v){if (length(v)>5) i=v[6] else i=0; i}),
                        resi7 = sapply(resi_list, function(v){if (length(v)>6) i=v[7] else i=0; i}),
                        taa1 = sapply(taa_list, "[", 1),
                        taa2 = sapply(taa_list, function(v){if (length(v)>1) i=v[2] else i=0; i}),
                        taa3 = sapply(taa_list, function(v){if (length(v)>2) i=v[3] else i=0; i}),
                        taa4 = sapply(taa_list, function(v){if (length(v)>3) i=v[4] else i=0; i}),
                        taa5 = sapply(taa_list, function(v){if (length(v)>4) i=v[5] else i=0; i}),
                        taa6 = sapply(taa_list, function(v){if (length(v)>5) i=v[6] else i=0; i}),
                        taa7 = sapply(taa_list, function(v){if (length(v)>6) i=v[7] else i=0; i}))
		
i = order(d$n_sub,
          sorting_df$resi1, sorting_df$taa1,
	  sorting_df$resi2, sorting_df$taa2,
	  sorting_df$resi3, sorting_df$taa3,
	  sorting_df$resi4, sorting_df$taa4,
	  sorting_df$resi5, sorting_df$taa5,
	  sorting_df$resi6, sorting_df$taa6,
	  sorting_df$resi7, sorting_df$taa7)

# data frame with variant counts
vc = d[i,]
rownames(vc) = NULL

# check that variant WT amino acid match expected WT sequence
subst_list = strsplit(vc$var, ":")
var_wt_list = lapply(subst_list, substr, start=1, stop=1)
resi_list = lapply(subst_list, function(v){  nv=as.numeric(substr(v,2,nchar(v)-1)); stopifnot(all(diff(nv) > 0)); nv })
if (any( unlist(var_wt_list[2:length(var_wt_list)]) != wt_vec[unlist(resi_list[2:length(resi_list)])] )) {
    # With variants have substitutions with WT different from expected WT
    wt_diff = sapply(seq(2,length(resi_list)), function(i){ any(wt_vec[resi_list[[i]]] != var_wt_list[[i]]) })
    print(sprintf("WARNING: %d variants have wild-type different from hardcoded given WT", sum(wt_diff)))    
    print(vc[which(wt_diff),])
}

# add source column names (scn) from source data frame (sdf) to target data frame (tdf) with column names (tcn)
add_data = function(tdf, tcn, sdf, scn) {
    stopifnot( length(tcn) == length(scn) )
    stopifnot( all( ! tcn %in% colnames(tdf) ) )
    tdf[,tcn] = sdf[match(tdf$var_orig,sdf$variant),scn]
    for ( cn in tcn ) {
        t_mask = is.na(tdf[,cn])
        tdf[which(t_mask),cn] = 0
    }

    # How many reads
    for ( cn in scn ) {
        # How many reads
        nc_wt = sdf[which(sdf$variant=="WT"),cn]
        if (length(nc_wt) != 1) {
	    nc_wt = 0
	    print(sprintf("    WARNING: File has %d variants named WT", length(nc_wt)))
	}
        nc = sum(sdf[,cn]) - nc_wt
        nv = sum(sdf[,cn] > 0) - 1
        print(sprintf("    Column %5s has %d read counts among %d variants (%.2f reads-per-var) excluding %d WT reads (%.2f%%)",
                      cn, nc, nv, nc/nv, nc_wt, nc_wt/(nc+nc_wt)*100))
    }

    # How many reads on multi-mutants >7
    for ( cn in tcn ) {
        nz_mask = tdf[,cn] > 0
        it = which(nz_mask  &  tdf$n_sub > 7)
        if (length(it) >0) {
            print(sprintf("    Found large multimut: Reads of %d variants with 8 to %d subst are %d (%.2f%% excl. WT)",
                          length(it), max(tdf[it,"n_sub"]), sum(tdf[it,cn]), sum(tdf[it,cn])/nc*100))
        }
	
        # How many substitutions
        subst_list = strsplit(tdf[which(nz_mask),"var"], ":")
        usub = unique(unlist(subst_list))
	n_usub_lib1 = sum(usub %in% expected_subst[which(expected_subst$lib1),"subst"])
	n_usub_lib2 = sum(usub %in% expected_subst[which(expected_subst$lib2),"subst"])
	# n_usub_lib3 = sum(usub %in% expected_subst[which(expected_subst$lib3),"subst"])
        print(sprintf("    Among %d unique subst, Theres are %d of %d from lib 1 (%.1f%%), %d of %d from lib2 (%.1f%%) and %d unexpected (%.1f%%)",
	              length(usub),
		      n_usub_lib1, sum(expected_subst$lib1), n_usub_lib1/sum(expected_subst$lib1)*100,
		      n_usub_lib2, sum(expected_subst$lib2), n_usub_lib2/sum(expected_subst$lib2)*100,
		      # n_usub_lib3, sum(expected_subst$lib3), n_usub_lib3/sum(expected_subst$lib3)*100,
		      sum(! usub %in% expected_subst$subst), sum(! usub %in% expected_subst$subst)/length(usub)*100))
		      
        tab = table(tdf[which(nz_mask),"n_sub"])
        print("N-mut distribution")
        print(tab)
    }

    return(tdf)
}

print("==== pyr1_dia_lib1 ====")
vc = add_data(vc, c("all1"), pyr1_dia_lib1_ref, c("sel_counts"))
vc = add_data(vc, c("bnd1"), pyr1_dia_lib1_bnd, c("sel_counts"))
print("==== pyr1_dia_lib2 ====")
vc = add_data(vc, c("all2"), pyr1_dia_lib2_ref, c("sel_counts"))
vc = add_data(vc, c("bnd2"), pyr1_dia_lib2_bnd, c("sel_counts"))
# print("==== pyr1_dia_lib3 ====")
# vc = add_data(vc, c("all3"), pyr1_dia_lib3_ref, c("sel_counts"))
# vc = add_data(vc, c("bnd3"), pyr1_dia_lib3_bnd, c("sel_counts"))


# Remove variants without reads in reference lib
for (ilib in seq(2)) {
    sel = paste0("bnd",ilib)
    ref = paste0("all",ilib)
    
    # zero read counts in all bound pools where no reads are observed in reference pool
    i = which(vc[,sel] > 0 & vc[,ref] == 0)
    print(sprintf("Removing %d counts (%.1f%% excl. WT) from %d varinats (%.1f%% excl. WT) in %s because these variants are not observed in %s",
                  sum(vc[i,sel]), sum(vc[i,sel])/sum(vc[2:nrow(vc),sel])*100, length(i), length(i)/(sum(vc[,sel]>0)-1)*100, sel, ref))
    vc[i,sel] = 0
}

# then remove variants that are not observed at all
# i_rm = which(vc$all1 == 0 & vc$all2 == 0 & vc$all3 == 0)
i_rm = which(vc$all1 == 0 & vc$all2 == 0)
print(sprintf("Removing %d of %d variants (%.1f%% excl. WT) that are now completely un-observed", length(i_rm), nrow(vc)-1, length(i_rm)/(nrow(vc)-1)*100))
vc = vc[-i_rm,]

# consider variant present in library if any reads are observed
vc$lib1 = apply(vc[,c("all1","bnd1")], MARGIN=1, function(v){ any(v>0) })
vc$lib2 = apply(vc[,c("all2","bnd2")], MARGIN=1, function(v){ any(v>0) })
# vc$lib3 = apply(vc[,c("all3","bnd3")], MARGIN=1, function(v){ any(v>0) })

vc$reads1 = apply(vc[,c("all1","bnd1")], MARGIN=1, sum)
vc$reads2 = apply(vc[,c("all2","bnd2")], MARGIN=1, sum)
# vc$reads3 = apply(vc[,c("all3","bnd3")], MARGIN=1, sum)


# look at complexity
# cns = c("all1","bnd1","all2","bnd2","all3","bnd3")
cns = c("all1","bnd1","all2","bnd2")
df = vc[2:nrow(vc),cns]
read_stats = data.frame(name=cns, reads=apply(df, MARGIN=2, sum))
# number of barcodes with non-zero number of reads
read_stats$nvar_nz = apply(df, MARGIN=2, function(v){ sum(v>0) })
# effective number of barcodes with reads, based on p*log(p)
read_stats$nvar_eff = apply(df, MARGIN=2, function(v){ p=v/sum(v)+1e-12; exp(-sum(p*log(p))) })
# average reads per variant
read_stats$reads_per_nz = read_stats$reads / read_stats$nvar_nz
read_stats$reads_per_eff = read_stats$reads / read_stats$nvar_eff
print(read_stats)


# # number of variants in the two libraries
# print(sprintf("Of %d variants, %d are in lib 1 (%.1f%%), %d are in lib 2 (%.1f%%), and %d are in lib 3 (%.1f%%)",
#               nrow(vc), sum(vc$lib1), sum(vc$lib1)/nrow(vc)*100, sum(vc$lib2), sum(vc$lib2)/nrow(vc)*100, sum(vc$lib3), sum(vc$lib3)/nrow(vc)*100))
print(sprintf("Of %d variants, %d are in lib 1 (%.1f%%) and %d are in lib 2 (%.1f%%)",
              nrow(vc), sum(vc$lib1), sum(vc$lib1)/nrow(vc)*100, sum(vc$lib2), sum(vc$lib2)/nrow(vc)*100))
nv_1or2 = sum(vc$lib1 | vc$lib2)
nv_1and2 = sum(vc$lib1 & vc$lib2)
print(sprintf("Of %d variants in libs 1 and 2, %d are in both (%.2f%% of all, %.2f%% of lib1 and %.2f%% of lib2). Multimut. distribution of common:",
              nv_1or2, nv_1and2, nv_1and2/nv_1or2*100, nv_1and2/sum(vc$lib1)*100, nv_1and2/sum(vc$lib2)*100))
print(table(vc[which(vc$lib1 & vc$lib2),"n_sub"]))

# nv_2or3 = sum(vc$lib2 | vc$lib3)
# nv_2and3 = sum(vc$lib2 & vc$lib3)
# print(sprintf("Of %d variants in libs 2 and 3, %d are in both (%.2f%% of all, %.2f%% of lib2 and %.2f%% of lib3). Multimut. distribution of common:",
#               nv_2or3, nv_2and3, nv_2and3/nv_2or3*100, nv_2and3/sum(vc$lib2)*100, nv_2and3/sum(vc$lib3)*100))
# print(table(vc[which(vc$lib2 & vc$lib3),"n_sub"]))


#############################################################################################
###  Calculate an enrichment score per variant
#############################################################################################
print("Calculate enrichment scores of combined library without applying thresholds")

# so, add pseudo counts and calc enrichment
pseudocounts = 1
# rpm = vc[,c("all1","bnd1","all2","bnd2","all3","bnd3")]
rpm = vc[,c("all1","bnd1","all2","bnd2")]

# normal pesudo count gives sequencing-depth dependent frequencies and non-zero enrichment
rpm = rpm + pseudocounts

# normalize to reads-per-million
rpm = as.data.frame(apply(rpm, MARGIN=2, function(v){ v/sum(v)*10^6 }))

# enrichment relative to WT - rpm normalization not necessary for this
calc_enrichment = function(input, output, i_wt=1) {
    if (is.na(i_wt)) {
        log2(output/input)
    } else {
        log2(output/input * input[i_wt]/output[i_wt])
    }
}

cut = 20
for (libi in seq(2)) {
    cn = paste0("enrich",libi)
    reads_cn = paste0("reads",libi)

    # normalizing to WT typically gives a good correlation between same variants in both libs but not necessarily similar distributions
    vc[,cn] = calc_enrichment(rpm[,paste0("all",libi)], rpm[,paste0("bnd",libi)], i_wt=1)

    # variants without reads
    # with pseudocout != 0, variants without counts have constant enrichment depending on the total read count
    vc[which(vc[,reads_cn] == 0), cn] = NA

    # variants with number of reads below threshold
    i = which(vc[,reads_cn] < cut & vc[,reads_cn] > 0)
    nreads = sum(vc[i,reads_cn])
    print(sprintf("Lib %d filtering: Removing %d variants (%.1f%%) with less than %d reads (all+bnd). Total reads filtered %d (%.1f%%)",
                  libi, length(i), length(i)/sum(vc[,reads_cn] > 0)*100, cut, nreads, nreads/sum(vc[,reads_cn])*100))
    vc[i,cn] = NA
}

# vc$enrich_all = apply(vc[,c("enrich1","enrich2","enrich3")], MARGIN=1, mean, na.rm=T)
vc$enrich_all = apply(vc[,c("enrich1","enrich2")], MARGIN=1, mean, na.rm=T)


# number of variants with scores
nall = sum(! is.na(vc$enrich_all))
nlib1 = sum(! is.na(vc$enrich1))
nlib2 = sum(! is.na(vc$enrich2))
nlib3 = sum(! is.na(vc$enrich3))
print(sprintf("Of %d enrichment scores, %d are in lib 1 (%.1f%%), %d are in lib 2 (%.1f%%), and %d are in lib 3 (%.1f%%)",
              nall, nlib1, nlib1/nall*100, nlib2, nlib2/nall*100, nlib3, nlib3/nall*100))
	      
nv_1or2 = sum((! is.na(vc$enrich1)) | (! is.na(vc$enrich2)))
nv_1and2 = sum((! is.na(vc$enrich1)) & (! is.na(vc$enrich2)))
print(sprintf("Of %d variants in libs 1 and 2, %d are in both (%.2f%% of common, %.2f%% of lib1 and %.2f%% of lib2). Multimut. distribution of common:",
              nv_1or2, nv_1and2, nv_1and2/nv_1or2*100, nv_1and2/sum(! is.na(vc$enrich1))*100, nv_1and2/sum(! is.na(vc$enrich2))*100))
print(table(vc[which((! is.na(vc$enrich1)) & (! is.na(vc$enrich2))),"n_sub"]))

nv_2or3 = sum((! is.na(vc$enrich2)) | (! is.na(vc$enrich3)))
nv_2and3 = sum((! is.na(vc$enrich2)) & (! is.na(vc$enrich3)))
print(sprintf("Of %d variants in libs 2 and 3, %d are in both (%.2f%% of common, %.2f%% of lib2 and %.2f%% of lib3). Multimut. distribution of common:",
              nv_2or3, nv_2and3, nv_2and3/nv_2or3*100, nv_2and3/sum(! is.na(vc$enrich2))*100, nv_2and3/sum(! is.na(vc$enrich3))*100))
print(table(vc[which((! is.na(vc$enrich2)) & (! is.na(vc$enrich3))),"n_sub"]))

print("Nonsence variants wuth reads above threshold")
i_nons = which(grepl('*', vc$var, fixed=T))
print(vc[intersect(i_nons,which(! is.na(vc$enrich_all))),])


# plot
i12 = which(vc$lib1 & vc$lib2)
# i23 = which(vc$lib2 & vc$lib3)

# cn = c("enrich1","enrich2","enrich3")
cn = c("enrich1","enrich2")
lim = c(min(vc[,cn], na.rm=T), max(vc[,cn], na.rm=T))
quartz(width=6, height=6)
plot(0,0,col=0, xlim=lim, ylim=lim, xlab="Enrichment", ylab="Enrichment")
abline(0,1)
points(vc[i12,"enrich1"], vc[i12,"enrich2"], col=2, pch=16)
# points(vc[i23,"enrich2"], vc[i23,"enrich3"], col=3, pch=16)
# legend("topleft", c("Lib 1 vs 2","Lib 2 vs 3"), pch=16, col=c(2,3))
legend("topleft", c("Lib 1 vs 2"), pch=16, col=c(2))
quartz.save("dia_Y81V_enrich_overlap_corr.png", type="png")

x = 10^0
cn = c("enrich1","enrich2")
lim = c(min(vc[,cn], na.rm=T), max(vc[,cn], na.rm=T))
breaks = seq(floor(lim[1]*x)/x, ceiling(lim[2]*x)/x, 1/x)
h1 = hist(vc[,"enrich1"], breaks=breaks, plot=F)
h2 = hist(vc[,"enrich2"], breaks=breaks, plot=F)
# h3 = hist(vc[,"enrich3"], breaks=breaks, plot=F)
ha = hist(vc[,"enrich_all"], breaks=breaks, plot=F)
hn = hist(vc[i_nons,"enrich_all"], breaks=breaks, plot=F)

quartz(width=8, height=5)
plot(0,0,col=0, xlim=lim, ylim=c(0,0.7), xlab="Enrichment", ylab="Density")
lines(ha$mids, ha$density, lwd=2, col=1)
lines(h1$mids, h1$density, lwd=2, col=2)
lines(h2$mids, h2$density, lwd=2, col=3)
# lines(h3$mids, h3$density, lwd=2, col=4)
# legend("topleft", c("Combined","Lib 1","Lib 2","Lib 3"), lty=1, lwd=2, col=c(1,2,3,4))
legend("topleft", c("Combined","Lib 1","Lib 2"), lty=1, lwd=2, col=c(1,2,3))
quartz.save("dia_Y81V_enrich_distributions.png", type="png")


#############################################################################################
###  Dump results
#############################################################################################

save(vc, wt, expected_subst, file="variant_counts_dia_Y81V.rda")

i1 = which(! is.na(vc$enrich1))
write.csv(vc[i1,c("var","all1","bnd1","enrich1")], row.names=F, file="pyr1_dia_Y81V_lib1.csv")

i2 = which(! is.na(vc$enrich2))
write.csv(vc[i2,c("var","all2","bnd2","enrich2")], row.names=F, file="pyr1_dia_Y81V_lib2.csv")

# i3 = which(! is.na(vc$enrich3))
# write.csv(vc[i3,c("var","all3","bnd3","enrich3")], row.names=F, file="pyr1_dia_Y81V_lib3.csv")

ia = which(! is.na(vc$enrich_all))
write.csv(vc[ia,c("var","reads1","reads2","enrich_all")], row.names=F, file="pyr1_dia_Y81V_lib12.csv")

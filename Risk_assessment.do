*set panel id
xtset id year

*trim data by farm type and country (otherwise a lot of data in one specific group are dropped due to its size, e.g., MLT)
foreach c in "BEL" "CYP" "CZE" "DAN" "DEU" "ELL" "ESP" "EST" "FRA" "HRV" "HUN" "IRE" "ITA" "LTU" "LUX" "LVA" "MLT" "NED" "OST" "POL" "POR" "ROU" "SUO" "SVE" "SVK" "SVN" {
foreach x in 1 5 {
foreach var of varlist land laborhour totalasset fert seed protection lsu veterinary feedgrazing {
summ `var' if type == `x' & country == "`c'", d
replace `var' = . if `var'>r(p99) & type == `x' & country == "`c'"
replace `var' = . if `var'<r(p1) & type == `x' & country == "`c'"
}
}
}

*drop the missed values after trimming the data
drop if missing(land,laborhour,totalasset,fert,seed,protection,lsu,veterinary,feedgrazing)

*ownership, current ratio, and irrigation rate in land
gen own = 1 - rentland/land
gen current = totalcurrent/totalasset
gen irr = irriland/land

sum income 
generate double std_income = income / r(sd)
sum land
gen std_land = land / r(sd)
sum laborhour
gen std_labor = laborhour / r(sd)
sum totalasset
gen std_asset = totalasset / r(sd)
sum lsu
gen std_lsu = lsu / r(sd)
sum own
gen std_own = own / r(sd)
sum current
gen std_curr = current / r(sd)
sum irr
gen std_irr = irr / r(sd)
sum fert
gen std_fert = fert/r(sd)
sum seed
gen std_seed = seed/r(sd)
sum feedgrazing
gen std_feed = feedgrazing/r(sd)
sum protection
gen std_prot = protection/r(sd)
sum veterinary
gen std_vet = veterinary/r(sd)

drop if missing(own,irr,current)

estpost tabstat income land laborhour seed fert protection feedgrazing own current irr totalasset lsu if (type == 1), by(country) statistic(mean sd N max min) column(stat)
esttab using dairy_trimmed.csv, replace cell((mean(fmt(3)) sd(fmt(3)) N))

global Vinput std_land std_own std_irr std_asset std_labor std_lsu std_curr std_fert std_seed std_feed std_prot std_vet c.std_labor#(c.std_labor c.std_fert c.std_seed c.std_feed c.std_prot c.std_vet) c.std_fert#(c.std_fert c.std_seed c.std_feed c.std_prot c.std_vet) c.std_seed#(c.std_seed c.std_feed c.std_prot c.std_vet) c.std_feed#(c.std_feed c.std_prot c.std_vet) c.std_prot#(c.std_prot c.std_vet) c.std_vet#c.std_vet

quietly reghdfe std_income $Vinput, a(id year) vce(cluster id year) residuals(e)
estimates store mean_mod

gen e2 = e*e
gen e3 = e*e*e

*calculate derivatives for the sureg
matrix a = e(b)
matrix list e(b)

scalar a_LB = a[1,5]
scalar a_FT = a[1,8]
scalar a_SD = a[1,9]
scalar a_FD = a[1,10]
scalar a_PR = a[1,11]
scalar a_VT = a[1,12]

scalar a_LBsq = a[1,13]
scalar a_FTsq = a[1,19]
scalar a_SDsq = a[1,24]
scalar a_FDsq = a[1,28]
scalar a_PRsq = a[1,31]
scalar a_VTsq = a[1,33]

*interaction terms
scalar a_LBFT = a[1,14]
scalar a_LBSD = a[1,15]
scalar a_LBFD = a[1,16]
scalar a_LBPR = a[1,17]
scalar a_LBVT = a[1,18]
scalar a_FTSD = a[1,20]
scalar a_FTFD = a[1,21]
scalar a_FTPR = a[1,22]
scalar a_FTVT = a[1,23]
scalar a_SDFD = a[1,25]
scalar a_SDPR = a[1,26]
scalar a_SDVT = a[1,27]
scalar a_FDPR = a[1,29]
scalar a_FDVT = a[1,30]
scalar a_PRVT = a[1,32]

scalar list a_VT a_VTsq a_LBVT

quietly gen mu1_LB = a_LB + 2*a_LBsq*std_labor + a_LBFT*std_fert + a_LBSD*std_seed + a_LBFD*std_feed + a_LBPR*std_prot + a_LBVT*std_vet
quietly gen mu1_FT = a_FT + 2*a_FTsq*std_fert + a_LBFT*std_labor + a_FTSD*std_seed + a_FTFD*std_feed + a_FTPR*std_prot + a_FTVT*std_vet
quietly gen mu1_SD = a_SD + 2*a_SDsq*std_seed + a_LBSD*std_labor + a_FTSD*std_fert + a_SDFD*std_feed + a_SDPR*std_prot + a_SDVT*std_vet
quietly gen mu1_FD = a_FD + 2*a_FDsq*std_feed + a_LBFD*std_labor + a_FTFD*std_fert + a_SDFD*std_seed + a_FDPR*std_prot + a_FDVT*std_vet
quietly gen mu1_PR = a_PR + 2*a_PRsq*std_prot + a_LBPR*std_labor + a_FTPR*std_fert + a_SDPR*std_seed + a_FDPR*std_feed + a_PRVT*std_vet
quietly gen mu1_VT = a_VT + 2*a_VTsq*std_vet + a_LBVT*std_labor + a_FTVT*std_fert + a_SDVT*std_seed + a_FDVT*std_feed + a_PRVT*std_prot

quietly reghdfe e2 $Vinput, a(id year) vce(cluster id year)
estimates store var_mod

matrix a2 = e(b)
matrix list e(b)

scalar a2_LB = a2[1,5]
scalar a2_FT = a2[1,8]
scalar a2_SD = a2[1,9]
scalar a2_FD = a2[1,10]
scalar a2_PR = a2[1,11]
scalar a2_VT = a2[1,12]

scalar a2_LBsq = a2[1,13]
scalar a2_FTsq = a2[1,19]
scalar a2_SDsq = a2[1,24]
scalar a2_FDsq = a2[1,28]
scalar a2_PRsq = a2[1,31]
scalar a2_VTsq = a2[1,33]

scalar a2_LBFT = a2[1,14]
scalar a2_LBSD = a2[1,15]
scalar a2_LBFD = a2[1,16]
scalar a2_LBPR = a2[1,17]
scalar a2_LBVT = a2[1,18]
scalar a2_FTSD = a2[1,20]
scalar a2_FTFD = a2[1,21]
scalar a2_FTPR = a2[1,22]
scalar a2_FTVT = a2[1,23]
scalar a2_SDFD = a2[1,25]
scalar a2_SDPR = a2[1,26]
scalar a2_SDVT = a2[1,27]
scalar a2_FDPR = a2[1,29]
scalar a2_FDVT = a2[1,30]
scalar a2_PRVT = a2[1,32]

quietly gen mu2_LB = a2_LB + 2*a2_LBsq*std_labor + a2_LBFT*std_fert + a2_LBSD*std_seed + a2_LBFD*std_feed + a2_LBPR*std_prot + a2_LBVT*std_vet
quietly gen mu2_FT = a2_FT + 2*a2_FTsq*std_fert + a2_LBFT*std_labor + a2_FTSD*std_seed + a2_FTFD*std_feed + a2_FTPR*std_prot + a2_FTVT*std_vet
quietly gen mu2_SD = a2_SD + 2*a2_SDsq*std_seed + a2_LBSD*std_labor + a2_FTSD*std_fert + a2_SDFD*std_feed + a2_SDPR*std_prot + a2_SDVT*std_vet
quietly gen mu2_FD = a2_FD + 2*a2_FDsq*std_feed + a2_LBFD*std_labor + a2_FTFD*std_fert + a2_SDFD*std_seed + a2_FDPR*std_prot + a2_FDVT*std_vet
quietly gen mu2_PR = a2_PR + 2*a2_PRsq*std_prot + a2_LBPR*std_labor + a2_FTPR*std_fert + a2_SDPR*std_seed + a2_FDPR*std_feed + a2_PRVT*std_vet
quietly gen mu2_VT = a2_VT + 2*a2_VTsq*std_vet + a2_LBVT*std_labor + a2_FTVT*std_fert + a2_SDVT*std_seed + a2_FDVT*std_feed + a2_PRVT*std_prot


quietly reghdfe e3 $Vinput, a(id year) vce(cluster id year)
estimates store skew_mod

matrix a3 = e(b)
matrix list e(b)

scalar a3_LB = a3[1,5]
scalar a3_FT = a3[1,8]
scalar a3_SD = a3[1,9]
scalar a3_FD = a3[1,10]
scalar a3_PR = a3[1,11]
scalar a3_VT = a3[1,12]

scalar a3_LBsq = a3[1,13]
scalar a3_FTsq = a3[1,19]
scalar a3_SDsq = a3[1,24]
scalar a3_FDsq = a3[1,28]
scalar a3_PRsq = a3[1,31]
scalar a3_VTsq = a3[1,33]

scalar a3_LBFT = a3[1,14]
scalar a3_LBSD = a3[1,15]
scalar a3_LBFD = a3[1,16]
scalar a3_LBPR = a3[1,17]
scalar a3_LBVT = a3[1,18]
scalar a3_FTSD = a3[1,20]
scalar a3_FTFD = a3[1,21]
scalar a3_FTPR = a3[1,22]
scalar a3_FTVT = a3[1,23]
scalar a3_SDFD = a3[1,25]
scalar a3_SDPR = a3[1,26]
scalar a3_SDVT = a3[1,27]
scalar a3_FDPR = a3[1,29]
scalar a3_FDVT = a3[1,30]
scalar a3_PRVT = a3[1,32]

quietly gen mu3_LB = a3_LB + 2*a3_LBsq*std_labor + a3_LBFT*std_fert + a3_LBSD*std_seed + a3_LBFD*std_feed + a3_LBPR*std_prot + a3_LBVT*std_vet
quietly gen mu3_FT = a3_FT + 2*a3_FTsq*std_fert + a3_LBFT*std_labor + a3_FTSD*std_seed + a3_FTFD*std_feed + a3_FTPR*std_prot + a3_FTVT*std_vet
quietly gen mu3_SD = a3_SD + 2*a3_SDsq*std_seed + a3_LBSD*std_labor + a3_FTSD*std_fert + a3_SDFD*std_feed + a3_SDPR*std_prot + a3_SDVT*std_vet
quietly gen mu3_FD = a3_FD + 2*a3_FDsq*std_feed + a3_LBFD*std_labor + a3_FTFD*std_fert + a3_SDFD*std_seed + a3_FDPR*std_prot + a3_FDVT*std_vet
quietly gen mu3_PR = a3_PR + 2*a3_PRsq*std_prot + a3_LBPR*std_labor + a3_FTPR*std_fert + a3_SDPR*std_seed + a3_FDPR*std_feed + a3_PRVT*std_vet
quietly gen mu3_VT = a3_VT + 2*a3_VTsq*std_vet + a3_LBVT*std_labor + a3_FTVT*std_fert + a3_SDVT*std_seed + a3_FDVT*std_feed + a3_PRVT*std_prot

*print regression results
esttab mean_mod var_mod skew_mod using regresults.csv, b(3) se(3) ar2 label interaction(" x ") star(* 0.1 ** 0.05 *** 0.01) nogaps

*Overall risk attitudes
constraint 1 mu2_LB = mu2_FT
constraint 2 mu2_FT = mu2_SD
constraint 3 mu2_FT = mu2_FD
constraint 4 mu2_FT = mu2_PR
constraint 5 mu2_FT = mu2_VT

constraint 6 mu3_LB = mu3_FT
constraint 7 mu3_FT = mu3_SD
constraint 8 mu3_FT = mu3_FD
constraint 9 mu3_FT = mu3_PR
constraint 10 mu3_FT = mu3_VT

*remove the observations that are not used for regression
drop if missing(e)

sureg (mu1_FT mu2_FT mu3_FT)(mu1_SD mu2_SD mu3_SD)(mu1_FD mu2_FD mu3_FD)(mu1_PR mu2_PR mu3_PR)(mu1_VT mu2_VT mu3_VT), constraint(2 3 4 5 7 8 9 10)
estimates store entire_sure

esttab entire_sure using sure_entire_trimmed.csv, b(3) se(3) ar2 star(* 0.1 ** 0.05 *** 0.01) nogaps

*Risk attitudes for different time ranges
gen time_range = 0
replace time_range = 1 if (year>=2008 & year<2011)
replace time_range = 2 if (year>=2011)

constraint 11 0b.time_range#c.mu2_LB = 0b.time_range#c.mu2_FT
constraint 12 0b.time_range#c.mu2_FT = 0b.time_range#c.mu2_SD
constraint 13 0b.time_range#c.mu2_FT = 0b.time_range#c.mu2_FD
constraint 14 0b.time_range#c.mu2_FT = 0b.time_range#c.mu2_PR
constraint 15 0b.time_range#c.mu2_FT = 0b.time_range#c.mu2_VT

constraint 16 0b.time_range#c.mu3_LB = 0b.time_range#c.mu3_FT
constraint 17 0b.time_range#c.mu3_FT = 0b.time_range#c.mu3_SD
constraint 18 0b.time_range#c.mu3_FT = 0b.time_range#c.mu3_FD
constraint 19 0b.time_range#c.mu3_FT = 0b.time_range#c.mu3_PR
constraint 20 0b.time_range#c.mu3_FT = 0b.time_range#c.mu3_VT

constraint 21 1.time_range#c.mu2_LB = 1.time_range#c.mu2_FT
constraint 22 1.time_range#c.mu2_FT = 1.time_range#c.mu2_SD
constraint 23 1.time_range#c.mu2_FT = 1.time_range#c.mu2_FD
constraint 24 1.time_range#c.mu2_FT = 1.time_range#c.mu2_PR
constraint 25 1.time_range#c.mu2_FT = 1.time_range#c.mu2_VT

constraint 26 1.time_range#c.mu3_LB = 1.time_range#c.mu3_FT
constraint 27 1.time_range#c.mu3_FT = 1.time_range#c.mu3_SD
constraint 28 1.time_range#c.mu3_FT = 1.time_range#c.mu3_FD
constraint 29 1.time_range#c.mu3_FT = 1.time_range#c.mu3_PR
constraint 30 1.time_range#c.mu3_FT = 1.time_range#c.mu3_VT

constraint 31 2.time_range#c.mu2_LB = 2.time_range#c.mu2_FT
constraint 32 2.time_range#c.mu2_FT = 2.time_range#c.mu2_SD
constraint 33 2.time_range#c.mu2_FT = 2.time_range#c.mu2_FD
constraint 34 2.time_range#c.mu2_FT = 2.time_range#c.mu2_PR
constraint 35 2.time_range#c.mu2_FT = 2.time_range#c.mu2_VT

constraint 36 2.time_range#c.mu3_LB = 2.time_range#c.mu3_FT
constraint 37 2.time_range#c.mu3_FT = 2.time_range#c.mu3_SD
constraint 38 2.time_range#c.mu3_FT = 2.time_range#c.mu3_FD
constraint 39 2.time_range#c.mu3_FT = 2.time_range#c.mu3_PR
constraint 40 2.time_range#c.mu3_FT = 2.time_range#c.mu3_VT

sureg (mu1_FT i.time_range#c.mu2_FT i.time_range#c.mu3_FT)(mu1_SD i.time_range#c.mu2_SD i.time_range#c.mu3_SD)(mu1_FD i.time_range#c.mu2_FD i.time_range#c.mu3_FD)(mu1_PR i.time_range#c.mu2_PR i.time_range#c.mu3_PR)(mu1_VT i.time_range#c.mu2_VT i.time_range#c.mu3_VT), constraint(12 13 14 15 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 37 38 39 40)
estimates store time_sure

esttab time_sure using sure_time_untrimmed.csv, b(3) se(3) ar2 star(* 0.1 ** 0.05 *** 0.01) nogaps

test _b[0b.time_range#c.mu2_FT] = _b[1.time_range#c.mu2_FT]
test _b[2.time_range#c.mu2_FT] = _b[1.time_range#c.mu2_FT]

test _b[0b.time_range#c.mu3_FT] = _b[1.time_range#c.mu3_FT]
test _b[2.time_range#c.mu3_FT] = _b[1.time_range#c.mu3_FT]

*Risk attitudes for different time range and different farm types
*For type 1 constraint time 0
constraint 41 0b.time_range#1b.type#c.mu2_LB = 0b.time_range#1b.type#c.mu2_FT
constraint 42 0b.time_range#1b.type#c.mu2_FT = 0b.time_range#1b.type#c.mu2_SD
constraint 43 0b.time_range#1b.type#c.mu2_FT = 0b.time_range#1b.type#c.mu2_FD
constraint 44 0b.time_range#1b.type#c.mu2_FT = 0b.time_range#1b.type#c.mu2_PR
constraint 45 0b.time_range#1b.type#c.mu2_FT = 0b.time_range#1b.type#c.mu2_VT

constraint 46 0b.time_range#1b.type#c.mu3_LB = 0b.time_range#1b.type#c.mu3_FT
constraint 47 0b.time_range#1b.type#c.mu3_FT = 0b.time_range#1b.type#c.mu3_SD
constraint 48 0b.time_range#1b.type#c.mu3_FT = 0b.time_range#1b.type#c.mu3_FD
constraint 49 0b.time_range#1b.type#c.mu3_FT = 0b.time_range#1b.type#c.mu3_PR
constraint 50 0b.time_range#1b.type#c.mu3_FT = 0b.time_range#1b.type#c.mu3_VT

*For type 1 constraint time 1
constraint 51 1.time_range#1b.type#c.mu2_LB = 1.time_range#1b.type#c.mu2_FT
constraint 52 1.time_range#1b.type#c.mu2_FT = 1.time_range#1b.type#c.mu2_SD
constraint 53 1.time_range#1b.type#c.mu2_FT = 1.time_range#1b.type#c.mu2_FD
constraint 54 1.time_range#1b.type#c.mu2_FT = 1.time_range#1b.type#c.mu2_PR
constraint 55 1.time_range#1b.type#c.mu2_FT = 1.time_range#1b.type#c.mu2_VT

constraint 56 1.time_range#1b.type#c.mu3_LB = 1.time_range#1b.type#c.mu3_FT
constraint 57 1.time_range#1b.type#c.mu3_FT = 1.time_range#1b.type#c.mu3_SD
constraint 58 1.time_range#1b.type#c.mu3_FT = 1.time_range#1b.type#c.mu3_FD
constraint 59 1.time_range#1b.type#c.mu3_FT = 1.time_range#1b.type#c.mu3_PR
constraint 60 1.time_range#1b.type#c.mu3_FT = 1.time_range#1b.type#c.mu3_VT

*For type 1 constraint time 2
constraint 61 2.time_range#1b.type#c.mu2_LB = 2.time_range#1b.type#c.mu2_FT
constraint 62 2.time_range#1b.type#c.mu2_FT = 2.time_range#1b.type#c.mu2_SD
constraint 63 2.time_range#1b.type#c.mu2_FT = 2.time_range#1b.type#c.mu2_FD
constraint 64 2.time_range#1b.type#c.mu2_FT = 2.time_range#1b.type#c.mu2_PR
constraint 65 2.time_range#1b.type#c.mu2_FT = 2.time_range#1b.type#c.mu2_VT

constraint 66 2.time_range#1b.type#c.mu3_LB = 2.time_range#1b.type#c.mu3_FT
constraint 67 2.time_range#1b.type#c.mu3_FT = 2.time_range#1b.type#c.mu3_SD
constraint 68 2.time_range#1b.type#c.mu3_FT = 2.time_range#1b.type#c.mu3_FD
constraint 69 2.time_range#1b.type#c.mu3_FT = 2.time_range#1b.type#c.mu3_PR
constraint 70 2.time_range#1b.type#c.mu3_FT = 2.time_range#1b.type#c.mu3_VT

*For type 5 constraint time 0
constraint 71 0b.time_range#5.type#c.mu2_LB = 0b.time_range#5.type#c.mu2_FT
constraint 72 0b.time_range#5.type#c.mu2_FT = 0b.time_range#5.type#c.mu2_SD
constraint 73 0b.time_range#5.type#c.mu2_FT = 0b.time_range#5.type#c.mu2_FD
constraint 74 0b.time_range#5.type#c.mu2_FT = 0b.time_range#5.type#c.mu2_PR
constraint 75 0b.time_range#5.type#c.mu2_FT = 0b.time_range#5.type#c.mu2_VT

constraint 76 0b.time_range#5.type#c.mu3_LB = 0b.time_range#5.type#c.mu3_FT
constraint 77 0b.time_range#5.type#c.mu3_FT = 0b.time_range#5.type#c.mu3_SD
constraint 78 0b.time_range#5.type#c.mu3_FT = 0b.time_range#5.type#c.mu3_FD
constraint 79 0b.time_range#5.type#c.mu3_FT = 0b.time_range#5.type#c.mu3_PR
constraint 80 0b.time_range#5.type#c.mu3_FT = 0b.time_range#5.type#c.mu3_VT

*For type 5 constraint time 1
constraint 81 1.time_range#5.type#c.mu2_LB = 1.time_range#5.type#c.mu2_FT
constraint 82 1.time_range#5.type#c.mu2_FT = 1.time_range#5.type#c.mu2_SD
constraint 83 1.time_range#5.type#c.mu2_FT = 1.time_range#5.type#c.mu2_FD
constraint 84 1.time_range#5.type#c.mu2_FT = 1.time_range#5.type#c.mu2_PR
constraint 85 1.time_range#5.type#c.mu2_FT = 1.time_range#5.type#c.mu2_VT

constraint 86 1.time_range#5.type#c.mu3_LB = 1.time_range#5.type#c.mu3_FT
constraint 87 1.time_range#5.type#c.mu3_FT = 1.time_range#5.type#c.mu3_SD
constraint 88 1.time_range#5.type#c.mu3_FT = 1.time_range#5.type#c.mu3_FD
constraint 89 1.time_range#5.type#c.mu3_FT = 1.time_range#5.type#c.mu3_PR
constraint 90 1.time_range#5.type#c.mu3_FT = 1.time_range#5.type#c.mu3_VT

*For type 5 constraint time 2
constraint 91 2.time_range#5.type#c.mu2_LB = 2.time_range#5.type#c.mu2_FT
constraint 92 2.time_range#5.type#c.mu2_FT = 2.time_range#5.type#c.mu2_SD
constraint 93 2.time_range#5.type#c.mu2_FT = 2.time_range#5.type#c.mu2_FD
constraint 94 2.time_range#5.type#c.mu2_FT = 2.time_range#5.type#c.mu2_PR
constraint 95 2.time_range#5.type#c.mu2_FT = 2.time_range#5.type#c.mu2_VT

constraint 96 2.time_range#5.type#c.mu3_LB = 2.time_range#5.type#c.mu3_FT
constraint 97 2.time_range#5.type#c.mu3_FT = 2.time_range#5.type#c.mu3_SD
constraint 98 2.time_range#5.type#c.mu3_FT = 2.time_range#5.type#c.mu3_FD
constraint 99 2.time_range#5.type#c.mu3_FT = 2.time_range#5.type#c.mu3_PR
constraint 100 2.time_range#5.type#c.mu3_FT = 2.time_range#5.type#c.mu3_VT

sureg (mu1_FT i.time_range#i.type#(c.mu2_FT c.mu3_FT))(mu1_SD i.time_range#i.type#(c.mu2_SD c.mu3_SD))(mu1_FD i.time_range#i.type#(c.mu2_FD c.mu3_FD))(mu1_PR i.time_range#i.type#(c.mu2_PR c.mu3_PR))(mu1_VT i.time_range#i.type#(c.mu2_VT c.mu3_VT)), constraint (42 43 44 45 47 48 49 50 52 53 54 55 57 58 59 60 62 63 64 65 67 68 69 70 72 73 74 75 77 78 79 80 82 83 84 85 87 88 89 90 92 93 94 95 97 98 99 100)
estimates store type_sure

esttab type_sure using sure_type_untrimmed.csv, b(3) se(3) ar2 star(* 0.1 ** 0.05 *** 0.01) nogaps

*How farmers change to a price shock
*H0: crop farmer' risk coefficients are the same as ones of dairy farmers
test _b[mu1_FT:0b.time_range#5.type#c.mu2_FT]  = _b[mu1_FT:0b.time_range#1b.type#c.mu2_FT]
test _b[mu1_FT:1.time_range#5.type#c.mu2_FT]  = _b[mu1_FT:1.time_range#1b.type#c.mu2_FT]
test _b[mu1_FT:2.time_range#5.type#c.mu2_FT]  = _b[mu1_FT:2.time_range#1b.type#c.mu2_FT]

test _b[mu1_FT:0b.time_range#5.type#c.mu3_FT]  = _b[mu1_FT:0b.time_range#1b.type#c.mu3_FT]
test _b[mu1_FT:1.time_range#5.type#c.mu3_FT]  = _b[mu1_FT:1.time_range#1b.type#c.mu3_FT]
test _b[mu1_FT:2.time_range#5.type#c.mu3_FT]  = _b[mu1_FT:2.time_range#1b.type#c.mu3_FT]

*H0: crop farmers' AP risk aversion changes the same as dairy farmers
test _b[mu1_FT:0b.time_range#5.type#c.mu2_FT] - _b[mu1_FT:1.time_range#5.type#c.mu2_FT] = _b[mu1_FT:0b.time_range#1b.type#c.mu2_FT] - _b[mu1_FT:1.time_range#1b.type#c.mu2_FT]
test _b[mu1_FT:1.time_range#5.type#c.mu2_FT] - _b[mu1_FT:2.time_range#5.type#c.mu2_FT] = _b[mu1_FT:1.time_range#1b.type#c.mu2_FT] - _b[mu1_FT:2.time_range#1b.type#c.mu2_FT]

*H0: crop farmers' DS risk aversion changes the same as dairy farmers
test _b[mu1_FT:0b.time_range#5.type#c.mu3_FT] - _b[mu1_FT:1.time_range#5.type#c.mu3_FT] = _b[mu1_FT:0b.time_range#1b.type#c.mu3_FT] - _b[mu1_FT:1.time_range#1b.type#c.mu3_FT]
test _b[mu1_FT:1.time_range#5.type#c.mu3_FT] - _b[mu1_FT:2.time_range#5.type#c.mu3_FT] = _b[mu1_FT:1.time_range#1b.type#c.mu3_FT] - _b[mu1_FT:2.time_range#1b.type#c.mu3_FT]

*Country specific risk attitudes
encode country, gen(ncountry)
label list ncountry

constraint 101 1b.ncountry#c.mu2_FT = 1b.ncountry#c.mu2_LB
constraint 102 1b.ncountry#c.mu2_FT = 1b.ncountry#c.mu2_SD
constraint 103 1b.ncountry#c.mu2_FT = 1b.ncountry#c.mu2_FD
constraint 104 1b.ncountry#c.mu2_FT = 1b.ncountry#c.mu2_PR
constraint 105 1b.ncountry#c.mu2_FT = 1b.ncountry#c.mu2_VT

constraint 106 1b.ncountry#c.mu3_FT = 1b.ncountry#c.mu3_LB
constraint 107 1b.ncountry#c.mu3_FT = 1b.ncountry#c.mu3_SD
constraint 108 1b.ncountry#c.mu3_FT = 1b.ncountry#c.mu3_FD
constraint 109 1b.ncountry#c.mu3_FT = 1b.ncountry#c.mu3_PR
constraint 110 1b.ncountry#c.mu3_FT = 1b.ncountry#c.mu3_VT

constraint 111 2.ncountry#c.mu2_FT = 2.ncountry#c.mu2_LB
constraint 112 2.ncountry#c.mu2_FT = 2.ncountry#c.mu2_SD
constraint 113 2.ncountry#c.mu2_FT = 2.ncountry#c.mu2_FD
constraint 114 2.ncountry#c.mu2_FT = 2.ncountry#c.mu2_PR
constraint 115 2.ncountry#c.mu2_FT = 2.ncountry#c.mu2_VT

constraint 116 2.ncountry#c.mu3_FT = 2.ncountry#c.mu3_LB
constraint 117 2.ncountry#c.mu3_FT = 2.ncountry#c.mu3_SD
constraint 118 2.ncountry#c.mu3_FT = 2.ncountry#c.mu3_FD
constraint 119 2.ncountry#c.mu3_FT = 2.ncountry#c.mu3_PR
constraint 120 2.ncountry#c.mu3_FT = 2.ncountry#c.mu3_VT

constraint 121 3.ncountry#c.mu2_FT = 3.ncountry#c.mu2_LB
constraint 122 3.ncountry#c.mu2_FT = 3.ncountry#c.mu2_SD
constraint 123 3.ncountry#c.mu2_FT = 3.ncountry#c.mu2_FD
constraint 124 3.ncountry#c.mu2_FT = 3.ncountry#c.mu2_PR
constraint 125 3.ncountry#c.mu2_FT = 3.ncountry#c.mu2_VT

constraint 126 3.ncountry#c.mu3_FT = 3.ncountry#c.mu3_LB
constraint 127 3.ncountry#c.mu3_FT = 3.ncountry#c.mu3_SD
constraint 128 3.ncountry#c.mu3_FT = 3.ncountry#c.mu3_FD
constraint 129 3.ncountry#c.mu3_FT = 3.ncountry#c.mu3_PR
constraint 130 3.ncountry#c.mu3_FT = 3.ncountry#c.mu3_VT

constraint 131 4.ncountry#c.mu2_FT = 4.ncountry#c.mu2_LB
constraint 132 4.ncountry#c.mu2_FT = 4.ncountry#c.mu2_SD
constraint 133 4.ncountry#c.mu2_FT = 4.ncountry#c.mu2_FD
constraint 134 4.ncountry#c.mu2_FT = 4.ncountry#c.mu2_PR
constraint 135 4.ncountry#c.mu2_FT = 4.ncountry#c.mu2_VT

constraint 136 4.ncountry#c.mu3_FT = 4.ncountry#c.mu3_LB
constraint 137 4.ncountry#c.mu3_FT = 4.ncountry#c.mu3_SD
constraint 138 4.ncountry#c.mu3_FT = 4.ncountry#c.mu3_FD
constraint 139 4.ncountry#c.mu3_FT = 4.ncountry#c.mu3_PR
constraint 140 4.ncountry#c.mu3_FT = 4.ncountry#c.mu3_VT

constraint 141 5.ncountry#c.mu2_FT = 5.ncountry#c.mu2_LB
constraint 142 5.ncountry#c.mu2_FT = 5.ncountry#c.mu2_SD
constraint 143 5.ncountry#c.mu2_FT = 5.ncountry#c.mu2_FD
constraint 144 5.ncountry#c.mu2_FT = 5.ncountry#c.mu2_PR
constraint 145 5.ncountry#c.mu2_FT = 5.ncountry#c.mu2_VT

constraint 146 5.ncountry#c.mu3_FT = 5.ncountry#c.mu3_LB
constraint 147 5.ncountry#c.mu3_FT = 5.ncountry#c.mu3_SD
constraint 148 5.ncountry#c.mu3_FT = 5.ncountry#c.mu3_FD
constraint 149 5.ncountry#c.mu3_FT = 5.ncountry#c.mu3_PR
constraint 150 5.ncountry#c.mu3_FT = 5.ncountry#c.mu3_VT

constraint 151 6.ncountry#c.mu2_FT = 6.ncountry#c.mu2_LB
constraint 152 6.ncountry#c.mu2_FT = 6.ncountry#c.mu2_SD
constraint 153 6.ncountry#c.mu2_FT = 6.ncountry#c.mu2_FD
constraint 154 6.ncountry#c.mu2_FT = 6.ncountry#c.mu2_PR
constraint 155 6.ncountry#c.mu2_FT = 6.ncountry#c.mu2_VT

constraint 156 6.ncountry#c.mu3_FT = 6.ncountry#c.mu3_LB
constraint 157 6.ncountry#c.mu3_FT = 6.ncountry#c.mu3_SD
constraint 158 6.ncountry#c.mu3_FT = 6.ncountry#c.mu3_FD
constraint 159 6.ncountry#c.mu3_FT = 6.ncountry#c.mu3_PR
constraint 160 6.ncountry#c.mu3_FT = 6.ncountry#c.mu3_VT

constraint 161 7.ncountry#c.mu2_FT = 7.ncountry#c.mu2_LB
constraint 162 7.ncountry#c.mu2_FT = 7.ncountry#c.mu2_SD
constraint 163 7.ncountry#c.mu2_FT = 7.ncountry#c.mu2_FD
constraint 164 7.ncountry#c.mu2_FT = 7.ncountry#c.mu2_PR
constraint 165 7.ncountry#c.mu2_FT = 7.ncountry#c.mu2_VT

constraint 166 7.ncountry#c.mu3_FT = 7.ncountry#c.mu3_LB
constraint 167 7.ncountry#c.mu3_FT = 7.ncountry#c.mu3_SD
constraint 168 7.ncountry#c.mu3_FT = 7.ncountry#c.mu3_FD
constraint 169 7.ncountry#c.mu3_FT = 7.ncountry#c.mu3_PR
constraint 170 7.ncountry#c.mu3_FT = 7.ncountry#c.mu3_VT

constraint 171 8.ncountry#c.mu2_FT = 8.ncountry#c.mu2_LB
constraint 172 8.ncountry#c.mu2_FT = 8.ncountry#c.mu2_SD
constraint 173 8.ncountry#c.mu2_FT = 8.ncountry#c.mu2_FD
constraint 174 8.ncountry#c.mu2_FT = 8.ncountry#c.mu2_PR
constraint 175 8.ncountry#c.mu2_FT = 8.ncountry#c.mu2_VT

constraint 176 8.ncountry#c.mu3_FT = 8.ncountry#c.mu3_LB
constraint 177 8.ncountry#c.mu3_FT = 8.ncountry#c.mu3_SD
constraint 178 8.ncountry#c.mu3_FT = 8.ncountry#c.mu3_FD
constraint 179 8.ncountry#c.mu3_FT = 8.ncountry#c.mu3_PR
constraint 180 8.ncountry#c.mu3_FT = 8.ncountry#c.mu3_VT

constraint 181 9.ncountry#c.mu2_FT = 9.ncountry#c.mu2_LB
constraint 182 9.ncountry#c.mu2_FT = 9.ncountry#c.mu2_SD
constraint 183 9.ncountry#c.mu2_FT = 9.ncountry#c.mu2_FD
constraint 184 9.ncountry#c.mu2_FT = 9.ncountry#c.mu2_PR
constraint 185 9.ncountry#c.mu2_FT = 9.ncountry#c.mu2_VT

constraint 186 9.ncountry#c.mu3_FT = 9.ncountry#c.mu3_LB
constraint 187 9.ncountry#c.mu3_FT = 9.ncountry#c.mu3_SD
constraint 188 9.ncountry#c.mu3_FT = 9.ncountry#c.mu3_FD
constraint 189 9.ncountry#c.mu3_FT = 9.ncountry#c.mu3_PR
constraint 190 9.ncountry#c.mu3_FT = 9.ncountry#c.mu3_VT

constraint 191 10.ncountry#c.mu2_FT = 10.ncountry#c.mu2_LB
constraint 192 10.ncountry#c.mu2_FT = 10.ncountry#c.mu2_SD
constraint 193 10.ncountry#c.mu2_FT = 10.ncountry#c.mu2_FD
constraint 194 10.ncountry#c.mu2_FT = 10.ncountry#c.mu2_PR
constraint 195 10.ncountry#c.mu2_FT = 10.ncountry#c.mu2_VT

constraint 196 10.ncountry#c.mu3_FT = 10.ncountry#c.mu3_LB
constraint 197 10.ncountry#c.mu3_FT = 10.ncountry#c.mu3_SD
constraint 198 10.ncountry#c.mu3_FT = 10.ncountry#c.mu3_FD
constraint 199 10.ncountry#c.mu3_FT = 10.ncountry#c.mu3_PR
constraint 200 10.ncountry#c.mu3_FT = 10.ncountry#c.mu3_VT

constraint 201 11.ncountry#c.mu2_FT = 11.ncountry#c.mu2_LB
constraint 202 11.ncountry#c.mu2_FT = 11.ncountry#c.mu2_SD
constraint 203 11.ncountry#c.mu2_FT = 11.ncountry#c.mu2_FD
constraint 204 11.ncountry#c.mu2_FT = 11.ncountry#c.mu2_PR
constraint 205 11.ncountry#c.mu2_FT = 11.ncountry#c.mu2_VT

constraint 206 11.ncountry#c.mu3_FT = 11.ncountry#c.mu3_LB
constraint 207 11.ncountry#c.mu3_FT = 11.ncountry#c.mu3_SD
constraint 208 11.ncountry#c.mu3_FT = 11.ncountry#c.mu3_FD
constraint 209 11.ncountry#c.mu3_FT = 11.ncountry#c.mu3_PR
constraint 210 11.ncountry#c.mu3_FT = 11.ncountry#c.mu3_VT

constraint 211 12.ncountry#c.mu2_FT = 12.ncountry#c.mu2_LB
constraint 212 12.ncountry#c.mu2_FT = 12.ncountry#c.mu2_SD
constraint 213 12.ncountry#c.mu2_FT = 12.ncountry#c.mu2_FD
constraint 214 12.ncountry#c.mu2_FT = 12.ncountry#c.mu2_PR
constraint 215 12.ncountry#c.mu2_FT = 12.ncountry#c.mu2_VT

constraint 216 12.ncountry#c.mu3_FT = 12.ncountry#c.mu3_LB
constraint 217 12.ncountry#c.mu3_FT = 12.ncountry#c.mu3_SD
constraint 218 12.ncountry#c.mu3_FT = 12.ncountry#c.mu3_FD
constraint 219 12.ncountry#c.mu3_FT = 12.ncountry#c.mu3_PR
constraint 220 12.ncountry#c.mu3_FT = 12.ncountry#c.mu3_VT

constraint 221 13.ncountry#c.mu2_FT = 13.ncountry#c.mu2_LB
constraint 222 13.ncountry#c.mu2_FT = 13.ncountry#c.mu2_SD
constraint 223 13.ncountry#c.mu2_FT = 13.ncountry#c.mu2_FD
constraint 224 13.ncountry#c.mu2_FT = 13.ncountry#c.mu2_PR
constraint 225 13.ncountry#c.mu2_FT = 13.ncountry#c.mu2_VT

constraint 226 13.ncountry#c.mu3_FT = 13.ncountry#c.mu3_LB
constraint 227 13.ncountry#c.mu3_FT = 13.ncountry#c.mu3_SD
constraint 228 13.ncountry#c.mu3_FT = 13.ncountry#c.mu3_FD
constraint 229 13.ncountry#c.mu3_FT = 13.ncountry#c.mu3_PR
constraint 230 13.ncountry#c.mu3_FT = 13.ncountry#c.mu3_VT

constraint 231 14.ncountry#c.mu2_FT = 14.ncountry#c.mu2_LB
constraint 232 14.ncountry#c.mu2_FT = 14.ncountry#c.mu2_SD
constraint 233 14.ncountry#c.mu2_FT = 14.ncountry#c.mu2_FD
constraint 234 14.ncountry#c.mu2_FT = 14.ncountry#c.mu2_PR
constraint 235 14.ncountry#c.mu2_FT = 14.ncountry#c.mu2_VT

constraint 236 14.ncountry#c.mu3_FT = 14.ncountry#c.mu3_LB
constraint 237 14.ncountry#c.mu3_FT = 14.ncountry#c.mu3_SD
constraint 238 14.ncountry#c.mu3_FT = 14.ncountry#c.mu3_FD
constraint 239 14.ncountry#c.mu3_FT = 14.ncountry#c.mu3_PR
constraint 240 14.ncountry#c.mu3_FT = 14.ncountry#c.mu3_VT

constraint 241 15.ncountry#c.mu2_FT = 15.ncountry#c.mu2_LB
constraint 242 15.ncountry#c.mu2_FT = 15.ncountry#c.mu2_SD
constraint 243 15.ncountry#c.mu2_FT = 15.ncountry#c.mu2_FD
constraint 244 15.ncountry#c.mu2_FT = 15.ncountry#c.mu2_PR
constraint 245 15.ncountry#c.mu2_FT = 15.ncountry#c.mu2_VT

constraint 246 15.ncountry#c.mu3_FT = 15.ncountry#c.mu3_LB
constraint 247 15.ncountry#c.mu3_FT = 15.ncountry#c.mu3_SD
constraint 248 15.ncountry#c.mu3_FT = 15.ncountry#c.mu3_FD
constraint 249 15.ncountry#c.mu3_FT = 15.ncountry#c.mu3_PR
constraint 250 15.ncountry#c.mu3_FT = 15.ncountry#c.mu3_VT

constraint 251 16.ncountry#c.mu2_FT = 16.ncountry#c.mu2_LB
constraint 252 16.ncountry#c.mu2_FT = 16.ncountry#c.mu2_SD
constraint 253 16.ncountry#c.mu2_FT = 16.ncountry#c.mu2_FD
constraint 254 16.ncountry#c.mu2_FT = 16.ncountry#c.mu2_PR
constraint 255 16.ncountry#c.mu2_FT = 16.ncountry#c.mu2_VT

constraint 256 16.ncountry#c.mu3_FT = 16.ncountry#c.mu3_LB
constraint 257 16.ncountry#c.mu3_FT = 16.ncountry#c.mu3_SD
constraint 258 16.ncountry#c.mu3_FT = 16.ncountry#c.mu3_FD
constraint 259 16.ncountry#c.mu3_FT = 16.ncountry#c.mu3_PR
constraint 260 16.ncountry#c.mu3_FT = 16.ncountry#c.mu3_VT

constraint 261 17.ncountry#c.mu2_FT = 17.ncountry#c.mu2_LB
constraint 262 17.ncountry#c.mu2_FT = 17.ncountry#c.mu2_SD
constraint 263 17.ncountry#c.mu2_FT = 17.ncountry#c.mu2_FD
constraint 264 17.ncountry#c.mu2_FT = 17.ncountry#c.mu2_PR
constraint 265 17.ncountry#c.mu2_FT = 17.ncountry#c.mu2_VT

constraint 266 17.ncountry#c.mu3_FT = 17.ncountry#c.mu3_LB
constraint 267 17.ncountry#c.mu3_FT = 17.ncountry#c.mu3_SD
constraint 268 17.ncountry#c.mu3_FT = 17.ncountry#c.mu3_FD
constraint 269 17.ncountry#c.mu3_FT = 17.ncountry#c.mu3_PR
constraint 270 17.ncountry#c.mu3_FT = 17.ncountry#c.mu3_VT

constraint 271 18.ncountry#c.mu2_FT = 18.ncountry#c.mu2_LB
constraint 272 18.ncountry#c.mu2_FT = 18.ncountry#c.mu2_SD
constraint 273 18.ncountry#c.mu2_FT = 18.ncountry#c.mu2_FD
constraint 274 18.ncountry#c.mu2_FT = 18.ncountry#c.mu2_PR
constraint 275 18.ncountry#c.mu2_FT = 18.ncountry#c.mu2_VT

constraint 276 18.ncountry#c.mu3_FT = 18.ncountry#c.mu3_LB
constraint 277 18.ncountry#c.mu3_FT = 18.ncountry#c.mu3_SD
constraint 278 18.ncountry#c.mu3_FT = 18.ncountry#c.mu3_FD
constraint 279 18.ncountry#c.mu3_FT = 18.ncountry#c.mu3_PR
constraint 280 18.ncountry#c.mu3_FT = 18.ncountry#c.mu3_VT

constraint 281 19.ncountry#c.mu2_FT = 19.ncountry#c.mu2_LB
constraint 282 19.ncountry#c.mu2_FT = 19.ncountry#c.mu2_SD
constraint 283 19.ncountry#c.mu2_FT = 19.ncountry#c.mu2_FD
constraint 284 19.ncountry#c.mu2_FT = 19.ncountry#c.mu2_PR
constraint 285 19.ncountry#c.mu2_FT = 19.ncountry#c.mu2_VT

constraint 286 19.ncountry#c.mu3_FT = 19.ncountry#c.mu3_LB
constraint 287 19.ncountry#c.mu3_FT = 19.ncountry#c.mu3_SD
constraint 288 19.ncountry#c.mu3_FT = 19.ncountry#c.mu3_FD
constraint 289 19.ncountry#c.mu3_FT = 19.ncountry#c.mu3_PR
constraint 290 19.ncountry#c.mu3_FT = 19.ncountry#c.mu3_VT

constraint 291 20.ncountry#c.mu2_FT = 20.ncountry#c.mu2_LB
constraint 292 20.ncountry#c.mu2_FT = 20.ncountry#c.mu2_SD
constraint 293 20.ncountry#c.mu2_FT = 20.ncountry#c.mu2_FD
constraint 294 20.ncountry#c.mu2_FT = 20.ncountry#c.mu2_PR
constraint 295 20.ncountry#c.mu2_FT = 20.ncountry#c.mu2_VT

constraint 296 20.ncountry#c.mu3_FT = 20.ncountry#c.mu3_LB
constraint 297 20.ncountry#c.mu3_FT = 20.ncountry#c.mu3_SD
constraint 298 20.ncountry#c.mu3_FT = 20.ncountry#c.mu3_FD
constraint 299 20.ncountry#c.mu3_FT = 20.ncountry#c.mu3_PR
constraint 300 20.ncountry#c.mu3_FT = 20.ncountry#c.mu3_VT

constraint 301 21.ncountry#c.mu2_FT = 21.ncountry#c.mu2_LB
constraint 302 21.ncountry#c.mu2_FT = 21.ncountry#c.mu2_SD
constraint 303 21.ncountry#c.mu2_FT = 21.ncountry#c.mu2_FD
constraint 304 21.ncountry#c.mu2_FT = 21.ncountry#c.mu2_PR
constraint 305 21.ncountry#c.mu2_FT = 21.ncountry#c.mu2_VT

constraint 306 21.ncountry#c.mu3_FT = 21.ncountry#c.mu3_LB
constraint 307 21.ncountry#c.mu3_FT = 21.ncountry#c.mu3_SD
constraint 308 21.ncountry#c.mu3_FT = 21.ncountry#c.mu3_FD
constraint 309 21.ncountry#c.mu3_FT = 21.ncountry#c.mu3_PR
constraint 310 21.ncountry#c.mu3_FT = 21.ncountry#c.mu3_VT

constraint 311 22.ncountry#c.mu2_FT = 22.ncountry#c.mu2_LB
constraint 312 22.ncountry#c.mu2_FT = 22.ncountry#c.mu2_SD
constraint 313 22.ncountry#c.mu2_FT = 22.ncountry#c.mu2_FD
constraint 314 22.ncountry#c.mu2_FT = 22.ncountry#c.mu2_PR
constraint 315 22.ncountry#c.mu2_FT = 22.ncountry#c.mu2_VT

constraint 316 22.ncountry#c.mu3_FT = 22.ncountry#c.mu3_LB
constraint 317 22.ncountry#c.mu3_FT = 22.ncountry#c.mu3_SD
constraint 318 22.ncountry#c.mu3_FT = 22.ncountry#c.mu3_FD
constraint 319 22.ncountry#c.mu3_FT = 22.ncountry#c.mu3_PR
constraint 320 22.ncountry#c.mu3_FT = 22.ncountry#c.mu3_VT

constraint 321 23.ncountry#c.mu2_FT = 23.ncountry#c.mu2_LB
constraint 322 23.ncountry#c.mu2_FT = 23.ncountry#c.mu2_SD
constraint 323 23.ncountry#c.mu2_FT = 23.ncountry#c.mu2_FD
constraint 324 23.ncountry#c.mu2_FT = 23.ncountry#c.mu2_PR
constraint 325 23.ncountry#c.mu2_FT = 23.ncountry#c.mu2_VT

constraint 326 23.ncountry#c.mu3_FT = 23.ncountry#c.mu3_LB
constraint 327 23.ncountry#c.mu3_FT = 23.ncountry#c.mu3_SD
constraint 328 23.ncountry#c.mu3_FT = 23.ncountry#c.mu3_FD
constraint 329 23.ncountry#c.mu3_FT = 23.ncountry#c.mu3_PR
constraint 330 23.ncountry#c.mu3_FT = 23.ncountry#c.mu3_VT

constraint 331 24.ncountry#c.mu2_FT = 24.ncountry#c.mu2_LB
constraint 332 24.ncountry#c.mu2_FT = 24.ncountry#c.mu2_SD
constraint 333 24.ncountry#c.mu2_FT = 24.ncountry#c.mu2_FD
constraint 334 24.ncountry#c.mu2_FT = 24.ncountry#c.mu2_PR
constraint 335 24.ncountry#c.mu2_FT = 24.ncountry#c.mu2_VT

constraint 336 24.ncountry#c.mu3_FT = 24.ncountry#c.mu3_LB
constraint 337 24.ncountry#c.mu3_FT = 24.ncountry#c.mu3_SD
constraint 338 24.ncountry#c.mu3_FT = 24.ncountry#c.mu3_FD
constraint 339 24.ncountry#c.mu3_FT = 24.ncountry#c.mu3_PR
constraint 340 24.ncountry#c.mu3_FT = 24.ncountry#c.mu3_VT

constraint 341 25.ncountry#c.mu2_FT = 25.ncountry#c.mu2_LB
constraint 342 25.ncountry#c.mu2_FT = 25.ncountry#c.mu2_SD
constraint 343 25.ncountry#c.mu2_FT = 25.ncountry#c.mu2_FD
constraint 344 25.ncountry#c.mu2_FT = 25.ncountry#c.mu2_PR
constraint 345 25.ncountry#c.mu2_FT = 25.ncountry#c.mu2_VT

constraint 346 25.ncountry#c.mu3_FT = 25.ncountry#c.mu3_LB
constraint 347 25.ncountry#c.mu3_FT = 25.ncountry#c.mu3_SD
constraint 348 25.ncountry#c.mu3_FT = 25.ncountry#c.mu3_FD
constraint 349 25.ncountry#c.mu3_FT = 25.ncountry#c.mu3_PR
constraint 350 25.ncountry#c.mu3_FT = 25.ncountry#c.mu3_VT

constraint 351 26.ncountry#c.mu2_FT = 26.ncountry#c.mu2_LB
constraint 352 26.ncountry#c.mu2_FT = 26.ncountry#c.mu2_SD
constraint 353 26.ncountry#c.mu2_FT = 26.ncountry#c.mu2_FD
constraint 354 26.ncountry#c.mu2_FT = 26.ncountry#c.mu2_PR
constraint 355 26.ncountry#c.mu2_FT = 26.ncountry#c.mu2_VT

constraint 356 26.ncountry#c.mu3_FT = 26.ncountry#c.mu3_LB
constraint 357 26.ncountry#c.mu3_FT = 26.ncountry#c.mu3_SD
constraint 358 26.ncountry#c.mu3_FT = 26.ncountry#c.mu3_FD
constraint 359 26.ncountry#c.mu3_FT = 26.ncountry#c.mu3_PR
constraint 360 26.ncountry#c.mu3_FT = 26.ncountry#c.mu3_VT

constraint 361 27.ncountry#c.mu2_FT = 27.ncountry#c.mu2_LB
constraint 362 27.ncountry#c.mu2_FT = 27.ncountry#c.mu2_SD
constraint 363 27.ncountry#c.mu2_FT = 27.ncountry#c.mu2_FD
constraint 364 27.ncountry#c.mu2_FT = 27.ncountry#c.mu2_PR
constraint 365 27.ncountry#c.mu2_FT = 27.ncountry#c.mu2_VT

constraint 366 27.ncountry#c.mu3_FT = 27.ncountry#c.mu3_LB
constraint 367 27.ncountry#c.mu3_FT = 27.ncountry#c.mu3_SD
constraint 368 27.ncountry#c.mu3_FT = 27.ncountry#c.mu3_FD
constraint 369 27.ncountry#c.mu3_FT = 27.ncountry#c.mu3_PR
constraint 370 27.ncountry#c.mu3_FT = 27.ncountry#c.mu3_VT

sureg (mu1_FT i.ncountry#(c.mu2_FT c.mu3_FT))(mu1_SD i.ncountry#(c.mu2_SD c.mu3_SD))(mu1_FD i.ncountry#(c.mu2_FD c.mu3_FD))(mu1_PR i.ncountry#(c.mu2_PR c.mu3_PR))(mu1_VT i.ncountry#(c.mu2_VT c.mu3_VT)), constraint(102 103 104 105 107 108 109 110 112 113 114 115 117 118 119 120 122 123 124 125 127 128 129 130 132 133 134 135 137 138 139 140 142 143 144 145 147 148 149 150 152 153 154 155 157 158 159 160 162 163 164 165 167 168 169 170 172 173 174 175 177 178 179 180 182 183 184 185 187 188 189 190 192 193 194 195 197 198 199 200 202 203 204 205 207 208 209 210 212 213 214 215 217 218 219 220 222 223 224 225 227 228 229 230 232 233 234 235 237 238 239 240 242 243 244 245 247 248 249 250 252 253 254 255 257 258 259 260 262 263 264 265 267 268 269 270 272 273 274 275 277 278 279 280 282 283 284 285 287 288 289 290 292 293 294 295 297 298 299 300 302 303 304 305 307 308 309 310 312 313 314 315 317 318 319 320 322 323 324 325 327 328 329 330 332 333 334 335 337 338 339 340 342 343 344 345 347 348 349 350 352 353 354 355 357 358 359 360 362 363 364 365 367 368 369 370) coefl

estimates store country_sure
esttab country_sure using sure_country_untrimmed.csv, b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) nogaps

foreach y in 1b 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 {
foreach x in 1b 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 {
test _b[mu1_FT:`y'.ncountry#c.mu2_FT]  = _b[mu1_FT:`x'.ncountry#c.mu2_FT]
}
}

foreach y in 1b 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 {
foreach x in 1b 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 {
test _b[mu1_FT:`y'.ncountry#c.mu3_FT]  = _b[mu1_FT:`x'.ncountry#c.mu3_FT]
}
}

foreach x in 27 {
test _b[mu1_FT:26.ncountry#c.mu3_FT]  = _b[mu1_FT:`x'.ncountry#c.mu3_FT]
}
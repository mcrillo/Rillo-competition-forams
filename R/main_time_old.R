#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################


# OLD 

##################
###### DATA ######
##################

# General data
species <- read.csv2(file="data/ssp_traps.csv", header=TRUE, stringsAsFactors = FALSE) # species names
traps_meta <- read.csv(file="data/traps_metadata.csv", header=TRUE, stringsAsFactors = FALSE) # species names
traps_names <- traps_meta$trap[-which(traps_meta$trap %in% c("WS1", "WS34"))] # removing metadata and trap with just one ssp (N_pac)
get_flux_info(traps_names, overwrite = FALSE) # Creating CSV files with additional information on data sampling period & year (data/traps_time/)


##############################
### ANALYSIS: Co-occurence ###
##############################

# Calculating species-pairs c-score within each sediment trap
list_c_score <- get_c_score(traps_names) # data.frame, string
# sister-species matrix (binary)
sisters_matrix <- get_sister_ssp(species) # data.frame, string
# get data.frame with ssp_pairs c-score and sister species info
sister_c_score <- get_sister_c_score(list_c_score, sisters_matrix) # data.frame, string
ggplot(sister_c_score, aes(c_score, fill = ssp_pair)) + geom_density(alpha = 0.2)
# subset for sediment traps with shorter time resolution
traps_short <- traps_meta[which(traps_meta[,6]<14),2]
sister_c_score_short <- sister_c_score[which(sister_c_score$trap %in% traps_short),]
ggplot(sister_c_score_short, aes(c_score, fill = ssp_pair)) + geom_density(alpha = 0.2)

# CHECK
# get_flux_overlap.R
# get_ranges_overlap.R


#########################################
### ANALYSIS: Time-series correlation ###
#########################################

# Calculating first differences of time-series (data/traps_first_diff/)
get_flux_first_diff(traps_names, days_closed = 10, overwrite = FALSE) 

# Plotting time-series, first diff series, and autocorrelation function (of first diffs)
plot_time_series(traps_names, species, overwrite = FALSE)

###########################################

# Colors for plot
cores21 <- c("#7c3f27","#7144ca","#c1dd44","#c949bd","#67ce59","#4b2d6b","#ccab3e","#6479c8","#d44430","#78dcc1","#ce3f71","#4d9060","#c588c4","#c3d088","#682b42","#7aa7c2","#cd7f41","#3d4342","#cf7a80","#5f692d","#ccb6a8")
cores10 <- c("#0037c6","#00b5ff","#004d41","#d0be00","#86442b","#db0011","#00b900","#8000b2","#ff60ff","#ff7314")
cores8 <- c("#00b900","#00b5ff","#0037c6","#8000b2","#86442b","#ff7314","#db0011","#ff60ff")

ssp_names <- read.csv("data/species_list.csv", header=TRUE)
data <- read.csv(paste("data/traps_time/GOM.csv",sep=""), header=TRUE)
trap_doy <- get_trap_doy(data, ssp_names[,"name_abbr"], trap_name = c("GOM"), overwrite = F)
trap_doy_na <- add_na_series(trap_doy, ssp_names[,"name_abbr"], trap_name = c("GOM"), overwrite = T)

################################################################################################
# main_traps

# Colors for plot
cores21 <- c("#7c3f27","#7144ca","#c1dd44","#c949bd","#67ce59","#4b2d6b","#ccab3e","#6479c8","#d44430","#78dcc1","#ce3f71","#4d9060","#c588c4","#c3d088","#682b42","#7aa7c2","#cd7f41","#3d4342","#cf7a80","#5f692d","#ccb6a8")
cores10 <- c("#0037c6","#00b5ff","#004d41","#d0be00","#86442b","#db0011","#00b900","#8000b2","#ff60ff","#ff7314")
cores8 <- c("#00b900","#00b5ff","#0037c6","#8000b2","#86442b","#ff7314","#db0011","#ff60ff")
plot(1:length(cores8), 1:length(cores8), col=cores8, pch=19, cex=3, xlab="", ylab="")


###############
### Max flux & zero flux for each ssp for each trap (for species trait matrix)
###############
max_mat <- get_flux_max(traps_names, species$abbr) 
max_mat <- data.frame(apply(max_mat, 2, function(x) max(x, na.rm=TRUE)))
zeros_mat <- get_flux_zeros(traps_names, species$abbr) 
zeros_mat <- data.frame(apply(zeros_mat, 2, function(x) median(x, na.rm=TRUE)))


###############
### Co-occurence matrix (considering time!)
###############
occur_mat <- get_flux_overlap(species$abbr, traps_names)
ssp_allop <- get_allopatry(occur_mat)
plot_flux_overlap(occur_mat, order_clust = c("original")) # "AOE", "FPC", "hclust", "alphabet"

# Ordering occur_mat by sampling, plot poster Evolution 2017
trait_mat_full <- read.csv("data/species_traits/traits_matrix.csv", header = TRUE)
range <- trait_mat_full[,c("abbr", "occupancy","tow_max","rel_abund_max","flux_max", "flux_zeros", "flux_sampling")]
srange <- range[which(range$abbr %in% row.names(occur_mat)),]
occur_mat2 <- occur_mat[order(srange$flux_zeros),order(srange$flux_zeros)]
occur_mat2 <- occur_mat[order(srange$flux_sampling, decreasing = TRUE),order(srange$flux_sampling, decreasing = TRUE)]
# plot_flux_overlap()

# no idea whats this! allo for allopatry?
allossp <- c("G_tru", "G_ada", "G_rubp", "S_deh", "H_dig", "T_par", "G_cra", "T_hum", "H_pel", "T_iot")
allorange <- range[which(range$abbr %in% allossp),]
allorange <- allorange[order(allorange$occupancy),]


###############
### Correlations (species pairs, first diff series)
###############
cormethod = c("kendall") # c("pearson", "kendall", "spearman")

# Calculating the first diff time-series correlation of each species pair of each sediment trap
get_flux_correlation(species, traps_names, cormethod) # creates array.RData of dim(ssp, ssp, traps)

# Plots of ssp X ssp correlation coeficients 
# dir.create(file.path("output_time/ssp_correlation/", paste(cormethod,"_plots", sep="")))
plot_flux_correlation(species, traps_names, cormethod, p = 0.05) # plots sspXssp correlation matrix for each sediment trap
corr_temp <- get_flux_corr_sst(traps_meta, traps_names, cormethod, p = 0.05) # gets significant correlations for each trap and its annual mean temperature (traps_metadata)
plot_flux_corr_sst(corr_temp) # plots: correlation ~ temp, prop_sig_corr ~ temp, richenss ~ prop_sig_corr

# Checking the consistency of ssp X ssp correlation across sediment traps (e.g. does ssp_i and ssp_j always correlate positively?)
corrsist <- run_flux_corrsistency(cormethod, p = 0.05)
plot_flux_corrsistency(corrsist)






########################################################################

# STILL WORKING (MARINA) - Cluster analysis
# Check what 'daisy' does exactly
corrsistx <- corrsist
corrsistx[which(is.na(corrsistx))] <- 0 # tranforming NAs in 0 for cluster analysis
diss <- daisy(corrsist, metric = "manhattan")
clust.res <- hclust(diss)
plot(clust.res)


# STILL WORKING (MARINA) - More auto-correlation stuff
# library(tseries)
# library(broom) # tidy
# x <- scale(data[,species])
# Kwiatkowski-Phillips-Schmidt-Shin test (non-parametric): H0 = stationarity
# tidy(kpss.test(x))
# Augmented Dickey–Fuller t-statistic: H0 = NOT stationary
# tidy(adf.test(x, alternative = "stationary"))
# Seasonal Decomposition time-series (email Renata)
# decomp <- decompose(stand(data[,7]))





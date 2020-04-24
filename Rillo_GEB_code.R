#====================================================
# PhD project competition-forams
#
# Code related to the paper:
# Rillo et al. "On the mismatch in the strength of competition among fossil and modern species of planktonic Foraminifera." 
# Global Ecol Biogeogr. 2019; 28: 1866– 1878. https://doi.org/10.1111/geb.13000
#
# GitHub: https://github.com/mcrillo/Rillo-competition-forams
#
# Author: Marina Costa Rillo
# Date: 30/01/2019
#====================================================

rm(list=ls())

# Working directory 
setwd("./Rillo-competition-forams")

# Output
if (!file.exists("output/")){ dir.create("output/")}


# Packages (with function that installs them if needed)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <-c('corrplot', # matrix correlation plot
             'dplyr',
             'GGally',
             'ggbeeswarm', # geom_quasirandom
             'ggplot2',
             'ggrepel', # geom_text_repel
             'ggforce', # geom_sina
             'igraph',
             'lubridate', # trap dates
             'mapproj', # map.grid
             'paleotree', # dropExtinct
             'phytools', # plotTree
             'plyr',
             'picante', # community phylogenetics
             'reshape', # function melt
             'tseries',
             'viridis') #  gradient color plot

ipak(packages)


#*************************************************************************************************
### DATA ***
#***********

## Relative abundance data (spatial) ------------------------------------------------------------------------

# Coretop data from ForCenS (Siccha & Kucera 2017) https://doi.org/10.1038/sdata.2017.109
# with added WOA SST https://www.nodc.noaa.gov/OC5/woa13/woa13data.html
d <- read.csv(file ="data/ForCenS_woa.csv", header=TRUE, stringsAsFactors = FALSE) 

# Remove dissolution-prone samples sites from >3,500m in the Pacific and Indian oceans and from >4,500m elsewhere were excluded (Tittensor et al. 2010 Nature)
d <- d[-which(d$Water_depth > 4500),] # 7 North Atlantic; 11 South Atlantic ; 513 Artic ; 1025 Mediterranean ; # 257 Southern
d <- d[-which(d$Water_depth > 3500 & d$Ocean %in% c(49,81,129, 2049)),] # 49 North Pacific; 81 South Pacific; 129 Indian; 2049 Red Sea

# Remove communities with conflicting species-level identification
remove <- c()
for(i in grep("_._", colnames(d))){
  if(any(d[,i]==0, na.rm = T)){ d[which(d[,i]==0),i] <- NA }
  remove <- c(remove,which(!is.na(d[,i])))
} # d[remove,grep("_._", colnames(d))]
d <- d[-remove,]
rm(remove)

# Fixing Globigerinoides tenellus to Globoturborotalita tenella
colnames(d)[grep("tenell",colnames(d))] <- "Globoturborotalita_tenella" 


## Sediment trap data (temporal) ------------------------------------------------------------------------
# Temporal: sediment trap data from 31 studies, see "data/traps_metadata.csv" or paper for full reference
traps <- sub(x=list.files(path = "data/traps_time/")[grep(list.files(path = "data/traps_time/"),pattern=".csv")], pattern=".csv", replacement = "")
traps_meta <- read.csv(file="data/traps_metadata.csv", header=TRUE, stringsAsFactors = FALSE)


## Phylogeny data (macroperforates) ---------------------------------------------------------------------

# Modified from Aze et al. 2011  https://doi.org/10.1111/j.1469-185X.2011.00178.x
phy <- read.tree("data/aLb_Rillo_modified.tre")
phy <- dropExtinct(phy)

png(file = "output/si_fig_phylogeny.png", width = 6, height = 5, unit = "in", res = 300)
 plotTree(phy)
dev.off()


## Shell diameter data  ------------------------------------------------------------------------

diam <- read.csv(file ="data/shell_diameter_data.csv", header=TRUE, stringsAsFactors = FALSE) 
sp_diam <- ddply(diam,~species,summarise,mean=mean(max_diam_mu),sd=sd(max_diam_mu), n_ind = length(max_diam_mu))
diam <- merge(diam, sp_diam)

size_boxplot <- ggplot(diam, aes(x=species_name2, y=max_diam_mu)) + 
  geom_boxplot() + 
  labs(y = expression(paste("Shell diameter (",mu,m,")", sep="")),
       x = "Species") +
  theme(axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12, face="italic", angle = 45, hjust = 1))

png(file = "output/si_fig_size_boxplot.png", width = 10, height = 7, unit = "in", res = 300)
 print(size_boxplot)
dev.off()

# Adding "Globigerinoides_white" as ruber to trait matrix (to match 'comm' matrix)
sp_diam <- join(sp_diam, diam[,c("species", "species_name2")], by = "species", match = "first") # G. elongatus and ruber merged
sp_diam <- rename(sp_diam, replace = c("species_name2" = "sspname"))
sp_diam[which(sp_diam$species == "ruber"),"sspname"] ="Globigerinoides_ruber" # pink ruber and white ruber same size
sp_diam <- rbind(sp_diam, c(sp_diam[which(sp_diam$species == "ruber"),c("species","mean","sd","n_ind")],sspname="Globigerinoides_white"))



#***************************************************************************************************
### Spatial analysis ***
#***********************

# Subseting data to species and lat, long and SST (woa)
comm <- as.matrix(d[,c(5,6,68,69,22:62)]) 
# Tranforming NA into zero, otherwise function in analysis returns NA
comm[which(is.na(comm))] <- 0 

# Community dispersion metrics to be analysed
ses_mntd <- data.frame()
ses_mpd <- data.frame()

# Calculating ecological distance matrices
phydist <- cophenetic(phy)
traitdist <- as.matrix(dist(data.frame(sp_diam$mean, row.names = sp_diam$sspname))) 

# Setting variables
null <- c("taxa","richness") # "sample.pool" 
nruns <- 5

# Running community phylogenetics analyses (picante) for both null models
print("------------------------------------------")
print("Spatial analysis: community phylogenetics")
print("------------------------------------------")

for(i in 1:length(null)){
  print(paste0("Null model: ", null[i]))
  
  # trait
  print("trait")
  result_mntd <- ses.mntd(samp = comm[,which(colnames(comm) %in% colnames(traitdist))], dis = traitdist, 
                          null.model = null[i], abundance.weighted=T, runs = nruns)
  result_mpd <- ses.mpd(samp = comm[,which(colnames(comm) %in% colnames(traitdist))], dis = traitdist, 
                        null.model = null[i], abundance.weighted=T, runs = nruns)

  result_mntd <- cbind(result_mntd, comm[,1:4], null = rep(null[i],nrow(result_mntd)), distance = rep("Shell size", nrow(result_mntd)))
  result_mpd <- cbind(result_mpd, comm[,1:4], null = rep(null[i],nrow(result_mpd)), distance = rep("Shell size", nrow(result_mpd)))
  
  ses_mntd <- rbind(ses_mntd,result_mntd)
  ses_mpd <- rbind(ses_mpd,result_mpd)
  
  rm(result_mntd)
  rm(result_mpd)
  
  
  # phylogeny
  print("phylogeny")
  result_mntd <- ses.mntd(samp = comm[,which(colnames(comm) %in% colnames(phydist))], dis = phydist, 
                          null.model = null[i], abundance.weighted=T, runs = nruns)

  result_mpd <- ses.mpd(samp = comm[,which(colnames(comm) %in% colnames(phydist))], dis = phydist, 
                        null.model = null[i], abundance.weighted=T, runs = nruns)
  
  result_mntd <- cbind(result_mntd, comm[,1:4], null = rep(null[i],nrow(result_mntd)), distance = rep("Phylogeny", nrow(result_mntd)))
  result_mpd <- cbind(result_mpd, comm[,1:4], null = rep(null[i],nrow(result_mpd)), distance = rep("Phylogeny", nrow(result_mpd)))
  
  ses_mntd <- rbind(ses_mntd,result_mntd)
  ses_mpd <- rbind(ses_mpd,result_mpd)
  
  rm(result_mntd)
  rm(result_mpd)
  
} # for


# Species not included in shell size analysis"
# colnames(comm)[which(colnames(comm) %!in% colnames(traitdist))][-c(1:4)]
# "Species not included in shell size analysis"
# colnames(comm)[which(colnames(comm) %!in% colnames(phydist))][-c(1:4)]

rownames(ses_mntd) <- c()
rownames(ses_mpd) <- c()

write.csv(ses_mntd, paste("output/spatial_ses_mntd_",nruns,".csv", sep = ""), row.names = F)
write.csv(ses_mpd, paste("output/spatial_ses_mpd_",nruns,".csv", sep = ""), row.names = F)

data_mntd <- ses_mntd
data_mpd <- ses_mpd

# Removing NAs
remove <- unique(c(which(is.na(data_mpd$mpd.obs.z)),which(is.na(data_mntd$mntd.obs.z))))
if(length(remove)>0){
  data_mntd <- data_mntd[-remove,]
  data_mpd <- data_mpd[-remove,]
}
all(row.names(data_mntd) == row.names(data_mpd)) # check

pmin <- 0.005
pmax <- 0.995

print("Spatial analysis: results ---------------------")
print(paste("Total number of communities analysed:", length(unique(data_mntd[,c("Latitude","Longitude")])[,1]))) # 3053

for(i in 1:length(null)){ 
  print("------------------------------------------")
  print(paste("NULL MODEL:", null[i]))
  
  print("----- Distance: Phylogeny -----")
  print("MNTD Phylogeny clustered NTI")
  print(length(which(data_mntd$mntd.obs.p < pmin & data_mntd$null == null[i] & data_mntd$distance == "Phylogeny" )))
  print("MNTD Phylogeny overdispersed NTI")
  print(length(which(data_mntd$mntd.obs.p > pmax & data_mntd$null == null[i] & data_mntd$distance == "Phylogeny" )))
  print("MPD Phylogeny clustered NRI")
  print(length(which(data_mpd$mpd.obs.p < pmin & data_mpd$null == null[i] & data_mpd$distance == "Phylogeny" )))
  print("MPD Phylogeny overdispersed NRI")
  print(length(which(data_mpd$mpd.obs.p > pmax & data_mpd$null == null[i] & data_mpd$distance == "Phylogeny" )))
  
  print("----- Distance: Shell size -----")
  print("MNTD Shell size clustered NTI")
  print(length(which(data_mntd$mntd.obs.p < pmin & data_mntd$null == null[i] & data_mntd$distance == "Shell size" )))
  print("MNTD Shell size overdispersed NTI")
  print(length(which(data_mntd$mntd.obs.p > pmax & data_mntd$null == null[i] & data_mntd$distance == "Shell size" )))
  print("MPD Shell size clustered NRI")
  print(length(which(data_mpd$mpd.obs.p < pmin & data_mpd$null == null[i] & data_mpd$distance == "Shell size" )))
  print("MPD Shell size overdispersed NRI")
  print(length(which(data_mpd$mpd.obs.p > pmax & data_mpd$null == null[i] & data_mpd$distance == "Shell size" )))
}  


#***************************************************************************************************
### Spatial plots ***
#********************
print("Generating spatial plots")

# Preparing to merge both datasets (for richness null model)
null_model <- 'richness'
data_mntd$metric <- 'NTI'
data_mpd$metric <- 'NRI'
data_mntd_null <-  data_mntd[which(data_mntd$null == null_model),]
data_mpd_null <-  data_mpd[which(data_mpd$null == null_model),]
colnames(data_mntd_null) <- sub("mntd.","",colnames(data_mntd_null))
colnames(data_mpd_null) <- sub("mpd.","",colnames(data_mpd_null))
all(colnames(data_mntd_null) == colnames(data_mpd_null))

# Merging:
data_plot <- rbind(data_mntd_null,data_mpd_null)

# double-check
unique(data_plot$null) == null_model
unique(data_plot$runs) == nruns
# Rounded columns to plot boxplot on top of the dots
data_plot$sst_round <- 1*round(data_plot$woa_tmn/1)


# SST - Note that we are plotting the NEGATION of SES (-obs.z), because this is the NRI and NTI metric
png(paste("output/fig_spatial_ses_", null_model, nruns,".png", sep =""), width = 8.5, height = 7, res = 400, units = "in")
ggplot(data = data_plot, aes(y = -obs.z, x = woa_tmn, group=sst_round)) + 
  geom_point(data = data_plot[which(data_plot$obs.p > pmin & data_plot$obs.p < pmax),],  size=1.5, color = "grey60", fill = alpha("grey80",.4),  shape = 21) +
  geom_point(data = data_plot[which(data_plot$obs.p < pmin),], color = "blue" , fill = alpha("blue",.4) , size=1.5, shape = 21, stroke =0.5) +
  geom_point(data = data_plot[which(data_plot$obs.p > pmax),], color = "red", fill = alpha("red",.4), size=1.5, shape = 21, stroke =0.8) +
  geom_boxplot(fill = alpha('white',0.7), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(x = expression("Annual mean sea-surface temperature ("*degree*"C)"), y = NULL) +
  theme_bw(base_size = 16) +
  scale_x_continuous(breaks = seq(0,30,5),limits=c(-2.5,31), expand = c(0, 0)) +
  #scale_y_continuous(breaks = seq(-4,4,2),limits=c(-5,5), expand = c(0, 0)) +
  facet_grid(metric~distance, scales = "fixed", switch = "y") +
  theme(axis.title = element_text(colour = 'black', size = 16),
        axis.text = element_text(colour = 'black', size = 16),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.5,0.2,0.2),"cm"),
        strip.placement.y = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black', size = 14))
dev.off()

# Calculating average shell size per commuity
row.names(sp_diam) <- sp_diam$sspname
size  <- sp_diam[,"mean", drop=F]
comm_size <- comm[,which(colnames(comm) %in% row.names(sp_diam))]
comm_size <- comm_size[,order(colnames(comm_size))]
size <- size[which(row.names(size) %in% colnames(comm_size)),,drop = F]
size <- size[order(row.names(size)),, drop=F]
# data.frame(colnames(comm_size),row.names(size)) # check 
d_size <- comm_size*size[col(comm_size)]
d_size[which(d_size == 0)] <- NA

d_size_mean <- cbind(comm[,c("Latitude","Longitude","woa_tmn","woa_tsd")], size_mean_0 = rowMeans(d_size))
d_size_mean <- cbind(d_size_mean, size_sum = rowSums(d_size, na.rm = T))

data_size <- merge(data_plot[which(data_plot$distance == "Shell size" ),], d_size_mean)
data_size <- merge(data_size, comm)
data_size$size_sum_round <- 25*round(data_size$size_sum/25)



png(paste("output/si_fig_ses_",null_model,nruns,"_comm_size.png", sep =""), width = 4, height = 5, unit = "in", res = 400)
ggplot(data_size, aes(x = size_sum, y = -obs.z, group = size_sum_round)) +
  geom_point(data = data_size[which(data_size$obs.p > pmin & data_size$obs.p < pmax),],  size=1.5, color = "grey60", fill = alpha("grey60",.0),  shape = 21) +
  geom_point(data = data_size[which(data_size$obs.p < pmin),], color = "blue" , fill = alpha("blue",.2) , size=1.5, shape = 21, stroke =0.5) +
  geom_point(data = data_size[which(data_size$obs.p > pmax),], color = "red", fill = alpha("red",.3), size=1.5, shape = 21, stroke =0.8) +
  geom_boxplot(fill = alpha('black',0.0), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  geom_hline(yintercept=0, linetype="dotted", lwd = 1) +
  theme_bw() +   
  labs(y = NULL, x = expression(paste("Abundance-weighted community size (",mu,m,")", sep=""))) +
  facet_grid(metric~., switch = "y") +
  theme(axis.text=element_text(size = 10, colour = "black"), 
        axis.title=element_text(size= 12, colour = "black"),
        strip.placement.y = "outside",
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 12)) +  
  scale_x_continuous(breaks = seq(0, 600, 100), limits=c(0,700), expand = c(0,0))
  # scale_y_continuous(breaks = seq(-4,2,2),limits=c(-5,4), expand = c(0, 0))
dev.off()

# for each species
metric = "NTI" # or NRI; community phylogenetic metric

if (!file.exists("output/spatial_species/")){ dir.create("output/spatial_species/")}

for (i in 19:59){
  ssp <- colnames(data_size)[i]
  # print(ssp)
  p <- ggplot() +
    geom_point(data = data_size[which(data_size$metric==metric),],  
               aes(x = data_size[which(data_size$metric==metric),i], y = -obs.z),
               size=1.5, color = "grey60", fill = alpha("grey60",.0),  shape = 21) +
    geom_hline(yintercept=0, linetype="dotted", lwd = 1) +
    theme_bw() +   
    theme(axis.text=element_text(size = 14, colour = "black"), 
          axis.title=element_text(size= 16, colour = "black")) +
    labs(x = ssp, y = metric)
    # scale_y_continuous(breaks = seq(-4,2,2),limits=c(-4.2,3), expand = c(0, 0))
  
  png(paste("output/spatial_species/si_fig_",metric,"_rel_abund_",ssp,".png", sep =""), width = 6, height = 4.5, unit = "in", res = 400)
     print(p)
  dev.off()
}


#***************************************************************************************************
### Temporal analysis & plots ***
#********************************
print("------------------------------------------")
print("------------------------------------------")
print("------------------------------------------")
print("Time series analysis")
print("------------------------------------------")

source("R/main_time.R")





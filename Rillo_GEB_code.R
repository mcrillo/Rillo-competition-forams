#########################
### ForCenS analysis  ###
#########################
# 11/Jun/2018

rm(list=ls())
setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./competition-forams")


if (!file.exists("output_spatial/")){ dir.create("output_spatial/")}

#  Libraries and functions
source("R/functions/functions_space.R")


##### DATA LOADING ################################################################################

### Spatial: coretop data from ForCenS (Siccha & Kucera 2017) https://doi.org/10.1038/sdata.2017.109
d <- read.csv(file ="data/coretop/ForCenS_woa.csv", header=TRUE, stringsAsFactors = FALSE) 

### Temporal: sediment trap data from 31 studies, see Table 



### Phylogeny data (Aze et al. 2011), just macroperforates
phy <- read.tree("data/phylogeny/aLb_Rillo_modified.tre")
phy <- dropExtinct(phy)
plotTree(phy)


### Shell size data
diam <- read.csv(file ="data/size/diameter_data_analysis.csv", header=TRUE, stringsAsFactors = FALSE) 
sp_diam <- ddply(diam,~species,summarise,mean=mean(max_diam_mu),sd=sd(max_diam_mu), n_ind = length(max_diam_mu))
sp_diam$sspname <- c("Dentigloborotalia_anfracta","Globigerina_bulloides","Globigerinella_calida","Turborotalita_clarkei",
                     "Globigerinoides_conglobatus","Globoquadrina_conglomerata","Globorotalia_crassaformis","Sphaeroidinella_dehiscens",
                     "Neogloboquadrina_dutertrei","Globigerina_falconensis","Tenuitella_fleisheri","Globigerinita_glutinata",
                     "Globorotaloides_hexagonus","Globorotalia_hirsuta","Turborotalita_humilis","Globoconella_inflata","Tenuitella_iota",
                     "Globorotalia_menardii","Globigerinita_minuta","Candeina_nitida","Pulleniatina_obliquiloculata","Neogloboquadrina_pachyderma",
                     "Tenuitella_parkerae","Berggrenia_pumilio","Turborotalita_quinqueloba","Orcardia_riedeli","Globigerinoides_ruber",
                     "Globoturborotalita_rubescens","Trilobatus_sacculifer","Globorotalia_scitula","Globigerinella_siphonifera",
                     "Globoturborotalita_tenella","Globorotalia_truncatulinoides","Globorotalia_tumida","Orbulina_universa","Globigerinita_uvula")

diam <- merge(diam, sp_diam[,c("species", "sspname")])

size_boxplot <- ggplot(diam, aes(x=sspname, y=max_diam_mu)) + 
  geom_boxplot() + 
  labs(y = expression(paste("Shell diameter (",mu,m,")", sep="")),
       x = "Species") +
  theme(axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12, face="italic", angle = 45, hjust = 1))

png(file = "si_fig_size_boxplot.png", width = 10, height = 7, unit = "in", res = 300)
 print(size_boxplot)
dev.off()


for (i in sp_diam$species){
  print("-----------------------")
  print(i)
  print(unique(diam[which(diam$species == i), "Reference"]))
  sp_diam[which(sp_diam$species == i), "Reference"] <-  paste(unique(diam[which(diam$species == i), "Reference"]), collapse = ";")
}

length(diam[which(diam$Reference == "Rillo et al. 2018"),1])
length(unique(diam[which(diam$Reference == "Rillo et al. 2018"),"species"]))

length(diam[which(diam$Reference == "Baranowski 2013 "),1])
length(unique(diam[which(diam$Reference == "Baranowski 2013 "),"species"]))

length(diam[which(diam$Reference == "Weinkauf et al. 2016"),1])
length(unique(diam[which(diam$Reference == "Weinkauf et al. 2016"),"species"]))

dim(diam)
length(unique(diam$species))

# Supp Info Table in LaTeX
sp_diam <- sp_diam[order(sp_diam$sspname),c("sspname","mean","n_ind", "Reference")]
latex_name <- paste("\textit{", gsub("_", " ", sp_diam$sspname), "}", sep="")
print(xtable(cbind(latex_name, round(sp_diam$mean,2), sp_diam$n_ind, sp_diam$Reference), include.rownames=F)) # LaTeX



##### DATA TRANSFORMATION ################################################################################

### Spatial (coretop data)

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
# Adding "Globigerinoides_white" as ruber to trait matrix (to match 'comm' matrix)
sp_diam <- rbind(sp_diam, c(sp_diam[which(sp_diam$species == "ruber"),c("species","mean","sd","n_ind")],sspname="Globigerinoides_white"))

# Subseting data to species and lat, long and SST (woa)
comm <- as.matrix(d[,c(5,6,68,69,22:62)]) 
# Tranforming NA into zero, otherwise function in analysis returns NA
comm[which(is.na(comm))] <- 0 


### Temporal: sediment trap


##### ANALYSIS ####################################################################################################

### Spatial data: community Phyogenetics (picante R package)

# Community dispersion metrics to be analysed
ses_mntd <- data.frame()
ses_mpd <- data.frame()

# Calculating ecological distance matrices
phydist <- cophenetic(phy)
traitdist <- as.matrix(dist(data.frame(sp_diam$mean, row.names = sp_diam$sspname))) 

# Setting variables
null <- c("taxa","richness") # "sample.pool" 
nruns <- 100

for(i in 1:length(null)){ # running for both null models
  
  print(null[i])
  
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

rownames(ses_mntd) <- c()
rownames(ses_mpd) <- c()

# Species not included in each shell size analysis
colnames(comm)[which(colnames(comm) %!in% colnames(traitdist))][-c(1:4)]
# Species not included in each phylogenetic analysis
colnames(comm)[which(colnames(comm) %!in% colnames(phydist))][-c(1:4)]


##### RESULTS ####################################################################################################

data_mntd <- ses_mntd
data_mpd <- ses_mpd

# REMOVING NAs
remove <- unique(c(which(is.na(data_mpd$mpd.obs.z)),which(is.na(data_mntd$mntd.obs.z))))
if(length(remove)>0){
  data_mntd <- data_mntd[-remove,]
  data_mpd <- data_mpd[-remove,]
}
all(row.names(data_mntd) == row.names(data_mpd)) # check


# Total number of communities analysed
pmin <- 0.005
pmax <- 0.995

print(paste("Total number of communities analysed:", length(unique(data_mntd[,c("Latitude","Longitude")])[,1]))) # 3053

for(i in 1:length(null)){ 
  print("***************************************")
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


##### PLOTS ####################################################################################################


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
png(paste("output_spatial/fig_ses_total_",null_model,nruns,"_sst.png", sep =""), width = 8.5, height = 7, res = 400, units = "in")
ggplot(data = data_plot, aes(y = -obs.z, x = woa_tmn, group=sst_round)) + 
  geom_point(data = data_plot[which(data_plot$obs.p > pmin & data_plot$obs.p < pmax),],  size=1.5, color = "grey60", fill = alpha("grey80",.4),  shape = 21) +
  geom_point(data = data_plot[which(data_plot$obs.p < pmin),], color = "blue" , fill = alpha("blue",.4) , size=1.5, shape = 21, stroke =0.5) +
  geom_point(data = data_plot[which(data_plot$obs.p > pmax),], color = "red", fill = alpha("red",.4), size=1.5, shape = 21, stroke =0.8) +
  geom_boxplot(fill = alpha('white',0.7), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(x = expression("Annual mean sea-surface temperature ("*degree*"C)"), y = NULL) +
  theme_bw(base_size = 16) +
  scale_x_continuous(breaks = seq(0,30,5),limits=c(-2.5,31), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-4,4,2),limits=c(-5,5), expand = c(0, 0)) +
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
trait <- sp_diam[,"mean", drop=F]
comm_size <- comm[,which(colnames(comm) %in% row.names(trait))]
comm_size <- comm_size[,order(colnames(comm_size))]
size <- trait[which(row.names(trait) %in% colnames(comm_size)),,drop = F]
size <- size[order(row.names(size)),, drop=F]
# data.frame(colnames(comm_size),row.names(size)) # check 
d_size <- comm_size*size[col(comm_size)]
d_size_mean <- cbind(comm[,c("Latitude","Longitude","woa_tmn","woa_tsd")], size_mean_0 = rowMeans(d_size))
d_size[which(d_size == 0)] <- NA
d_size_mean <- cbind(d_size_mean, size_mean_no0 = rowMeans(d_size, na.rm = T))
d_size_mean <- cbind(d_size_mean, size_sum = rowSums(d_size, na.rm = T))
d_size_mean <- as.data.frame(d_size_mean)
dim(d_size_mean)

data_size <- merge(data_plot[which(data_plot$distance == "Shell size" ),], d_size_mean)
data_size$size_mean_norm <- data_size$size_mean_0 / max(data_size$size_mean_0)
data_size$size_mean_no0_norm <- data_size$size_mean_no0 / max(data_size$size_mean_no0)
data_size$size_sum_norm <- data_size$size_sum / max(data_size$size_sum)

data_size$size_mean_norm_round <- 0.05*round(data_size$size_mean_norm/0.05)
data_size$size_mean_no0_norm_round <- 0.05*round(data_size$size_mean_no0_norm/0.05)
data_size$size_mean_norm_round <- 0.05*round(data_size$size_mean_norm/0.05)
data_size$size_sum_norm_round <- 0.05*round(data_size$size_sum_norm/0.05)
data_size$size_sum_round <- 25*round(data_size$size_sum/25)

data_size <- merge(data_size, comm)


png(paste("output_spatial/si_fig_ses_total_",null_model,nruns,"_comm_size.png", sep =""), width = 4, height = 5, unit = "in", res = 400)
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
  scale_x_continuous(breaks = seq(0, 600, 100), limits=c(0,700), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-4,2,2),limits=c(-5,4), expand = c(0, 0))
dev.off()

# for each species
if (!file.exists("output_spatial/species/")){ dir.create("output_spatial/species/")}
 metric = "NTI"
for (i in 27:ncol(data_size)){
  ssp <- colnames(data_size)[i]
  
  p <- ggplot() +
    geom_point(data = data_size[which(data_size$metric==metric & data_size$obs.p > pmin & data_size$obs.p < pmax),],  
               aes(x = data_size[which(data_size$metric==metric  & data_size$obs.p > pmin & data_size$obs.p < pmax),i], y = -obs.z),
               size=1.5, color = "grey60", fill = alpha("grey60",.0),  shape = 21) +
    geom_point(data = data_size[which(data_size$metric==metric & data_size$obs.p < pmin),], 
               aes(x = data_size[which(data_size$metric==metric & data_size$obs.p < pmin),i], y = -obs.z),
               color = "blue" , fill = alpha("blue",.2) , size=1.5, shape = 21, stroke =0.5) +
    geom_point(data = data_size[which(data_size$metric==metric & data_size$obs.p > pmax),], 
               aes(x = data_size[which(data_size$metric==metric & data_size$obs.p > pmax),i], y = -obs.z),
               color = "red", fill = alpha("red",.3), size=1.5, shape = 21, stroke =0.8) +
    geom_hline(yintercept=0, linetype="dotted", lwd = 1) +
    theme_bw() +   
    theme(axis.text=element_text(size = 14, colour = "black"), 
          axis.title=element_text(size= 16, colour = "black")) +
    labs(x = ssp, y = metric) +
    scale_y_continuous(breaks = seq(-4,2,2),limits=c(-4.2,3), expand = c(0, 0))
  
  png(paste("output_spatial/species/si_fig_",metric,"_",ssp,".png", sep =""), width = 6, height = 4.5, unit = "in", res = 400)
     print(p)
  dev.off()
  
}




#########################
### ForCenS analysis  ###
#########################
# 11/Jun/2018

rm(list=ls())

setwd("/Users/marinacostarillo/Dropbox/PhD/projects")
setwd("./competition-forams")
if (!file.exists("output_spatial/")){ dir.create("output_spatial/")}

#  Libraries and functions
source("R/functions/functions_space.R")


############################################################################################################################
### DATA ###################################################################################################################
############################################################################################################################

###
### Community relative abundance data (ForCenS, Siccha & Kucera 2017)
###

d <- read.csv(file ="data/coretop/ForCenS_woa.csv", header=TRUE, stringsAsFactors = FALSE) 

# Finding out which Ocean is which
# ocean = 2049
# world <- rnaturalearth::countries110
# plot(world)
# points(x=d[which(d$Ocean==ocean), "Longitude"], y=d[which(d$Ocean==ocean), "Latitude"], col = "blue")

# Remove dissolution-prone samples
# Tittensor et al. 2010 Nature: sites from >3,500m in the Pacific and Indian oceans and from >4,500m elsewhere were excluded
d <- d[-which(d$Water_depth > 4500),] # 7 North Atlantic; 11 South Atlantic ; 513 Artic ; 1025 Mediterranean ; # 257 Southern
d <- d[-which(d$Water_depth > 3500 & d$Ocean %in% c(49,81,129, 2049)),] 
# 49 North Pacific; 81 South Pacific; 129 Indian; 2049 Red Sea

# Remove communities with two species identification
remove <- c()
for( i in grep("_._", colnames(d))){
  if(any(d[,i]==0, na.rm = T)){ d[which(d[,i]==0),i] <- NA }
  remove <- c(remove,which(!is.na(d[,i])))
} # d[remove,grep("_._", colnames(d))]
d <- d[-remove,]
rm(remove)

colnames(d)[grep("tenell",colnames(d))] <- "Globoturborotalita_tenella" 
colnames(d)[grep("Globigerinoides",colnames(d))]
d[,which(colnames(d)=="Globigerinoides_ruber")] <- rowSums(d[,which(colnames(d)%in%c("Globigerinoides_ruber", "Globigerinoides_white","Globigerinoides_ruber_._Globigerinoides_white"))], na.rm = T)

# Function that assigns names to each grid
d <- function_name_grid(data = d)
# MEAN: average abundance data for each grid
d_mean <- data.frame(do.call("rbind", by(d[,c(5,6,22:62,68,69)], d$grid_name, function(x) colMeans(x, na.rm = T))))
# RANDOM SAMPLE: randomly selecting a sample within each named grid 
d_random <- do.call('rbind', lapply(split(d[,c(5,6,22:62,68,69)], d$grid_name), function(x) x[sample(1:nrow(x), 1),]))

# Subseting data to species and lat, long
comm <- as.matrix(d[,c(5,6,22:62,68,69)]) 
comm <- comm[,c(1,2,44,45,3:43)] # simply re-ordering columns
# Tranforming NA into zero,  otherwise function in analysis returns NA
comm[which(is.na(comm))] <- 0 
 
 

### Phylogeny data (based on Aze et al. 2011)

# Macroperforates
phy <- read.tree("data/phylogeny/aLb_Rillo_modified.tre")
phy <- dropExtinct(phy)
plotTree(phy)

# Microperforates
phy <- read.tree("data/phylogeny/aLb_Rillo_modified_plus_micro.tre")
plot(phy, show.tip.label = T); axisPhylo()

# Prune the phylogeny to include only the species that occurred in the ForCenS communities
prunedphy <- prune.sample(comm, phy)
plotTree(prunedphy)
length(prunedphy$tip.label)

# Calculate distance matrix
phydist <- cophenetic(prunedphy)

### Shell size data

diam <- read.csv(file ="data/size/diameter_data_analysis.csv", header=TRUE, stringsAsFactors = FALSE) 
diam$species <- as.factor(diam$species)
sp_diam <- ddply(diam, .(species_name2), summarize , mean=mean(max_diam_mu, na.rm =TRUE))
sp_diam <- merge(sp_diam, diam[,grep("species",colnames(diam))])

size_boxplot <- ggplot(diam, aes(x =  species, y=max_diam_mu)) + 
  geom_boxplot() + 
  labs(y = expression(paste("Shell diameter (",mu,m,")", sep="")),
       x = "Species") +
  theme(axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 13, 
                                   face="italic", angle = 45, hjust = 1))
pdf(file = "si_fig_size_boxplot.pdf", width=10, height=5, paper = "special")
 print(size_boxplot)
dev.off()

traitdist <- as.matrix(dist(data.frame(sp_diam$mean, row.names = sp_diam$species))) 


################################################################################################################################
### ANALYSIS ###################################################################################################################
################################################################################################################################

###
### Community Phyogenetics
### (picante package)


# Phylogenetic signal
# trait <- as.data.frame(sp_diam[,2], row.names = sp_diam$species_name2)
# oi <- data.frame(trait[phy$tip.label,])
# row.names(oi) <- phy$tip.label
# multiPhylosignal(oi, prunedphy)

ses_mntd <- data.frame()
ses_mpd <- data.frame()
null <- c("taxa","richness") # "sample.pool" 
nruns <- 100

for(i in 1:length(null)){
  
   print(null[i])
  
   # trait
   print("trait")
   result_mntd <- ses.mntd(samp = comm[,which(colnames(comm) %in% colnames(traitdist))], dis = traitdist, 
                      null.model = null[i], abundance.weighted=T, runs = nruns)
   print(dim(result_mntd))
   result_mpd <- ses.mpd(samp = comm[,which(colnames(comm) %in% colnames(traitdist))], dis = traitdist, 
                      null.model = null[i], abundance.weighted=T, runs = nruns)
   print(dim(result_mpd))
   
   if(nrow(result_mntd) < 400) { # d_mean or d_random (both have nrow == 307)
    result_mntd <- cbind(result_mntd[order(rownames(result_mntd)),], grid = row.names(result_mntd)[order(rownames(result_mntd))], 
                as.data.frame(comm[order(rownames(comm)),1:4]), null = rep(null[i],nrow(result_mntd)),
               distance = rep("Shell size", nrow(result_mntd)))
    result_mpd <- cbind(result_mpd[order(rownames(result_mpd)),], grid = row.names(result_mpd)[order(rownames(result_mpd))], 
                         as.data.frame(comm[order(rownames(comm)),1:4]), null = rep(null[i],nrow(result_mpd)),
                         distance = rep("Shell size", nrow(result_mpd)))
   }else{ # d_total
     result_mntd <- cbind(result_mntd, comm[,1:4], null = rep(null[i],nrow(result_mntd)), distance = rep("Shell size", nrow(result_mntd)))
     result_mpd <- cbind(result_mpd, comm[,1:4], null = rep(null[i],nrow(result_mpd)), distance = rep("Shell size", nrow(result_mpd)))
   }
   
   ses_mntd <- rbind(ses_mntd,result_mntd)
   ses_mpd <- rbind(ses_mpd,result_mpd)
  
   rm(result_mntd)
   rm(result_mpd)
   
   
   # phylogeny
   print("phylogeny")
   result_mntd <- ses.mntd(samp = comm[,which(colnames(comm) %in% colnames(phydist))], dis = phydist, 
                           null.model = null[i], abundance.weighted=T, runs = nruns, iterations = 1000)
   print(dim(result_mntd))
   result_mpd <- ses.mpd(samp = comm[,which(colnames(comm) %in% colnames(phydist))], dis = phydist, 
                         null.model = null[i], abundance.weighted=T, runs = nruns, iterations = 1000)
   print(dim(result_mpd))
   
   if(nrow(result_mntd) < 400) { # d_mean or d_random (both have nrow == 307)
     result_mntd <- cbind(result_mntd[order(rownames(result_mntd)),], grid = row.names(result_mntd)[order(rownames(result_mntd))], 
                          as.data.frame(comm[order(rownames(comm)),1:4]), null = rep(null[i],nrow(result_mntd)),
                          distance = rep("Phylogeny", nrow(result_mntd)))
     result_mpd <- cbind(result_mpd[order(rownames(result_mpd)),], grid = row.names(result_mpd)[order(rownames(result_mpd))], 
                         as.data.frame(comm[order(rownames(comm)),1:4]), null = rep(null[i],nrow(result_mpd)),
                         distance = rep("Phylogeny", nrow(result_mpd)))
     
   }else{ # d_total
     result_mntd <- cbind(result_mntd, comm[,1:4], null = rep(null[i],nrow(result_mntd)), distance = rep("Phylogeny", nrow(result_mntd)))
     result_mpd <- cbind(result_mpd, comm[,1:4], null = rep(null[i],nrow(result_mpd)), distance = rep("Phylogeny", nrow(result_mpd)))
   }
   
   ses_mntd <- rbind(ses_mntd,result_mntd)
   ses_mpd <- rbind(ses_mpd,result_mpd)
   
   rm(result_mntd)
   rm(result_mpd)

} # for

rownames(ses_mntd) <- c()
rownames(ses_mpd) <- c()

write.csv(ses_mntd, paste("output_spatial/ses_mntd_dtotal_",nruns,".csv", sep = ""), row.names = F)
write.csv(ses_mpd, paste("output_spatial/ses_mpd_dtotal_",nruns,".csv", sep = ""), row.names = F)

# Species not included in each analysis (phylogenetic and shell size)
colnames(comm)[which(colnames(comm) %!in% colnames(traitdist))][-c(1:4)]
colnames(comm)[which(colnames(comm) %!in% colnames(phydist))][-c(1:4)]


###
### Community Turnover
###


# Subseting data to species and lat, long
d_mean <- as.matrix(d_mean[,c(1,2,44,45,3:43)]) # simply re-ordering columns
# Tranforming NA into zero,  otherwise function in analysis returns NA
d_mean[which(is.na(d_mean))] <- 0 
dim(d_mean)

turnover <- data.frame()
pairs <- combn(nrow(d_mean),2)
length(pairs)

if(!file.exists("output_spatial/turnover_dmean.csv")){
 for(i in 1:ncol(pairs)){
   turnover[i,"sst_diff"] <- as.numeric(dist(rbind(d_mean[pairs[,i][1],"woa_tmn"],d_mean[pairs[,i][2],"woa_tmn"])))
   turnover[i,"bray"]     <- as.numeric(vegdist(x=d_mean[pairs[,i],5:ncol(d_mean)], method = "bray"))
   turnover[i,"jaccard"]  <- as.numeric(vegdist(x=d_mean[pairs[,i],5:ncol(d_mean)], method = "jaccard"))
   turnover[i,"horn"]     <- as.numeric(vegdist(x=d_mean[pairs[,i],5:ncol(d_mean)], method = "horn"))
   turnover[i,"chao"]     <- suppressWarnings(as.numeric(vegdist(x=d_mean[pairs[,i],5:ncol(d_mean)], method = "chao")))
 }
 turnover[,"sst_diff_round"] <- 1*(round(turnover[,"sst_diff"]/1))
 write.csv(turnover, "output_spatial/turnover_dmean.csv", row.names = F)
}

pdf("output_spatial/fig_turnover_dmean_chao.pdf", paper = "special", width = 7, height = 5)
ggplot(data = turnover, aes(y = chao)) + 
  geom_point(aes(x = sst_diff), size=1.5, col = alpha('seagreen',0.5)) +
  geom_boxplot(aes(x = sst_diff_round, group=sst_diff_round), fill = alpha('black',0.2), col = 'black', outlier.shape = NA) +
  labs(x = expression("Pair-wise difference in temperature ("*degree*"C)"), 
       y = "Pair-wise community dissimilarity (Chao)") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = seq(0,30,5),limits=c(-2,32), expand = c(0, 0)) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) 
dev.off()


###
### Community Evenness
###

richness <- specnumber(x = comm[,5:ncol(comm)])
div_shan <- diversity(x = comm[,5:ncol(comm)], index = "shannon")
div_simp <- diversity(x = comm[,5:ncol(comm)], index = "simpson")
div_invsimp <- diversity(x = comm[,5:ncol(comm)], index = "invsimpson")
div_simp_even <- div_simp/richness

comm_div <- cbind(data.frame(richness,div_shan, div_simp, div_invsimp,div_simp_even), comm)
write.csv(comm_div, "output_spatial/comm_diversity_dtotal.csv", row.names = F)

comm_div$sst_round <- 1*round(comm_div$woa_tmn/1)

pdf("output_spatial/fig_diversity_shannon_total.pdf", paper = "special", width = 7, height = 5)
ggplot(data = comm_div, aes(x = woa_tmn,y = div_shan, group = sst_round)) + 
  geom_point(size=1.5, col = alpha('seagreen',0.5)) +
  geom_boxplot(fill = alpha('black',0.0), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  labs(x = expression("Annual mean sea surface temperature ("*degree*"C)"), 
       y = "Shannon diversity") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = seq(0,30,5),limits=c(-2,32), expand = c(0, 0)) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) 
dev.off()

pdf("output_spatial/fig_diversity_shannon_tsd_total.pdf", paper = "special", width = 7, height = 5)
ggplot(data = comm_div, aes(x = woa_tsd,y = div_shan)) + 
  geom_point(size=1.5, col = alpha('seagreen',0.5)) +
  labs(x = expression("Sea surface temperature standard deviation ("*degree*"C)"), 
       y = "Shannon diversity") +
  theme_bw(base_size = 14) +
  # scale_x_continuous(breaks = seq(0,30,5),limits=c(-2,32), expand = c(0, 0)) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) 
dev.off()

################################################################################################################################
### RESULTS ####################################################################################################################
################################################################################################################################
# rm(list=ls())
null <- c("taxa","richness")
nruns <- 100

data_mntd <- read.csv(paste("output_spatial/ses_mntd_dtotal_",nruns,".csv", sep = ""), header=T)
data_mpd <- read.csv(paste("output_spatial/ses_mpd_dtotal_",nruns,".csv", sep = ""), header=T)
comm_div <- read.csv("output_spatial/comm_diversity_dtotal.csv", header=T)

data_mntd <- merge(data_mntd, comm_div)
data_mpd <- merge(data_mpd, comm_div)


# REMOVING NAs
# unique(data_mntd[which(is.na(data_mntd$mntd.obs.z)),"ntaxa"])
# unique(data_mpd[which(is.na(data_mpd$mpd.obs.z)),"ntaxa"])
remove <- unique(c(which(is.na(data_mpd$mpd.obs.z)),which(is.na(data_mntd$mntd.obs.z))))
if(length(remove)>0){
  data_mntd <- data_mntd[-remove,]
  data_mpd <- data_mpd[-remove,]
}
all(row.names(data_mntd) == row.names(data_mpd)) # check


# Preparing to merge both datasets
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

colMeans(data_mntd[which(data_mntd$woa_tmn<0),c(10,20:60)], na.rm = TRUE)


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


### MANUSCRIPT PLOT ****************************************************************************************************************************

# Rounded columns to plot boxplot on top of the dots
data_plot$lat_round <- 3*round(data_plot$Latitude/3)
data_plot$sst_round <- 1*round(data_plot$woa_tmn/1)
data_plot$sstsd_round <- 0.1*round(data_plot$woa_tsd/0.1)
data_plot$shan_round <- 0.1*round(data_plot$div_shan/0.1)

# If SES_metric is negative: clustered (i.e., smaller distance between species pairs)

# SST - Note that we are plotting the NEGATION of SES (-obs.z), because this is the NRI and NTI metric
png(paste("output_spatial/fig_ses_total_",null_model,nruns,"_sst_v2019_p05.png", sep =""), width = 8.5, height = 7, res = 400, units = "in")
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


# Shannon diversity - Note that we are plotting the NEGATION of SES (-obs.z)
png(paste("output_spatial/fig_ses_total_",null_model,nruns,"_shannon.png", sep =""), width = 8.5, height = 6, res = 400,  units = "in")
 ggplot(data = data_plot, aes(y = -obs.z, x = div_shan, group = shan_round)) + 
   geom_point(data = data_plot[which(data_plot$obs.p > pmin & data_plot$obs.p < pmax),],  size=1.5, color = "grey60", fill = alpha("grey60",.0),  shape = 21) +
   geom_point(data = data_plot[which(data_plot$obs.p < pmin),], color = "blue" , fill = alpha("blue",.2) , size=1.5, shape = 21, stroke =0.5) +
   geom_point(data = data_plot[which(data_plot$obs.p > pmax),], color = "red", fill = alpha("red",.3), size=1.5, shape = 21, stroke =0.8) +
   geom_boxplot(fill = alpha('black',0.0), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
   geom_hline(yintercept=0, linetype="dotted") +
   labs(x = "Shannon diversity", y = NULL) +
  facet_grid(metric~distance, scales = "fixed", switch = "y") +
  theme_bw() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 14),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.5,0.2,0.2),"cm"),
        strip.placement.y = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black', size = 14)) +
   scale_y_continuous(breaks = seq(-4,3,2),limits=c(-5,4), expand = c(0, 0))
dev.off()

# Simpson diversity - Note that we are plotting the NEGATION of SES (-obs.z)
png(paste("output_spatial/fig_ses_total_",null_model,nruns,"_simpson.png", sep =""), width = 8.5, height = 6, res = 400,  units = "in")
ggplot(data = data_plot, aes(y = -obs.z, x = div_simp)) + 
  geom_hline(yintercept = 0, linetype="dashed", lwd = 1) +
  geom_point(data = data_plot[which(data_plot$obs.p > pmin & data_plot$obs.p < pmax),],  size=1.5, color = "grey60", fill = alpha("grey60",.0),  shape = 21) +
  geom_point(data = data_plot[which(data_plot$obs.p < pmin),], color = "blue" , fill = alpha("blue",.2) , size=1.5, shape = 21, stroke =0.5) +
  geom_point(data = data_plot[which(data_plot$obs.p > pmax),], color = "red", fill = alpha("red",.3), size=1.5, shape = 21, stroke =0.8) +
  labs(x = "Simpson diversity", y = NULL) +
  facet_grid(metric~distance, scales = "fixed", switch = "y") +
  theme_bw() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 14),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.5,0.2,0.2),"cm"),
        strip.placement.y = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black', size = 14)) +
  scale_y_continuous(breaks = seq(-4,3,2),limits=c(-5,4), expand = c(0, 0))
dev.off()


# Richness - Note that we are plotting the NEGATION of SES (-obs.z)
png(paste("output_spatial/fig_ses_total_",null_model,nruns,"_richness.png", sep =""), width = 8.5, height = 6, res = 400,  units = "in")
ggplot(data = data_plot, aes(y = -obs.z, x = ntaxa, group = ntaxa)) + 
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_point(data = data_plot[which(data_plot$obs.p > pmin & data_plot$obs.p < pmax),],  size=1.5, color = "grey60", fill = alpha("grey60",.0),  shape = 21) +
  geom_point(data = data_plot[which(data_plot$obs.p < pmin),], color = "blue" , fill = alpha("blue",.2) , size=1.5, shape = 21, stroke =0.5) +
  geom_point(data = data_plot[which(data_plot$obs.p > pmax),], color = "red", fill = alpha("red",.3), size=1.5, shape = 21, stroke =0.8) +
  geom_boxplot(fill = alpha('black',0.0), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  labs(x = "Species richness", y = NULL) +
  facet_grid(metric~distance, scales = "fixed", switch = "y") +
  scale_x_continuous(breaks = seq(2,26,2),limits=c(1,27), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 14),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.5,0.2,0.2),"cm"),
        strip.placement.y = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black', size = 14)) +
  scale_y_continuous(breaks = seq(-4,3,2),limits=c(-5,4), expand = c(0, 0))
dev.off()




### ANALYSING COMMUNITY SIZE ****************************************************************************************************************************
# Overdispersed
data_plot[which(data_plot$obs.p > pmax & data_plot$distance == "Shell size" ),]

# Clustered
data_plot[which(data_plot$obs.p< pmin & data_plot$distance == "Shell size" ),]


# Calculating average shell size per commuity
trait <- unique(sp_diam[,c("species_name2","mean")])
row.names(trait) <- trait$species_name2
comm_size <- comm[,which(colnames(comm) %in% row.names(trait))]
comm_size <- comm_size[,order(colnames(comm_size))]
size <- trait[which(row.names(trait) %in% colnames(comm_size)),"mean",drop = F]
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
  scale_x_continuous(breaks = seq(0, 600, 100), limits=c(0,650), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-4,3,2),limits=c(-5,4), expand = c(0, 0))
dev.off()


plot(data_size$obs.z ~ data_size$size_mean_0)
plot(data_size$obs.z ~ data_size$div_shan)
plot(data_size$size_mean_0 ~ data_size$div_shan)
plot(data_size$size_mean_0 ~ data_size$woa_tmn)


# for each species, save plot
if (!file.exists("output_spatial/species/")){ dir.create("output_spatial/species/")}

metric = "NTI"
for (i in 20:60){
  ssp <- colnames(data_size)[i]
  print(ssp)
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
         scale_y_continuous(breaks = seq(-4,4,1),limits=c(-4.2,4.2), expand = c(0, 0))
  
  png(paste("output_spatial/species/si_fig_",metric,"_",ssp,".png", sep =""), width = 6, height = 4.5, unit = "in", res = 400)
   print(p)
  dev.off()
}



# fig_shannon_comm_size.pdf
ggplot(data_size, aes(x = size_mean_norm, y = div_shan)) +
  geom_point(aes(color = woa_tmn)) +theme_bw() +   
  labs(y = "Shannon diversity", x = "Community average shell size") +
  theme(axis.text=element_text(size = 14, colour = "black"), 
        axis.title=element_text(size= 16, colour = "black")) +  
  scale_y_continuous(breaks = seq(0, 2.5, 0.5), limits=c(-0.1, 3), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits=c(0, 1.05), expand = c(0,0)) +
  scale_color_viridis(name= c("SST"), option="plasma")



##############################



# Size as a function of relative abundance 
ggplot(traits_all, aes(x = abund, y = mean_area_sq )) +
  geom_point(pch = 21, color ="black", size = 4, stroke = 1) + 
  theme_bw() +
  geom_text_repel(data = trait[which(trait$abund > 0.019),], aes(label = species), fontface = 'italic',
                  segment.color = 'black', size =5, nudge_x= 0.05, direction = "y")  +    
  labs(y = expression(paste("Species' average shell size  (",mu,m,")", sep="")), 
       x = expression("Mean relative abundance (%)")) +
  theme(axis.text=element_text(size = 14, colour = "black"), 
        axis.title=element_text(size= 16, colour = "black"),
        plot.margin=unit(c(0.5,0.8,0.2,0.2),"cm")) +
  scale_y_continuous(breaks = seq(100, 800, 100), limits=c(100, 800), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0, 0.20, 0.025), limits=c(-0.01, 0.20), expand = c(0,0))




# Dispersion metrics as function of relative abundance of specific species

ggplot(subset(data_mntd, null == "richness" & distance == "Shell size"), 
       aes(x=Globorotalia_menardii, y=mntd.obs.z)) +
  geom_hline(yintercept=0, linetype="dashed", lwd = 1) +
  geom_point(size=2, pch = 21, stroke = 0.7, fill = alpha('grey',0.4)) +
  geom_point(data = subset(data_mntd, null == "richness" & distance ==  "Shell size" & mntd.obs.p < pmin), 
             aes(x=Globorotalia_menardii, y=mntd.obs.z), size=3, pch = 21, stroke = 0.7, fill = alpha("blue", 0.6)) +
  geom_point(data = subset(data_mntd, null == "richness" & distance ==  "Shell size" & mntd.obs.p > pmax), 
             aes(x=Globorotalia_menardii, y=mntd.obs.z), size=3, pch = 21, stroke = 0.7, fill = alpha("red", 0.6)) +
  #labs(y = "MNTD value (shell size)", x = "Shannon diversity") +
  theme_bw() +
  theme(axis.text=element_text(size=14, colour = "black"), 
        axis.title=element_text(size=16, colour = "black"))



comm_mntd <- merge(data_mntd[which(data_mntd$null == "richness" & data_mntd$distance == "Shell size" ),], comm,
                   by = c("Latitude","Longitude","woa_tmn","woa_tsd"))
comm_mntd[which(comm_mntd$mntd.obs.p>pmax),]

size_mntd <- merge(d_size_mean, data_mntd[which(data_mntd$null == "richness" & data_mntd$distance == "Shell size" ),], 
                   by = c("Latitude","Longitude","woa_tmn","woa_tsd"))

size_mntd$sst_round <- 1*round(size_mntd$woa_tmn/1)

ggplot(data = size_mntd, aes(y = mean_no0, x = woa_tmn, group = sst_round)) + 
  geom_point(size=1.5, col = "grey60") +
  geom_boxplot(fill = alpha('black',0.0), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  geom_point(data = size_mntd[which(size_mntd$mntd.obs.z<(-2)),], size=1.5, pch = 21, col = "blue", fill = alpha("blue", 0.2)) +
  geom_point(data = size_mntd[which(size_mntd$mntd.obs.z>3),], size=1.5, pch = 21, col = "red", fill = alpha("red", 0.2)) +
  labs(x = expression("Annual mean sea surface temperature ("*degree*"C)"), 
       y = "Assemblage average shell size") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) 

ggplot(data = size_mntd, aes(y = ntaxa, x = woa_tmn, group = sst_round)) + 
  geom_point(size=1.5, col = "grey60") +
  geom_boxplot(fill = alpha('black',0.0), col = alpha('black',0.8), outlier.shape = NA, lwd=0.3) +
  geom_point(data = size_mntd[which(size_mntd$mntd.obs.p<pmin),], size=1.5, col = "blue") +
  geom_point(data = size_mntd[which(size_mntd$mntd.obs.p>pmax),], size=1.5, col = "red") +
  labs(x = expression("Annual mean sea surface temperature ("*degree*"C)"), 
       y = "Sample richness") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) 



### ANALYSING COMMUNITY PHYLOGENETICS ************************************************************************************************************
# MNTD Phylogeny overdispersed - both null models!
comm_phy <- data_mntd[which(data_mntd$mntd.obs.p > pmax & data_mntd$null == "richness" & data_mntd$distance == "Phylogeny" ),]
comm_phy <- rbind(comm_phy, data_mntd[which(data_mntd$mntd.obs.p > pmax & data_mntd$null == "taxa" & data_mntd$distance == "Phylogeny" ),])
# All the 9 assemblages in the richness model are included in the taxa null model:
comm_phy[order(comm_phy[,c("Latitude", "Longitude")]),]
comm_phy <-comm_phy[-duplicated(comm_phy[,c("Latitude", "Longitude")]),]
comm_phy <- merge(comm_phy, comm, by=c("Latitude", "Longitude","woa_tmn"))
comm_phy <- comm_phy[,!apply(comm_phy == 0 , 2, all)]
comm_phy # Discussion: species-poor assemblages with bulloides and incompta/pachyderma.
mean(comm_phy$ntaxa)
median(comm_phy$ntaxa)

# MPD Phylogeny clustered 
comm_phy <- data_mpd[which(data_mpd$mpd.obs.z < (-2.75) & data_mpd$null == "richness" & data_mpd$distance == "Phylogeny" ),]
comm_phy <- rbind(comm_phy, data_mpd[which(data_mpd$mpd.obs.z < (-2.75) & data_mpd$null == "taxa" & data_mpd$distance == "Phylogeny" ),])
# All the 9 assemblages in the richness model are included in the taxa null model:
comm_phy[order(comm_phy[,c("Latitude", "Longitude")]),]
comm_phy <-comm_phy[-duplicated(comm_phy[,c("Latitude", "Longitude")]),]
comm_phy <- merge(comm_phy, comm, by=c("Latitude", "Longitude","woa_tmn"))
comm_phy <- comm_phy[,!apply(comm_phy == 0 , 2, all)]
comm_phy # Discussion: species-poor assemblages with bulloides and incompta/pachyderma.
mean(comm_phy$ntaxa)
median(comm_phy$ntaxa)


length(which(data_mpd$mpd.obs.p < pmin & data_mpd$null == "richness" & data_mpd$distance == "Phylogeny" ))



# NHM 17 Oct - plots
comm_phy <- data_mntd[which(data_mntd$null == "richness" & data_mntd$distance == "Phylogeny" ),]

ggplot(comm_phy, aes(x = ntaxa, y = mntd.obs.z )) +
  geom_point(pch = 21, color ="grey60", size = 3, stroke = 1) +
  geom_point(data=comm_phy[which(comm_phy$mntd.obs.p < 0.05),], pch = 21, color ="blue", size = 3, stroke = 1) + 
  geom_point(data=comm_phy[which(comm_phy$mntd.obs.p > 0.95),], pch = 21, color ="red", size = 3, stroke = 1) + 
  theme_bw() 




################################################################################################################################
# check: https://rstudio-pubs-static.s3.amazonaws.com/7433_4537ea5073dc4162950abb715f513469.html

###
### World map
###
coords_clust <- data[which(data$sign == -1 & data$null == "richness"),c("Latitude", "Longitude")]
mappoints <- mapplot + geom_point(data = coords_clust, aes(x = Longitude, y = Latitude, group = row.names(coords_clust)),  
                                  shape = 21, color = "navyblue", fill = alpha("blue",0.2), size = 3, stroke = 1) # add points
# geom_label_repel(aes(x = lonE, y = latN, group = trap, label = trap), data = traps_meta ) 

png(file = "output_spatial/map_mntd_total.png", width = 12, height = 7, units = "in", res = 300)
 print(mappoints)
dev.off()


###
### Boxplots of SES
###

# Significance as factor
data[,"sign"] <- 0
data[which(data$mntd.obs.p > pmax),"sign"] <- 1
data[which(data$mntd.obs.p < pmin),"sign"] <- -1
data[,"sign"] <- as.factor(data[,"sign"])

### d_mean and d_random

p <- ggplot(data = data, aes(x = null, y = round(mntd.obs.z/0.1)*0.1)) + facet_grid(. ~ distance) +
  geom_beeswarm(aes(color = sign), cex=1.7, size = 2.5) + 
  geom_boxplot(alpha = 0, width=0.2, size = 1) +
  scale_color_manual(values=c("blue", "grey70", "red")) +
  theme_bw() + ylim(-2,2)

p + labs(x = "Null Model", y = "SES Mean Nearest Distance") +
  theme(strip.text = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size=16),
        axis.text = element_text(size = 14, colour = "black"),
        legend.position = "none")

# Poster
png(file = "output_spatial/boxplot_dmean_001sign.png", width=14, height=8, units = "in", res = 500)
print(
  p +  labs(x = "Null Model", y = "SES Mean Nearest Distance") +
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 30, margin = margin(0.5,0,0.5,0, "cm")),
          axis.title.x = element_text(size = 30, margin = margin(0.5,0,0,0, "cm")),
          axis.title.y = element_text(size = 30, margin = margin(0,0.5,0,0, "cm")),
          axis.text = element_text(size = 24, colour = "black"),
          legend.position = "none")
)
dev.off()

### d_total

p <- ggplot(data = data, aes(x = null, y = mntd.obs.z)) + facet_grid(. ~ distance) +
  geom_beeswarm(aes(color = sign), cex=0.5, size = 1.5) + 
  geom_boxplot(alpha = 0, width=0.2, size = 1) +
  scale_color_manual(values=c("blue", "grey70", "red")) +
  theme_bw() + ylim(-2,3)

# Poster
png(file = "output_spatial/boxplot_dtotal_001sign2.png", width=16, height=8, units = "in", res = 500)
print(
  p +  labs(x = "Null Model", y = "SES Mean Nearest Distance") +
    scale_x_discrete(labels = c('Richness','Tree Tips')) +
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 30, margin = margin(0.5,0,0.5,0, "cm")),
          axis.title.x = element_text(size = 30, margin = margin(0.5,0,0,0, "cm")),
          axis.title.y = element_text(size = 30, margin = margin(0,0.5,0,0, "cm")),
          axis.text = element_text(size = 24, colour = "black"),
          legend.position = "none")
)
dev.off()




################################################################################################################################
### Readings
# Tutorial from: http://picante.r-forge.r-project.org/picante-intro.pdf
# https://daijiang.name/en/2014/05/29/null-model-for-phylogenetic-functional-diversity/
# http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

# names(ses_mntd_result)
  # ntaxa Number of taxa in community
  # mpd.obs Observed mpd in community
  # mpd.rand.mean Mean mpd in null communities
  # mpd.rand.sd Standard deviation of mpd in null communities
  # mpd.obs.rank Rank of observed mpd vs. null communities
  # mpd.obs.z Standardized effect size of mpd vs. null communities 
  # mpd.obs.p P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank/ runs + 1)
  # runs Number of randomizations
  
# MNTD (mean nearest taxon/species distance)
# SES_metric = mntd.obs.z
# SES_metric = (MNTD_obs - mean(MNTD_null))/sd(MNTD_null)
# SES means standardized effect size of MNTD vs. 999 null communities
# If SES_metric is negative, this means that the distance between species of the mean(MNTD_null) is larger than the observed distance between species in the community
# So, if SES_metric is negative, the observed community is phylogenetically clustered (i.e., smaller distance between species pairs)

# MPD (mean pairwise distance)
# Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) indicate phylogenetic evenness, 
# or a greater phylogenetic distance among co-occurring species than expected (overdispersion). 
# Negative SES values and low quantiles (mpd.obs.p < 0.05) indicate phylogenetic clustering, 
# or small phylogenetic distances among co-occurring species than expected.

# Map plot with grid
# map("world", interior = FALSE) 
# map.grid(c(-180, 180, -90, 90), nx=36, ny=18, labels=FALSE, pretty=FALSE) +
#  points(d$Longitude, d$Latitude, xlim = c(0, 360), ylim = c(0, 180), col= "red", cex= 0.5)


### DAMOCLES package




### PEZ package




################################################################################################################################################
# OLD - 2017

# abund_raw <- read.csv(file ="data/coretop/ForCenS_woa.csv", header=TRUE, stringsAsFactors = FALSE) 


# Plotting relative abundance histogram for each species
# plot_rel_abund_histo(abund_raw, bin = 0.01) # = 1% relative abundance
# warnings due to NAs in the rel abund data


### Max relative abundance
# max <- get_max_relabund(abund_raw)
# max <- max[order(max[,"species"]),]


### Co-occurence matrix (considering time!)

# Reshaping data: rows are species, columns are sites
# abund_mat <- get_abund_matrix(abund_raw) # ASSUMPTION: 'NA' entry -> 0 , see plot output/coretop_ranges/range_overlap_NA.pdf
# abund_mat[1:10,1:10]


# Creating binary presence-absence matrix
# occur_mat <- abund_mat
# occur_mat[which(occur_mat>0)] <- 1
# occur_mat[1:10,1:10]


# Range overlap
# ranges_overlap <- get_ranges_overlap(occur_mat)
# plot_ranges_overlap(ranges_overlap, method="color", order_clust = c("alphabet")) # "AOE", "FPC", "hclust", "alphabet


# Allopatric species pairs
# presence <- as.matrix(ranges_overlap)
# presence[which(presence>0)] <- 1 # presence in the sample (=1)
# ssp_allop <- get_allopatry(presence)


# Plot range overlap clustering dendogram
# STILL WORKING
# diss <- daisy(ranges_overlap, metric = "manhattan")
# clust.res <- hclust(diss)
# plot(clust.res)


# C-score (checkboard)
# c_score(occur_mat)


# Defining rare X abundant species
# rare_abund <- get_rare_abund(ranges_overlap, occur_mat, species)
# plot_rare_abund(rare_abund)


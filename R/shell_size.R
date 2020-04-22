

library(plyr)

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./competition-forams")
diam <- read.csv(file ="data/size/diameter_data_analysis.csv", header=TRUE, stringsAsFactors = FALSE) 

sp_diam <- ddply(diam,~species,summarise,mean=mean(max_diam_mu),sd=sd(max_diam_mu), n_ind = length(max_diam_mu))

sp_diam$sspname <- c("Dentigloborotalia anfracta","Globigerina bulloides","Globigerinella calida","Turborotalita clarkei",
                     "Globigerinoides conglobatus","Globoquadrina conglomerata","Globorotalia crassaformis","Sphaeroidinella dehiscens",
                     "Neogloboquadrina dutertrei","Globigerina falconensis","Tenuitella fleisheri","Globigerinita glutinata",
                     "Globorotaloides hexagonus","Globorotalia hirsuta","Turborotalita humilis","Globoconella inflata","Tenuitella iota",
                     "Globorotalia menardii","Globigerinita minuta","Candeina nitida","Pulleniatina obliquiloculata","Neogloboquadrina pachyderma",
                     "Tenuitella parkerae","Berggrenia pumilio","Turborotalita quinqueloba","Orcardia riedeli","Globigerinoides ruber",
                     "Globoturborotalita rubescens","Trilobatus sacculifer","Globorotalia scitula","Globigerinella siphonifera",
                     "Globoturborotalita tenella","Globorotalia truncatulinoides","Globorotalia tumida","Orbulina universa","Globigerinita uvula")

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



### OLD 


traits_all <- read.csv(file ="data/size_species_all.csv", header=TRUE, stringsAsFactors = FALSE) 
row.names(traits_all)<- traits_all$species

trait <- data.frame(mean_area_sq = sqrt(rowMeans(traits_all[,c("mean_area","book_scan_average_area_m2")], na.rm = T)))
if(any(is.na(trait[,1]))){ trait <- trait[-which(is.na(trait[,1])),, drop=FALSE]}
if(length(grep("clarkei", colnames(comm))) == 0){ trait <- trait[-grep("clarkei", row.names(trait)),, drop=F]}
# Calculate distance matrix
traitdist <- as.matrix(dist(trait)) 

size_measur <- read.csv(file ="data/bias_size_ind_resamples.csv", header=TRUE, stringsAsFactors = FALSE) 
size_measur$area_sq <-  sqrt(size_measur$area)


# Species with more than 20 individuals measured
unique_measur <- unique(size_measur[,c("total_ind", "species", "sample_no")])
ssp_hist <- names(which(by(unique_measur$total_ind, unique_measur$species, function(x) sum(x))>20))

# Unimodality - not...
# by(size_measur[which(size_measur$species %in% ssp_hist),"area_sq"], 
#   size_measur[which(size_measur$species %in% ssp_hist),"species"], 
#   function(x) summary(Mclust(x), parameters = F))

if (!file.exists("si_fig_size_histo.pdf")){
  hist_ssp <- ggplot(size_measur[which(size_measur$species %in% ssp_hist),], aes(x=area_sq)) + 
    geom_histogram(binwidth = 50) + 
    facet_grid(sspname ~ ., scales = "free") + # horiz.: facet_grid(. ~ species)
    labs(x = expression(paste("Square-root shell area (",mu,m,")", sep="")),
         y = "Number of individuals measured") +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18),
          strip.text = element_text(size = 13, face="italic"))
  
  pdf(file = "si_fig_size_histo.pdf", width=5, height=20, paper = "special")
  print(hist_ssp)
  dev.off()
}

if (!file.exists("si_fig_size_boxplot.pdf")){
  hist_ssp <- ggplot(size_measur[which(size_measur$species %in% ssp_hist),], aes(x =  sspname, y=area_sq)) + 
    geom_boxplot() + 
    labs(y = expression(paste("Shell area (",mu,m,")", sep="")),
         x = "Species") +
    theme(axis.text.y=element_text(size=16),
          axis.title.y=element_text(size=16),
          axis.title.x=element_blank(),
          axis.text.x = element_text(size = 13, 
                                     face="italic", angle = 45, hjust = 1))
  
  pdf(file = "si_fig_size_boxplot.pdf", width=7, height=4, paper = "special")
  print(hist_ssp)
  dev.off()
}

rm(list=ls())

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./competition-forams")


#  Libraries and functions
source("R/functions/functions_time.R")
source("R/functions/functions_space.R")

if (!file.exists("relative_absolute/")){ dir.create("relative_absolute/")}

###
### Data
###

# ForCenS
d <- read.csv(file ="data/coretop/ForCenS_woa.csv", header=TRUE, stringsAsFactors = FALSE) 

# Time-series
ssp_abbr <- c("B_dig","B_pum","C_nit","D_anf","G_bul","G_fal","G_ada","G_cal","G_sip","G_glu","G_min","G_uvu","G_con","G_rubp","G_rubw","G_inf","G_coa","G_cav","G_cra","G_hir","G_men","G_sci","G_the","G_tru","G_tum","G_ung","G_hex","G_rus","G_ten","H_pel","H_dig","N_dut","N_inc","N_pac","O_uni","P_obl", "S_deh","T_fle","T_iot","T_par","G_sac","T_cla","T_hum","T_qui")
traps <- sub(x=list.files(path = "data/traps_time/")[grep(list.files(path = "data/traps_time/"),pattern=".csv")], pattern=".csv", replacement = "")
traps_meta <- read.csv(file="data/traps_metadata.csv", header=TRUE, stringsAsFactors = FALSE)
traps_meta <- traps_meta[which(traps_meta$trap %in% traps),which(colnames(traps_meta) %in% c("trap","sample_no","latN","lonE","woa13_SSTmn", "mean_open_period_days", 
                                                                                             "length_days","length_years","length_time.series", "from", "to"))] 





### ABSOLUTE ABUNDANCE

abund_mat <- data.frame(matrix(rep(NaN,length(ssp_abbr)*length(traps)), 
                           nrow = length(ssp_abbr),
                           ncol = length(traps),
                           dimnames = list(ssp_abbr,traps)))

for(i in traps){ # i = traps[2]
  data <- read.csv(paste("data/traps_raw/",i,".csv",sep=""), header=TRUE)
  data <- data[, -c(1,2)]
  
  fluxsum <- data.frame(colSums(data))
  abund_mat[match(rownames(fluxsum), rownames(abund_mat)), i] <- fluxsum
}

# Species present in at least 10 traps
ntraps <- data.frame(ntraps = rowSums(!is.na(abund_mat)))
sort(ntraps, decreasing = T)
species <- sort(rownames(ntraps[which(ntraps$ntraps>10),,drop = F]))
species <- cbind(species, c("Globigerina_bulloides","Globigerinella_calida","Globigerinita_glutinata","Globoconella_inflata",
                                 "Globigerinoides_white","Trilobatus_sacculifer","Globorotalia_scitula","Globigerinella_siphonifera",
                                 "Globorotalia_truncatulinoides","Neogloboquadrina_dutertrei","Neogloboquadrina_incompta",
                                 "Neogloboquadrina_pachyderma","Orbulina_universa","Turborotalita_quinqueloba"))
colnames(species) <- c("trap", "forcens")
species <- as.data.frame(species)

trap_abund <- t(abund_mat[which(rownames(abund_mat) %in% species$trap),])
trap_abund <- round(trap_abund/traps_meta$length_years)

### RELATIVE ABUNDANCE

coords <- traps_meta[,c("lonE","latN")]
rownames(coords) <- traps_meta$trap
coords[,c("row_findin","distance")] <- NA
coords[,c("row_findin","distance")] <- ldply(1:nrow(coords), .fun = function(x) find_neighbours(coords[x,], findin = d[,c("Longitude", "Latitude")], distance = 0))
coords <- cbind(coords, d[coords$row_findin,which(colnames(d) %in% species$forcens)])
write.csv(coords, "relative_absolute/traps_forcens_neighbours.csv")


# Merging

abunds <- merge(trap_abund, coords, by = "row.names")
colnames(abunds)[1] <- "trap"
abunds <- merge(abunds, traps_meta[,c("trap","length_years")], by = "trap")

species <- as.data.frame(lapply(species, as.character), stringsAsFactors = F)

for (i in 1:nrow(species)){

   abundssp <- abunds[, which(colnames(abunds) %in% paste(c(species[i,],"length_years")))]
   names(abundssp) <- c("absolute", "relative", "Years")
   abundssp$absolute <- log(abundssp$absolute)
     
   fit1 <- lm(relative ~ absolute, data = abundssp)

   l = ifelse(coef(summary(fit1))[2,4] < 0.05, 1 , 2)
   c = ifelse(coef(summary(fit1))[2,4] < 0.05, "red" , "grey60")
   
   plot <- ggplot(abundssp, aes(y=relative, x=absolute, fill = Years)) +
     ggtitle(sub("_"," ",species[i,2])) + 
     geom_smooth(method=lm, se=FALSE, linetype = l, color = c, size = 0.8) + theme_bw() +
     geom_point(size = 2.5, pch = 21, color = "black") + 
     labs(y = "Relative abundance", x = "Absolute abundance (log)") +
     theme(axis.text=element_text(size=12, colour = "black"), 
           axis.title=element_text(size=14, colour = "black"),
           legend.position = "none",
           plot.margin = margin(0.5, 0.5,0.5,0.5, "cm"),
           plot.title = element_text(face = "italic", size = 14, vjust = 0.5, hjust = 0.5)) +
     scale_fill_viridis(alpha = 0.3, option = "inferno")
  
   png(file = paste("relative_absolute/", species[i,2],"_log.png", sep = ""), width = 5.75, height = 4, units = "in", res = 300)
     print(plot)
   dev.off()
}



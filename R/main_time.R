
#  Libraries and functions
source("R/functions_time.R")

if (!file.exists("output/traps")){ dir.create("output/traps")}

cormethod = c("kendall")
ssp_abbr <- c("B_dig","B_pum","C_nit","D_anf","G_bul","G_fal","G_ada","G_cal","G_sip","G_glu","G_min","G_uvu","G_con","G_rubp","G_rubw","G_inf","G_coa","G_cav","G_cra","G_hir","G_men","G_sci","G_the","G_tru","G_tum","G_ung","G_hex","G_rus","G_ten","H_pel","H_dig","N_dut","N_inc","N_pac","O_uni","P_obl", "S_deh","T_fle","T_iot","T_par","G_sac","T_cla","T_hum","T_qui")
traps_meta <- traps_meta[which(traps_meta$trap %in% traps),which(colnames(traps_meta) %in% c("trap","sample_no","latN","lonE","woa13_SSTmn", "mean_open_period_days", "length_days","length_years","length_time.series", "from", "to"))] 

# Creating 3D matrices
array_initial <- array(rep(NaN,length(ssp_abbr)*length(ssp_abbr)*length(traps)), 
                      dim = c(length(ssp_abbr), length(ssp_abbr), length(traps)),
                      dimnames = list(ssp_abbr,ssp_abbr,traps))

# 'Raw' correlation
array_corr <- array_initial
array_p <- array_initial

# First differences
array_fd_corr <- array_initial
array_fd_p <- array_initial

# Relative abundance
array_rel_corr <- array_initial
array_rel_p <- array_initial

# Seasonality - surrogates
array_surr_ses <- array_initial
array_surr_p <- array_initial

residuals_list <- list()
resid_review <- data.frame()

print("Running analysis for each sediment trap:")

for (trap_name in traps){
    # trap_name <- "AB"
  
    print(trap_name) 
    # Creating output folder for this sediment trap
    if (!file.exists(paste("output/traps/",trap_name,sep=""))){ dir.create(paste("output/traps/",trap_name,"/",sep=""))}
    
    raw_data <- read.csv(paste("data/traps_time/",trap_name,".csv", sep = ""), header = T)
    
  
    if(trap_name != "JAM"){
      raw_data$open <- ymd(raw_data$open)
      raw_data$close <- ymd(raw_data$close)
    }else{ # trap_name == "JAM"
      raw_data$open <- dmy(raw_data$open)
      raw_data$close <- dmy(raw_data$close)
    }
    
    
    data <- cbind(open = ymd(raw_data$open), raw_data[,(which(colnames(raw_data)=="cum_days")+1):ncol(raw_data)]) 
    # head(data)
    
    # Removing time-series with just one value (i.e. species present in just one sample)
    data <- data[,sapply(data, function(x) length(which(x!=0)) > 2)]
    
    # Total richness and diversity per trap
    # colSums(data[,-1])
    traps_meta[which(traps_meta$trap == trap_name), "richness"] <- ncol(data) - 1
    traps_meta[which(traps_meta$trap == trap_name), "div_shannon"] <- vegan::diversity(colSums(data[,-1]),index = "shannon") 
    traps_meta[which(traps_meta$trap == trap_name), "div_simpson"] <- vegan::diversity(colSums(data[,-1]),index = "simpson")
    traps_meta[which(traps_meta$trap == trap_name), "total_ind"] <- sum(colSums(data[,-1], na.rm = T), na.rm = T)
    
    
    #############################
    # review: analysis residuals
    # rescaling time-series so that residuals are comparable among sediment traps and species
    data_norm <- data.frame(open = data[,1], apply(data[,2:ncol(data)],2, function(x) (x-min(x))/(max(x)-min(x)) ))
    # calculating splines
    splines_norm <- seasonal_splines(DataSeries = data_norm, DateFormat='%d/%m/%y', SavePlots = F, overwrite = F)
    # calculating residuals
    residuals_norm <- data_norm[,2:ncol(data_norm)] - splines_norm[,2:ncol(data_norm)]
    # merging data
    resid_review <- rbind(resid_review, data.frame(
      species = names(data_norm[,2:ncol(data_norm)]),
      resid_mean_abs = apply(abs(residuals_norm),2, function(x) mean(x, na.rm = T)), 
      resid_median_abs = apply(abs(residuals_norm),2, function(x) median(x, na.rm = T)), 
      resid_mean = apply(residuals_norm,2, function(x) mean(x, na.rm = T)), 
      resid_median = apply(residuals_norm,2, function(x) median(x, na.rm = T)), 
      resid_sd = apply(residuals_norm,2, function(x) sd(x)), 
      n_ind = round(apply(data[2:ncol(data)],2, function(x) sum(x))),
      traps_meta[which(traps_meta$trap == trap_name), 1:14]))
    #############################
    
    
    # Getting columns and rows names of 3D matrix, to fill it with correlations
    cols <- colnames(array_initial[,,which(traps==trap_name)])[colnames(array_initial[,,which(traps==trap_name)]) %in% names(data)]
    rows <- rownames(array_initial[,,which(traps==trap_name)])[rownames(array_initial[,,which(traps==trap_name)]) %in% names(data)]
    
    # Plot scatterplot of everything
    if(!file.exists(paste("output/traps/",trap_name,"/ggpairs_",trap_name,".png",sep=""))){
      png(paste("output/traps/",trap_name,"/ggpairs_",trap_name,".png",sep=""), width = 15, height = 15, units = "in",res = 300)
        print(ggpairs(data[,2:ncol(data)]))
      dev.off()
    }
    
    # Testing if time-series are stationary (see output/traps/trap_name/stationary_test.csv)
    stationary <- test_stationary(data, trap_name, overwrite = F)
    
    # Cross-Correlation (lags)
    lags_correlation <- cross_correlation(data, trap_name, overwrite = F)
    
    
    #######################
    ### Raw correlation ###
    #######################
    
    corr_list <- cormatrix(data[,-1], cormethod, conf.level = 0.95)
    colnames(corr_list$corr) <- rownames(corr_list$corr)<- colnames(data)[-1]
    colnames(corr_list$p_value) <- rownames(corr_list$p_value) <- rownames(corr_list$corr)
    
    corr_vector <- data.frame(corrs = corr_list$corr[lower.tri(corr_list$corr, diag = FALSE)])
    traps_meta[which(traps_meta$trap == trap_name), "corr_mean"] <- mean(corr_vector$corrs, na.rm = T)
    traps_meta[which(traps_meta$trap == trap_name), "corr_sd"] <- sd(corr_vector$corrs, na.rm = T)
    traps_meta[which(traps_meta$trap == trap_name), "corr_median"] <- median(corr_vector$corrs, na.rm = T)
    
    # Filling 3D matrix
    array_corr[rows, cols,which(traps==trap_name)] <- corr_list$corr[rows, cols]
    array_p[rows, cols, which(traps==trap_name)] <- corr_list$p_value[rows, cols]
    # Tranforming NaN into NA
    array_corr[,,which(traps==trap_name)][which(is.nan(array_corr[,,which(traps==trap_name)]))] <- NA
    array_p[,,which(traps==trap_name)][which(is.nan(array_p[,,which(traps==trap_name)]))] <- NA
    
    
    # Plot histogram correlation 
    if(!file.exists(paste("output/traps/",trap_name,"/corrs_hist_",trap_name,".png",sep=""))){
      p <- ggplot(corr_vector, aes(x=corrs)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.05)+
        geom_density(alpha=.2, fill="#FF6666") + labs(x = "Pair-wise correlations") +
        xlim(-1,1)
      png(file =  paste("output/traps/",trap_name,"/corrs_hist_",trap_name,".png",sep=""), width=12, height=10, units = "in", res = 350)
      print(p)
      dev.off()
    }
    
    # Plot matrix correlation 
    if(!file.exists(paste("output/traps/",trap_name,"/corrplot_",trap_name,".png",sep=""))){
      png(file =  paste("output/traps/",trap_name,"/corrplot_",trap_name,".png",sep=""), width=12, height=12, units = "in", res = 350)
      #corrplot.mixed(corr_list$corr)
      corrplot.mixed(corr_list$corr, p.mat = corr_list$p_value, insig = "pch", cl.cex = 1.5, tl.cex = 1.1, tl.col = "black")
      dev.off()
    }
    
    rm(corr_list)
    rm(corr_vector)
    
    
    #########################
    ### First differences ###
    #########################
    if (!file.exists("output/traps_first_diff")){ dir.create("output/traps_first_diff/")}
    
    datafd <- get_first_diff(raw_data, trap_name, days_closed = 10, overwrite = F)
    # days_closed : maximum number of days inbetween samples,if that gap between two samples is bigger than 10 days (next_open > 10 days), then the difference between these two samples is not included
    
    corr_list <- cormatrix(datafd[,-1], cormethod, conf.level = 0.95)
    colnames(corr_list$corr) <- rownames(corr_list$corr)<- colnames(datafd)[-1]
    colnames(corr_list$p_value) <- rownames(corr_list$p_value) <- rownames(corr_list$corr)
    
    corr_vector <- data.frame(corrs = corr_list$corr[lower.tri(corr_list$corr, diag = FALSE)])
    traps_meta[which(traps_meta$trap == trap_name), "corr_fd_mean"] <- mean(corr_vector$corrs, na.rm = T)
    traps_meta[which(traps_meta$trap == trap_name), "corr_fd_sd"] <- sd(corr_vector$corrs, na.rm = T)
    traps_meta[which(traps_meta$trap == trap_name), "corr_fd_median"] <- median(corr_vector$corrs, na.rm = T)
    
    # Filling 3D matrix
    array_fd_corr[rows, cols,which(traps==trap_name)] <- corr_list$corr[rows, cols]
    array_fd_p[rows, cols, which(traps==trap_name)] <- corr_list$p_value[rows, cols]
    # Tranforming NaN into NA
    array_fd_corr[,,which(traps==trap_name)][which(is.nan(array_fd_corr[,,which(traps==trap_name)]))] <- NA
    array_fd_p[,,which(traps==trap_name)][which(is.nan(array_fd_p[,,which(traps==trap_name)]))] <- NA
    
    
    # Plot histogram correlation 
    if(!file.exists(paste("output/traps/",trap_name,"/corrs_hist_",trap_name,"_fd.png",sep=""))){
      p <- ggplot(corr_vector, aes(x=corrs)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.05)+
        geom_density(alpha=.2, fill="#FF6666") + labs(x = "Pair-wise correlations")+
        xlim(-1,1)
      png(file =  paste("output/traps/",trap_name,"/corrs_hist_",trap_name,"_fd.png",sep=""), width=12, height=10, units = "in", res = 350)
      print(p)
      dev.off()
    }
    
    # Plot matrix correlation 
    if(!file.exists(paste("output/traps/",trap_name,"/corrplot_",trap_name,"_fd.png",sep=""))){
      png(file =  paste("output/traps/",trap_name,"/corrplot_",trap_name,"_fd.png",sep=""), width=12, height=12, units = "in", res = 350)
        corrplot.mixed(corr_list$corr, p.mat = corr_list$p_value, insig = "pch", cl.cex = 1.5, tl.cex = 1.1, tl.col = "black")
      dev.off()
    }
    
    rm(corr_list)
    rm(corr_vector)
   
    
    
    ##################
    ### Surrogates ###
    ##################
    # Creating null model based on species' yearly seasonality
  
    # Splines and species seasonality
    splines <- seasonal_splines(DataSeries = data, DateFormat='%d/%m/%y', SavePlots = T, overwrite = T)

    # Residuals 
    residuals <- data[,2:ncol(data)] - splines[,2:ncol(data)]
    residuals_list[[trap_name]]<- residuals
      
    # Correlation and Surrogates
    corr_surrog <- corr_surrogates(data, splines, trap_name, nreps=100, cormethod, overwrite = T)   
    # Surrogates: randomization of residuals, summed with splines to generate null series (nreps = 500 null series, columns V1 - V500)
    # corr_surrog[1:50,1:15]
    
    # Plotting box-plots: correlations with surrogate distribution for each focal variable
    corr_surrogates_boxplots(corr_surrog, trap_name, overwrite=T)
    
    # Symmetrical matrix
    # Observed correlation Standardized Effect Size (SES) related to surrogate null distribution
    surr_ses_matrix <- graph.data.frame(corr_surrog[,c("var1","var2","ses_surr")], directed=FALSE)
    surr_ses_matrix <- get.adjacency(surr_ses_matrix, attr="ses_surr", sparse=FALSE)
    # Observed correlation p-value related to surrogate null distribution
    surr_p_matrix <- graph.data.frame(corr_surrog[,c("var1","var2","surr_p")], directed=FALSE)
    surr_p_matrix <- get.adjacency(surr_p_matrix, attr="surr_p", sparse=FALSE)
    # Transforming two tailed p-value into one tailed
    surr_p_matrix[which(surr_p_matrix > 0.975)] <- 1 - surr_p_matrix[which(surr_p_matrix > 0.975)]
    surr_p_matrix[which(surr_p_matrix > 0.025 & surr_p_matrix < 0.05)] <- surr_p_matrix[which(surr_p_matrix > 0.025 & surr_p_matrix < 0.05)] + 0.025
    
    diag(surr_ses_matrix) <- max(surr_ses_matrix)
    
    # Plot matrix SES 
    if(!file.exists(paste("output/traps/",trap_name,"/corrplot_",trap_name,"_surr.png",sep=""))){
      png(file = paste("output/traps/",trap_name,"/corrplot_",trap_name,"_surr.png",sep=""), width=12, height=12, units = "in", res = 350)
        corrplot.mixed(surr_ses_matrix, p.mat = surr_p_matrix, is.corr = FALSE, insig = "pch", cl.cex = 1.5, tl.cex = 1.1, tl.col = "black")
      dev.off()
    }
    
    # Vector for histogram
    ses_vector <- data.frame(corrs = surr_ses_matrix[lower.tri(surr_ses_matrix, diag = FALSE)])
    traps_meta[which(traps_meta$trap == trap_name), "ses_surg_mean"] <- mean(ses_vector$corrs, na.rm = T)
    traps_meta[which(traps_meta$trap == trap_name), "ses_surg_sd"] <- sd(ses_vector$corrs, na.rm = T)
    traps_meta[which(traps_meta$trap == trap_name), "ses_surg_median"] <- median(ses_vector$corrs, na.rm = T)
    
    # Plot histogram SES
    if(!file.exists(paste("output/traps/",trap_name,"/corrs_hist_",trap_name,"_surr.png",sep=""))){
      p <- ggplot(ses_vector, aes(x=corrs)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.2)+
        geom_density(alpha=.2, fill="#FF6666") + labs(x = "SES observed correlation vs. surrogate")+
        geom_vline(xintercept = 0, linetype = "dashed")
      png(file =  paste("output/traps/",trap_name,"/corrs_hist_",trap_name,"_surr.png",sep=""), width=12, height=10, units = "in", res = 350)
       print(p)
      dev.off()
    }
    
    # Filling 3D matrix
    # Seasonality - surrogates
    array_surr_ses[rows, cols,which(traps==trap_name)] <- surr_ses_matrix[rows, cols]
    array_surr_p[rows, cols, which(traps==trap_name)] <- surr_p_matrix[rows, cols]
    # Tranforming NaN into NA
    array_surr_ses[,,which(traps==trap_name)][which(is.nan(array_surr_ses[,,which(traps==trap_name)]))] <- NA
    array_surr_p[,,which(traps==trap_name)][which(is.nan(array_surr_p[,,which(traps==trap_name)]))] <- NA
    
    rm(surr_ses_matrix)
    rm(surr_p_matrix)
    rm(ses_vector)
    
} # for traps

# Saving trap_meta with average correlations
write.csv(traps_meta, "output/traps_meta.csv", row.names = F)

# Paper review: analysis residuals
write.csv(resid_review,  "output/resid_review.csv", row.names = F)


# Saving 3D matrices as R files
if(!file.exists("output/array_corr.RData")){
  # Raw correlation
  save(array_corr, file =  "output/time_corr.RData")
  save(array_p, file =  "output/time_corr_p.RData")
  # First differences
  save(array_fd_corr, file =  "output/time_fd_corr.RData")
  save(array_fd_p, file =  "output/time_fd_p.RData")
  # Seasonality - surrogates
  save(array_surr_ses, file =  "output/time_surr_ses.RData")
  save(array_surr_p, file =  "output/time_surr_p.RData")
  # Residuals list
  save(residuals_list, file =  "output/residuals_list.RData")
}


###
### Results
###

pval <- 0.01

### Correlation
ses <- melt(array_surr_ses)
p <- melt(array_surr_p)[,4]
colnames(ses) <- c("ssp1","ssp2","trap","ses")
ses <- cbind(ses, p)
# Remove diagonal
ses <- ses[-which(ses$ssp1 == ses$ssp2),] 
# Remove symmetricals
ses <- ses[!duplicated(t(apply(ses[,1:3], 1, sort))),]
# Remove NAs
ses <- ses[complete.cases(ses),]

# Length significant 
print("-----------------------")
print("Time series analysis results")
print(paste0("Total pair-wise comparisons: ", length(ses[,1]) ))
print(paste0("p-value = ", pval))
print(paste0("Significant comparisons: ", length(which(ses$p <= pval)) ))
print(paste0("Significant positive: ", length(which(ses$p <= pval & ses$ses > 0)) ))
print(paste0("Significant negative: ", length(which(ses$p <= pval & ses$ses < 0))))

# Negative interactions
# ses[which(ses$p < pval & ses$ses < 0),]

ses_meta <- merge(ses, traps_meta)
names(ses_meta)


###
### Plots
###
print("-----------------------")
print("Generating temporal plots")

# Theme for phylogeny, trait and symbiosis plots
th <- theme_bw() + 
  theme(axis.title = element_text(colour = 'black', size = 28),
        axis.text = element_text(colour = 'black', size = 26),
        plot.margin=unit(c(0.2,0.2,0.2,0.5),"cm")) 

# Phylogeny
phy <- read.tree("data/aLb_Rillo_modified.tre")
phy <- dropExtinct(phy)
phy$tip.label <-  c("G_coa","G_ung","G_tum","G_men","G_cav","G_tru","G_cra","G_sci","G_the","G_hir","G_inf","P_obl",
                    "N_dut","N_pac","N_inc","G_hex","G_rubw","G_rubp","G_con","G_ten","G_rus","S_deh","G_sac","O_uni",
                    "T_hum","T_cla","T_qui","G_ada","G_cal","G_sip","B_dig","G_bul","G_fal")
phydist <- cophenetic(phy)
phydist <- melt(phydist)
names(phydist) <- c("ssp1","ssp2", "phydist")

ses_phy <- merge(ses_meta, phydist)
ses_phy$div_time <- ses_phy$phydist/2
ses_phy$div_time <- round(ses_phy$div_time,2)

p_phy <- ggplot(data = ses_phy, aes(y=ses, x = div_time, group = div_time)) +  
  geom_jitter(data = ses_phy[which(ses_phy$p > pval),],  color = "grey60", fill = alpha("grey80",.4),  shape = 21, size = 3, stroke = 1) +
  geom_jitter(data = ses_phy[which(ses_phy$p <= pval & ses_phy$ses > 0),], color = "blue", fill = alpha("blue",.4),  shape = 21, size = 3, stroke = 1) +
  geom_jitter(data = ses_phy[which(ses_phy$p <= pval & ses_phy$ses < 0),],  color = "red", fill = alpha("red",.4),  shape = 21, size = 3, stroke = 1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  geom_boxplot(data = ses_phy[which(ses_phy$div_time %in% names(which(table(ses_phy$div_time)> 40))),],
               fill = alpha('white',0.7), col =  alpha('black',1), outlier.shape = NA, lwd=0.7, width = 1.5) +
  xlab("Divergence time (My)") + ylab("SES (time-series)") +
  scale_x_continuous(breaks = seq(0,65,15),limits=c(-2,70), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-4,12,4),limits=c(-5,13.5), expand = c(0, 0)) +
  th

png(file = "output/fig_time_ses_phy.png", width=10, height=7, units = "in", res = 400)
 print(p_phy)
dev.off()



# Shell size data
diam <- read.csv(file ="data/shell_diameter_data.csv", header=TRUE, stringsAsFactors = FALSE) 
sp_diam <- ddply(diam,~species,summarise,mean=mean(max_diam_mu),sd=sd(max_diam_mu), n_ind = length(max_diam_mu))
sp_diam <- rbind(sp_diam, c(sp_diam[which(sp_diam$species == "ruber"),c("species","mean","sd","n_ind")]))
sp_diam <- sp_diam[-which(sp_diam$species == "riedeli"),]
sp_diam$sspname <- row.names(sp_diam) <- c("D_anf","G_bul","G_cal","T_cla","G_con","G_coa","G_cra","S_deh","N_dut","G_fal","T_fle","G_glu",
                     "G_hex","G_hir","T_hum","G_inf","T_iot","G_men","G_min","C_nit","P_obl","N_pac","T_par","B_pum",
                     "T_qui","G_rubw","G_rus","G_sac","G_sci","G_sip","G_ten","G_tru",
                     "G_tum","O_uni","G_uvu", "G_rubp")
# not included "B_dig" "G_ada" "G_cav"  "G_rus" "G_the" "G_ung" "H_dig" "H_pel" "N_inc"                   

traitdist <- as.matrix(dist(sp_diam[,"mean",drop=F])) 
# traitdist <- traitdist/max(traitdist) # normalized
traitdist  <- melt(traitdist)
names(traitdist) <- c("ssp1","ssp2", "traitdist")
traitdist$round <- 50*round(traitdist$traitdist/50)

ses_trait <- merge(ses_meta, traitdist)

p_size <- ggplot(data = ses_trait, aes(y=ses, x = traitdist, group = round)) +  
  geom_point(data = ses_trait[which(ses_trait$p > pval),],  color = "grey60", fill = alpha("grey80",.4),  shape = 21, size = 3, stroke = 1) +
  geom_point(data = ses_trait[which(ses_trait$p <= pval & ses_trait$ses > 0),], color = "blue", fill = alpha("blue",.4),  shape = 21, size = 3, stroke = 1) +
  geom_point(data = ses_trait[which(ses_trait$p <= pval & ses_trait$ses < 0),],  color = "red", fill = alpha("red",.4),  shape = 21, size = 3, stroke = 1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  geom_boxplot(fill = alpha('white',0.7), col =  alpha('black',1), outlier.shape = NA, lwd=0.7, width = 1) +
  xlab(expression(paste("Size difference (",mu,m,")", sep=""))) + ylab("SES (time-series)") +
  scale_x_continuous(breaks = seq(0,800,100),limits=c(-25,850), expand = c(0, 0)) +
  #scale_y_continuous(breaks = seq(-4,12,4),limits=c(-5,13.5), expand = c(0, 0)) +
  th

png(file = "output/fig_time_ses_size.png", width=10, height=7, units = "in", res = 400)
 print(p_size)
dev.off()


# Symbiosis correlation
ssp_symb <- read.csv(file="data/symbiosis.csv", header=TRUE, stringsAsFactors = FALSE)
ssp_symb <- ssp_symb[which(ssp_symb$abbr %in% ssp_abbr),]
row.names(ssp_symb) <- ssp_symb$abbr
ssp_symb <- ssp_symb[!is.na(ssp_symb$binary),3, drop = FALSE]
ssp_symb_mat <- as.matrix(ssp_symb)%*%as.matrix(t(ssp_symb))
symbpair <- melt(ssp_symb_mat)
names(symbpair) <- c("ssp1","ssp2", "symbpair")


ses_symb <- merge(ses_meta, symbpair)
ses_symb$symbpair <- as.factor(ses_symb$symbpair)
levels(ses_symb$symbpair) <- list(nonsymb_nonsymb="1", symb_nonsymb="2", symb_symb="4")
ses_symb$sign_color <- NA
ses_symb[which(ses_symb$p > pval),"sign_color"] <- "grey60"
ses_symb[which(ses_symb$p <= pval & ses_symb$ses > 0),"sign_color"] <- "blue"
ses_symb[which(ses_symb$p <= pval & ses_symb$ses < 0),"sign_color"] <- "red"
ses_symb$sign_color <- as.factor(ses_symb$sign_color)

p_symb <-  ggplot(data = ses_symb, aes(y=ses, x = symbpair, fill = sign_color, color = sign_color)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  geom_quasirandom(varwidth = TRUE, bandwidth = 0.5, width = 0.5, pch = 21,size = 3, stroke = 1) +
  geom_boxplot(fill = alpha('white',0.7), col =  alpha('black',1), outlier.shape = NA, lwd=0.7, width = 0.2) +
  ylab("SES (time-series)") + 
  scale_fill_manual(values = alpha(c("blue","grey80","red"),0.5)) +
  scale_color_manual(values = c("blue","grey60","red")) +
  scale_x_discrete(labels = c("Heterotrophs", "Heterotroph \n & Phototroph", "Phototrophs")) +
  scale_y_continuous(breaks = seq(-4,12,4),limits=c(-5,13.5), expand = c(0, 0)) +
  th +
  theme(legend.position="none", axis.title.x = element_blank())

png(file = "output/fig_time_ses_symb.png", width=10, height=7, units = "in", res = 400)
  print(p_symb)
dev.off()


# OR: https://stackoverflow.com/questions/38477616/overlapping-points-when-using-fill-aesthetic-in-ggplot2-geom-dotplot-in-r/38479399#38479399

# Plot per trap

load(file =  "output/array_surr_ses.RData")
load(file =  "output/array_surr_p.RData")

dim(array_surr_ses)
str(array_surr_ses)

trap_ses <- data.frame(matrix(NA, ncol = 3, nrow = 1))
colnames(trap_ses) <- c("trap_name", "surr_ses", "surr_p")

for (i in 1:35) {
  t <- as.matrix(array_surr_ses[,,i])
  p <- as.matrix(array_surr_p[,,i])
  t <- t[lower.tri(t)]
  p <- p[lower.tri(p)]
  p <- p[!is.na(t)]
  t <- t[!is.na(t)]
  trap_ses <- rbind(trap_ses, data.frame(trap_name = dimnames(array_surr_ses)[[3]][i], surr_ses = t, surr_p = p))
}
trap_ses <- trap_ses[-1,]
trap_ses$traps_col <- ifelse(trap_ses$surr_ses < 0, "red", "blue")
trap_ses[which(trap_ses$surr_p > pval),"traps_col"] <- "grey80"
head(trap_ses)


png(file = "output/fig_time_ses_traps.png", width=13, height=6, units = "in", res = 400)
ggplot(data = trap_ses , aes(y = surr_ses, x = trap_name, group = trap_name)) +   theme_bw() +
  geom_quasirandom(color = as.character(trap_ses[,"traps_col"]), fill = alpha(as.character(trap_ses[,"traps_col"]),.3),  shape = 21, stroke = 0.7,varwidth=TRUE) +
  geom_boxplot(fill = alpha('white',0), col = 'black', outlier.shape = NA, lwd=0.7, width = 0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  ylab("SES (time-series)") +
  theme(axis.title.y = element_text(colour = 'black', size = 20),
        axis.text.y = element_text(colour = 'black', size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 14, angle = 45, hjust = 1),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  scale_y_continuous(breaks = seq(-4,12,2),limits=c(-5,13), expand = c(0, 0)) 
dev.off()



# Plot first differences per trap

load(file =  "output/array_fd_corr.RData")
load(file =  "output/array_fd_p.RData")

trap_fd <- data.frame(matrix(NA, ncol = 3, nrow = 1))
colnames(trap_fd) <- c("trap_name", "fd_corr", "fd_p")

for (i in 1:35) {
  t <- as.matrix(array_fd_corr[,,i])
  p <- as.matrix(array_fd_p[,,i])
  t <- t[lower.tri(t)]
  p <- p[lower.tri(p)]
  p <- p[!is.na(t)]
  t <- t[!is.na(t)]
  trap_fd <- rbind(trap_fd, data.frame(trap_name = dimnames(array_fd_corr)[[3]][i], fd_corr = t, fd_p = p))
}
trap_fd <- trap_fd[-1,]
trap_fd$traps_col <- ifelse(trap_fd$fd_corr < 0, "red", "blue")
trap_fd[which(trap_fd$fd_p > pval),"traps_col"] <- "grey60"
head(trap_fd)

png(file = "output/si_fig_time_ses_fd.png", width=13, height=6, units = "in", res = 500)
ggplot(data = trap_fd, aes(y = fd_corr, x = trap_name, group = trap_name)) +   theme_bw() +
  geom_quasirandom(color = as.character(trap_fd[,"traps_col"]), fill = alpha(as.character(trap_fd[,"traps_col"]),.3),  shape = 21, stroke = 0.7,varwidth=TRUE) +
  geom_boxplot(fill = alpha('white',0), col = 'black', outlier.shape = NA, lwd=0.7, width = 0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  ylab("Correlation (first differences)") +
  theme(axis.title.y = element_text(colour = 'black', size = 18),
        axis.text.y = element_text(colour = 'black', size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 13, angle = 45, hjust = 1),
        plot.margin=unit(c(0.5,0.2,0.2,0.2),"cm")) +
  scale_y_continuous(breaks = seq(-1,1,0.2),limits=c(-1,1), expand = c(0, 0)) 
dev.off()



# Plot Interannual variation (residuals)

load(file =  "output/residuals_list.RData")
mean_res_cor <- c()

for(i in 1 : length(residuals_list)){
  print(names(residuals_list)[i])
  cormat <- cor(residuals_list[[i]], method = cormethod)
  mean_res_cor <- c(mean_res_cor, mean(cormat))
}

traps_meta$mean_res_cor <- mean_res_cor

png(file = "output/si_fig_interannual_residuals.png", width=8, height=6, units = "in", res = 500)
ggplot(data = traps_meta[which(traps_meta$length_years>2),], 
       aes(x = mean_res_cor, y = length_years)) + 
  geom_point(pch = 21, fill = 'orangered3', color = "black", size = 4, stroke = 1) + theme_bw() +
  geom_text_repel(aes(label = trap), point.padding = 0.2, force = 2, box.padding = 0.3, size = 4.3, fontface = "bold" ) +
  labs(y = "Length sediment trap sampling (years)", x = "Mean correlation of species' residuals series") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size=16, colour = "black")) +
  scale_x_continuous(breaks = seq(0.1, 0.6, 0.1), limits=c(0.05, 0.65), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(2, 12, 1), limits=c(1, 13), expand = c(0, 0))
dev.off()


##################################
# Paper review: analysis residuals

resid_review <- read.csv("output/resid_review.csv", header = T)
resid_review$resolution <- round(resid_review$mean_open_period_days, 0)

label_data <- unique(resid_review[,c("trap","resolution", "length_years" , "length_time.series")])
max_res <- as.data.frame(cbind(max_res = by(resid_review$resid_sd, resid_review$trap, max)))
max_res$trap <- row.names(max_res)
label_data <- merge(label_data, max_res)

# unique(resid_review[order(resid_review$length_days), c("trap", "length_days")])

# Resolution
png(file = "output/si_fig_resid_mn_res.png", width=12, height=7, units = "in", res = 500)
ggplot(data = resid_review, 
       aes(x = as.factor(resolution), # resolution , length_years , length_time.series
           y = resid_mean_abs, group = trap, fill = as.factor(trap))) + 
  geom_boxplot(fill = alpha('white',0.2), col = 'black', lwd=0.5, width = 0.5) +
  labs(x = "Resolution time series (days per sample)", # Resolution time series (days per sample)" "Duration of sampling (years)" "Length time series (number of samples)",
       y = "Average residuals (norm.)") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size=12, colour = "black"),
        legend.position = "none") +
  #scale_x_continuous(breaks = seq(0, 12, 0.5), limits=c(0, 12.5), expand = c(0, 0)) +
  scale_fill_discrete(drop=FALSE) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits=c(-0.1, 0.4), expand = c(0, 0))
dev.off()

# Years
png(file = "output/si_fig_resid_mn_years.png", width=12, height=7, units = "in", res = 500)
ggplot(data = resid_review, 
       aes(x = as.factor(length_years), 
           y = resid_mean_abs, group = trap, fill = as.factor(trap))) + 
  geom_boxplot(fill = alpha('white',0.2), col = 'black', lwd=0.5, width = 0.5) +
  labs(x = "Duration of sampling (years)",
       y = "Average residuals (norm.)") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size=12, colour = "black"),
        legend.position = "none") +
  #scale_x_continuous(breaks = seq(0, 12, 0.5), limits=c(0, 12.5), expand = c(0, 0)) +
  scale_fill_discrete(drop=FALSE) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits=c(-0.1, 0.4), expand = c(0, 0))
dev.off()


# Length
png(file = "output/si_fig_resid_mn_length.png", width=12, height=7, units = "in", res = 500)
ggplot(data = resid_review, 
       aes(x = as.factor(length_time.series),
           y = resid_mean_abs, group = trap, fill = as.factor(trap))) + 
  geom_boxplot(fill = alpha('white',0.2), col = 'black', lwd=0.5, width = 0.5) +
  labs(x = "Length time series (number of samples)",
       y = "Average residuals (norm.)") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size=12, colour = "black"),
        legend.position = "none") +
  #scale_x_continuous(breaks = seq(0, 12, 0.5), limits=c(0, 12.5), expand = c(0, 0)) +
  scale_fill_discrete(drop=FALSE) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits=c(-0.1, 0.4), expand = c(0, 0))
dev.off()







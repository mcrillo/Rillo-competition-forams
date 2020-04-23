
########################################################################
### Functions main_time-series (sediment trap data)
########################################################################

get_first_diff <- function(data, trap_name, days_closed, overwrite){ # data.frame, string
  # days_closed : maximum number of days inbetween samples,if that gap between two samples is bigger than 10 days (next_open > 10 days), then the difference between these two samples is not included
  
  if(overwrite == TRUE | !file.exists(paste("output/traps_first_diff/",trap_name,"_first_diff.csv",sep=""))){
    
    # Differentiating data
    datafd <- as.data.frame(apply(data[,(which(colnames(data)=="cum_days")+1):ncol(data)], 2, diff))
    datafd <- cbind(open=ymd(data[-1,"open"]), datafd)
    
    # Excluding first difference between samples that had a gap inbetween bigger than days_closed = 10 days (next_open > 10 days)
    if (any(na.omit(data$next_open) > days_closed)){
      rows <- which(data$next_open > days_closed)
      datafd <- datafd[-rows,] 
    }
    
    write.csv(datafd, paste("output/traps_first_diff/",trap_name,"_first_diff.csv",sep=""), row.names = FALSE)
    return (datafd)
    
  }else{
    
    datafd <- read.csv(paste("output/traps_first_diff/",trap_name,"_first_diff.csv",sep=""), header = TRUE)
    return (datafd)
  }
  
}

########################################################################

cormatrix <- function(mat, cormethod, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  cor.mat <- p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(cor.mat) <- 1
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level, 
                      method=cormethod, alternative="two.sided")
      cor.mat[i,j] <- cor.mat[j,i] <- tmp$estimate
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      
      if(cormethod=="pearson"){
        lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
        uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
      }
    }
  }
  return(list(corr = cor.mat, p_value = p.mat, lower_ci = lowCI.mat, upper_ci = uppCI.mat))
}

########################################################################

test_stationary <- function(data, trap_name, overwrite){ 
  
  
  if(overwrite == TRUE | !file.exists(paste("output/traps/",trap_name,"/stationary_test_",trap_name,".csv",sep=""))){
    
    
    statio <- data.frame()
    
    for (i in 2:ncol(data)){
      # print(i)
      statio[i, "variable"] <- names(data[i])
      
      # ADF
      adf_test <- adf.test(data[,i])
      statio[i,"adf_method"] <- adf_test$method
      statio[i,"adf_stats"] <- adf_test$statistic
      statio[i,"adf_lag_order"] <- adf_test$parameter
      statio[i,"alternative_hyp"] <- adf_test$alternative
      statio[i,"adf_p_value"] <- adf_test$p.value
      
      # KPSS
      kpss_test <- kpss.test(data[,i])
      statio[i,"kpss_method"] <- kpss_test$method
      statio[i,"kpss_level"] <- kpss_test$statistic
      statio[i,"kpss_trunc_lag"] <- kpss_test$parameter
      statio[i,"kpss_p_value"] <- kpss_test$p.value
      
    }
    
    statio <- statio[-1,]
    
    write.csv(statio, paste("output/traps/",trap_name,"/stationary_test_",trap_name,".csv",sep=""), row.names = FALSE)
    
    return(statio)
    
    #paste("data/",trap_name,"/",trap_name,"_first_diff.csv",sep="")
    
    
  }else{
    
    statio <- read.csv(paste("output/traps/",trap_name,"/stationary_test_",trap_name,".csv",sep=""), header = TRUE)
    return (statio)
  }
  
}

########################################################################

cross_correlation <- function(data, trap_name, overwrite){ 
  
  
  if(overwrite == TRUE | !file.exists(paste("output/traps/",trap_name,"/cross_correlation_",trap_name,".csv",sep=""))){
    
    
    ccf_all <- data.frame()
    
    # combination of all pairs of variables to test for cross-correlation and lags
    comb_two <- combn(2:ncol(data),2)
    
    for(i in 1:length(comb_two[1,])){ # for each pair of variables
      
      cross_corr <- ccf(data[,comb_two[1,i]], data[,comb_two[2,i]], plot=FALSE)
      cross_df <- data.frame(lag = cross_corr$lag, cr_corr = cross_corr$acf, cr_corr_abs = abs(cross_corr$acf))
      
      ccf_all <- rbind(ccf_all, cbind(var1=names(data)[comb_two[1,i]],
                                      var2=names(data)[comb_two[2,i]],
                                      cross_df[which.max(cross_df$cr_corr_abs),],
                                      cross_df[which.max(cross_df$cr_corr_abs[cross_df$cr_corr_abs!=max(cross_df$cr_corr_abs)]),],
                                      cross_df[which(cross_df$lag == 0),]
      ) ) 
      
      
    } # for
    
    # Regression (think about negative lag values and who is y
    
    #lag1 <- cross_df[which.max(cross_df$cr_corr_abs),"lag"]
    #lag2 <- cross_df[which.max(cross_df$cr_corr_abs[cross_df$cr_corr_abs!=max(cross_df$cr_corr_abs)]),"lag"]
    
    # Lag 2 (-39)
    # yplus39=tail(diffnumberpeople,26025)
    # reg2<-lm(yplus39~difftemperature[1:26025])
    # summary(reg2)
    
    # reg2 <- lm(y ~ x)
    # summary(reg2)$adj.r.squared
    # summary(reg2)$r.squared
    # summary(reg2)$coefficients[2,1] # Estimate slope regression 
    # summary(reg2)$coefficients[2,4] # p-value
    # library(lmtest)
    # dwtest(reg)
    
    ccf_all[,grep(names(ccf_all), pattern=c("corr"))] <- round(ccf_all[,grep(names(ccf_all), pattern=c("corr"))],2)
    
    write.csv(ccf_all, paste("output/traps/",trap_name,"/cross_correlation_",trap_name,".csv",sep=""), row.names = FALSE)
    return(ccf_all)
    
  }else{
    
    ccf_all <- read.csv(paste("output/traps/",trap_name,"/cross_correlation_",trap_name,".csv",sep=""), header = TRUE)
    return (ccf_all)
  }
  
}
# code from http://www.michaeljgrogan.com/cross-correlation-r/

########################################################################

seasonal_splines <- function(DataSeries,DateFormat, SavePlots, overwrite){
  #---
  # Author: Brenno Cabella (SeasonalSplines.R modified by Marina Rillo)
  # This code extracts the sazonality from time series using splines.
  # Inputs:   DataSeries - Time-series with only one date ("open") column
  #           DateFormat - Format of date (usually '%d/%m/%Y')
  #           SavePlots - Save sazonality results as figures (T or F)
  # Outputs:  Figure if SavePlots=T
  #           csv file with the seasonality for all variables in DataFile   
  #---
  
  if(overwrite == TRUE | !file.exists(paste("output/traps/",trap_name,"/splines_",trap_name,".csv",sep=''))){
    
    # read data
    DateCol=1 # "open" column
    # names variables
    VariablesNames<-names(DataSeries)
    # attribute date format to date column
    DataSeries[,DateCol]<-as.Date(DataSeries[,DateCol], DateFormat)
    # save original date format for plot
    OriginalDate<-DataSeries[,DateCol]
    # change date format to track days and months only
    AuxDate<-format(DataSeries[,DateCol],"%m/%d")
    DataSeries[,DateCol]<-format(as.Date(DataSeries[,DateCol], DateFormat), format='%m/%d')
    # create output matrix
    Seasonal<-matrix(NA,nrow(DataSeries),ncol(DataSeries))
    # agregate same month-day, taking the mean of species' flux in the same month-day over the years
    AgDados<-aggregate(DataSeries[,-DateCol], by=list(DataSeries[,DateCol]), FUN=mean, na.rm=TRUE)
    # weight for the splines considering the repetitions of each month-day
    WDados <- table(DataSeries[,DateCol])
    # match the dates in the original data with the aggregated ones
    MatchDate<-match(AuxDate,AgDados[,1])
    
    # loop for variables
    for (i in 2:ncol(DataSeries)){
      # replicate 3 times weights and smooth to avoid discontinuities
      # 3 times so that the middle series (the one we are using) does not have weird discontinuities in the beginning and the end
      AgDadosRep<-rep(AgDados[,i],3)
      spline <- smooth.spline(AgDadosRep,w=rep(WDados,3))
      
      # str(spline)
      # all(spline$data$y == AgDadosRep) # spline$data is the original data
      SeasonalDataRep<-spline$y # estimated seasonality
      
      # take the middle part of the spline result (all 3 parts are identical)
      SeasonalDataAux<-SeasonalDataRep[(floor(length(SeasonalDataRep)/3)+1):floor((2*length(SeasonalDataRep)/3))]
      # store result in output matrix
      Seasonal[,i]<-SeasonalDataAux[MatchDate]
      
      # test of how seasonal a species is: correlation of raw data with spline
      # SspSeasonal <- rbind(SspSeasonal, data.frame(variable = VariablesNames[i], correlation = cor(spline$data$y, spline$y)))
      # write.csv(SspSeasonal,"output/splines/seasonal_species_corr.csv",row.names = F)
      # cor(DataSeries[,i],SeasonalDataAux[MatchDate])
      
      # save plot
      if (SavePlots==T){
        if (!file.exists(paste("output/traps/",trap_name,"/splines_plots",sep=""))){ dir.create(paste("output/traps/",trap_name,"/splines_plots/",sep=""))}
        setEPS()
        postscript(paste("output/traps/",trap_name,"/splines_plots/",VariablesNames[i],'_',trap_name,'.eps',sep=''),width=8.5,height=6)
        plot(OriginalDate,DataSeries[,i],ylab=VariablesNames[i],xlab='date')
        lines(OriginalDate,Seasonal[,i],col='blue',lwd=2)
        dev.off()
      } # if
    } # for
    
    # store original date in output
    Seasonal[,1]<-format(as.Date(OriginalDate,'%y-%m-%d'),'%Y-%m-%d')
    # insert columns names
    colnames(Seasonal)<-VariablesNames
    # write csv output file
    
    write.csv(Seasonal,paste("output/traps/",trap_name,"/splines_",trap_name,".csv",sep=''),row.names = F)
    
    splines <- as.data.frame(Seasonal, stringsAsFactors = FALSE)
    splines[,2:ncol(splines)]<- sapply(splines[,2:ncol(splines)], as.numeric)
    splines$open <- ymd(splines$open) # transforming to date format
    # data.frame(names(columns), names(splines))
    
    
    # Seasonality test: correlation between a time-series and its spline
    seasonality_species <- data.frame(species = rep(NA,ncol(splines)), corr_spline = rep(NA,ncol(splines)))
    for (i in 2:ncol(splines)){
      seasonality_species[i,"species"] <- names(splines)[i]
      seasonality_species[i,"corr_spline"] <- cor(splines[,i], data[,i])
    }
    seasonality_species <- seasonality_species[-1,]
    write.csv(seasonality_species, paste("output/traps/",trap_name,"/seasonality_",trap_name,".csv",sep=''),row.names = F)
    
    
    return(splines)
    
  }else{
    
    splines <- read.csv(paste("output/traps/",trap_name,"/splines_",trap_name,".csv",sep=''), header = TRUE)
    return (splines)
  }
  
}

########################################################################

corr_surrogates <- function(data, splines, trap_name, nreps, corr_method, overwrite){     
  
  
  if(overwrite == TRUE | !file.exists(paste("output/traps/",trap_name,"/corr_surrogates_",trap_name,".csv",sep=""))){
    
    # (1) # Randomization of residuals of (data - splines)
    # data.frame(data = data[,i], spline = splines[,i], resid = data[,i] - splines[,i])
    random_resid <- list()
    
    for(i in (2:ncol(splines))){
      set.seed(i)
      random_resid[[i-1]] <- replicate(nreps, sample(c(data[,i] - splines[,i]), size = length(splines[,i]), replace = FALSE, prob = NULL))
    }
    
    # OLD: random_resid <- lapply(2:ncol(splines), function(i) replicate(nreps, sample(c(data[,i] - splines[,i]), size = length(splines[,i]), replace = FALSE, prob = NULL)))
    # str(random_resid) # Each element of the list is a matrix with 500 randomized residual series
    # var(random_resid[[i]][,90]) # double check if variance is always the same, i.e., if vectors are all identical just with values in diferent order (randomized)
    
    random_resid <- lapply(1:length(random_resid), function(i) random_resid[[i]][,order(sample(1:nreps))])
    names(random_resid) <- names(splines)[2:ncol(splines)]
    # Shuffle residuals order among species (so that they were not all sampled in the same order, just in case)
    sample(1:nreps, size = length(random_resid))
    
    # (2) # Generating surrogates for each species: adding random residuals to splines (i.e., seasonal time-series)
    # data.frame(resid = random_resid[[i]][,500], spline = splines[,i],  surrogate = random_resid[[i]][,500]+ splines[,i])
    surrogates <- lapply(1:length(random_resid), function(i){random_resid[[i]] + splines[,i+1]}) # i+1 because first column of splines is "open"
    names(surrogates) <- names(random_resid)
    # str(surrogates) # Each element of the list is a matrix with 500 null series (= spline + random residual)
    # all(round(surrogates[[i]][,500] - random_resid[[i]][,500] - splines[i+1],4)==0) # double check: this must be all 0 for any value of i in [1,15] and j[1,500]
    
    # Calculate correlations column by columns of each variable pair
    corr_pairs <- data.frame()
    comb_two <- combn(1:length(surrogates),2)
    
    for(i in 1:length(comb_two[1,])){ # i = 1
      
      col_var1 <- which(names(data)==names(surrogates)[(comb_two[1,i])])
      col_var2 <- which(names(data)==names(surrogates)[(comb_two[2,i])])
      
      # correlation of each variable with its own spline (measure of seasonality of the time-series)
      season_var1 <- cor(data[,col_var1], splines[,col_var1])
      season_var2 <-  cor(data[,col_var2], splines[,col_var2])
      
      # correlations
      corr_obs <-  cor.test(data[,col_var1], data[,col_var2], method = corr_method)$estimate
      corr_obs_p <- cor.test(data[,col_var1], data[,col_var2], method = corr_method)$p.value
      
      corr_splines <- cor.test(splines[,col_var1], splines[,col_var2], method = corr_method)$estimate
      corr_splines_p <- cor.test(splines[,col_var1], splines[,col_var2], method = corr_method)$p.value
      
      corr_resid <- cor.test((data[,col_var1]-splines[,col_var1]),(data[,col_var2]-splines[,col_var2]),method = corr_method)$estimate 
      corr_resid_p <- cor.test((data[,col_var1]-splines[,col_var1]),(data[,col_var2]-splines[,col_var2]),method = corr_method)$p.value
      
      # correlations surrogates (null model)
      corr_surr <- mapply(function(x, y) cor.test(x, y, method = corr_method)$estimate, as.data.frame(surrogates[[(comb_two[1,i])]]), as.data.frame(surrogates[[(comb_two[2,i])]]))
      # correlates column i of the matrices x with column i of matrix y (ncol = nreps = null series)
      # corr_surr is then a vector with 'nreps' correlations
      corr_surr_ses = (corr_obs - mean(corr_surr))/sd(corr_surr)
      
      
      # function for  p-value significance of time-series correlation with respect to surrogate null distribution
      ecdf_fun <- function(x,perc){ecdf(x)(perc)} # estimate quantile of value in a vector x
      
      # binding all variables pairs in one big data.frame
      corr_pairs <- rbind.data.frame(corr_pairs, cbind.data.frame(
        var1 = names(surrogates)[(comb_two[1,i])],
        var2 = names(surrogates)[(comb_two[2,i])],
        # Seasonality of each variables (i.e. correlation of series ssp i with spline of ssp i)
        season_var1 = round(season_var1,5),
        season_var2 = round(season_var2,5),
        # Correlation of series of ssp i with series of ssp j)
        corr_obs = round(corr_obs,5), 
        corr_obs_p = round(corr_obs_p,5), 
        # Correlation of spline of ssp i with series of ssp j)
        corr_splines = round(corr_splines,5),
        corr_splines_p = round(corr_splines_p,5), 
        # Correlation of residuals of splines of ssp i with series of ssp j)
        corr_resid = round(corr_resid,5), 
        corr_resid_p = round(corr_resid_p,5), 
        # Standardized size effect of observed correlation and surrogate correlation (null distribution)
        ses_surr = round(corr_surr_ses,5),
        surr_p = ecdf_fun(corr_surr,corr_obs), # p-value
        
        round(t(corr_surr),5)), # add surrogates to data.frame
        make.row.names = FALSE) 
    }
    
    # corr_pairs[1:50,1:10]
    # str(corr_pairs)
    
    corr_pairs[,1:2] <- data.frame(lapply(corr_pairs[,1:2], as.character), stringsAsFactors=FALSE)
    
    write.csv(corr_pairs, paste("output/traps/",trap_name,"/corr_surrogates_",trap_name,".csv",sep=""), row.names = F)
    return(corr_pairs)
    
  }else{
    
    corr_pairs <- read.csv(paste("output/traps/",trap_name,"/corr_surrogates_",trap_name,".csv",sep=""), header = TRUE, stringsAsFactors = FALSE)
    return (corr_pairs)
  }
  
}

########################################################################

corr_surrogates_boxplots <- function(corr_surrog, trap_name, overwrite){     
  
  if(overwrite == TRUE | !file.exists(paste("output/traps/",trap_name,"/corr_surrogates_boxplots",sep=""))){
    
    if (!file.exists(paste("output/traps/",trap_name,"/corr_surrogates_boxplots",sep=""))){ dir.create(paste("output/traps/",trap_name,"/corr_surrogates_boxplots",sep="")) }
    
    for(i in unique(c(corr_surrog$var1,corr_surrog$var2))){ # i = "G_cal"
      
      corr_subset <- corr_surrog[c(which(corr_surrog$var1 == i), which(corr_surrog$var2 == i)),]
      # corr_subset[1:14,1:10]
      
      # Generating vector with variables names that are being compared with the focus variable i
      corr_subset$group <- c(corr_subset$var2[which(corr_subset$var2 != i)],corr_subset$var1[which(corr_subset$var1 != i)])
      # corr_subset[,c("var1","var2", "group")]
      
      # Plot
      subset_melt <- melt(corr_subset[, (which(names(corr_subset)=="surr_p")+1):(ncol(corr_subset))], id = "group")
      # (which(names(corr_subset)=="corr_resid_p")+1) is the same as which(names(corr_subset)=="V1"), but the name "V1" can be "V1.tau" depending on the correlation test used, so this 'which()+1'is to avoid this problem when plotting
     
      p <- ggplot() + 
        geom_boxplot(data = subset_melt, aes(factor(group), value)) +
        labs(y = "Correlation", x = element_blank()) +
        theme_bw() + ylim(-1.1, 1.1) +
        geom_point(data = corr_subset, aes(x = factor(group), y = corr_obs), color = 'red',size = 2) + # red dot (observed correlation)
        geom_point(data = corr_subset, aes(x = factor(group), y = corr_splines), color = 'blue',shape = 115,size = 2) + # shape: s
        geom_point(data = corr_subset, aes(x = factor(group), y = corr_resid), color = 'blue',shape = 114,size = 2) + # shape: r
        theme(axis.title=element_text(size=18, face="bold"))
      
      if(any(corr_subset$surr_p<0.025)){
        p <- p + geom_point(data = corr_subset[which(corr_subset$surr_p<0.05),], aes(x = factor(group), y = 1.05), color = 'red',shape = 8, size = 3)
      }
      
      if(any(corr_subset$surr_p>0.975)){
        p <- p + geom_point(data = corr_subset[which(corr_subset$surr_p>0.975),], aes(x = factor(group), y = 1.05), color = 'red',shape = 8, size = 3)
      }
      
      p <- p + geom_text(data = corr_subset, aes(x = factor(group), y = -1.05, label=round(ses_surr, 2)), color = 'red', size = 3)
      
      pdf(file =  paste("output/traps/",trap_name,"/corr_surrogates_boxplots/",i,"_",trap_name,"2.pdf",sep=""), width=length(corr_subset$group), height=5, paper = "special")
        print(p)
      dev.off()  
      
      rm(corr_subset)
      rm(subset_melt)
      
    }
  }# if
}

########################################################################

find_neighbours <- function(point, findin, distance) { # vector, data.frame, numeric

  ## point: vector of two numbers (longitude, latitude) 
  ## findin: a matrix of 2 columns (first one is longitude, second is latitude) 
  ## distance: if 0 finds nearest neighbour, if a positive value (in meters) finds neighbours within a circle with the value radius
  
  ## Returns a data.frame: 
  ## "row" = the row number of the neighbour in data.frame
  ## "distance" = the distance between the points
  dist_data <- apply(findin, 1, function(x) distCosine(point, x)) # Matrix
  
  if(distance>0) { # find neighbours within radius of distance
    neighb <- data.frame(row_findin = which(dist_data<=distance), distance = dist_data[which(dist_data<=distance)])
    if(length(neighb[,1])==0) distance = 0 
  }  
  
  if(distance==0) { # find nearest neighbour
    neighb <- data.frame(row_findin = which.min(dist_data), distance = min(dist_data))
  }
  
  return(neighb)
  
}

########################################################################


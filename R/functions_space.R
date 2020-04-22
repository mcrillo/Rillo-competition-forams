
########################################################################
### Functions main_spatial (ForCenS)
########################################################################

########################################################################
# Functions

'%!in%' <- function(x,y)!('%in%'(x,y))

########################################################################

function_calculate_props <- function(data, all_comb){ 
### Description
### Arguments
# data
# all_comb
  
  # New rows with non-available (NA) data which can be filled up with data.
  prop1 <- rep(NA, length(all_comb[1,]))
  prop2 <- rep(NA, length(all_comb[1,]))
  
  species1 <- rep(NA, length(all_comb[1,]))
  species2 <- rep(NA, length(all_comb[1,]))
  
  
  for (i in 1:length(all_comb[1,])) {   
    # columns of d to compare: all_comb[,i]
    # for all names of the combinations---------print(names(d[,all_comb[,i]]))
    d_pair <- data[,(all_comb[,i])]
    # the data of relative abundances for both species in a pair for each i
    cooccur <- ifelse(d_pair[,1] != 0 & d_pair[,2] != 0, yes = 1, no = 0)
    # Comparison of data of relative abundance
    # If both are different than 0 then the outcome will be 1 if one = 0 then outcome will be 0
    
    cooccur <- cooccur[-which(is.na(cooccur))]
    # deletes all the NA's in comparison data.
    
    species1[i] <- names(data)[(all_comb[,i])][1]
    species2[i] <- names(data)[(all_comb[,i])][2]
    
    
    # Calculating proportions
    
    prop1[i] <- sum(cooccur)/length(cooccur)
    # shows the sum of amounts of overlap in relative abundances per sample divided by the amount of data points
    # which in this case can be either 0 or higher than 0 (but NA were excluded)
    
    
    range_ssp1 <- row.names(data[which(d_pair[,1] != 0),])
    range_ssp2 <- row.names(data[which(d_pair[,2] != 0),])
    prop2[i] <- sum(cooccur)/length(unique(c(range_ssp1,range_ssp2)))
    if(is.nan(prop2[i])){prop2[i] <- 0}
    # shows the sum of amounts of overlap in relative abundances per sample
    # divided by the maximum amount of times either one of the species has a value differing from 0.
  }
  
  all_comb_new <- data.frame(species1, species2, prop1, prop2)
  # binding in the extra rows in the matrix which contain the proportional relative abundances.
  
  return(all_comb_new)
  
}

########################################################################

function_name_grid <- function(data){ 
### Description
# This function receives the ForCenS data, and returns the ForCenS data with three extra columns: 'factorLatitude', 'factorLongitude' & 'grid_name'
### Arguments
# data
  
  # transform the latitudes and longitudes to positive numbers
  data$factorLatitude <- (data$Latitude + 90)
  data$factorLongitude <- (data$Longitude + 180)
  # print(data$factorLatitude)
  # print(data$factorLongitude)
  
  # creating names for each longitude interval (numbers 1 to 36)
  long_min <- seq(0, 350 , by=10)
  long_max <- seq(10, 360, by=10)
  long_name <- seq(1:36)
  long_data <- data.frame(long_min, long_max, long_name)
  
  # creating names for each latitude interval (letters)
  lat_min <- seq(0, 170, by=10)
  lat_max <- seq(10, 180, by= 10)
  lat_name <-  letters [1:18]
  lat_data <- data.frame(lat_min, lat_max, lat_name)
  
  # creating new column to add square/grid name
  data[,"grid_name"] <- NA 
  
  for (i in 1:nrow(data)){ # i goes through the rows of the ForCenS data set
    # i = 1
    
    lat_smaller <- which(lat_data$lat_min < data[i,"factorLatitude"]) # latitudes smaller than latitude from data
    lat_larger <- which(lat_data$lat_max >= data[i,"factorLatitude"]) # latitudes larger than latitude from data
    lat_line <- intersect(lat_smaller,lat_larger) # row of the lat_data that you want: intersect between all smaller and all larger gives the one you want
    # Double check:
    # print(paste("ForCenS lat:", data[i,"factorLatitude"])) 
    # print(lat_data[lat_line,])
    
    long_smaller <- which(long_data$long_min < data[i,"factorLongitude"]) # longitudes smaller than longitudes from data
    long_larger <- which(long_data$long_max >= data[i,"factorLongitude"]) # longitudes larger than longitudes from data
    long_line <- intersect(long_smaller,long_larger) # row of long_data that you want
    # Double check:
    # print(paste("ForCenS long:", data[i,"factorLongitude"]))
    # print(long_data[long_line,]) 
    
    data[i,"grid_name"] <- paste(lat_data[lat_line,"lat_name"], long_data[long_line,"long_name"], sep="")
    
  }
  
  # grid names as factors
  data$grid_name <- as.factor(data$grid_name)
  
  return(data)
  
}


########################################################################
# mapplot
##### 

#mapplot <- ggplot(map_data("world"), mapping = aes(x = long, y = lat, group = group)) +
#  geom_polygon(fill = "grey60", colour = "grey60") +
#  theme_bw() + coord_fixed(1.3) +
#  theme(axis.text=element_text(size=16, colour = "black"), 
#        axis.title=element_blank(),
#        axis.ticks = element_blank()) +
#                     labels=c(expression("180"*degree*W), 
#  scale_x_continuous(breaks = seq(-180, +180, 40), position = "top", limits=c(-180, 180), expand = c(0.02, 0.02),
#                              expression("140"*degree*W),
#                              expression("100"*degree*W),
#                              expression("60"*degree*W),
#                              expression("20"*degree*E),
#                              expression("60"*degree*E),
#                              expression("100"*degree*E),
#                              expression("140"*degree*E),
#                              expression("180"*degree*E))
#                             expression("20"*degree*W),
#  ) +
#  scale_y_continuous(breaks = seq(-80,80,20),limits=c(-90,90), expand = c(0, 0),
#                     labels=c(expression("80"*degree*S), 
#                              expression("60"*degree*S),
#                              expression("40"*degree*S),
#                              expression("20"*degree*S),
#                              expression("0"*degree),
#                              expression("20"*degree*N),
#                              expression("40"*degree*N),
#                              expression("60"*degree*N),
#                              expression("80"*degree*N))
#  )
#####
########################################################################


########################################################################
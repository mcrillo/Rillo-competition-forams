# From : https://bhaskarvk.github.io/user2017.geodataviz/notebooks/02-Static-Maps.nb.html#plotting_using_ggplot2

rm(list=ls())

setwd("/Users/marinacostarillo/Google Drive USP/PhD/projects")
setwd("./competition-forams")



###
library(maps)
library(rnaturalearth)
library(sp)
library(scales)

# traps data
traps <- sub(x=list.files(path = "data/traps_time/")[grep(list.files(path = "data/traps_time/"),pattern=".csv")], pattern=".csv", replacement = "")
traps_meta <- read.csv(file="data/traps_metadata.csv", header=TRUE, stringsAsFactors = FALSE)
traps_meta <- traps_meta[which(traps_meta$trap %in% traps),which(colnames(traps_meta) %in% c("trap","sample_no","latN","lonE","woa13_SSTmn", "mean_open_period_days", "length_days","length_years","length_time.series", "from", "to"))] 
traps_coords <- traps_meta[,c(1,3,4)]
names(traps_coords) <- c("name", "lat","long")

# forcens data
forcens <- d
forcens_coords <- unique(forcens[,c("Latitude", "Longitude")])
names(forcens_coords) <- c("lat","long")


# Guides
fdata <- forcens_coords
tdata <- traps_coords
fdata$Dataset <- "Core sample"
tdata$Dataset <- "Sediment trap"
gdata <- rbind(fdata, tdata[,c("lat","long","Dataset")])
ggplot(gdata, aes(x=lat, y=long, group = Dataset, shape = Dataset)) + 
  geom_point(aes(color =  Dataset, fill = Dataset), size = 3, stroke = 1.5) +
  scale_shape_manual(values=c(21,25)) +
  scale_color_manual(values=c('seagreen',alpha('orangered4',0.8))) +
  scale_fill_manual(values=c(alpha('seagreen',0.3),alpha('orangered',0.8)))


# transforming points data into spatial
coordinates(traps_coords) <- ~long+lat 
proj4string(traps_coords) <- '+init=epsg:4326'
traps_coords <- spTransform(traps_coords, CRS("+proj=wintri"))

coordinates(forcens_coords) <- ~long+lat 
proj4string(forcens_coords) <- '+init=epsg:4326'
forcens_coords <- spTransform(forcens_coords, CRS("+proj=wintri"))



# buckley data
buckley <- read.csv(file ="data/buckley.csv", header=TRUE, stringsAsFactors = FALSE) 
buckley_coords <- unique(buckley[,c("Lat.decimal", "Long.decimal")])
names(buckley_coords) <- c("lat","long")
buckley_coords <- buckley_coords[complete.cases(buckley_coords),]
coordinates(buckley_coords) <- ~long+lat 
proj4string(buckley_coords) <- '+init=epsg:4326'
buckley_coords <- spTransform(buckley_coords, CRS("+proj=wintri"))


# world data
world <- rnaturalearth::countries110
# world$name 

# grid lines
grid.lines.mj <- gridlines(world,easts = seq(-180,180,by=30), norths = seq(-90,90,by=30))
grid.lines.mi <- gridlines(world,easts = seq(-165,195,by=15), norths = seq(-90,90,by=15))

# transform all to Winkel Tripel projection
traps_coords <- spTransform(traps_coords, CRS("+proj=wintri"))
forcens_coords <- spTransform(forcens_coords, CRS("+proj=wintri"))

world <- spTransform(world, CRS("+proj=wintri"))
grid.lines.mj <- spTransform(grid.lines.mj,CRS("+proj=wintri"))
grid.lines.mi <- spTransform(grid.lines.mi,CRS("+proj=wintri"))


# plot
# dev.off()

pdf("fig_map_ldg/map_traps.pdf", paper = "special",  width = 15, height = 8.8)
#pdf("fig_map_Buckley.pdf", paper = "special",  width = 15, height = 8.8)

 par(mar = c(0, 0, 0, 0))
 plot(methods::as(world, 'Spatial'), expandBB=c(0,0,0,0))

 plot(grid.lines.mi, col=grey(0.9), add=T)
 plot(grid.lines.mj, col=grey(0.4), add=T)

 plot(world, add=TRUE, border=grey(0.6), col=grey(0.6))
 #plot(buckley_coords, add=TRUE, col='blue', bg=alpha('blue',0.3), pch=21, cex = 2, lwd = 3)
 plot(forcens_coords, add=TRUE, col='seagreen', bg=alpha('seagreen',0.3), pch=21, cex = 1, lwd = 3)
 #plot(traps_coords, add=TRUE, col=alpha('orangered4',0.8), bg=alpha('orangered',0.8), pch=25, cex = 2.5, lwd = 4)
 text(labels(grid.lines.mj, side=2, labelCRS = CRS("+init=epsg:4326")), col = grey(.4), offset=0.5, cex = 2)
dev.off()

# Richness LDG
ldg_raw <- forcens[,c(5,68, 22:62)]
names(ldg_raw)

ldg <- cbind(lat = ldg_raw[,c(1,2)], richness = apply(ldg_raw[,3:ncol(ldg_raw)], 1, FUN = function(x) sum(x!=0, na.rm = T)))
ldg <- as.data.frame(ldg)
names(ldg) <- c("Latitude","SST","Richness")
ldg$lat_round <- 5*round(ldg$Latitude/5)

pdf("fig_map_ldg/ldg.pdf", paper = "special",  width = 4, height = 5)
ggplot(ldg, aes(x=Latitude, y=Richness, group = lat_round)) + 
  geom_boxplot(fill = alpha('seagreen',0.7), outlier.size = 0.5) +
  labs(x = "Latitude", y = "Species richness") +
  theme_classic(base_size = 14) + 
  scale_y_continuous(breaks = seq(0,30,5),limits=c(0,30), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-90,90,30),limits=c(-75,90), expand = c(0, 0), position = "top",
                     labels=c(expression("90"*degree*S), 
                              expression("60"*degree*S),
                              expression("30"*degree*S),
                              expression("0"*degree),
                              expression("30"*degree*N),
                              expression("60"*degree*N),
                              expression("90"*degree*N))) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = grey(0.4)),
        axis.line = element_line(colour = grey(0.4)),
        axis.ticks = element_line(colour = grey(0.4)),
        panel.grid.major = element_line(colour = grey(0.8))) +
  coord_flip()
dev.off()



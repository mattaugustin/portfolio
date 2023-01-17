################################################################################################################################################
#
#    Creation of a regional seismic hazard map using various amplitude smoothing techniques based on Peak Ground Acceleration (PGA) values 
#    seismic hazard data based on Neo-Deterministic Seismic Hazard Analysis (NDSHA) methodology
#    Developed by Matt A.  ,   December 2021
#    
################################################################################################################################################


library(rspatial) ; library(stringr) ;
library(sp)    ; library(sf)
library(rgdal) ; library(dismo)
library(gstat) ; library(fields)
library(dplyr) ; library(readr) ; 
library(reshape2); library(ggplot2) ;


### Load hazard data and explore data
#
regio_pga_df <- read_csv("portfolio/kriging/Kriging_data/regiomap_pga_uk_example.csv") %>% as.data.frame()
regio_pga_df$index <- seq(1, nrow(regio_pga_df), 1)
head( regio_pga_df , 20)


###  Define geocoordinates and reference system for the dataframe before transforming into a Spatial Polygon suited to local coordinate reference system (here British National Grid)
#
coordinates(regio_pga_df) <- ~ longitude + latitude
proj4string(regio_pga_df) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")   
BNG       <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=km +no_defs")
regio_pga_sp <- spTransform(regio_pga_df, BNG)
plot(regio_pga_sp)

regio_pga_df <- regio_pga_sp %>% data.frame() %>% left_join( . , regio_pga_df[, -c(1,2)] %>% data.frame(), 
                                                             by = c("index") ) %>% 
  select( . , c("pga_cms2","pga_g","index","longitude.x", "latitude.x","longitude.y","latitude.y")) %>% 
  rename( . , c("longitude_easting"="longitude.x", "latitude_northing"="latitude.x", 
                "longitude"="longitude.y","latitude"="latitude.y"))



###   Load and Create a SpatialPolygonDataframe of the UK land suited to BNG and create a corresponding raster layer with a higher resolution
###
ukmap_spdf      <-  rgdal::readOGR( dsn =  "portfolio/kriging/Kriging_data/ukmap.shp" , stringsAsFactors = F , ) 
ukmap_spdf_bng  <-  spTransform (ukmap_spdf, BNG)
ukmap_ras <-  raster (ukmap_spdf_bng)
res (ukmap_ras) <-  1     #  The lower the selected number as compared to original value, the higher the resolution ; here resolution of 1 km
ukmap_grid <- as(ukmap_ras, 'SpatialGrid')
plot(ukmap_grid)


###   Load  relevant maps of the UK and meighbouring countries
###
ukmap_BNG_df      <-  read.csv("portfolio/kriging/Kriging_data/UKmap_BNG_df.csv" , as.is = 7 , stringsAsFactors =  TRUE) %>% as.data.frame() 
IsleOfMan_BNG_df  <-  read.csv("portfolio/kriging/Kriging_data/IsleOfMan_df.csv" , as.is = 7 ) %>% as.data.frame()
class(IsleOfMan_BNG_df)

test_ukmap_BNG_df$group <- test_ukmap_BNG_df$group %>% as.character() %>% factor( )




###   create the variogram required for kriging interpolation
### 
regio_pga_sp_stat   <-  gstat( formula = pga_g ~ 1 , locations = regio_pga_sp )
regio_pga_vario     <- variogram ( regio_pga_sp_stat  , width = 10)     ##  consecutive bin size in km for semivariance estimates
plot(regio_pga_vario)

regio_pga_fitvar    <-  gstat::fit.variogram( regio_pga_vario, vgm(c("Exp","Sph","Mat")) )  ## fit.variogram will select the best model to fit the data between spherical, exponential and Matern models
regio_pga_fitvar
plot(regio_pga_vario , regio_pga_fitvar ) 



###   Interpolate hazard data using the variogram-fitting model, based on locations where hazard was defined, and predict across the previously-defined UK map grid ###
###   NOTE : this will take several minutes
startime <- Sys.time()
regio_pga_krig      <- gstat(formula= pga_g ~ 1 , locations = regio_pga_sp , model = regio_pga_fitvar )
hazard_uk_krig_pred <- predict( regio_pga_krig ,  ukmap_grid )
endtime <- Sys.time()
paste("time difference is " , difftime(endtime, startime, units = "mins") , "mins")


### Convert the SpatialGridDataFrame of predictions into a multi-layer raster and mask with UK grid ###
###
regio_pga_brick <- brick(hazard_uk_krig_pred)
regio_pga_brick <- mask(regio_pga_brick, ukmap_spdf_bng)
names(regio_pga_brick) <- c("prediction","variance")
plot(regio_pga_brick )


### Explore range of hazard PGA values to decide on a suitable bin size
###
regio_pga_brick_spdf <- as(  regio_pga_brick , "SpatialPixelsDataFrame")
regio_pga_brick_df   <- regio_pga_brick_spdf %>% data.frame()

min( regio_pga_brick_df$prediction)
max( regio_pga_brick_df$prediction)

ggplot(  regio_pga_brick_df  ) +
  geom_vline( xintercept = 0, size = 2 )  +
  geom_histogram(   aes ( x = prediction )   , binwidth =  0.02  , boundary = 0  , alpha = 0.3   , col = "blue"  )


### Define bins and colour code
###
regio_pga_df$pga_g_bin <- case_when(
  regio_pga_df$pga_g < 0.02 ~ "0.00 - 0.02",
  regio_pga_df$pga_g < 0.04 ~ "0.02 - 0.04",
  regio_pga_df$pga_g < 0.06 ~ "0.04 - 0.06",
  regio_pga_df$pga_g < 0.08 ~ "0.06 - 0.08",
  regio_pga_df$pga_g < 0.10 ~ "0.08 - 0.10",
  regio_pga_df$pga_g < 0.12 ~ "0.10 - 0.12",
  regio_pga_df$pga_g < 0.14 ~ "0.12 - 0.14",
  TRUE ~ "0.14 and above"
)



regio_pga_brick_df$prediction_bin <- case_when(
  regio_pga_brick_df$prediction < 0.02 ~ "0.00 - 0.02",
  regio_pga_brick_df$prediction < 0.04 ~ "0.02 - 0.04",
  regio_pga_brick_df$prediction < 0.06 ~ "0.04 - 0.06",
  regio_pga_brick_df$prediction < 0.08 ~ "0.06 - 0.08",
  regio_pga_brick_df$prediction < 0.10 ~ "0.08 - 0.10",
  regio_pga_brick_df$prediction < 0.12 ~ "0.10 - 0.12",
  regio_pga_brick_df$prediction < 0.14 ~ "0.12 - 0.14",
  TRUE ~ "0.14 and above"
)

new_colours11 <- colorRampPalette(c("#329cd7","#5DB7CB","#87d1bf","#C3E2BE","#FFEEA8","#F6CE6A","#FE9E62","#e76a68"))
chosencolours <- new_colours11
code_colours  <- c("0.00 - 0.02" = new_colours11(7)[1],
                   "0.02 - 0.04" = new_colours11(7)[2],
                   "0.04 - 0.06" = new_colours11(7)[3],
                   "0.06 - 0.08" = new_colours11(7)[4],
                   "0.08 - 0.10" = new_colours11(7)[5],
                   "0.10 - 0.12" = new_colours11(7)[6],
                   "0.12 - 0.14" = new_colours11(7)[7]
) 


###   CUSTOMIZED PLOTTING STYLE DEFINED AT END OF SCRIPT


### Plot the original map distribution and save to an external file for publication
###
ggplot() + 
  portfolio_map_style() +
  coord_equal( xlim = c( -30 ,  720 ) ,  ylim = c( -70 , 1080) ) +
  #
  #geom_polygon(test_ukmap_BNG_df[test_ukmap_BNG_df$FID == "UK", ] , mapping = aes(long, lat, group = group_char), colour = "black", fill = "white", size = .2) +
  geom_polygon(ukmap_BNG_df    , mapping = aes(long, lat, group = group_char), colour = "black"  ,   fill = "white", size = .2) +
  geom_polygon(IsleOfMan_BNG_df                         , mapping = aes(long, lat, group = group), colour = "black", fill = "white" , size = .2) +
  geom_point(regio_pga_df , mapping = aes( x = longitude_easting, y = latitude_northing, fill = pga_g_bin) , col = "white", shape = 21) +
  geom_polygon(ukmap_BNG_df[ukmap_BNG_df$FID == "UK", ] , mapping = aes(long, lat, group = group_char), colour = "black", fill = NA, size = .2) +
  #
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(values = code_colours,
                    name = str_wrap(string = "PGA (g)", width = 15)) +
  #
  guides(fill = guide_legend(title.position = "top", ncol = 3)) +
  ggtitle("Seismic Hazard Map for the United Kingdom (NDSHA)") +
  #
  theme(     plot.title = ggplot2::element_text(family=font_bold,
                                                size=13,
                                                face = "bold",
                                                color="#3e3e3e") ,
             legend.box.margin = margin(3,0,0,0) )


ggsave(  "portfolio/kriging/Kriging_data/regiomap_pga_original.png" , 
          plot = ggplot2::last_plot() , width = 21   , height = 21 , units = "cm" , type = "cairo") 
  


### Plot the kriged map and save to an external file for publication
###
ggplot() + 
  portfolio_map_style() +
  coord_equal( xlim = c( -30 ,  720 ) ,  ylim = c( -70 , 1080) ) +
  #
  geom_tile(regio_pga_brick_df , mapping = aes(x, y, fill = prediction_bin)) +
  # geom_contour(regio_pga_brick_df  ,     # NOTE : add contours if desired
  #              mapping = aes(x, y, z = prediction),  colour = "#444444", size = 0.1, linetype ="dotted"  , breaks = seq(0, 0.12 , 0.02 )) +
  #
  geom_polygon(ukmap_BNG_df[ukmap_BNG_df$FID == "UK", ] , mapping = aes(long, lat, group = group_char), colour = "black", fill = NA      , size = .2) +
  geom_polygon(ukmap_BNG_df[ukmap_BNG_df$FID != "UK", ] , mapping = aes(long, lat, group = group_char), colour = "black", fill = "white ", size = .2) +
  geom_polygon(IsleOfMan_BNG_df                         , mapping = aes(long, lat, group = group), colour = "black", fill = "white" , size = .2) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(values = code_colours,
                    name = str_wrap(string = "PGA (g)", width = 15)) +
  #
  guides(fill = guide_legend(title.position = "top", ncol = 3)) +
  ggtitle("Seismic Hazard Map for the United Kingdom (NDSHA)") +
  #
  theme(     plot.title = ggplot2::element_text(family=font_bold,
                                                size=13,
                                                face = "bold",
                                                color="#3e3e3e") ,
             legend.box.margin = margin(3,0,0,0) )


ggsave(  "portfolio/kriging/Kriging_data/regiomap_pga_kriged.png"  , 
          plot = ggplot2::last_plot() , width = 21   , height = 21 , units = "cm" , type = "cairo") 



### Define preferred style for plotting the map
###
portfolio_map_style <- function() {
  font <- "Be Vietnam Light"
  font_bold  <- "Be Vietnam ExtraBold"
  #
  extrafont::loadfonts(device = "win", quiet = TRUE)
  
  new_colours = colorRampPalette(c("#0f3057", "#00587a", "#008891", "#88d2d7", "#e7e7de", "#ffcb8e", "#e97171", "#cf1b1b", "#900d0d"))
  greens    = colorRampPalette(c("white", new_colours(10)[5:1]))
  YlGreens  = colorRampPalette(new_colours(11)[7:1])
  YlReds    = colorRampPalette(new_colours(11)[5:11])
  
  water2.5  = colorRampPalette(c("#E2F0F6", "#E5F5FD"))(4)[3]
  water4    = colorRampPalette(c(water2.5, "white"))(3)[2]
  
  ggplot2::theme(
    #
    #Text format:
    #
    plot.title = ggplot2::element_blank(),              
    #
    # plot.subtitle = ggplot2::element_text(family=font,
    #                                       size=22,
    #                                       margin=ggplot2::margin(9,0,9,0)),
    plot.caption = ggplot2::element_blank(),
    #
    #  Legend format
    # 
    legend.position = "bottom",
    legend.text.align = 0,
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(family=font_bold,
                                         size=13,
                                         face = "bold",
                                         color="#3e3e3e"),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family=font,
                                        size=11,
                                        color="#3e3e3e"),                
    
    # Axis format
    #
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),                
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    #
    #Grid lines
    # 
    panel.grid = ggplot2::element_blank(),
    panel.border = element_rect(fill = NA) ,
    #Blank background
    #
    panel.background = ggplot2::element_rect(fill = water4, colour = "black"),
    #
    #Strip background 
    strip.background = ggplot2::element_rect(fill="white"),
    strip.text = ggplot2::element_text(size  = 11,  
                                       # hjust = 0, 
                                       face = "bold",
                                       family = font_bold)
  )
}





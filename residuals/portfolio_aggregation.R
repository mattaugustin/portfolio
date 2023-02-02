################################################################################################################################################
#
#    Creation of a regional seismic hazard map using various amplitude smoothing techniques based on Peak Ground Acceleration (PGA) values 
#    seismic hazard data based on Neo-Deterministic Seismic Hazard Analysis (NDSHA) methodology
#    Developed by Matt A.  ,   December 2021
#    
################################################################################################################################################


library(dplyr)   ; library(stringr) ;
library(readr)   ; library(purrr)   ; library(tidyr)
library(reshape2); library(ggplot2) ; library(geosphere) ;

parameter_table      <-  read_csv("portfolio/aggregation/NDSHA_parameters_unique.csv", trim_ws = FALSE)  %>% as.data.frame()
head(  parameter_table , 4)
unique(  parameter_table$directivity ) ; unique(  parameter_table$event ) ; unique(  parameter_table$component ) ; 
unique(  parameter_table$groundmotion_bis ) ; unique(  parameter_table$real )        



#####
#####################################################################
###         AGGREGATE  RESULTS OF SIMULATIONS                     ###
#####################################################################


## Aggregate based on "real" and "directivity"
mean_parameter_table <-   aggregate(  parameter_table[,16]  , by = list( parameter_table$event   , parameter_table$directivity , parameter_table$epidistance ,
                                                                         parameter_table$channel , parameter_table$station     , parameter_table$groundmotion_bis , parameter_table$component ,
                                                                         parameter_table$station_longitude, parameter_table$station_latitude)     , FUN = mean )

names(mean_parameter_table)  <-  c("event","directivity","epidistance","channel","station","groundmotion_bis","component","station_longitude","station_latitude","mean_amplitude")
#
#
median_parameter_table <-    aggregate(  parameter_table[,16] , by = list( parameter_table$event   , parameter_table$directivity , parameter_table$epidistance ,
                                                                         parameter_table$channel , parameter_table$station     , parameter_table$groundmotion_bis , parameter_table$component ,
                                                                         parameter_table$station_longitude, parameter_table$station_latitude)     , FUN = median )
#
names(median_parameter_table) <- c("event","directivity","epidistance","channel","station","groundmotion_bis","component","station_longitude","station_latitude","median_amplitude")

#
parameter_table_residual <- left_join( mean_parameter_table , median_parameter_table )  %>%
                            left_join( . , parameter_table %>% select(., -amplitude,-real, ) %>% unique )


parameter_table_residual <- parameter_table_residual[, c("event","event_longitude","event_latitude","directivity","channel","epidistance",
                                                         "station","station_id","station_longitude","station_latitude","structure","Vs_30","site_type",
                                                         "groundmotion","groundmotion_bis","component","mean_amplitude","median_amplitude")]



### Aggregate based on "real" only
mean_parameter_table <-  aggregate(  parameter_table[,16]  , by = list( parameter_table$event   , parameter_table$epidistance ,
                                                                        parameter_table$channel , parameter_table$station     , parameter_table$groundmotion_bis , parameter_table$component ,
                                                                        parameter_table$station_longitude, parameter_table$station_latitude)     , FUN = mean )

names(mean_parameter_table) <- c("event","epidistance","channel","station","groundmotion_bis","component","station_longitude","station_latitude","mean_amplitude")
#
#
median_parameter_table <-  aggregate(  parameter_table[,16] , by = list( parameter_table$event   , parameter_table$epidistance ,
                                                                         parameter_table$channel , parameter_table$station     , parameter_table$groundmotion_bis , parameter_table$component ,
                                                                         parameter_table$station_longitude, parameter_table$station_latitude)     , FUN = median )
#
names(median_parameter_table)  <- c("event","epidistance","channel","station","groundmotion_bis","component","station_longitude","station_latitude","median_amplitude")
#
parameter_table_residual_directivity <- left_join( mean_parameter_table , median_parameter_table ) %>%
                                        left_join( . , parameter_table %>% select(., -amplitude,-real,-directivity ) %>% unique ) %>%
                                        mutate(. , directivity = "combined" )

parameter_table_residual_directivity <- parameter_table_residual_directivity[ , c("event","event_longitude","event_latitude","directivity","channel","epidistance",
                                                                                  "station","station_id","station_longitude","station_latitude","structure","Vs_30","site_type",
                                                                                  "groundmotion","groundmotion_bis","component","mean_amplitude","median_amplitude") ]


parameter_table_residual_directivity <- rbind ( parameter_table_residual_directivity , parameter_table_residual  ) 



write.csv(parameter_table_residual_directivity , "C:/ndsha_maps/gmpe_1/portfolio/aggregation/results/NDSHA_predictions_aggregated.csv", row.names = FALSE)

# rm ( mean_parameter_table, median_parameter_table , parameter_table_residual_directivity , parameter_table_residual )



#####
###############################################################################
###         CALCULATE RESIDUALS BETWEEN OBSERVATIONS AND SIMULATIONS        ###
###############################################################################


SAC_table   <- read_csv("portfolio/aggregation/SAC_parameters_unique.csv", trim_ws = FALSE)  %>% as.data.frame()
head(  SAC_table , 4) ; names( SAC_table )
unique( SAC_table$event )  ; unique( SAC_table$component ) ;   unique( SAC_table$groundmotion_bis ) ; unique( SAC_table$filtering_bis ) ;
SAC_table   <-  SAC_table %>% filter( . , component %in% c("maximum","geometric mean","resultant") )


NDSHA_agg <-  read_csv("portfolio/aggregation/results/NDSHA_predictions_aggregated.csv", trim_ws = FALSE)  %>% as.data.frame()
unique( NDSHA_agg$event )  ; unique( NDSHA_agg$component ) ;   unique( NDSHA_agg$groundmotion_bis ) ;


#####
###############################################################################
###         VISUAL COMPARISON OBSERVATIONS vs PREDICTIONS                   ###
###############################################################################


NDSHAagg_SAC_long <- left_join( NDSHA_agg , SAC_table[,-c(5,9)] , 
                                  by = c("event","channel","station","station_longitude","station_latitude","Vs_30"="Vs_30_value","groundmotion","groundmotion_bis","component") ) %>%
                      rename(., "observations_amplitude"="amplitude") %>%
                      pivot_longer(  . , cols = c("mean_amplitude", "median_amplitude","observations_amplitude") ,
                                         names_to = "data" , values_to = "amplitude"       )  %>%
                      mutate( data_bis = case_when(.$data == "mean_amplitude"        ~ "predictions" ,    
                                                   .$data == "median_amplitude"       ~ "predictions" ,
                                                   .$data == "observations_amplitude" ~ "observations" ))


## DO NOT FORGET TO SELECT GROUND MOTION, FILTERING LEVEL & EARTHQUAKE RUPTURE DIRECTIVITY LEVEL
ggplot( NDSHAagg_SAC_long[ NDSHAagg_SAC_long$groundmotion_bis == "PGA" &
                           NDSHAagg_SAC_long$directivity == "combined" &
                          !str_detect( NDSHAagg_SAC_long$data , "median") &
                           str_detect( NDSHAagg_SAC_long$filtering_bis , "35 Hz") ,]  )  +
  #
  martina_new_style () +
  theme(panel.grid.minor.x = element_line(colour = "grey80") , 
        panel.grid.minor.y = element_line(colour = "grey90") , 
        legend.box = "vertical" ,
        legend.spacing.y = unit( 0.02 , "mm"  ) ) +
  #
  geom_point(  mapping = aes(x= epidistance , y = amplitude , shape = site_type , size = site_type , col = data_bis ) , alpha = 1 ) +
  #
  annotation_logticks(sides = "lr" ,size = .001 , alpha = 0.2) +
  scale_y_log10( limits = c(0.01 , 100) , breaks = c(0.01 , 0.1 , 1, 10 ,100 ), minor_breaks =  rep(1:9, 21)*(10^rep(-5:5, each = 9))   ) +    #
  #
  scale_x_continuous( limits = c(0, 350 ) , breaks = seq(0,300,50)  ) +
  scale_fill_manual(   values = rep(NA, length(unique(NDSHAagg_SAC_long$station)))  ) +
  scale_colour_manual( values = new_colours11(8)[c(1,6,8)]   )   + 
  scale_shape_manual( values = c( "rock"= 18, "soil" = 21   )  ) +
  scale_size_manual( values = c( "rock"= 2.6  , "soil" = 2.2  )) +
  ylab( bquote('PGA (cm/'~s^2~ ')') ) + xlab("epidistance (km)") +
  facet_wrap( ~ event , labeller = labeller(  event = c("dudley" ="Dudley", "folkestone"="Folkestone","rasen"="Market Rasen","swansea"="Swansea") ) ,  ncol = 2 ) +
  guides(fill = "none", 
         colour = guide_legend(title = "Data" , override.aes = c( alpha=1 , size = 2 )  ),
         shape  = guide_legend(title = "Observed values on"),
         size   = guide_legend(title = "Observed values on")  ) +
  coord_cartesian(  xlim = c(0,350 ) , expand = c(0,0 ))



ggsave("portfolio/aggregation/img/NDSHA_vs_SAC_median.png", width = 16, height = 18, units = "cm", device = "png")



#####
###############################################################################
###         RESIDUAL CALCULATION & PREPARATION                              ###
###############################################################################


residual_table  <-   left_join( NDSHA_agg , SAC_table[,-c(5,9)] , 
                                by = c("event","channel","station","station_longitude","station_latitude","Vs_30"="Vs_30_value","groundmotion","groundmotion_bis","component")   ) %>%
                                relocate(., filtering_bis , .before = amplitude )  %>%  rename(., "observations_amplitude"="amplitude")

residual_table$residual_mean   <- log10(residual_table$mean_amplitude   / residual_table$observations_amplitude )
residual_table$residual_median <- log10(residual_table$median_amplitude / residual_table$observations_amplitude )

residual_table$distance_bin <- case_when(
  residual_table$epidistance <= 50  ~ "[0-50]" ,
  residual_table$epidistance <= 100 ~ "(50-100]" ,
  residual_table$epidistance <= 150 ~ "(100-150]" ,
  residual_table$epidistance <= 200 ~ "(150-200]" ,
  residual_table$epidistance <= 250 ~ "(200-250]" ,
  residual_table$epidistance <= 300 ~ "(250-300]" ,
  residual_table$epidistance <= 350 ~ "(300-350]"
)
residual_table$distance_bin_centre <- case_when(
  residual_table$distance_bin == "[0-50]"     ~ 25  ,
  residual_table$distance_bin == "(50-100]"   ~ 75  ,
  residual_table$distance_bin == "(100-150]"  ~ 125  ,
  residual_table$distance_bin == "(150-200]"  ~ 175  ,
  residual_table$distance_bin == "(200-250]"  ~ 225  ,
  residual_table$distance_bin == "(250-300]"  ~ 275  ,
  residual_table$distance_bin == "(300-350]"  ~ 325
)




## CALCULATE MEAN RESIDUAL FOR ALL DISTANCE BINS, BY SEPARATING BETWEEN ROCK DATA ONLY & ALL GEOLOGIES, 
residual_table_bin  <-  residual_table   %>%    group_by( . , distance_bin , groundmotion_bis , component , directivity , filtering_bis , ) %>% 
                                                summarise(. , "mean_residual_bin" = mean( residual_mean, na.rm = T ) , "median_residual_bin" = mean( residual_median , na.rm = T)   )          %>% 
                                                #
                                                left_join(. , 
                                                              residual_table %>% filter(. , site_type == "rock")  %>% 
                                                              group_by( . , distance_bin , groundmotion_bis , component , directivity , filtering_bis ) %>% 
                                                              summarise(. , "rock_mean_residual_bin" = mean( residual_mean , na.rm = T  ) , "rock_median_residual_bin" = mean( residual_median, na.rm = T))   )  %>%
                                                #
                                                left_join( residual_table , .)
                             

residual_NDSHA_SAC <- 
residual_table_bin  %>%   select(. , event , directivity , channel , station, station_id, station_longitude, station_latitude, groundmotion_bis, component,  mean_amplitude , median_amplitude , filtering_bis ) %>%
                          pivot_longer( ., cols = c("mean_amplitude","median_amplitude") , 
                                           names_to  = "prediction" ,
                                           values_to = "prediction_amplitude" ) %>%
                          mutate( quantile = case_when( .$prediction == "mean_amplitude"   ~ "mean" , 
                                                        .$prediction == "median_amplitude" ~ "median" )) %>% 
                          #
                          left_join( . , 

residual_table_bin  %>%   select(. , event , directivity , channel , station, station_id, station_longitude, station_latitude, groundmotion_bis, component,  residual_mean , residual_median , filtering_bis ) %>%
                          pivot_longer( ., cols = c("residual_mean","residual_median") , 
                          names_to  = "residual" ,
                          values_to = "residual_amplitude" ) %>%
                          mutate( quantile = case_when( .$residual == "residual_mean"   ~ "mean" , 
                                                        .$residual == "residual_median" ~ "median" ))  )  %>%
                          #
                          left_join( . , 
                          #
residual_table_bin  %>%   select(. , event , directivity , channel , station, station_id, station_longitude, station_latitude, groundmotion_bis, component, rock_mean_residual_bin , mean_residual_bin , filtering_bis ) %>%
                          pivot_longer( ., cols = c("rock_mean_residual_bin","mean_residual_bin") , 
                          names_to  = "residual_bin" ,
                          values_to = "residual_bin_amplitude" ) %>%
                          mutate( geology_type = case_when( .$residual_bin == "rock_mean_residual_bin"  ~ "Rock" , 
                                                            .$residual_bin == "mean_residual_bin"       ~ "All geologies" ))  %>% 
                          mutate( quantile = "mean")  %>% 
                          #
                          rbind( . , 
                          #
residual_table_bin  %>%   select(. , event , directivity , channel , station, station_id, station_longitude, station_latitude, groundmotion_bis, component, rock_median_residual_bin , median_residual_bin , filtering_bis ) %>%
                          pivot_longer( ., cols = c("rock_median_residual_bin","median_residual_bin") , 
                          names_to  = "residual_bin" ,
                          values_to = "residual_bin_amplitude" ) %>%
                          mutate( geology_type = case_when( .$residual_bin == "rock_median_residual_bin"  ~ "Rock" , 
                                                            .$residual_bin == "median_residual_bin"       ~ "All geologies" ))  %>% 
                          mutate( quantile = "median")  )  )  %>% 
                          #
                          left_join( residual_table_bin %>% 
                                          select(., event, event_longitude, event_latitude, directivity, channel, epidistance, station, station_id, station_longitude, station_latitude, structure, Vs_30,
                                                    site_type, groundmotion, groundmotion_bis,component, filtering, filtering_bis, observations_amplitude,distance_bin,distance_bin_centre )  %>% unique , .)



write.csv(residual_NDSHA_SAC  , "portfolio/aggregation/results/residuals_bin_NDSHA_SAC.csv", row.names = FALSE )

rm ( residual_table, residual_table_bin , residual_NDSHA_SAC , SAC_table , NDSHA_agg , NDSHAagg_SAC_long  )


#####
###############################################################################
###         RESIDUAL PLOTTING                                               ###
###############################################################################


residual_NDSHA_SAC  <-  read_csv("portfolio/aggregation/results/residuals_bin_NDSHA_SAC.csv", trim_ws = FALSE)  %>% as.data.frame()

residual_NDSHA_SAC$groundmotion_bis <- factor(residual_NDSHA_SAC$groundmotion_bis  ,
                                               levels = c("PGA",  "Sa(T=0.05s)" ,"Sa(T=0.1s)","Sa(T=0.2s)","Sa(T=0.3s)"  ,"Sa(T=0.4s)"  ,"Sa(T=0.5s)",  "Sa(T=0.6s)" , "Sa(T=0.7s)" ,
                                                          "Sa(T=0.8s)",  "Sa(T=0.9s)"  ,"Sa(T=1s)","Sa(T=2s)","Sa(T=3s)","PGV","PGD" ) )


## The table contains results depending on different settings, aka directivity, component, filtering, geologies and quantile.
## Assessing the influence of each of these settings is done by fixing all settings BUT the ones of interest and plotting them


## 1) VISUALIZE INFLUENCE OF GEOLOGIES (aka FIX "quantile","directivity", "component" and "filtering")
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                                   directivity == "combined", 
                                   component == "geometric mean" ,
                                   str_detect(filtering_bis , "35") ) %>%
ggplot() + 
  #
  martina_new_style() +
  #
  geom_hline(yintercept = 0 , size = 0.5, linetype = "longdash") +
  #
  geom_point(mapping = aes(x= epidistance, y= residual_amplitude ,  fill = "Residual data" )  , size = 0.9 , alpha = 0.20  , col = "grey20" ) +
  #
  geom_point(mapping = aes(x= distance_bin_centre , y= residual_bin_amplitude , colour = geology_type, shape = geology_type, size = geology_type  ), alpha =1 ) + # size = 2.3
  #
  ggtitle ( "Predictions VS Observations | 35Hz | Median component | Rock vs All Geologies" ) +
  xlab ("epicentral distance (km)") + ylab (bquote(Residuals~~~epsilon)) +
  coord_cartesian(ylim = c( -2,2 ) , xlim = c(0,300 ) , expand = c(0,0)) +
  guides(fill   = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = ""),
         shape  = guide_legend(order = 2, title = ""),
         size   = guide_legend(order = 2, title = "")) +
  #
  facet_wrap( ~ groundmotion_bis )  +
  #
  theme( plot.title = element_text()   ) +
  #
  scale_colour_manual(values = c("All geologies" = new_colours11(8)[1] , "Rock"= "grey40" )) +
  scale_shape_manual( values = c("All geologies" = 16  , "Rock"= 18    )) +
  scale_size_manual(  values = c("All geologies" = 2.2 , "Rock"= 2.6   ))


ggsave("portfolio/aggregation/img/Allgeo_vs_Rock.png", width = 16, height = 13, units = "cm", device = "png")


## 2) VISUALIZE INFLUENCE OF FILTERING (aka FIX "quantile","directivity", "component" and "geology_type")
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                                   directivity == "combined", 
                                   component == "geometric mean" ,
                                   geology_type == "All geologies"  ) %>%
  ggplot() + 
  #
  martina_new_style() +
  #
  geom_hline(yintercept = 0 , size = 0.5, linetype = "longdash") +
  #
  geom_point(mapping = aes(x= epidistance, y= residual_amplitude ,  fill = "Residual data" )  , size = 0.9 , alpha = 0.20  , col = "grey20" ) +
  #
  geom_point(mapping = aes(x= distance_bin_centre , y= residual_bin_amplitude , colour = filtering_bis, shape = filtering_bis, size = filtering_bis  ), alpha =1 ) + # size = 2.3
  #
  ggtitle ( "Predictions VS Observations | Median component | All Geologies | 10 Hz vs 35 Hz" ) +
  xlab ("epicentral distance (km)") + ylab (bquote(Residuals~~~epsilon)) +
  coord_cartesian(ylim = c( -2,2 ) , xlim = c(0,300 ) , expand = c(0,0)) +
  guides(fill = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = ""),
         shape = guide_legend(order = 2, title = ""),
         size = guide_legend(order = 2, title = "")) +
  #
  facet_wrap( ~ groundmotion_bis )  +
  #
  theme(plot.title = element_text() ) +
  #
  scale_colour_manual(values = c("low frequency filtering (10 Hz)" = new_colours11(8)[1] , "high frequency filtering (35 Hz)" = new_colours11(8)[8] )) +   #  new_colours11(8)[7]
  scale_shape_manual( values = c("low frequency filtering (10 Hz)" = 16  , "high frequency filtering (35 Hz)"  = 18   )) +
  scale_size_manual(  values = c("low frequency filtering (10 Hz)" = 2.2 , "high frequency filtering (35 Hz)"  = 2.6   ))

ggsave("portfolio/aggregation/img/10_vs_35Hz.png", width = 16, height = 13, units = "cm", device = "png")



## 3) VISUALIZE INFLUENCE OF DIRECTIVITY (aka FIX "quantile","filtering", "component" and "geology_type")
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                                   str_detect(filtering_bis , "35"),
                                   component   == "geometric mean"    ,
                                   directivity != "combined"    ,
                                   geology_type == "All geologies"  )  %>%
  ggplot() + 
  #
  martina_new_style() +
  #
  geom_hline(yintercept = 0 , size = 0.5, linetype = "longdash") +
  #
  geom_point(mapping = aes(x= epidistance, y= residual_amplitude ,  fill = "Residual data" )  , size = 0.9 , alpha = 0.10  , col = "grey20" ) +
  #
  geom_point(mapping = aes(x= distance_bin_centre , y= residual_bin_amplitude , colour = directivity , shape = directivity , size = directivity   ), alpha =1 ) +
  #
  ggtitle ( "Predictions VS Observations | Median component | All Geologies | Directivity impact" ) +
  xlab ("epicentral distance (km)") + ylab (bquote(Residuals~~~epsilon)) +
  coord_cartesian(ylim = c( -2,2 ) , xlim = c(0,300 ) , expand = c(0,0)) +
  guides(fill = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = "Directivity:"),
         shape = guide_legend(order = 2 , title = "Directivity:"),
         size = guide_legend(order = 2  , title = "Directivity:")) +
  #
  facet_wrap( ~ groundmotion_bis )  +
  #
  theme(plot.title = element_text() ) +
  #
  scale_colour_manual(values = c("bil_dir000" = new_colours11(8)[1] , "bil_dir090" = new_colours11(8)[6], "bil_dir180" = new_colours11(8)[8] )) +   
  scale_shape_manual( values = c("bil_dir000" = 15 , "bil_dir090" = 17, "bil_dir180" = 19 )) +   
  scale_size_manual( values = c("bil_dir000" = 2.2 , "bil_dir090" = 2.4, "bil_dir180" = 2.2  ))

ggsave("portfolio/aggregation/img/directivity_impact.png", width = 16, height = 13, units = "cm", device = "png")


#####
############################################################################################
###       PLOTTING OF HISTOGRAMS ASSOCIATED TO PREVIOUS RESIDUAL  PLOTS                  ###
############################################################################################


residual_NDSHA_SAC  <-  read_csv("portfolio/aggregation/results/residuals_bin_NDSHA_SAC.csv", trim_ws = FALSE)  %>% as.data.frame()

residual_NDSHA_SAC$groundmotion_bis <- factor(residual_NDSHA_SAC$groundmotion_bis  ,
                                              levels = c("PGA",  "Sa(T=0.05s)" ,"Sa(T=0.1s)","Sa(T=0.2s)","Sa(T=0.3s)"  ,"Sa(T=0.4s)"  ,"Sa(T=0.5s)",  "Sa(T=0.6s)" , "Sa(T=0.7s)" ,
                                                         "Sa(T=0.8s)",  "Sa(T=0.9s)"  ,"Sa(T=1s)","Sa(T=2s)","Sa(T=3s)","PGV","PGD" ) )

names(  residual_NDSHA_SAC)

## The table contains results depending on different settings, aka directivity, component, filtering, geologies and quantile.
## Assessing the influence of each of these settings is done by fixing all settings BUT the ones of interest and plotting them

## 1bis) HISTOGRAM ASSOCIATED INFLUENCE OF GEOLOGIES
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                              directivity == "combined", 
                              component == "geometric mean" ,
                              geology_type  == "All geologies" ,    # this is related to "residual bin average" and needs being fixed to avoid repetitions as there exists only 55 different stations
                              str_detect(filtering_bis , "35") ) %>%

ggplot( ) +
  #
  martina_new_style() +
  #
  geom_hline( yintercept = 0 , linetype = "longdash" , size = 0.5  ) +
  #
  geom_histogram(aes(y = residual_amplitude , fill = site_type ), binwidth = 0.5, boundary = -4.5, alpha = 0.5, position = "identity") +
  scale_x_continuous( breaks = seq(0,40,4) , limits = c(0, 20) , expand = c(0,0) ) +  #, limits = c(0, 26)
  scale_y_continuous( breaks = seq(-4,4,1) , limits = c(-2,2,1), expand = c(0,0) ) +
  facet_wrap(~ groundmotion_bis) +
  guides(shape  = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = "")) +
  #coord_cartesian(ylim = c(-2 , 1.5)   ) + 
  xlab("Residual counts") +  ylab (bquote(Residuals~~~epsilon)) +
  ggtitle ( "Residual Distribution | 35 Hz | Median component | Influence of geologies" ) +
  theme(plot.title = element_text() ) +
  #theme( panel.background = element_rect(fill = "white", colour = "white"),
  #       strip.text.x = element_text(size = 15), 
  #       legend.text = element_text(size=13),
  #       legend.title = element_text(size=16),
  #       panel.grid.major = element_line(colour = "grey70",linetype = "longdash"), 
  #       panel.grid.minor = element_line(colour = "grey90",linetype = "dashed"),
  #       axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #       axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
  #       axis.title.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = 0, face = "plain"),
  #       axis.title.y = element_text(color = "grey20", size = 13, angle = 90, hjust = .5, vjust = .5, face = "plain") )    +
  guides(fill = guide_legend (title = "Geology") )    +
  scale_fill_manual( values = c( "rock"= "grey40"  , "soil" = new_colours11(8)[6]   )  )   # "grey40"                      # red and brown  (rock) coloured
#scale_fill_manual( values = c( "filtered at 35Hz"= "#C62828" , "filtered at 10Hz" = "#2E7D32")  )   # red and green (10Hz) coloured


ggsave("portfolio/aggregation/img/All_geo_histogram1.png", width = 18, height = 13, units = "cm", device = "png")


## 2bis) HISTOGRAM ASSOCIATED INFLUENCE OF FILTERING
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                              directivity == "combined", 
                              component == "geometric mean" ,
                              geology_type == "All geologies"  ) %>% 
  
  ggplot( ) +
  #
  martina_new_style() +
  #
  geom_hline( yintercept = 0 , linetype = "longdash" , size = 0.5  ) +
  #
  geom_histogram(aes(y = residual_amplitude , fill = filtering_bis ), binwidth = 0.5, boundary = -4.5, alpha = 0.5, position = "identity") +
  scale_x_continuous( breaks = seq(0,40,4) , limits = c(0, 28) , expand = c(0,0) ) + 
  scale_y_continuous( breaks = seq(-4,4,1) , limits = c(-2,2,1), expand = c(0,0) ) +
  facet_wrap(~ groundmotion_bis) +
  guides(shape  = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = "")) +
  xlab("Residual counts") +  ylab (bquote(Residuals~~~epsilon)) +
  ggtitle ( "Residual Distribution | All directivities | Median component | Influence of filtering" ) +
  theme(plot.title = element_text() ) +
  #theme( panel.background = element_rect(fill = "white", colour = "white"),
  #       strip.text.x = element_text(size = 15), 
  #       legend.text = element_text(size=13),
  #       legend.title = element_text(size=16),
  #       panel.grid.major = element_line(colour = "grey70",linetype = "longdash"), 
  #       panel.grid.minor = element_line(colour = "grey90",linetype = "dashed"),
  #       axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #       axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
  #       axis.title.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = 0, face = "plain"),
  #       axis.title.y = element_text(color = "grey20", size = 13, angle = 90, hjust = .5, vjust = .5, face = "plain") )    +
  guides(fill = guide_legend (title = "") )    +
  scale_fill_manual( values = c( "low frequency filtering (10 Hz)"= "grey40"  , "high frequency filtering (35 Hz)" = new_colours11(8)[6]   )  )


ggsave("portfolio/aggregation/img/Filtering_histogram.png", width = 18, height = 13, units = "cm", device = "png")


## 3bis) HISTOGRAM ASSOCIATED INFLUENCE OF DIRECTIVITY
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                              directivity != "combined" ,
                              filtering_bis == "high frequency filtering (35 Hz)", 
                              component == "geometric mean" ,
                              geology_type == "All geologies"  ) %>% 
                       mutate( . , directivity_plot = case_when(.$directivity == "bil_dir000" ~ "0 °"   , 
                                                                .$directivity == "bil_dir090" ~ "90 °"  ,
                                                                .$directivity == "bil_dir180" ~ "180 °" )) %>% 
ggplot( ) +
  #
  martina_new_style() +
  #
  geom_hline( yintercept = 0 , linetype = "longdash" , size = 0.5  ) +
  #
  geom_histogram(aes(y = residual_amplitude , fill = directivity_plot ), binwidth = 0.5, boundary = -4.5, alpha = 0.5, position = "identity") +
  scale_x_continuous( breaks = seq(0,80,4) , minor_breaks = seq(0,80,4) , limits = c(0, 28) , expand = c(0,0) ) +
  scale_y_continuous( breaks = seq(-4,4,1)  , limits = c(-2,2,1), expand = c(0,0) ) +
  facet_wrap(~ groundmotion_bis) +
  guides(shape  = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = "")) +
  xlab("Residual counts") +  ylab (bquote(Residuals~~~epsilon)) +
  ggtitle ( "Residual Distribution | 35 Hz | Median component | Influence of directivity" ) +
  theme(plot.title = element_text() , 
        panel.grid.minor.x = element_line(colour = "grey90",linetype = "dashed") ) +
  #
  guides(fill = guide_legend (title = "Directivity") )    +
  scale_fill_manual( values = c( "0 °"= "grey60"  , "90 °" = new_colours11(8)[6] , "180 °" =   "red"  ))  

ggsave("portfolio/aggregation/img/Directivity_histogram1.png", width = 18, height = 13, units = "cm", device = "png")

# This last histogram is not easy to read with 3 coloured directivity variables, 
# so let us try other ways to redraw it, using stacked counts and split bars on directivity



## 3bis) STACKED : HISTOGRAM ASSOCIATED INFLUENCE OF DIRECTIVITY
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                              directivity != "combined" ,
                              filtering_bis == "high frequency filtering (35 Hz)", 
                              component == "geometric mean" ,
                              geology_type == "All geologies"  ) %>% 
  mutate( . , directivity_plot = case_when(.$directivity == "bil_dir000" ~ "0 °"   , 
                                           .$directivity == "bil_dir090" ~ "90 °"  ,
                                           .$directivity == "bil_dir180" ~ "180 °" )) %>% 
  ggplot( ) +
  #
  martina_new_style() +
  #
  geom_hline( yintercept = 0 , linetype = "longdash" , size = 0.5  ) +
  #
  geom_histogram(aes(y = residual_amplitude , fill = directivity_plot ), binwidth = 0.5, boundary = -4.5, alpha = 0.5, position = "stack") +
  scale_x_continuous( breaks = seq(0,80,10) , minor_breaks = seq(0,80,5) , limits = c(0, 75) , expand = c(0,0) ) +
  scale_y_continuous( breaks = seq(-4,4,1) , limits = c(-2,2,1), expand = c(0,0) ) +
  facet_wrap(~ groundmotion_bis) +
  guides(shape  = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = "")) +
  xlab("Residual counts") +  ylab (bquote(Residuals~~~epsilon)) +
  ggtitle ( "Residual Distribution | 35 Hz | Median component | Influence of directivity" ) +
  theme(plot.title = element_text() , 
        panel.grid.minor.x = element_line(colour = "grey90",linetype = "dashed") ) +
  #
  guides(fill = guide_legend (title = "Directivity") )    +
  scale_fill_manual( values = c( "0 °"= "grey60"  , "90 °" = new_colours11(8)[6] , "180 °" =   "red"  ))  

ggsave("portfolio/aggregation/img/Directivity_histogram2.png", width = 18, height = 13, units = "cm", device = "png")


## 3bis) SPLIT BARS : HISTOGRAM ASSOCIATED INFLUENCE OF DIRECTIVITY
residual_NDSHA_SAC %>% filter(. ,  quantile == "median" , 
                              directivity != "combined" ,
                              filtering_bis == "high frequency filtering (35 Hz)", 
                              component == "geometric mean" ,
                              geology_type == "All geologies"  ) %>% 
  mutate( . , directivity_plot = case_when(.$directivity == "bil_dir000" ~ "0 °"   , 
                                           .$directivity == "bil_dir090" ~ "90 °"  ,
                                           .$directivity == "bil_dir180" ~ "180 °" )) %>%     # select(., residual_amplitude) %>% min( ., na.rm =  T)
  mutate(., residual_bin_manual = cut( residual_amplitude , breaks = seq(-2,2,0.5) , right = T ))  %>%   #  select(. , residual_bin_manual ) %>% unlist
  group_by( . , groundmotion_bis , directivity_plot , residual_bin_manual ) %>% 
  summarise(. , "count_residual_bin_manual" = length( residual_bin_manual  )    )       %>%
  filter(. , !is.na(residual_bin_manual)  ) %>%
  #
  ggplot( ) +
  #
  martina_new_style() +
  #
  geom_col(aes(x = count_residual_bin_manual , y = residual_bin_manual  , fill = directivity_plot , col = directivity_plot ), width = 0.6 , alpha = 1 , position = position_dodge2(
    width = 0.5,
    preserve = "single",
    padding = -0.2,
    reverse = F
  )) +
  scale_x_continuous( breaks = seq(0,80,5) , minor_breaks = seq(0,80,5) , limits = c(0, 30) , expand = c(0,0)  ) + 
  scale_y_discrete( drop = FALSE) +
  facet_wrap(~ groundmotion_bis) +
  guides(shape  = guide_legend(order = 1, title = ""),
         colour = guide_legend(order = 2, title = "")) +
  xlab("Residual counts") +  ylab (bquote(Residuals~~~epsilon)) +
  ggtitle ( "Residual Distribution | 35 Hz | Median component | Influence of directivity" ) +
  theme(plot.title = element_text() , 
        panel.grid.minor.x = element_line(colour = "grey90",linetype = "dashed") ) +
  #
  guides(fill = guide_legend (title = "Directivity") ,
         colour = guide_legend ( title = "Directivity") )    +
  scale_fill_manual( values = c( "0 °"= "grey60"  , "90 °" = new_colours11(8)[6] , "180 °" =   new_colours11(8)[8]  )  ,  drop = FALSE ) +
  scale_colour_manual( values = c( "0 °"= "grey60"  , "90 °" = new_colours11(8)[6] , "180 °" =  new_colours11(8)[8]  )  ,  drop = FALSE )


ggsave("portfolio/aggregation/img/Directivity_histogram3.png", width = 18, height = 18, units = "cm", device = "png")



#####
############################################################################################
###         PREPARE DATA FOR PLOTTING  OBSERVATIONS vs DISTRIBUTION OF PREDICTIONS       ###
############################################################################################


# histogram missing

SAC_table <-  read_csv("portfolio/aggregation/SAC_parameters_unique.csv", trim_ws = FALSE)  %>% as.data.frame()
unique( SAC_table$filtering_bis ) ; unique( SAC_table$component ) ;  unique( SAC_table$groundmotion_bis ) ; 

NDSHA_table   <-   read_csv("portfolio/aggregation/NDSHA_parameters_unique.csv", trim_ws = FALSE)  %>% as.data.frame()
unique( NDSHA_table$directivity ) ; unique( NDSHA_table$groundmotion_bis ) ; unique( NDSHA_table$component ) ; 

SAC_table     <-   SAC_table[   SAC_table$component   == "geometric mean" & !str_detect(SAC_table$filtering_bis , "11 Hz")  , ] 
NDSHA_table   <-   NDSHA_table[ NDSHA_table$component == "geometric mean"   , ] 


#### PRepare table mixing obs and simulations distribution for plotting
NDSHA_SAC_table   <-  left_join( NDSHA_table  ,  SAC_table %>% select(., -epidistance,-site_type_binary) , 
                       by = c("event","station_longitude","station_latitude","groundmotion_bis","groundmotion","channel","station","component","Vs_30"="Vs_30_value") ) %>% 
                   rename( ., "prediction_amplitude"="amplitude.x", "observation_amplitude"="amplitude.y") %>%
                   relocate(. , prediction_amplitude, observation_amplitude, .after = filtering_bis)


NDSHA_SAC_table_calc <- group_by( NDSHA_SAC_table  , event , station , groundmotion_bis , filtering_bis ) %>% summarize( . , maxpred = max( prediction_amplitude) , minpred = min(prediction_amplitude) , meanpred = mean(prediction_amplitude ) )
NDSHA_SAC_table      <- left_join( NDSHA_SAC_table , NDSHA_SAC_table_calc )
NDSHA_SAC_table$performance  <- case_when(
  NDSHA_SAC_table$observation_amplitude > NDSHA_SAC_table$maxpred ~ "underestimating observations" ,
  NDSHA_SAC_table$observation_amplitude < NDSHA_SAC_table$minpred ~ "overestimating" ,
 TRUE ~ "fitting" 
)


## SELECT A GROUND MOTION AND A FILTERING LEVEL & PLOT
#

NDSHA_SAC_table[ NDSHA_SAC_table$groundmotion_bis == "PGA" &
                   str_detect( NDSHA_SAC_table$filtering_bis , "35 Hz") ,] %>% View

ggplot( NDSHA_SAC_table[ NDSHA_SAC_table$groundmotion_bis == "PGA" &
                        str_detect( NDSHA_SAC_table$filtering_bis , "35 Hz") ,]  ) +
  #
  martina_new_style () +
  theme(panel.grid.minor.x = element_line(colour = "grey80") , 
        panel.grid.minor.y = element_line(colour = "grey90") , 
        legend.box = "vertical" ,
        legend.spacing.y = unit( 0.02 , "mm"  ) ) +
  #
  ggbeeswarm::geom_quasirandom( mapping = aes(x= epidistance , y = prediction_amplitude , colour = performance ), alpha = .09, width = 5 , size = 0.2 ) +
  #
  geom_point(  mapping = aes(x= epidistance , y = observation_amplitude , shape = site_type , size = site_type  ) , alpha = 1 ) +
  #
  annotation_logticks(sides = "lr" ,size = .001 , alpha = 0.2) +
  scale_y_log10( limits = c(0.01 , 100) , breaks = c(0.01 , 0.1 , 1, 10 ,100 ), minor_breaks =  rep(1:9, 21)*(10^rep(-5:5, each = 9))   ) +    #
  #
  # scale_y_log10( limits = c(0.01, 200) ,  breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)), minor)   +
  scale_x_continuous( limits = c(0, 350 ) , breaks = seq(0,300,50)  ) +
  scale_fill_manual(   values = rep(NA, length(unique(NDSHA_SAC_table$station)))  ) +
  scale_colour_manual( values = new_colours11(8)[c(1,6,8)]   )   + 
  scale_shape_manual( values = c( "rock"= 18, "soil" = 21   )  ) +
  scale_size_manual( values = c( "rock"= 2.6  , "soil" = 2.2  )) +
  ylab( bquote('PGA (cm/'~s^2~ ')') ) + xlab("epidistance (km)") +
  facet_wrap( ~ event , labeller = labeller(  event = c("dudley" ="Dudley", "folkestone"="Folkestone","rasen"="Market Rasen","swansea"="Swansea") ) ,  ncol =1 ) +
  guides(fill = "none", 
         colour = guide_legend(title = "simulations" , override.aes = c(alpha=0.5 , size = 2 )  ),
         shape  = guide_legend(title = "Observed values on"),
         size   = guide_legend(title = "Observed values on")  ) +
  coord_cartesian(  xlim = c(0,350 ) , expand = c(0,0 ))


ggsave("portfolio/aggregation/img/NDSHA_vs_SAC_swarm.png", width = 16, height = 18, units = "cm", device = "png")




rm ( parameter_table_residual_directivity ,parameter_table_residual , station_soil, station_soil_bis , station_soil_ter , median_parameter_table , mean_parameter_table,
     parameter_table, SAC_table , NDSHA_table, NDSHA_SAC_table, residual_NDSHA_SAC  )




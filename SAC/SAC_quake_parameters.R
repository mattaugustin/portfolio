###################################################################################################################
#                                                                                                                 #
#    Processing earthquake time histories (acceleration) and extracting relevant parameters:                      #
#       PGA, spectral accelerations and Arias-intensity related parameters                                        #
#    seismic hazard data from files saved in SAC (Seismic Analysis Code), associated to earthquake observations   #
#    Developed by Matt A.  ,   December 2021                                                                      #
#                                                                                                                 #
###################################################################################################################

###################################
##    LOAD LIBRARIES             ##
###################################

library(dplyr)   ; library(stringr)  ;
library(readr)   ; library(reshape2) ;

library(purrr)    ; library(geosphere)
library(eseis)    ; library(expm)
library(scales)   ; library(lubridate)



## accel1 and accel2 defined in cm/s2
arias_parameters   <- function(step , accel_1 , accel_2){
  #
  time_vec <- data.frame(seq(0,step*(length(accel_1)-1), step))
  arias_df <- cbind (  time_vec , data.frame(accel_1) , data.frame(accel_2))
  names(arias_df) <- c("time","acc_ew","acc_ns")
  #
  arias_ew_g <- c() ; arias_ns_g <- c()  ;arias_ew_m <- c() ; arias_ns_m <- c() ; arias_stats <- data.frame() ;
  #
  arias_ew_g    <- cumsum((arias_df[,2]/980.665 )^2)*pi*step/2/G_const
  arias_ns_g    <- cumsum((arias_df[,3]/980.665 )^2)*pi*step/2/G_const
  arias_ew_m    <- cumsum((arias_df[,2]/100     )^2)*pi*step/2/G_const
  arias_ns_m    <- cumsum((arias_df[,3]/100     )^2)*pi*step/2/G_const      
  #
  arias_df <- cbind( arias_df , data.frame(arias_ew_g) , data.frame(arias_ns_g) , data.frame(arias_ew_m),data.frame(arias_ns_m)  )
  ew_g_5pct  <- 0.05*max( arias_df$arias_ew_g  )
  ew_g_time_5 <-  arias_df[abs(arias_df$arias_ew_g - ew_g_5pct) == min(abs(arias_df$arias_ew_g - ew_g_5pct)),"time"]
  ew_g_95pct  <-  0.95*max( arias_df$arias_ew_g  )
  ew_g_time_95 <-  arias_df[abs(arias_df$arias_ew_g - ew_g_95pct) == min(abs(arias_df$arias_ew_g - ew_g_95pct)),"time"]
  ew_g_time_diff <- ew_g_time_95 - ew_g_time_5
  #
  ns_g_5pct   <-  0.05*max( arias_df$arias_ns_g  )
  ns_g_time_5 <-  arias_df[abs(arias_df$arias_ns_g - ns_g_5pct) == min(abs(arias_df$arias_ns_g - ns_g_5pct)),"time"]
  ns_g_95pct  <-  0.95*max( arias_df$arias_ns_g  )
  ns_g_time_95 <-  arias_df[abs(arias_df$arias_ns_g - ns_g_95pct) == min(abs(arias_df$arias_ns_g - ns_g_95pct)),"time"]
  ns_g_time_diff <- ns_g_time_95 - ns_g_time_5
  #
  ew_m_5pct   <-  0.05*max( arias_df$arias_ew_m  )
  ew_m_time_5 <-  arias_df[abs(arias_df$arias_ew_m - ew_m_5pct) == min(abs(arias_df$arias_ew_m - ew_m_5pct)),"time"]
  ew_m_95pct  <-  0.95*max( arias_df$arias_ew_m  )
  ew_m_time_95 <-  arias_df[abs(arias_df$arias_ew_m - ew_m_95pct) == min(abs(arias_df$arias_ew_m - ew_m_95pct)),"time"]
  ew_m_time_diff <- ew_m_time_95 - ew_m_time_5
  #
  ns_m_5pct   <-  0.05*max( arias_df$arias_ns_m  )
  ns_m_time_5 <-  arias_df[abs(arias_df$arias_ns_m - ns_m_5pct) == min(abs(arias_df$arias_ns_m - ns_m_5pct)),"time"]
  ns_m_95pct  <-  0.95*max( arias_df$arias_ns_m  )
  ns_m_time_95 <-  arias_df[abs(arias_df$arias_ns_m - ns_m_95pct) == min(abs(arias_df$arias_ns_m - ns_m_95pct)),"time"]
  ns_m_time_diff <- ns_m_time_95 - ns_m_time_5
  #
  parameters <- list ( max(arias_ew_g),max(arias_ns_g),max(arias_ew_m), max(arias_ns_m) , ew_m_time_diff, ns_m_time_diff , ew_m_time_5 ,ew_m_time_95 , ns_m_time_5 ,ns_m_time_95   )
  names(parameters) <- c("arias_ew_g","arias_ns_g","arias_ew_m","arias_ns_m","ew_m_time_diff","ns_m_time_diff","ew_m_time_5" ,"ew_m_time_95" , "ns_m_time_5" ,"ns_m_time_95")
  return(parameters)
}


resp_spectra_harsh <- function( acceleration , duration , nStep, per ){ 
  #
  #
  acceleration <- c(0, acceleration)
  # per <- seq(0.01, 10, by = 0.1)
  #
  freq <- 1/per
  nPeriod <- length(per)
  #
  dw <- 2*pi/duration
  w <- seq(0, nStep*dw , by = dw)
  #
  u <- c()
  utime <- c()
  umax <- c()
  vmax <- c()
  amax <- c()
  H <- c()
  #
  Afft <- fft(acceleration)
  #
  k <- 1000
  damp <- 0.05
  #
  for (j in 1:nPeriod) {
    #  
    m <- (per[j]/(2*pi))^2  * k
    c <- 2*damp*(k*m)^(0.5)
    #
    if(schoolmath::is.odd(nStep) == TRUE ){ nStep = nStep - 1  }
    for (l in 1:(nStep/2+1)){
      H[l] <- 1/(-m*w[l]*w[l] + sqrt(as.complex(-1))*c*w[l] + k) 
      H[nStep +3 -l] <- Conj(H[l])
    }
    #
    #
    Qfft <- -m * Afft
    for(l in 1:(nStep+1)){
      u[l]<- H[l] * Qfft[l]
    }
    #
    #
    utime <- Re(signal::ifft(u))
    umax[j] = max (abs(utime))
    vmax[j] = (2*pi/per[j])*umax[j]
    amax[j] = (2*pi/per[j])*vmax[j]
  }
  #
  #
  resultsharsh <- list ( amax)
  names(resultsharsh) <- c("PSA")
  return(resultsharsh)
}


#####
###################################
##    CONSRUCT  FILE  STRUCTURE  ##
###################################


earthquakes_info <- read_csv("portfolio/observations/earthquakes_info.csv", trim_ws = FALSE)  %>% as.data.frame()
head(earthquakes_info)


data_directory   <- "portfolio/observations/"
channel          <- c("Bchan","Hchan")
filtering        <- c(10,25,35)
motion           <- c("accel","vel","disp")
orientation      <- c("HHE","HHN","BHE","BHN")


file_list_station       <- expand.grid(  earthquakes_info$event , channel  )
file_list_station$path  <- paste0(  data_directory , file_list_station$Var1, "/Base_" , file_list_station$Var2 , "/PZ/" )
list_station            <-  list.files( path = file_list_station$path , recursive = F) %>%
                            str_split_fixed(  . , "\\." , n =3 ) %>% as.data.frame %>% select( . , 2  ) %>% unlist() %>% unique() %>% as.character()


#####
#####################################################################
###   READ  POLE-ZERO FILES  FOR ALL EVENTS  AND STATIONS         ###
#####################################################################

## using "apply" method
list_pz <- apply ( file_list_station , 1 , function(x){ if( length(list.files( path = x[3] , recursive = F))>0  ){ 
                               #
                               list.files( path = x[3] , recursive = F) %>% data.frame() %>% 
                               cbind ( . , x %>% t %>% as.data.frame() ) }}   
           ) %>% do.call("rbind", .) %>% rename(. , "PZ"="." ,"event"="Var1","channel"="Var2","stem_path"="path")



##  using "loop for" method
list_pz <- list()
for(i in 1:nrow(file_list_station)){
  if( length( list.files( path = file_list_station$path[i] , recursive = F)) > 0  ){
    list_pz[[i]] <-  list.files( path = file_list_station$path[i] , recursive = F) %>% data.frame() %>% cbind ( . , file_list_station[i,])
    }
}
list_pz <- do.call("rbind", list_pz) %>% rename(. , "PZ"="." ,"event"="Var1","channel"="Var2","stem_path"="path")


head(list_pz)
list_pz$PZ       <-  list_pz$PZ %>% as.character()
list_pz$station  <- str_sub(  list_pz$PZ , 7 , nchar(list_pz$PZ %>% as.character()) - 4 )
list_pz$path_pz  <- paste0( list_pz$stem_path , list_pz$PZ )
list_pz$exp_name <- paste( list_pz$event , list_pz$channel ,list_pz$station , 
                           str_sub( list_pz$PZ , nchar(list_pz$PZ) - 2 ,  nchar(list_pz$PZ)  ) ,sep = "_" )

list_pz$exists <- case_when(  file.exists( list_pz$path_pz) ~ TRUE,   TRUE ~ FALSE )
list_pz <- list_pz[list_pz$exists == TRUE,]


pzdata <- list()
for (i in 1:nrow(list_pz)){
  pzdata[[i]] <- read.csv(list_pz$path_pz[[i]], stringsAsFactors = FALSE)
}
names(pzdata) <- list_pz$exp_name
View(pzdata)



#####
#####################################################################
###   READ  SAC-DATA  FROM FOLDERS                                ###
#####################################################################


file_list_sac       <- expand.grid(  earthquakes_info$event %>% unique , channel , filtering , motion ) 
file_list_sac$path  <- paste0(  data_directory , file_list_sac$Var1, "/", file_list_sac$Var2 , "_" , 
                                file_list_sac$Var3 , "Hz/" , file_list_sac$Var4 )

## using "apply" method
list_sac <- apply ( file_list_sac , 1 , function(x){ if( length(list.files( path = x[5] , recursive = F))>0  ){ 
            #
            list.files( path = x[5] , recursive = F) %>% data.frame() %>% 
            cbind ( . , x %>% t %>% as.data.frame() ) }}   
)           %>%  
            do.call("rbind", .) %>% rename(. , "SAC"="." ,"event"="Var1","channel"="Var2","filtering"="Var3",
                                               "motion"="Var4","stem_path"="path")

list_sac$SAC       <- list_sac$SAC %>% as.character()
list_sac$station   <- str_sub(  list_sac$SAC ,  1 ,  nchar(list_sac$SAC %>% as.character()) - 4 )
list_sac$path_sac  <- paste0( list_sac$stem_path , "/" , list_sac$SAC )
list_sac$exp_name  <- paste0( list_sac$event   , "_"  , list_sac$channel , "_" , list_sac$filtering , "Hz_" , list_sac$motion , "_" , 
                              list_sac$station , "_"  , str_sub( list_sac$SAC , nchar(list_sac$SAC) - 2 ,  nchar(list_sac$SAC)  )  )

list_sac$exists <- case_when(  file.exists( list_sac$path_sac) ~ TRUE,   TRUE ~ FALSE )
list_sac <- list_sac[list_sac$exists == TRUE,]


sacdata <- list()
start_sac <- Sys.time()
for (i in 1:nrow(list_sac)){
  sacdata[[i]] <- read_sac(list_sac$path_sac[[i]])
}
names(sacdata) <- list_sac$exp_name
end_sac <- Sys.time()
paste("time difference is " , difftime(end_sac, start_sac, units = "mins") )


###   Detailed characteristics for each event with SAC-data available stations and coordinates
sac_features <- str_split_fixed( string = names(sacdata) , pattern = "_" , n = 6 ) %>% as.data.frame() %>% 
  rename( . , "event"="V1","channel"="V2","filtering"="V3","motion"="V4","station"="V5","component"="V6") 

earthquakes_info <- left_join(earthquakes_info , sac_features ) %>% relocate( . , c(depth,magnitude, Strike, Dip, Rake) , .before = event )



#####
###########################################################################################
###                     VIEW  SIGNAL FOR A GIVEN SAC FILE                               ###
###########################################################################################

## choose an i index between 1 and length(sacdata) to plot the signals
i <- 1  
main_title_sac <-  paste0( str_split_fixed( names(sacdata)[i] , "_", n=6)[,1] %>% toupper(), "  " ,
                           str_split_fixed( names(sacdata)[i] , pattern = "_" , n =6)[,4] , 
                           "  signal recorded at station  " ,  sacdata[[i]]$meta$station )
plot_signal( data = sacdata[[i]]  , col = "darkred" ,  lwd = 0.5 ,
             main  = main_title_sac  ,
             xlab = "time of event" , ylab = "amplitude"  )
# abline(h=seq(-1,1,0.25),col="grey80" , lty=2, lwd = 0.2)
abline(v=seq(sacdata[[i]]$meta$starttime, 
             sacdata[[i]]$meta$starttime + 60*60, 5*60)  , col="grey80", lty=2 , lwd = 0.2)



#####
###########################################################################################
###  ASSIGN  OCCURRENCE TIMES TO SAC SIGNALS  DUE TO INCOMPLETENESS OF SAC META DATA    ###
###########################################################################################


###  Some SAC files have wrong metadata start times in spite of correct time registration and earthquake-consistent signals on BGS website. 
###  These erroneous metadata  cause "SAC clipping" operation to generate too large SAC files.
###  Hence, find out which SAC files includes inconsistent "start times", aka "start times" more than 1 day before or after the earthquake occurrence date and time. 
df_test_time <- as.POSIXct(x = "2000-01-01 10:00:00" , tz = "UTC" )
df_test_name <- as.character()
for(i in 1:length(sacdata)){
  #
  df_test_time[i]  <-   sacdata[[i]]$meta$starttime %>% as.POSIXlt()
  df_test_name[i]  <-   names(sacdata)[i]
}
df_test   <-   data.frame(  time_sac = df_test_time %>% as.POSIXlt()  , name =  df_test_name )
df_test$event <- str_split_fixed(  df_test$name , fixed("_") , n = 6)[,1]
df_test <- left_join( df_test , earthquakes_info[ , c("event","date_bis") ] %>% unique()  ) %>% rename( . , "date_event"="date_bis")
df_test$time_difference <-  difftime(  df_test$time_sac  ,  df_test$date_event , units = "mins") %>% as.numeric() %>% round ( . , digits = 2)
head(df_test , 4 )
max( abs(df_test$time_difference)  ) < 24*60  ## 1 day converted in minutes


##  For these time-inconsistent SAC files, redefine start-times shortly before earthquake occurrence
for(i in 1:length(sacdata)){
  #
  sac_event <- str_split_fixed(names(sacdata)[i],pattern = "_",n = 6)[,1]
  sacdata[[i]]$meta$starttime    <- case_when(
      abs( difftime(   sacdata[[i]]$meta$starttime , 
                       earthquakes_info[earthquakes_info$event == sac_event, ]$date_bis %>% unique() ,
                       units = "min" ) ) > 24*60  ~  earthquakes_info[earthquakes_info$event == sac_event, ]$date_bis %>% unique() - 60 ,
      TRUE ~ sacdata[[i]]$meta$starttime
    )
}



##   Verify using time difference values between SAC metadata and earthquake occurrence times that NO SAC metadata is containing errors.
df_test_time <- as.POSIXct(x = "2000-01-01 10:00:00" , tz = "UTC" )
df_test_name <- as.character()
for(i in 1:length(sacdata)){
  #
    df_test_time[i]  <-   sacdata[[i]]$meta$starttime %>% as.POSIXlt()
    df_test_name[i]  <-   names(sacdata)[i]
  }
df_test_new <- data.frame(  time_sac = df_test_time %>% as.POSIXlt()  , name =  df_test_name )
df_test_new$event <- str_split_fixed(  df_test_new$name , fixed("_") , n = 6)[,1]
df_test_new <- left_join( df_test_new , earthquakes_info[ , c("event","date_bis") ] %>% unique() ) %>% rename( . , "date_event"="date_bis")
df_test_new$time_difference <-  difftime(  df_test_new$time_sac  ,  df_test_new$date_event , units = "mins") %>% as.numeric() %>% round ( . , digits = 2)
head(df_test_new , 4 )
max( abs(df_test_new$time_difference)  ) < 24*60  ## 1 day in minute




#####
###########################################################################################
###      CLIP SAC SIGNALS TO SMALLER SIZE IN ORDER TO APPLY RESPONSE SPECTRA LATER      ###
###########################################################################################



## clip signals to 15 seconds before and 5 minutes after earthquake occurrence time
sacdata_clip <- list()
for( i in 1:length(sacdata) ){
  #
  sac_event          <- str_split_fixed(names(sacdata)[i],pattern = "_",n = 6)[,1]
  sacdata_clip[[i]]  <- signal_clip( data  = sacdata[[i]] , 
                                     limits = c( earthquakes_info[earthquakes_info$event == sac_event, c("event","date_bis")]$date_bis %>% unique() - 15 , 
                                                 earthquakes_info[earthquakes_info$event == sac_event, c("event","date_bis")]$date_bis %>% unique() + 5*60 ) )
}
names(sacdata_clip) <- names(sacdata)



#####
###########################################################################################
###      PLOT SAC AND CLIPPED SAC SIGNALS  TO  VISUALIZE CLIPPING EFFECT                ###
###########################################################################################

## choose an i index between 1 and length(sacdata) to plot the signal and the corresponding clipped version
i <- 5  
main_title_sac     <-  paste0( str_split_fixed( names(sacdata)[i] , pattern = "_" , n =6)[,1]  %>% toupper(), "  " ,
                               str_split_fixed( names(sacdata)[i] , pattern = "_" , n =6)[,4] , 
                               "  signal recorded at station  "           ,  sacdata[[i]]$meta$station )
main_title_sacclip <-  paste0( str_split_fixed( names(sacdata_clip)[i] , pattern = "_" , n =6)[,1]  %>% toupper(), "  " ,
                               str_split_fixed( names(sacdata_clip)[i] , pattern = "_" , n =6)[,4] , 
                               "  signal clipped & recorded at station  " ,  sacdata_clip[[i]]$meta$station )

par(mfrow=c(1,2))
plot_signal( data = sacdata[[i]]  , col = "black" ,  lwd  = 0.3 ,
             main  = main_title_sac  ,
             xlim = c( sacdata[[i]]$meta$starttime , 
                       sacdata[[i]]$meta$starttime + sacdata[[i]]$meta$dt * sacdata[[i]]$meta$n   ) ,
             xlab = "time of event" , ylab = "amplitude"  )
abline(h=0,col="grey80" , lty=1, lwd = 0.2)
abline(v=seq(sacdata[[i]]$meta$starttime, 
             sacdata[[i]]$meta$starttime + 120*60, 2*60)  , col="grey80", lty=2 , lwd = 0.2 )  ## minor plot lines every 2 min
plot_signal( data = sacdata_clip[[i]]  , col = "darkred" ,  lwd = 0.3 ,
             main  = main_title_sacclip  ,
             xlim = c( sacdata[[i]]$meta$starttime , 
                       sacdata[[i]]$meta$starttime + sacdata[[i]]$meta$dt * sacdata[[i]]$meta$n   ) ,
             xlab = "time of event" , ylab = "amplitude"  )
# abline(h=seq(-1,1,0.25),col="grey80" , lty=2, lwd = 0.2)
abline(v=seq(sacdata[[i]]$meta$starttime, 
             sacdata[[i]]$meta$starttime + 120*60, 2*60)  , col="grey80", lty=2 , lwd = 0.2)  ## minor plot lines every 2 min



#####
###########################################################################################
###    RETRIEVE METADATA FROM POLE ZERO FILES DUE TO INCOMPLETENESS OF SAC META         ###
###########################################################################################


station_name = c() ; station_latitude = c() ; station_longitude = numeric() ; event = c() ; sampling = c()
for (i in 1:length(pzdata)){
  station_name[i] <- str_split_fixed(names(pzdata)[[i]] , pattern = "_" , n = 4)[,3]
  linlat <- grep( "LATITUDE" , pzdata[[i]][,] , fixed = TRUE)
  linlon <- linlat + 1
  station_latitude[i]    <-    str_split_fixed( string =   pzdata[[i]][linlat,1] , pattern = ":" , n =2  )[,2]   %>% trimws( . , which = c("both")) %>% as.numeric()
  station_longitude[i]   <-    str_split_fixed( string =   pzdata[[i]][linlon,1] , pattern = ":" , n =2  )[,2]   %>% trimws( . , which = c("both")) %>% as.numeric()
  event[i]    <-    str_split_fixed( string =   names(pzdata)[[i]] , pattern = "_" , n = 4  )[,1] 
  sampling[i] <-    str_split_fixed( string =   names(pzdata)[[i]] , pattern = "_" , n = 4  )[,2]
}   
station_info <- cbind.data.frame ( station_name , station_latitude , station_longitude , event , sampling  , stringsAsFactors = FALSE ) %>% unique()

head(station_info, 4)



#####
###########################################################################################
#####    ASSOCIATE POLEZERO METADATA TO "CLIPPED"  SAC FILES                        #######
###########################################################################################


##  Metadata in SAC files are sometimes incomplete, hence the need  to assign POLEZERO (PZ) metadata to SAC metadata
for(i in 1:length(sacdata_clip) ){
    #
    eventi  <- str_split_fixed(  names(sacdata_clip)[[i]] , "_" , n = 6 )[,1]
    chani   <- str_split_fixed(  names(sacdata_clip)[[i]] , "_" , n = 6 )[,2]
    stanami <- str_split_fixed(  names(sacdata_clip)[[i]] , "_" , n = 6 )[,5]
    #
    sacdata_clip[[i]]$meta$longitude <- station_info[station_info$station_name == stanami & station_info$event == eventi & station_info$sampling == chani,]$station_longitude
    sacdata_clip[[i]]$meta$latitude  <- station_info[station_info$station_name == stanami & station_info$event == eventi & station_info$sampling == chani,]$station_latitude
}



#####
###########################################################################################
#####    ORGANIZE SAC DATA IN A TREE STRUCTURE AS STACKED LISTS                     #######
###########################################################################################



# create list of acceleration and velocity
stacked_list <-  rep( list( list() ) , levels(earthquakes_info$motion) %>% length  ) %>% setNames( . , levels(earthquakes_info$motion ))  
# fill each item of the list with event names
for (i in 1:length(stacked_list)) {
  stacked_list[[i]] <- rep( list( list() ) , unique(earthquakes_info$event ) %>% length  ) %>% setNames( . , unique(earthquakes_info$event) )
  #
  # fill each event with broadband type (channel)
  for (j in 1:length(stacked_list[[i]])) {
    stacked_list[[i]][[j]] <- rep( list( list() ) , levels(earthquakes_info$channel) %>% length  ) %>% setNames( . , levels(earthquakes_info$channel) )
    #
    # fill each channel with filtering in Hz
    for (k in 1:length(stacked_list[[i]][[j]])) {
      #
      stacked_list[[i]][[j]][[k]] <- rep( list( list() ) , levels(earthquakes_info$filtering ) %>% length  ) %>% 
                                             setNames( . , levels(earthquakes_info$filtering ) )
      #
      # fill each filtering with the data
      for (l in 1:length(stacked_list[[i]][[j]][[k]])) {
        searchname <- paste0( names(stacked_list[[i]])[j], "_", names(stacked_list[[i]][[j]])[k], "_", names(stacked_list[[i]][[j]][[k]])[l],
                             "_", names(stacked_list)[i])
        sac_element <- grep( searchname  , names(sacdata_clip)  , fixed = TRUE)
        stacked_list[[i]][[j]][[k]][[l]] <- sacdata_clip[sac_element]
        
      }
    }
  }
}

# View(stacked_list)


#####
###########################################################################################
#####    CALCULATE EARTHQUAKE PARAMETERS ON SAC FILES POPULATING THE TREE           #######
###########################################################################################


stacked_tables    <-  rep( list( list() ) , levels(earthquakes_info$motion) %>% length  ) %>% setNames( . , levels(earthquakes_info$motion ))     
parameters_table  <-  rep( list( list() ) , levels(earthquakes_info$motion) %>% length  ) %>% setNames( . , levels(earthquakes_info$motion ))       
#
xi <- 0.05 ;                                                                          #  damping factor for spectral acceleration (SA) calculation
sPeriod <- c(0.05 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1 , 2 , 3)  # structural periods at which to calculate SA 

startime <- Sys.time()
#
# loop on motion level
for (i in 1:length(stacked_list)) {
  stacked_tables[[i]]   <- rep( list( list() ) , unique(earthquakes_info$event ) %>% length  ) %>% setNames( . , unique(earthquakes_info$event) )
  parameters_table[[i]] <- rep( list( list() ) , unique(earthquakes_info$event ) %>% length  ) %>% setNames( . , unique(earthquakes_info$event) )
  #
  # loop on event level
  for (j in 1:length(stacked_list[[i]])) {
    stacked_tables[[i]][[j]]   <- rep( list( list() ) , levels(earthquakes_info$channel) %>% length  ) %>% setNames( . , levels(earthquakes_info$channel) )
    parameters_table[[i]][[j]] <- rep( list( list() ) , levels(earthquakes_info$channel) %>% length  ) %>% setNames( . , levels(earthquakes_info$channel) )
    #
    # loop on  CHANNEL (B or H) level
    for (k in 1:length(stacked_list[[i]][[j]])) {
      #
      stacked_tables[[i]][[j]][[k]]   <- rep( list( list() ) , levels(earthquakes_info$filtering ) %>% length  ) %>% setNames( . , levels(earthquakes_info$filtering ) )
      parameters_table[[i]][[j]][[k]] <- rep( list( list() ) , levels(earthquakes_info$filtering ) %>% length  ) %>% setNames( . , levels(earthquakes_info$filtering ) )
      #
      # loop on Hz filtering level
      for (l in 1:length(stacked_list[[i]][[j]][[k]]) ) {
        #
        #
        if( length(stacked_list[[i]][[j]][[k]][[l]]) != 0 ) {
          maxHe <- c() ; maxHn <- c() ; maxGM <- c() ; maxEN <- c()  ; maxRes <- c() 
          station <- c() ; stalon <- c() ; stalat <- c() ; epidistance <- c() ; spacEW <- list() ; spacNS <- list() ;
          Sa_max <- list () ; Sa_gm <- list() ; Sa_res <- list() ; spac_max <- numeric() ; spac_gm <- numeric() ; spac_res <- numeric() ;
          parameters <- list() ; Sa_max <- list () ; Sa_gm <- list() ; Sa_res <- list()  ; 
          arias_ew_g  <- c()  ; arias_ns_g <- c()   ; arias_ew_m <- c()  ; arias_ns_m <- c()   ; ew_std_m <- c() ; ns_std_m <- c() ; 
          ew_m_time_5 <- c() ; ew_m_time_95 <- c() ; ns_m_time_5 <- c() ; ns_m_time_95 <- c() ;
          #
          for (m in 0:(length(stacked_list[[i]][[j]][[k]][[l]])/2-1) ){
            # 
            #
            stacked_list_subset_ew   <- stacked_list[[i]][[j]][[k]][[l]][[2*m+1]]
            stacked_list_subset_ns   <- stacked_list[[i]][[j]][[k]][[l]][[2*m+2]]
            time_ew                  <- seq(  0  , stacked_list_subset_ew$meta$n * stacked_list_subset_ew$meta$dt , stacked_list_subset_ew$m$dt ) %>% max()
            srate                    <- stacked_list_subset_ew$meta$dt
            #
            quake_tag <- earthquakes_info[ earthquakes_info$event == names(stacked_tables[[i]])[j] , ]
            station[m+1] <- str_split_fixed( names(stacked_list[[i]][[j]][[k]][[l]])[2*m+1] , "_", n = 6)[,5]    # "jump" retrieval of station associated to HHN as its same as HHE
            # station[m+1] <- stacked_list[[i]][[j]][[k]][[l]][[2*m+1]]$meta$station                             # station can also be retrieved from SAC meta data but possible mismatch owing to data web origin
            stalat[m+1]  <- stacked_list[[i]][[j]][[k]][[l]][[2*m+1]]$meta$latitude  %>% round( . , digits = 4)
            stalon[m+1]  <- stacked_list[[i]][[j]][[k]][[l]][[2*m+1]]$meta$longitude %>% round( . , digits = 4)
            epidistance[m+1] <- round( distGeo( unlist(quake_tag[,c("longitude" , "latitude")]  %>% unique ) , c(stalon[m+1], stalat[m+1]), a=6378.137, f=1/298.257223563 ), 3)
            #
            # in the following,   convert signal samplitudes from mm/s2 to cm/s2 as per earthquake metric standards
            maxHe[m+1] <- round( max(abs( stacked_list_subset_ew$signal), na.rm = T ) , 3 ) / 10  
            maxHn[m+1] <- round( max(abs( stacked_list_subset_ns$signal), na.rm = T ) , 3 ) / 10
            maxEN[m+1] <- max (maxHe[m+1], maxHn[m+1])
            maxGM[m+1] <- round( max( ( abs(stacked_list_subset_ew$signal) * abs(stacked_list_subset_ns$signal) ) ^.5 ), 3 ) / 10
            maxRes[m+1]<- round( max( ( (stacked_list_subset_ew$signal)^2 + (stacked_list_subset_ns$signal)^2 ) ^.5 ), 3 )   / 10
            #
            #
            if ( names(stacked_tables)[i] == "accel"){
              #
              #  signals need to be input in m/s2 as per definition of "resp_spectra_harsh"
              spacEW <- resp_spectra_harsh( stacked_list_subset_ew$signal / 1000 , time_ew[length(time_ew)] , length(stacked_list_subset_ew$signal), sPeriod )
              spacNS <- resp_spectra_harsh( stacked_list_subset_ns$signal / 1000 , time_ew[length(time_ew)] , length(stacked_list_subset_ns$signal), sPeriod )
              #
              spac_max <- 100 * pmax(spacEW$PSA, spacNS$PSA)
              spac_gm  <- 100 * sqrt( abs(  (100*spacEW$PSA) * (100*spacNS$PSA)  ))
              spac_res <- 100 * sqrt( (100*spacEW$PSA)^2 + (100*spacNS$PSA)^2     )
              Sa_max[[m+1]] <- as.data.frame(t(as.data.frame(spac_max)))
              Sa_gm[[m+1]]  <- as.data.frame(t(as.data.frame(spac_gm )))
              Sa_res[[m+1]] <- as.data.frame(t(as.data.frame(spac_res)))   
              #
              parameters         <- arias_parameters( srate, stacked_list_subset_ew$signal /10 , stacked_list_subset_ns$signal / 10 )
              arias_ew_g[m+1]    <- parameters[[1]]
              arias_ns_g[m+1]    <- parameters[[2]]
              arias_ew_m[m+1]    <- parameters[[3]]
              arias_ns_m[m+1]    <- parameters[[4]]
              ew_std_m[m+1]      <- parameters[[5]]
              ns_std_m[m+1]      <- parameters[[6]]
              ew_m_time_5[m+1]   <- parameters[[7]]
              ew_m_time_95[m+1]  <- parameters[[8]]
              ns_m_time_5[m+1]   <- parameters[[9]]
              ns_m_time_95[m+1]  <- parameters[[10]]
              }
            # 
          }
          if ( names(stacked_tables)[i] == "accel"){
            # create Sa table
            Sa_max_table <- do.call("rbind", Sa_max)
            for (n in 1:length(sPeriod)) {
              names(Sa_max_table)[n] <- paste0("max_Sa_", sPeriod[n], "_s")   }
            Sa_gm_table <- do.call("rbind", Sa_gm)
            for (n in 1:length(sPeriod)) {
              names(Sa_gm_table)[n] <- paste0("GM_Sa_", sPeriod[n], "_s")   }
            Sa_res_table <- do.call("rbind", Sa_res)
            for (n in 1:length(sPeriod)) {
              names(Sa_res_table)[n] <- paste0("Res_Sa_", sPeriod[n], "_s")   }
          }
          #
          # fill them in a table
          #
          stacked_tables[[i]][[j]][[k]][[l]] <-  if ( names(stacked_tables)[i] == "accel"){ cbind.data.frame(  station, epidistance , stalat, stalon , maxHe , maxHn , maxEN , maxGM , maxRes, 
                                                                                                               Sa_max_table , Sa_gm_table , Sa_res_table,
                                                                                                               arias_ew_g , arias_ns_g , arias_ew_m , arias_ns_m , ew_std_m , ns_std_m ,
                                                                                                               ew_m_time_5, ew_m_time_95, ns_m_time_5, ns_m_time_95,
                                                                                            stringsAsFactors = FALSE  )   } else {
                                                                                            cbind.data.frame(  station, epidistance , stalat, stalon , maxHe , maxHn , maxEN , maxGM , maxRes,
                                                                                                      stringsAsFactors = FALSE  )                    
                                                                                 }
          stacked_tables[[i]][[j]][[k]][[l]]$filt   <- names(stacked_tables[[i]][[j]][[k]])[l]
          stacked_tables[[i]][[j]][[k]][[l]]$chan   <- names(stacked_tables[[i]][[j]])[k]
          stacked_tables[[i]][[j]][[k]][[l]]$event  <- names(stacked_tables[[i]])[j]
          stacked_tables[[i]][[j]][[k]][[l]]$motion <- names(stacked_tables)[i]
        }
      }
      parameters_table[[i]][[j]][[k]] <- do.call("rbind", stacked_tables[[i]][[j]][[k]])
    }
    parameters_table[[i]][[j]] <- do.call("rbind", parameters_table[[i]][[j]])
  }
  parameters_table[[i]] <- do.call("rbind", parameters_table[[i]])
}
parameters_table <- left_join(  parameters_table[[1]] , parameters_table[[2]] ,   by = c("station","stalat","stalon", "epidistance", "event" , "filt" , "chan" ) ) %>%
                         left_join( . , parameters_table[[3]] ,by = c("station","stalat","stalon", "epidistance", "event" , "filt" , "chan" ) ) 
endtime <- Sys.time()
paste("time difference is " , difftime(endtime,startime, units="mins") )



#####
###########################################################################################
#####    RE-ORGANIZE FINAL TABLE AND SAVE                                           #######
###########################################################################################


names(parameters_table)  <- c("station"   ,  "epidistance", "station_latitude"    ,  "station_longitude"   ,  
                               "Acc_Max_EW"   ,  "Acc_Max_NS" ,    "Acc_Max"   ,  "Acc_GM_max"   ,  "Acc_Res_max"  ,
                               "Max_Sa_0.05_s" ,  "Max_Sa_0.1_s", "Max_Sa_0.2_s" , "Max_Sa_0.3_s" , "Max_Sa_0.4_s" , "Max_Sa_0.5_s",
                               "Max_Sa_0.6_s",  "Max_Sa_0.7_s", "Max_Sa_0.8_s",  "Max_Sa_0.9_s",  "Max_Sa_1_s",  "Max_Sa_2_s" , "Max_Sa_3_s", 
                               "GM_Sa_0.05_s" ,  "GM_Sa_0.1_s", "GM_Sa_0.2_s" , "GM_Sa_0.3_s" , "GM_Sa_0.4_s" , "GM_Sa_0.5_s", 
                               "GM_Sa_0.6_s",  "GM_Sa_0.7_s", "GM_Sa_0.8_s",  "GM_Sa_0.9_s",  "GM_Sa_1_s",  "GM_Sa_2_s" , "GM_Sa_3_s", 
                               "Res_Sa_0.05_s" ,  "Res_Sa_0.1_s", "Res_Sa_0.2_s" , "Res_Sa_0.3_s" , "Res_Sa_0.4_s" , "Res_Sa_0.5_s", 
                               "Res_Sa_0.6_s",  "Res_Sa_0.7_s", "Res_Sa_0.8_s",  "Res_Sa_0.9_s",  "Res_Sa_1_s",  "Res_Sa_2_s" , "Res_Sa_3_s",
                               "arias_ew_g","arias_ns_g","arias_ew_m","arias_ns_m","arias_ew_std","arias_ns_std", "ew_m_time_5","ew_m_time_95","ns_m_time_5"   ,"ns_m_time_95",
                               "filtering" ,"channel"   , "event"    , "motionacc",
                               "Dis_Max_EW"  , "Dis_Max_NS" , "Dis_Max"    ,"Dis_GM_max"  ,  "Dis_Res_max" ,  "motiondisp",
                               "Vel_Max_EW"  , "Vel_Max_NS" , "Vel_Max"    ,"Vel_GM_max"  ,  "Vel_Res_max" ,  "motionvel" ) 

 
 
parameters_table <- parameters_table[,c("event","channel" ,"filtering", "station" ,"epidistance", "station_latitude","station_longitude",
                             "Dis_Max_EW"  , "Dis_Max_NS" , "Dis_Max",  "Vel_Max_EW"  , "Vel_Max_NS" , "Vel_Max"  ,   "Acc_Max_EW" ,   "Acc_Max_NS" , "Acc_Max"   ,
                             #
                             "Max_Sa_0.05_s" ,  "Max_Sa_0.1_s", "Max_Sa_0.2_s" , "Max_Sa_0.3_s" , "Max_Sa_0.4_s" , "Max_Sa_0.5_s", 
                             "Max_Sa_0.6_s",  "Max_Sa_0.7_s", "Max_Sa_0.8_s",  "Max_Sa_0.9_s",  "Max_Sa_1_s",  "Max_Sa_2_s" , "Max_Sa_3_s",
                             #
                             "Dis_GM_max"  , "Vel_GM_max" , "Acc_GM_max"  ,
                             "GM_Sa_0.05_s" ,  "GM_Sa_0.1_s", "GM_Sa_0.2_s" , "GM_Sa_0.3_s" , "GM_Sa_0.4_s" , "GM_Sa_0.5_s", 
                             "GM_Sa_0.6_s",  "GM_Sa_0.7_s", "GM_Sa_0.8_s",  "GM_Sa_0.9_s",  "GM_Sa_1_s",  "GM_Sa_2_s" , "GM_Sa_3_s",
                             #
                             "Dis_Res_max" ,"Vel_Res_max" ,"Acc_Res_max"  ,
                             "Res_Sa_0.05_s" ,  "Res_Sa_0.1_s", "Res_Sa_0.2_s" , "Res_Sa_0.3_s" , "Res_Sa_0.4_s" , "Res_Sa_0.5_s", 
                             "Res_Sa_0.6_s",  "Res_Sa_0.7_s", "Res_Sa_0.8_s",  "Res_Sa_0.9_s",  "Res_Sa_1_s",  "Res_Sa_2_s" , "Res_Sa_3_s",
                             #
                             "arias_ew_g","arias_ns_g","arias_ew_m","arias_ns_m","arias_ew_std","arias_ns_std", "ew_m_time_5","ew_m_time_95","ns_m_time_5"   ,"ns_m_time_95")   ]


head(parameters_table, 4 )

##  combine original earthquake infos with station information (coordinates...)
earthquakes_info             <- left_join( earthquakes_info, station_info , by = c( "event","station"="station_name","channel"="sampling")) #%>% 
earthquakes_info$epidistance <- round( distGeo( earthquakes_info[,c("longitude","latitude")]   , earthquakes_info[,c("station_longitude","station_latitude")]   , a = 6378.137, f=1/298.257223563 ), 3)
earthquakes_info             <- earthquakes_info %>% relocate( . , component, .before = filtering) %>%  relocate( . , epidistance, .before = station_latitude)

write.csv(parameters_table    , "portfolio/observations/parameters_table_sample.csv"            , row.names = FALSE)
write.csv(earthquakes_info ,    "portfolio/observations/earthquake_stations_characteristics.csv", row.names = FALSE)


# remove size-chunk data from environment
# rm( sacdata_clip_reduced ,  earthquakes_info_bis , earthquakes_info , sacdata_clip , sacdata , parameters , parameters_cum , parameters_table,
#   stacked_tables , stacked_list , df_test , df_test_new , df_test_name , df_test_time , list_sac, list_pz , pzdata )

#####
###########################################################################################
#####     PLOT ENGINEERING PARAMETERS  ACROSS EPIDISTANCE                           #######
###########################################################################################


## Assign soil conditions to stations, re-organize table prior to plotting
parameters_table  <- read_csv("portfolio/observations/parameters_table_sample.csv",  trim_ws = FALSE)  %>% as.data.frame()
station_soil      <- read_csv("portfolio/observations/station_soil_info.csv"      ,  trim_ws = FALSE)  %>% as.data.frame()

head( parameters_table , 4  )
head( station_soil     , 4  )

parameters_table       <-    select( parameters_table , which ( !str_detect( names(parameters_table) , "EW|NS|time_5|time_95"  ) ) ) %>% 
                             left_join( .  , station_soil , by = c("station","station_latitude", "station_longitude")  ) %>% 
                             relocate(., c( Vs_30_value, site_type_binary)  , .before = Dis_Max)  

parameters_table_long  <- parameters_table   %>% melt(.  ,  
                                                     id.vars = c("event","channel","filtering","station","epidistance",
                                                                 "station_latitude" , "station_longitude" ,"Vs_30_value","site_type_binary") , 
                                                     measure.vars = c(  "Dis_Max","Vel_Max","Acc_Max","Max_Sa_0.05_s","Max_Sa_0.1_s","Max_Sa_0.2_s","Max_Sa_0.3_s","Max_Sa_0.4_s","Max_Sa_0.5_s",
                                                                        "Max_Sa_0.6_s","Max_Sa_0.7_s","Max_Sa_0.8_s","Max_Sa_0.9_s","Max_Sa_1_s","Max_Sa_2_s","Max_Sa_3_s",
                                                                        "Dis_GM_max","Vel_GM_max","Acc_GM_max","GM_Sa_0.05_s","GM_Sa_0.1_s","GM_Sa_0.2_s","GM_Sa_0.3_s","GM_Sa_0.4_s",
                                                                        "GM_Sa_0.5_s","GM_Sa_0.6_s","GM_Sa_0.7_s","GM_Sa_0.8_s","GM_Sa_0.9_s","GM_Sa_1_s","GM_Sa_2_s","GM_Sa_3_s",
                                                                        "Dis_Res_max","Vel_Res_max","Acc_Res_max","Res_Sa_0.05_s","Res_Sa_0.1_s","Res_Sa_0.2_s","Res_Sa_0.3_s","Res_Sa_0.4_s",
                                                                        "Res_Sa_0.5_s","Res_Sa_0.6_s","Res_Sa_0.7_s","Res_Sa_0.8_s","Res_Sa_0.9_s","Res_Sa_1_s","Res_Sa_2_s","Res_Sa_3_s",
                                                                        "arias_ew_g","arias_ns_g","arias_ew_m","arias_ns_m","arias_ew_std","arias_ns_std"), 
                                                     variable.name = "groundmotion"  ,value.name = "amplitude"  )


parameters_table_long$groundmotion <- parameters_table_long$groundmotion %>% as.character()
parameters_table_long$component <- case_when(
  str_detect( parameters_table_long$groundmotion , "res|Res" )   ~ "resultant" ,
  str_detect( parameters_table_long$groundmotion , "GM"  )       ~ "geometric mean" ,
  str_detect( parameters_table_long$groundmotion , "max|Max" )   ~ "maximum" ,
  str_detect( parameters_table_long$groundmotion , "ew" )        ~ "East-West" ,
  str_detect( parameters_table_long$groundmotion , "ns" )        ~ "North-South" ,
  TRUE ~ "none"
)

parameters_table_long$groundmotion_bis <- case_when(
  str_detect( parameters_table_long$groundmotion , "Dis" )       ~ "PGD" ,
  str_detect( parameters_table_long$groundmotion , "Vel" )       ~ "PGV" ,
  str_detect( parameters_table_long$groundmotion , "Acc" )       ~ "PGA" ,
  str_detect( parameters_table_long$groundmotion , "std" )       ~ "Arias STD" ,
  str_detect( parameters_table_long$groundmotion , "_g" )        ~ "Arias intensity (g)" ,
  str_detect( parameters_table_long$groundmotion , "ew_m|ns_m" ) ~ "Arias intensity (m)" ,
  #
  str_detect( parameters_table_long$groundmotion , "Sa_0.05")    ~ "Sa(T=0.05 s)",
  str_detect( parameters_table_long$groundmotion , "Sa_0.1" )    ~ "Sa(T=0.1s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.2" )    ~ "Sa(T=0.2s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.3" )    ~ "Sa(T=0.3s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.4" )    ~ "Sa(T=0.4s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.5" )    ~ "Sa(T=0.5s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.6" )    ~ "Sa(T=0.6s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.7" )    ~ "Sa(T=0.7s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.8" )    ~ "Sa(T=0.8s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_0.9" )    ~ "Sa(T=0.9s)" ,
  str_detect( parameters_table_long$groundmotion , "Sa_1" )      ~ "Sa(T=1s)"   ,
  str_detect( parameters_table_long$groundmotion , "Sa_2" )      ~ "Sa(T=2s)"   ,
  str_detect( parameters_table_long$groundmotion , "Sa_3" )      ~ "Sa(T=3s)"   ,
  TRUE ~ "none"
)

parameters_table_long$filtering_bis <- case_when(
  str_detect( parameters_table_long$filtering , "10Hz" )  ~ "low frequency filtering (10 Hz)" ,
  TRUE ~ "high frequency filtering (35 Hz)"
)

# View(parameters_table)

## NB : the custom theme is provided at end of the script

####   PLOT IMPACT OF FILTERING
ggplot() +
  #
  portfolio_parameters_style () +
  #
  geom_point(parameters_table_long %>% filter( . ,  component == "resultant" , groundmotion_bis == "PGA" )  ,
             mapping = aes(x = epidistance, y = amplitude , shape = site_type_binary , size = site_type_binary , col = filtering_bis ) , alpha = 1 ) +
  #
  scale_y_log10( limits = c(0.01 , 50) , breaks = c(0.01 , 0.1 , 1, 10 ,100 ), minor_breaks =  rep(1:9, 21)*(10^rep(-10:10, each = 9)) ) +    #
  annotation_logticks(sides = "lr" ,size = .001 ,alpha = 0.2) +
  # scale_x_log10( limits = c(1, 400 )   , breaks = c(1,10, 100 , 1000 ), minor_breaks =  rep(1:9, 21)*(10^rep(-10:10, each = 9))) +
  scale_x_continuous( limits = c(0, 350 ) , breaks = seq(0,300,50)  ) +
  scale_colour_manual( values = new_colours11(8)[c(1,6,8)]   )        + 
  scale_shape_manual ( values = c( "Rock"= 18, "Soil" = 21   )  )     +
  scale_size_manual  ( values = c( "Rock"= 2.6  , "Soil" = 2.2   )  ) +
  ylab( bquote('PGA (cm/'~s^2~ ')') ) + xlab("epidistance (km)")      +
  facet_wrap( ~ event , labeller = labeller(  event = c("folkestone"="Folkestone","rasen"="Market Rasen","swansea"="Swansea") ) ,  ncol =1 ) +
  guides(fill = "none", 
         colour = guide_legend(title = "Filtering" , override.aes = c(alpha=0.5 , size = 2 )  ),
         shape = guide_legend(title = "Observed values on"),
         size  = guide_legend(title = "Observed values on")  ) +
  coord_cartesian(  xlim = c(0,350 ) , expand = c(0,0 ))      

ggsave("portfolio/observations/ground_motion_PGA.png", type = "cairo", width = 10, height = 10, units = "cm")





###  Customized theme for plotting
portfolio_parameters_style <- function() {
  font <- "Be Vietnam Light"
  font_bold  <- "Be Vietnam ExtraBold"
  
  extrafont::loadfonts(device = "win", quiet = TRUE)
  
  new_colours11 = colorRampPalette(c("#329cd7","#5DB7CB","#87d1bf","#C3E2BE","#FFEEA8","#F6CE6A","#FE9E62","#e76a68"))
  
  ggplot2::theme(
    
    #Text format:
    #This sets the font, size, type and colour of text for the chart's title
    plot.title = ggplot2::element_blank(),
    #This sets the font, size, type and colour of text for the chart's subtitle, as well as setting a margin between the title and the subtitle
    plot.subtitle = ggplot2::element_text(family=font,
                                          size=10,
    ),
    plot.caption = ggplot2::element_blank(),
    #This leaves the caption text element empty, because it is set elsewhere in the finalise plot function
    
    #Legend format
    #This sets the position and alignment of the legend, removes a title and backround for it and sets the requirements for any text within the legend. The legend may often need some more manual tweaking when it comes to its exact position based on the plot coordinates.
    legend.position = "top",
    legend.text.align = 0,
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(family=font_bold,
                                         face = "bold",
                                         size=8,
                                         color="#3e3e3e"),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family=font,
                                        size=7.5,
                                        color="#3e3e3e"),         
    legend.margin = margin(0.1,0,0,0, unit = "cm"),
    legend.box = "vertical" ,
    legend.spacing.y = unit( 0.02 , "mm"  )
    
    #Axis format
    #This sets the text font, size and colour for the axis test, as well as setting the margins and removes lines and ticks. In some cases, axis lines and axis ticks are things we would want to have in the chart - the cookbook shows examples of how to do so.
    axis.title = ggplot2::element_text(family=font_bold,
                                       size=8,
                                       face = "bold",
                                       color="#3e3e3e"),
    axis.title.y = element_text(hjust= .5, vjust = 3, angle = 90),
    axis.text = ggplot2::element_text(family=font,
                                      size=7.5,
                                      color="#3e3e3e"),
    axis.ticks = ggplot2::element_blank(),
        
    #Grid lines
    #This removes all minor gridlines and adds major y gridlines. In many cases you will want to change this to remove y gridlines and add x gridlines. The cookbook shows you examples for doing so
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey30"),          
    #
    panel.grid.minor.y = element_line(colour = "grey80",linetype = "dotted"  ) ,
    panel.grid.major.y = element_line(colour = "grey80",linetype = "dotted"  ) ,
    panel.grid = ggplot2::element_line(linetype = "dotted", size=0.15), # 0.2
    
    
    #Blank background
    #This sets the panel background as blank, removing the standard grey ggplot background colour from the plot
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, colour = "black", size = .2),  
    
    #  Space between faceted wrap figures to accomodate for legend on x axis
    plot.margin = margin(1, 3 , 1 , 3 , "mm")  , 
    
    #Strip background (#This sets the panel background for facet-wrapped plots to white, removing the standard grey ggplot background colour and sets the title size of the facet-wrap title to font size 22)
    strip.background = ggplot2::element_rect(fill="white", colour = NA),
    strip.text = ggplot2::element_text(size  = 8,  
                                       #hjust = 0, 
                                       vjust = 1,
                                       family = font_bold,
                                       face = "bold",
                                       margin = margin(1,0, 1,0, unit = "mm")),
    panel.spacing.x = unit(4, "mm") ,
   
  )
}


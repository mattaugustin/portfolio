#################################################################################################################
#                                                                                                               #
#    Processing earthquake time histories (acceleration) and extracting relevant parameters:         #
#       PGA, spectral accelerations and Arias-intensity related parameters                                 #
#    seismic hazard data based on Neo-Deterministic Seismic Hazard Analysis (NDSHA) methodology                 #
#    Developed by Matt A.  ,   December 2021                                                                    #
#                                                                                                               #
#################################################################################################################

##    LOAD LIBRARIES
##
library(dplyr)   ; library(stringr)  ;
library(readr)   ; library(purrr)   ; library(plyr) ;
library(reshape2); library(ggplot2) ;


##########################################################
##  CREATE  FUNCTIONS REQUIRED FOR LATER OPERATIONS     ##
##########################################################

##  reorganize a pair of dataframes combining EW and NS earthquake acceleration signals
prepare_df <- function(table_a, table_b){
  df_merge <- dplyr::left_join( table_a , table_b , by = c("time") )
  names(df_merge) <- c("time","acceleration_ew","acceleration_ns")
  return(df_merge)
}

  
##  Calculate spectral accelerations from earthquake signals, using as inputs a set of spectral periods and 
##  accelerations signals measured in m/s2 units
##
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

 

##  Calculate parameters related to Arias intensity (Arias intensities, significant time durations) 
##  using as inputs pair of accelerations signals on EW/NS components measured in cm/s2 units and corresponding time step

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


####
###########################################################################
##    CONSTRUCT FOLDER AND FILE ARCHITECTURE FOR SUBSEQUENT READING      ##
###########################################################################


data_directory = "portfolio/signals/"
directivity <-sprintf("%03d",c(0,90,180))
real_nb    <- sprintf("%02d",seq(1,99,1))
station_nb <- sprintf("%02d",seq(1,25,1))
struct_nb  <- sprintf("%04d",seq(1,9,1))


file_list_signal  <- list.files( path="portfolio/signals/bil000/" , recursive = F) 
file_list_signal  <- file_list_signal[str_detect( file_list_signal , "quake_") %>% which()]
head(file_list_signal)


path_signal <- expand.grid( directivity  , file_list_signal , real_nb , struct_nb , station_nb  )  %>% 
      dplyr::rename(. , "directivity" = "Var1", "event_path"="Var2"  , "real"="Var3"  ,  "structure"="Var4" , "stn"="Var5" )

path_signal_bis <- expand.grid( directivity  , file_list_signal , real_nb , struct_nb , station_nb  ) %>% 
      dplyr::rename(. , "directivity"="Var1", "event_path"="Var2"  , "real"="Var3"  ,  "structure"="Var4" , "stn"="Var5" )
# head( path_signal , 3) ; head( path_signal_bis , 3)


#  examples of file path for reading an earthquake signal:
#          signals/"directivity"/"event"/"real_nb"/ukr"structure"f2f2.sew."station_nb".plt
#         "signals/bil000/quake_du/real01/signal/ukr0001f2f2.sew.000001.plt"     if quake signal generated by MS  mode
#         "signals/bil090/quake_mh/real54/signal/ukr0010puf2f2.sns.000018.plt"   if quake signal generated by DWN mode
# 

## quake files generated by MS mode
path_signal$accel_ew     <- paste(  data_directory, "bil", path_signal$directivity , "/" , path_signal$event_path  , "/real" , path_signal$real , "/signal/ukr",
                                    path_signal$structure,"f2f2.sew.0000",path_signal$stn ,".plt",  sep="") 
path_signal$accel_ns     <- paste(  data_directory, "bil", path_signal$directivity , "/"  ,path_signal$event_path  , "/real" , path_signal$real , "/signal/ukr",
                                    path_signal$structure,"f2f2.sns.0000",path_signal$stn ,".plt",  sep="")

## quake files generated by DWN mode
path_signal_bis$accel_ew <- paste(  data_directory, "bil", path_signal_bis$directivity , "/" , path_signal_bis$event_path  , "/real" , path_signal_bis$real , "/signal/ukr",
                                    path_signal_bis$structure , "puf2f2.sew.0000" , path_signal_bis$stn , ".plt",  sep="")
path_signal_bis$accel_ns <- paste(  data_directory, "bil", path_signal_bis$directivity , "/" , path_signal_bis$event_path  , "/real" , path_signal_bis$real , "/signal/ukr",
                                    path_signal_bis$structure , "puf2f2.sns.0000" , path_signal_bis$stn , ".plt",  sep="")



path_signal$exists_acc_ew <- case_when(
  file.exists( path_signal$accel_ew) ~ TRUE ,
  TRUE ~ FALSE )
path_signal$exists_acc_ns <- case_when(
  file.exists( path_signal$accel_ns) ~ TRUE ,
  TRUE ~ FALSE )

path_signal_bis$exists_acc_ew <- case_when(
  file.exists( path_signal_bis$accel_ew) ~ TRUE ,
  TRUE ~ FALSE )
path_signal_bis$exists_acc_ns <- case_when(
  file.exists( path_signal_bis$accel_ns) ~ TRUE ,
  TRUE ~ FALSE )

path_signal     <- path_signal[path_signal$exists_acc_ew == TRUE   , ]
path_signal_bis <- path_signal_bis[path_signal_bis$exists_acc_ew == TRUE   , ]
path_signal     <- rbind ( path_signal , path_signal_bis )
# View(path_signal)


path_signal$event   <- str_split_fixed( path_signal$event_path , "_" , 3)[,2]
path_signal$channel <- str_sub( path_signal$event_path , start = nchar(path_signal$event_path %>% as.character()) - 4  , end = nchar(path_signal$event_path %>% as.character())    )

path_signal$exp_name_ew  <-  paste0(  path_signal$event , "_", path_signal$chan ,  "_dir", path_signal$directivity , "_real" ,  path_signal$real , 
                                     "_stn" , path_signal$stn , "_" , as.numeric(path_signal$struct),"ew") 
path_signal$exp_name_ns  <-  paste0(  path_signal$event , "_", path_signal$chan ,  "_dir", path_signal$directivity , "_real" ,  path_signal$real , 
                                    "_stn" , path_signal$stn , "_" , as.numeric(path_signal$struct),"ns") 
#
path_signal <- path_signal[ , c("directivity" ,  "event_path" ,  "event",    "channel" ,  "real" , "structure"  ,   "stn"  , 
                                "accel_ew" , "accel_ns"   ,   "exists_acc_ew", "exists_acc_ns", "exp_name_ew" ,  "exp_name_ns"  ) ]
#

ndsha_ew_ns <- rbind( path_signal[, c("accel_ew","exp_name_ew")]  %>% dplyr::rename(. , "path_file"="accel_ew"  , "exp_name"="exp_name_ew") ,  
                      path_signal[, c("accel_ns","exp_name_ns")]  %>% dplyr::rename(. , "path_file"="accel_ns"  , "exp_name"="exp_name_ns")) %>%
                      arrange( . , exp_name )


#####
###########################################
##        LOADING SIGNAL FILES           ##
###########################################

signal_ndsha_acc <- list()
startime <- Sys.time()
for(i in 1:nrow(ndsha_ew_ns)){
  signal_ndsha_acc[[i]]        <- as.data.frame(readr::read_table2(  ndsha_ew_ns$path_file[[i]] , col_names = FALSE, skip = 1  ))
  names(signal_ndsha_acc[[i]]) <- c("time","acceleration")
}
names(signal_ndsha_acc) <- ndsha_ew_ns$exp_name
endtime <- Sys.time()
paste("time difference is " , difftime(endtime, startime, units = "mins") , "mins")


##  Join corresponding EW and NS files into one dataframe for reducing size

signal_ndsha_acc_ew <- signal_ndsha_acc[str_ends( names(signal_ndsha_acc) , pattern = "ew"  )  ==TRUE  ]
signal_ndsha_acc_ns <- signal_ndsha_acc[str_ends( names(signal_ndsha_acc) , pattern = "ns"  ) ==TRUE  ]
# rm(signal_ndsha_acc_ns, signal_ndsha_acc_ew) 

signal_ndsha_acc_ewns <- list()
startime <- Sys.time()
signal_ndsha_acc_ewns <- mapply( prepare_df , signal_ndsha_acc_ew , signal_ndsha_acc_ns , SIMPLIFY = FALSE )
endtime <- Sys.time()
paste("time difference is " , difftime(endtime, startime, units = "mins") , "mins")
names( signal_ndsha_acc_ewns) <- str_sub ( names( signal_ndsha_acc_ewns) ,  1 ,   nchar(names(signal_ndsha_acc_ewns)) - 2  )
# View( signal_ndsha_acc_ewns )


#####
######################################################
##        ORGANIZE FILES INTO A TREE STRUCTURE      ##
######################################################


stacked_list_ndsha_acc <-  rep( list( list() ) , levels(path_signal$directivity) %>% length  ) %>% setNames( . , paste0( "dir" , levels(path_signal$directivity) ))  
# fill each directivity with event
for (i in 1:length(stacked_list_ndsha_acc)) {
  stacked_list_ndsha_acc[[i]] <- rep( list( list() ) , unique(path_signal$event ) %>% length  ) %>% setNames( . , unique(path_signal$event) )
  #
  # fill each event with sampling rate (channel)
  for (j in 1:length(stacked_list_ndsha_acc[[i]])) {
    stacked_list_ndsha_acc[[i]][[j]] <- rep( list( list() ) , unique(path_signal$chan ) %>% length  ) %>% setNames( . , unique(path_signal$chan) )
    #
    for (k in 1:length(stacked_list_ndsha_acc[[i]][[j]])) {
      real <- list()
      stacked_list_ndsha_acc[[i]][[j]][[k]]        <- rep( list(real), length(unique(path_signal$real))  )
      names(stacked_list_ndsha_acc[[i]][[j]][[k]]) <- c(paste0("real",unique(path_signal$real)))
      #  
      # fill each filtering with the data
      for (l in 1:length(stacked_list_ndsha_acc[[i]][[j]][[k]])) {
        searchname  <- paste(names(stacked_list_ndsha_acc[[i]])[j], str_sub(names(stacked_list_ndsha_acc[[i]][[j]])[k], 1, 5), 
                            names(stacked_list_ndsha_acc)[i] , names(stacked_list_ndsha_acc[[i]][[j]][[k]])[l]   , sep = "_")
        sac_element <- grep( searchname  , names(signal_ndsha_acc_ewns)  , fixed = TRUE)
        stacked_list_ndsha_acc[[i]][[j]][[k]][[l]] <- signal_ndsha_acc_ewns[sac_element]
      }
    }
  }
}

# View( stacked_list_ndsha_acc )


#####
###################################################################################################
##     EXTRACT AND DERIVE PARAMETERS FROM EARTHQUAKE SIGNALS CLASSIFIED INTO TREE STRUCTURE      ##
###################################################################################################


stacked_ndsha_acc_tables <-  rep( list( list() ) , levels(path_signal$directivity) %>% length  ) %>% setNames( . , paste0( "dir" , levels(path_signal$directivity) ))  
parameter_table          <-  rep( list( list() ) , levels(path_signal$directivity) %>% length  ) %>% setNames( . , paste0( "dir" , levels(path_signal$directivity) ))  
xi <- 0.05     #  damping factor for spectral acceleration
sPeriod <- c(0.05 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1 , 2 , 3)
startime <- Sys.time()
#
for (i in 1:length(stacked_list_ndsha_acc)) {
  #
  # define EVENT at j level
  stacked_ndsha_acc_tables[[i]] <- rep( list( list() ) , unique(path_signal$event ) %>% length  ) %>% setNames( . , unique(path_signal$event) )
  parameter_table[[i]]          <- rep( list( list() ) , unique(path_signal$event ) %>% length  ) %>% setNames( . , unique(path_signal$event) )
  #
  for (j in 1:length(stacked_list_ndsha_acc[[i]])) {
    #
    # define CHANNEL at k level
    stacked_ndsha_acc_tables[[i]][[j]] <-  rep( list( list() ) , unique(path_signal$chan ) %>% length  ) %>% setNames( . , unique(path_signal$chan) )
    parameter_table[[i]][[j]]          <-  rep( list( list() ) , unique(path_signal$chan ) %>% length  ) %>% setNames( . , unique(path_signal$chan) )
    #
    for (k in 1:length(stacked_list_ndsha_acc[[i]][[j]])) {
      #  
      # define REALISATION at l level
      real <- list()
      stacked_ndsha_acc_tables[[i]][[j]][[k]] <- rep( list(real), length(unique(path_signal$real))  )
      parameter_table[[i]][[j]][[k]]          <- rep( list(real), length(unique(path_signal$real))  )
      names(stacked_ndsha_acc_tables[[i]][[j]][[k]]) <- c(paste0("real",unique(path_signal$real)))
      names( parameter_table[[i]][[j]][[k]])         <- c(paste0("real",unique(path_signal$real)))
      #
      for (l in 1:length(stacked_list_ndsha_acc[[i]][[j]][[k]])) {
        #
        #  # verify whether list is empty at REALISATION level before processing
        if( length(stacked_list_ndsha_acc[[i]][[j]][[k]][[l]]) != 0 ) {  
          #
          stacked_list_subset_ew <- c () ; stacked_list_subset_ns <- c() ; srate <- c() ; time_ew <- c() ; station_id <- c() ;
          maxEW <- c() ; maxNS <- c() ; maxAccGM <- c() ; maxAccRes <- c() ; maxAcc <- c() ; 
          spacEW <- list() ; spacNS <- list() ; spac_max <- list () ; spac_gm <- list ()   ; spac_res <- list () ;
          parameters <- list() ; arias_ew_g <- c() ; arias_ns_g <- c() ; arias_ew_m <- c() ; arias_ns_m <- c() ; ew_std_m <- c() ; ns_std_m <- c() ;
          ew_m_time_5 <- c() ;  ew_m_time_95 <- c() ;  ns_m_time_5 <- c() ;  ns_m_time_95 <- c() ;
          Sa_max <- list () ; Sa_gm <- list() ; Sa_res <- list() ; Samax_table <- list() ; Sageom_table <- list() ; Sares_table <- list() ;
          #
          for (m in 1:(length(stacked_list_ndsha_acc[[i]][[j]][[k]][[l]])) ){
            #
            time_ew                  <- stacked_list_ndsha_acc[[i]][[j]][[k]][[l]][[m]][,1]
            stacked_list_subset_ew   <- stacked_list_ndsha_acc[[i]][[j]][[k]][[l]][[m]][,2]
            stacked_list_subset_ns   <- stacked_list_ndsha_acc[[i]][[j]][[k]][[l]][[m]][,3]
            srate                    <- time_ew[2]
            #
            maxEW[m] <- round( max(abs(stacked_list_subset_ew )) , 3 )
            maxNS[m] <- round( max(abs(stacked_list_subset_ns )) , 3 )
            maxAcc[m] <- max (maxEW[m] , maxNS[m] )
            #
            station_id[m] <- str_sub( names(stacked_list_ndsha_acc[[i]][[j]][[k]][[l]])[m] , nchar(names(stacked_list_ndsha_acc[[i]][[j]][[k]][[l]])[m]  ) - 6 ,
                                   nchar(names(stacked_list_ndsha_acc[[i]][[j]][[k]][[l]])[m]  ) - 2   )
            #
            maxAccGM[m]    <- max (   sqrt( abs(   stacked_list_subset_ew  *   stacked_list_subset_ns )  ) , na.rm = TRUE  )  
            maxAccRes[m]   <- max (   sqrt(     (stacked_list_subset_ew)^2 +  (stacked_list_subset_ns)^2 ) , na.rm = TRUE  ) 
            #
            #  signals currently defined in cm/s2 need to be input in m/s2 as per definition of "resp_spectra_harsh"
            spacEW <- resp_spectra_harsh( stacked_list_subset_ew / 100 , time_ew[length(time_ew)] , length(stacked_list_subset_ew), sPeriod )
            spacNS <- resp_spectra_harsh( stacked_list_subset_ns / 100 , time_ew[length(time_ew)] , length(stacked_list_subset_ns), sPeriod )
            #
            spac_max <- 100* pmax( spacEW$PSA , spacNS$PSA )
            spac_gm  <- sqrt(abs(( 100*spacEW$PSA) * (100*spacNS$PSA) )  )
            spac_res <- sqrt( (100*spacEW$PSA)^2 + (100*spacNS$PSA)^2   )  
            Sa_max[[m]] <- data.frame(t(data.frame(spac_max, row.names = NULL)), row.names =NULL)
            Sa_gm[[m]]  <- data.frame(t(data.frame(spac_gm , row.names = NULL)), row.names =NULL)
            Sa_res[[m]] <- data.frame(t(data.frame(spac_res, row.names = NULL)), row.names =NULL)
            #  
            parameters       <- arias_parameters( srate, stacked_list_subset_ew , stacked_list_subset_ns )
            arias_ew_g[m]    <- parameters[[1]]
            arias_ns_g[m]    <- parameters[[2]]
            arias_ew_m[m]    <- parameters[[3]]
            arias_ns_m[m]    <- parameters[[4]]
            ew_std_m[m]      <- parameters[[5]]
            ns_std_m[m]      <- parameters[[6]]
            ew_m_time_5[m]   <- parameters[[7]]
            ew_m_time_95[m]  <- parameters[[8]]
            ns_m_time_5[m]   <- parameters[[9]]
            ns_m_time_95[m]  <- parameters[[10]]
            #
          }
          #  
          Samax_table <- do.call("rbind", Sa_max)
          for (p in 1:length(sPeriod)) {
            names(Samax_table)[p] <- paste0("Max_Sa_", sPeriod[p], "_s")   }
          Sageom_table <- do.call("rbind", Sa_gm)
          for (p in 1:length(sPeriod)) {
            names(Sageom_table)[p] <- paste0("GM_Sa_", sPeriod[p], "_s")   }
          Sares_table <- do.call("rbind", Sa_res)
          for (p in 1:length(sPeriod)) {
            names(Sares_table)[p] <- paste0("Res_Sa_", sPeriod[p], "_s")   }
          #
          #  define and stack data information together for all stations before processing next realisation 
          stacked_ndsha_acc_tables[[i]][[j]][[k]][[l]] <-  cbind.data.frame(  station_id, maxEW , maxNS , maxAcc , maxAccGM , maxAccRes ,
                                                                              Samax_table , Sageom_table , Sares_table , 
                                                                              arias_ew_g , arias_ns_g , arias_ew_m , arias_ns_m , ew_std_m , ns_std_m ,
                                                                              ew_m_time_5, ew_m_time_95, ns_m_time_5, ns_m_time_95, 
                                                                              stringsAsFactors = FALSE )  
          # 
          stacked_ndsha_acc_tables[[i]][[j]][[k]][[l]]$directivity <- names(stacked_ndsha_acc_tables)[i]
          stacked_ndsha_acc_tables[[i]][[j]][[k]][[l]]$event       <- names(stacked_ndsha_acc_tables[[i]])[j]
          stacked_ndsha_acc_tables[[i]][[j]][[k]][[l]]$channel     <- names(stacked_ndsha_acc_tables[[i]][[j]])[k]
          stacked_ndsha_acc_tables[[i]][[j]][[k]][[l]]$real        <- names(stacked_ndsha_acc_tables[[i]][[j]][[k]])[l]
          #
        }    
      }    
      #
      # stack all "l" levels (realisations) together
      parameter_table[[i]][[j]][[k]] <- do.call("rbind", stacked_ndsha_acc_tables[[i]][[j]][[k]])
      #
    }
    #   stack all " k " levels (channel) together   
    parameter_table[[i]][[j]] <- do.call("rbind", parameter_table[[i]][[j]])
  }
  #   stack all " j " levels (event) together 
  parameter_table[[i]] <- do.call("rbind", parameter_table[[i]])
}
#   stack all " i " levels (directivity) together 
parameter_table <- do.call("rbind", parameter_table)
endtime <- Sys.time()
paste("time difference is " , difftime(endtime, startime, units = "mins") , "mins")


#####
########################################
##     REORGANIZE AND SAVE DATA      ###
########################################


names(parameter_table)  <- c( "station_id" ,
                              "Acc_Max_EW"   ,  "Acc_Max_NS" ,    "Acc_Max"   ,  "Acc_GM_max"   ,  "Acc_Res_max"  ,
                              "Max_Sa_0.05_s" ,  "Max_Sa_0.1_s", "Max_Sa_0.2_s" , "Max_Sa_0.3_s" , "Max_Sa_0.4_s" , "Max_Sa_0.5_s",
                              "Max_Sa_0.6_s",  "Max_Sa_0.7_s", "Max_Sa_0.8_s",  "Max_Sa_0.9_s",  "Max_Sa_1_s",  "Max_Sa_2_s" , "Max_Sa_3_s",
                              #
                              "GM_Sa_0.05_s" ,  "GM_Sa_0.1_s", "GM_Sa_0.2_s" , "GM_Sa_0.3_s" , "GM_Sa_0.4_s" , "GM_Sa_0.5_s",
                              "GM_Sa_0.6_s",  "GM_Sa_0.7_s", "GM_Sa_0.8_s",  "GM_Sa_0.9_s",  "GM_Sa_1_s",  "GM_Sa_2_s" , "GM_Sa_3_s",
                              #
                              "Res_Sa_0.05_s" ,  "Res_Sa_0.1_s", "Res_Sa_0.2_s" , "Res_Sa_0.3_s" , "Res_Sa_0.4_s" , "Res_Sa_0.5_s",
                              "Res_Sa_0.6_s",  "Res_Sa_0.7_s", "Res_Sa_0.8_s",  "Res_Sa_0.9_s",  "Res_Sa_1_s",  "Res_Sa_2_s" , "Res_Sa_3_s",
                              #
                              "arias_ew_g","arias_ns_g","arias_ew_m","arias_ns_m","arias_ew_std","arias_ns_std", 
                              "ew_m_time_5","ew_m_time_95","ns_m_time_5"   ,"ns_m_time_95",
                              "directivity"   ,  "event", "channel", "real" 
                              )


parameter_table  <-  parameter_table %>% relocate( . , c(directivity,event,channel,real) , .before = station_id  ) %>%
                                         relocate( . , Acc_GM_max  , .before = GM_Sa_0.05_s  ) %>%
                                         relocate( . , Acc_Res_max , .before = Res_Sa_0.05_s )

write.csv(parameter_table, "portfolio/signals/quake_parameters.csv", row.names = FALSE)



# Clean data environment
rm(  parameter_table ,  stacked_list_ndsha_acc  ,stacked_ndsha_acc_tables, ndsha_ew_ns , signal_ndsha_acc_ew ,
     signal_ndsha_acc_ns ,path_signal , path_signal_bis ,  signal_ndsha_acc_ewns)


#####
#############################################################
##   CONNECTING EARTHQUAKE PARAMETERS TO STATION DATA      ##
#############################################################

## In the following sections, instead of the example dataset, the full dataset including all 99 realisations for each directivity angle and all monitoring stations 
## is connected to station information and re-arranged prior to plotting

parameter_table  <-  read_csv("portfolio/signals/quake_parameters_fullset.csv"    , trim_ws = FALSE)  %>% as.data.frame()
head( parameter_table )
  
station_metadata <-  read_csv("portfolio/signals/station_instru_info.csv" , trim_ws = FALSE)  %>% as.data.frame()


parameter_table  <- left_join( parameter_table , station_metadata  ) 
parameter_table$epidistance  <-  distGeo( parameter_table %>% select(., event_longitude   , event_latitude)   , 
                                          parameter_table %>% select(., station_longitude , station_latitude) ,
                                            a=6378.137, f=1/298.257223563   )  %>% round( . , digits = 3)

parameter_table  <- parameter_table  %>%  relocate(. , c(structure,station_longitude,station_latitude, epidistance,station,Vs_30, site_type), .before = Acc_Max_EW ) %>%
                                          relocate(. , c(event_longitude, event_latitude) , .after = event)


parameter_table_long  <- parameter_table   %>% melt(.  ,  
                                                      id.vars = c("directivity","event","event_longitude", "event_latitude",
                                                                  "channel","real","station","station_id","structure",
                                                                  "epidistance","station_latitude" , "station_longitude" ,
                                                                  "Vs_30","site_type") , 
                                                      measure.vars = c(  "Dis_Max"  ,"Vel_Max","Acc_Max","Max_Sa_0.05_s","Max_Sa_0.1_s","Max_Sa_0.2_s","Max_Sa_0.3_s","Max_Sa_0.4_s",
                                                                         "Max_Sa_0.5_s","Max_Sa_0.6_s","Max_Sa_0.7_s","Max_Sa_0.8_s","Max_Sa_0.9_s","Max_Sa_1_s","Max_Sa_2_s","Max_Sa_3_s",
                                                                         #
                                                                         "Dis_GM_max","Vel_GM_max","Acc_GM_max","GM_Sa_0.05_s","GM_Sa_0.1_s","GM_Sa_0.2_s","GM_Sa_0.3_s","GM_Sa_0.4_s",
                                                                         "GM_Sa_0.5_s","GM_Sa_0.6_s","GM_Sa_0.7_s","GM_Sa_0.8_s","GM_Sa_0.9_s","GM_Sa_1_s","GM_Sa_2_s","GM_Sa_3_s",
                                                                         #
                                                                         "Dis_Res_max","Vel_Res_max","Acc_Res_max","Res_Sa_0.05_s","Res_Sa_0.1_s","Res_Sa_0.2_s","Res_Sa_0.3_s","Res_Sa_0.4_s",
                                                                         "Res_Sa_0.5_s","Res_Sa_0.6_s","Res_Sa_0.7_s","Res_Sa_0.8_s","Res_Sa_0.9_s","Res_Sa_1_s","Res_Sa_2_s","Res_Sa_3_s",
                                                                         "arias_ew_g","arias_ns_g","arias_ew_m","arias_ns_m","arias_ew_std","arias_ns_std"), 
                                                      variable.name = "groundmotion"  ,value.name = "amplitude"  )


parameter_table_long$groundmotion <- parameter_table_long$groundmotion %>% as.character()
parameter_table_long$component <- case_when(
  str_detect( parameter_table_long$groundmotion , "res|Res" )   ~ "resultant" ,
  str_detect( parameter_table_long$groundmotion , "GM"  )       ~ "geometric mean" ,
  str_detect( parameter_table_long$groundmotion , "max|Max" )   ~ "maximum" ,
  str_detect( parameter_table_long$groundmotion , "ew" )        ~ "East-West" ,
  str_detect( parameter_table_long$groundmotion , "ns" )        ~ "North-South" ,
  TRUE ~ "none"
)

parameter_table_long$groundmotion_bis <- case_when(
  str_detect( parameter_table_long$groundmotion , "Dis" )       ~ "PGD" ,
  str_detect( parameter_table_long$groundmotion , "Vel" )       ~ "PGV" ,
  str_detect( parameter_table_long$groundmotion , "Acc" )       ~ "PGA" ,
  str_detect( parameter_table_long$groundmotion , "std" )       ~ "Arias STD" ,
  str_detect( parameter_table_long$groundmotion , "_g" )        ~ "Arias intensity (g)" ,
  str_detect( parameter_table_long$groundmotion , "ew_m|ns_m" ) ~ "Arias intensity (m)" ,
  #
  str_detect( parameter_table_long$groundmotion , "Sa_0.05")    ~ "Sa(T=0.05 s)",
  str_detect( parameter_table_long$groundmotion , "Sa_0.1" )    ~ "Sa(T=0.1s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.2" )    ~ "Sa(T=0.2s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.3" )    ~ "Sa(T=0.3s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.4" )    ~ "Sa(T=0.4s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.5" )    ~ "Sa(T=0.5s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.6" )    ~ "Sa(T=0.6s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.7" )    ~ "Sa(T=0.7s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.8" )    ~ "Sa(T=0.8s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_0.9" )    ~ "Sa(T=0.9s)" ,
  str_detect( parameter_table_long$groundmotion , "Sa_1" )      ~ "Sa(T=1s)"   ,
  str_detect( parameter_table_long$groundmotion , "Sa_2" )      ~ "Sa(T=2s)"   ,
  str_detect( parameter_table_long$groundmotion , "Sa_3" )      ~ "Sa(T=3s)"   ,
  TRUE ~ "none"
)



write.csv(  parameter_table_long %>% filter( . , !str_detect(  groundmotion_bis , "Arias")) , "portfolio/signals/groundmotion_PGX_SA.csv", row.names = FALSE)
write.csv(  parameter_table_long %>% filter( . ,  str_detect(  groundmotion_bis , "Arias")) , "portfolio/signals/groundmotion_Arias.csv", row.names = FALSE)


#####
####################################################
##   PLOTTING SELECTED EARTHQUAKE PARAMETERS      ##
####################################################


parameter_table  <-  read_csv("portfolio/signals/groundmotion_PGX_SA.csv", trim_ws = FALSE)  %>% as.data.frame()
head(parameter_table , 4)


parameter_table$event_format <- recode(parameter_table$event,"dudley"="Dudley","folkestone"="Folkestone","rasen"="Market Rasen","swansea"="Swansea") %>%
  factor(levels = c("Dudley","Folkestone", "Market Rasen" ,"Swansea"))


###  PLOT  GEOLOGY DISTRIBUTION VS EPICENTRAL DISTANCE BY SELECTING THE UNIQUE SET OF SEISMIC MONITORING STATIONS
### 
ggplot(parameter_table %>% select( . , event, event_format, epidistance,site_type, channel) %>% unique ) + 
    #
  geom_quasirandom(aes(epidistance, factor(site_type, levels = c("soil", "rock")), fill = site_type , shape = site_type ), 
                   size = 2.5, alpha = .6, stroke = .1, groupOnX = F, width = .5) + #, shape = 21
  xlab("Distance (km)") +
  scale_fill_manual(values = new_colours11(2), name = "Site type") +
  scale_shape_manual(values = c("rock" = 23, "soil" = 21), name = "Site type") +
  guides() +
  facet_wrap(~ event_format , ncol = 1) +
  martina_new_style() + theme(axis.title.y = element_blank(), legend.margin = margin(0,0,-.5,0,"cm"), panel.spacing = unit(.2, "lines"))


ggsave(filename = "portfolio/signals/events_geology_vs_distance.png", type = "cairo", width = 15, height = 9, units = "cm")




###  PLOT  AMPLITUDE DISTRIBUTION DUE TO DIRECTIVITY ANGLE ACROSS EPICENTRAL DISTANCE BY SELECTING AN EARTHQUAKE
###
ggplot() +
  #
  portfolio_parameters_style () +
  #
  ggbeeswarm::geom_quasirandom(  parameter_table[parameter_table$groundmotion_bis == "PGA" & parameter_table$component == "resultant" & parameter_table$event == "swansea" ,] ,
                                 mapping = aes(x= epidistance , y = amplitude, col = directivity  ),  #
                                alpha = 0.8 , width = 5 , size = 0.2  ) +
  #
  annotation_logticks(sides = "lr" ,size = .001 ,alpha = 0.2) +
  #
  scale_y_log10( limits = c(0.01, 200)    ,  breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)) ,
                                             minor_breaks =  rep(1:9, 21)*(10^rep(-10:10, each = 9))  )   +
  scale_x_continuous( limits = c(0, 350 ) , breaks = seq(0,300,50)  ) +
  scale_fill_manual(   values = rep(NA, length(unique(parameter_table$station)))  ) +
  scale_colour_manual( values = new_colours11(8)[c(1,6,8)] ,
                       labels = c("bil_dir000"="0°","bil_dir090"="90°","bil_dir180"="180°") ) +
  #
  ylab( bquote('PGA (cm/'~s^2~ ')') ) + xlab("epidistance (km)") +  ggtitle("Amplitude distribution due to directivity angle (Swansea earthquake)") +
  # facet_wrap( ~ event_format , labeller = labeller(  event = c("dudley" ="Dudley", "folkestone"="Folkestone","mkrasen"="Market Rasen","swansea"="Swansea") ) ,  ncol =1 ) +
  facet_wrap( ~ directivity , ncol = 3 , labeller = labeller ( directivity = c("bil_dir000"="","bil_dir090"="","bil_dir180"="") ) ) +
  guides(fill = "none", 
         colour = guide_legend(title = "Directivity" , override.aes = c(alpha=0.5 , size = 3 )  )) +
  coord_cartesian(  xlim = c(0,350 ) , expand = c(0,0 ))


ggsave(filename = "portfolio/signals/PGA_amp_directivity_distr.png", type = "cairo", width = 25, height = 14, units = "cm")


## CLEAN DATA ENVIRONMENT                  
rm(  parameter_table ,  stacked_list_ndsha_acc  ,stacked_ndsha_acc_tables, ndsha_ew_ns , signal_ndsha_acc_ew ,
     signal_ndsha_acc_ns ,path_signal , path_signal_bis ,  signal_ndsha_acc_ewns)

                
                

 ###   CUSTOM THEME FOR PLOTTING        
 ###                
portfolio_parameters_style <- function() {
  font <- "Be Vietnam Light"
  font_bold  <- "Be Vietnam ExtraBold"
  
  extrafont::loadfonts(device = "win", quiet = TRUE)
  
  new_colours11 = colorRampPalette(c("#329cd7","#5DB7CB","#87d1bf","#C3E2BE","#FFEEA8","#F6CE6A","#FE9E62","#e76a68"))
  
  ggplot2::theme(
    
    #Text format:
    #This sets the font, size, type and colour of text for the chart's title
    plot.title = ggplot2::element_text(),
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
    panel.grid.minor.x = element_line(colour = "grey80",linetype = "dotted"  ) ,
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
    strip.text = element_blank() ,
    panel.spacing.x = unit(4, "mm") ,
   
  )
}
                
                
 

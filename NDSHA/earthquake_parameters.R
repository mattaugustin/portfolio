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


parameter_table  <-  parameter_table %>% relocate( . , c(directivity,event,channel,real) , .before = station_id  )
names(parameter_table) ; View(parameter_table)
write.csv(parameter_table, "portfolio/signals/quake_parameters.csv", row.names = FALSE)


# Clean data environment
rm(  parameter_table ,  stacked_list_ndsha_acc  ,stacked_ndsha_acc_tables, ndsha_ew_ns , signal_ndsha_acc_ew ,
     signal_ndsha_acc_ns ,path_signal , path_signal_bis ,  signal_ndsha_acc_ewns)


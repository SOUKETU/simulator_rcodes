
# example spatial points
library(sp);library(INLA);library(RandomFields)
library(magick); library(maptools)
ocean<-readShapeSpatial("/Users/jiecao/Desktop/snow_crab/simulator/50m_land.shp")
setwd('/Users/jiecao/Desktop/snow_crab/simulator')
SimFile = ('/Users/jiecao/Desktop/snow_crab/simulator/results')
SimFile = paste0(SimFile,'/',format(Sys.time(), "%d-%b-%Y %H.%M"))
dir.create(SimFile)
#source('Calc_Kmeans.r'); source('Calc_Anisotropic_Mesh.r'); 
source('/Users/jiecao/Desktop/snow_crab/simulator/simulator_rcodes/grow_tran.r')

# domain
# Specify the grid 
# angles <- c(0, 90)         # Orientations (in degrees) of easting and northing 
# length <- c(0.5, 0.5)        # Grid spacings (degree), east-west and north-south
# origin <- c(-179, 55)        # Grid origin coordinates (degree, projected)
# nrows <- 20                  # Number of west-east strips
# ncols <- 20                  # Number of north-south strips
# 
# # Create the points on the grid.
# basis <- rbind(cos(angles * pi/180), sin(angles * pi/180)) %*% diag(length)
# x <- as.vector(outer(1:ncols, 1:nrows, FUN=function(x,y) x-1))
# y <- as.vector(outer(1:ncols, 1:nrows, FUN=function(x,y) y-1))
# grid <- t(basis %*% rbind(x, y) + origin)
# coords = grid
# P1 = SpatialPoints(coords)
# plot(P1, axes = TRUE, cex=0.1)

LP = FALSE # Lognormal-Poisson for generating correlated counts

library(SpatialDeltaGLMM)
strata.limits <- data.frame('STRATA'="All_areas")
Region = "Eastern_Bering_Sea"
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
Data_Extrap = Extrapolation_List$Data_Extrap[Extrapolation_List$Data_Extrap$Area_in_survey_km2>0,]
coords = cbind(Data_Extrap$Lon,Data_Extrap$Lat)
coords = coords[coords[,1]<0,]

P1 = SpatialPoints(coords)

#plot(P1, axes = TRUE)

# Display the grid.
#plot(grid, asp=1, pch=19, col="black", cex=2 * (nrows * ncols)^(-1/4))

#data <- read.csv('snow_crab_survey_2016.csv',header = T)
#coords = cbind(data$LONGITUDE,data$LATITUDE)
#P1 = SpatialPoints(coords)
#plot(P1, axes = TRUE)

# data input 

  sp.points = P1
  n_s = nrow(sp.points@coords)
  loc_x = P1@coords

  # define population structure
    #n_p = 20 # numnber of size classes  
    n_sex = 2 # number of sexes
    n_g = 4  # 1-immature; 2-mature; 3-new shell; 4-old shell
    #n_v = 2  # number of shell condition/stages
    
    low_s = 25  # lower boundary of size range
    up_s  = 125 # upper boundary of size range
    bin_width = 20 
    
    n_p = (up_s-low_s)/bin_width # number of bins
    s_mid_points = seq(low_s+bin_width/2, up_s-bin_width/2, bin_width)
    
  # define time priod
    n_t = 10
    
    mat_at_size_female = c(0,0.1,0.3,0.5,0.9)#,0.3,0.4,0.5,0.8,0.9)#,0.5,0.7,0.9,0.9,0.9,0.9,0.95)
    mat_at_size_male = c(0,0.1,0.3,0.5,0.9)#,0.5,0.7,0.9,0.9,0.9,0.9,0.95)
    
    prob_old_female = c(0,0.1,0.3,0.5,0.9)#,0.5,0.7,0.9,0.9,0.9,0.9,0.95) # probability of old shell for female  
    prob_old_male = c(0,0.1,0.3,0.5,0.9)#,0.5,0.7,0.9,0.9,0.9,0.9,0.95) 
    
  # R
    R_mean_v = rpois(n_t,lambda=1e6/n_s); # per grid 
    if (LP==TRUE){
      R_mean_v = log(R_mean_v)
    }
    SD_R_v = rep(0.2,n_t) # spatial variation for each year 
    R_size_pert_v = c(1,0,rep(0,n_p-2)) # percentages for each size bin
    R_sex_ratio = 1 # male/total
    R_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t))
    
    for (i in 1:n_t){
      spde_model_R <- RMgauss(var=SD_R_v[i]^2, scale=0.5)
      R_kmeans = R_mean_v[i] * exp(RFsimulate(model = spde_model_R, x=loc_x[,1], y=loc_x[,2])@data[,1])
      if (LP==TRUE){
        R_kmeans = exp(R_kmeans)
        for (ii in 1:(length(R_kmeans))){
          R_kmeans[ii] = ifelse(R_kmeans[ii]<0.00001, 0, rpois(1, R_kmeans[ii]))
        }
      }
      R_at_size_s[,,i] = as.matrix(R_kmeans)%*%R_size_pert_v
    }

    R_male_at_size_s = R_sex_ratio*R_at_size_s
    R_female_at_size_s = (1-R_sex_ratio)*R_at_size_s
    
  # specify life history parameters
    # M
    M_spatial = FALSE # spatial or non-spatial
    M_mean_v = rep(0.23,n_t); 
    SD_M_v = rep(0.1,n_t) # spatial variation of M for each year
    M_size_v = rep(1,n_p) # size-specific M
  
    M_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t))
    if (!M_spatial) SD_M_v = rep(0,n_t)
    for(i in 1:n_t){
      spde_model_M <- RMgauss(var=SD_M_v[i]^2, scale=2)
      M_kmeans = M_mean_v[i] * exp(RFsimulate(model = spde_model_M, x=loc_x[,1], y=loc_x[,2])@data[,1])
      M_at_size_s[,,i] = as.matrix(M_kmeans)%*%M_size_v
    }
    
    # F
    
    # selectivity logistic
    
    ############# selectivities ###############
    # cal_sels <- function(pars_v,sizem_v)
    # {
    #     sels <- 1.0/(1.0+exp((pars_v[1]- sizem_v)* pars_v[2]));
    #     sels <- sels/max(sels)
    #   return(sels)
    #  }
    
    n_f = 3 # total catch; retained catch; bycatch
    pars_male = c(70,0.05);
    #pars_female = c(130,1);
    pars_male_re = c(109.59,0.09);
    pars_bycatch = c(109.59,0.09);
    
    sels_male <- 1.0/(1.0+exp((pars_male[1]- s_mid_points)* pars_male[2]));
    #sels_male <- sels_male/max(sels_male)
    
    #sels_female <- 1.0/(1.0+exp((pars_female[1]- s_mid_points)* pars_female[2]));
    #sels_female <- sels_female/max(sels_female)
    
    #sels_male_re <- 1.0/(1.0+exp((pars_male_re[1]- s_mid_points)* pars_male_re[2]));
    #sels_male_re <- sels_male_re/max(sels_male_re)
    
    #sels_bycatch <- 1.0/(1.0+exp((pars_bycatch[1]- s_mid_points)* pars_bycatch[2]));
    #sels_bycatch <- sels_bycatch/max(sels_bycatch)
    
    #sels_male = c(0.05,0.1,0.5,0.8,1)
    sels_female = c(0,0,0,0,0)
    sels_male_re = c(1,1,1,1,1)
    sels_bycatch = c(0,0,0,0,0)
    
    F_mean_v =  rnorm(n_t,mean=0.5,sd=0.1); 
    F_spatial = TRUE # spatial or non-spatial
    SD_F_v = rep(0.3,n_t) # spatial variation by year
    if (!F_spatial)  SD_F_v = rep(0.0,n_t)
    
    F_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t,n_sex))
    F_kmeans = matrix(NA, nrow=n_s, ncol=n_t)
    
    for(i in 1:n_t){
      spde_model_F <- RMgauss(var=SD_F_v[i]^2, scale=1)
      F_kmeans[,i] = F_mean_v[i] * exp(RFsimulate(model = spde_model_F, x=loc_x[,1], y=loc_x[,2])@data[,1])
      F_at_size_s[,,i,1] = as.matrix(F_kmeans[,i])%*%sels_male
      F_at_size_s[,,i,2] = as.matrix(F_kmeans[,i])%*%sels_female
    }
    
    # Growth
    G_spatial = FALSE # spatial or non-spatial
    
    n_grpar = 3
    int.male = rep(1,n_t); int.female = rep(1,n_t)
    slope.male = rep(1.5,n_t); slope.female = rep(1.5,n_t)
    beta.male = rep(0.5,n_t); beta.female = rep(0.5,n_t)
    
    sd.int.male = rep(0.2,n_t); sd.int.female = rep(0.2,n_t)
    sd.slope.male = rep(0.2,n_t); sd.slope.female = rep(0.2,n_t)
    sd.beta.male = rep(0.2,n_t); sd.beta.female = rep(0.2,n_t)
    
    if (!G_spatial)
      sd.int.male=sd.int.female=sd.slope.male=sd.slope.female=sd.beta.male=sd.beta.female=rep(0,n_t)
    
    growth_par = array (data = NA, dim = c(n_s,n_grpar,n_t,n_sex)) # 1-male; 2-female
    
    for(i in 1:n_t){
      
      spde_model_G_int_male        <- RMgauss(var=sd.int.male[i]^2, scale=2)
      spde_model_G_slope_male      <- RMgauss(var=sd.slope.male[i]^2, scale=2)
      spde_model_G_beta_male       <- RMgauss(var=sd.beta.male[i]^2, scale=2)
      int_kmeans_male = int.male[i] + RFsimulate(model = spde_model_G_int_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
      slope_kmeans_male = slope.male[i] + RFsimulate(model = spde_model_G_slope_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
      beta_kmeans_male = beta.male[i] + RFsimulate(model = spde_model_G_beta_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
      growth_par[,1,i,1] = int_kmeans_male
      growth_par[,2,i,1] = slope_kmeans_male
      growth_par[,3,i,1] = beta_kmeans_male
      
      spde_model_G_int_female        <- RMgauss(var=sd.int.female[i]^2, scale=2)
      spde_model_G_slope_female      <- RMgauss(var=sd.slope.female[i]^2, scale=2)
      spde_model_G_beta_female       <- RMgauss(var=sd.beta.female[i]^2, scale=2)
      int_kmeans_female = int.female[i] + RFsimulate(model = spde_model_G_int_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
      slope_kmeans_female = slope.female[i] + RFsimulate(model = spde_model_G_slope_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
      beta_kmeans_female = beta.female[i] + RFsimulate(model = spde_model_G_beta_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
      growth_par[,1,i,2] = int_kmeans_female
      growth_par[,2,i,2] = slope_kmeans_female
      growth_par[,3,i,2] = beta_kmeans_female
    }
    
    # n_grpar = 5
    # L_inf_v_male = rep(105,n_t); L_inf_v_female = rep(105,n_t);
    # SD_Linf_v_male = rep(0.1,n_t); SD_Linf_v_female = rep(0.1,n_t)
    # 
    # K_v_male = rep(0.1,n_t); K_v_female = rep(0.1,n_t)
    # SD_K_v_male = rep(0.01,n_t); SD_K_v_female = rep(0.01,n_t)
    # 
    # L_inf_sd_v_male = rep(0.5,n_t); L_inf_sd_v_female = rep(0.5,n_t)
    # SD_Linf_sd_v_male = rep(0.05,n_t); SD_Linf_sd_v_female = rep(0.05,n_t)
    # 
    # K_sd_v_male = rep(0.01,n_t); K_sd_v_female = rep(0.01,n_t)
    # SD_K_sd_v_male = rep(0.05,n_t); SD_K_sd_v_female = rep(0.05,n_t)
    # 
    # rho_v_male = rep(0.9,n_t); rho_v_female = rep(0.9,n_t)
    # SD_rho_v_male = rep(0.05,n_t); SD_rho_v_female = rep(0.05,n_t)
    # 
    # if (!G_spatial)  
    #   SD_Linf_v_male = SD_Linf_v_female = SD_K_v_male = SD_K_v_female = SD_Linf_sd_v_male = SD_Linf_sd_v_female = SD_K_sd_v_male = SD_K_sd_v_female = SD_rho_v_male = SD_rho_v_female = rep(0.0,n_t)
    # 
    # 
    # growth_par = array (data = NA, dim = c(n_s,n_grpar,n_t,n_sex)) # 1-male; 2-female
    # 
    # for(i in 1:n_t){
    #   spde_model_G_Linf_male   <- RMgauss(var=SD_Linf_v_male[i]^2, scale=0.1)
    #   spde_model_G_K_male      <- RMgauss(var=SD_K_v_male[i]^2, scale=0.1)
    #   spde_model_G_Linfsd_male <- RMgauss(var=SD_Linf_sd_v_male[i]^2, scale=0.1)
    #   spde_model_G_Ksd_male    <- RMgauss(var=SD_K_sd_v_male[i]^2, scale=0.1)
    #   spde_model_G_rho_male    <- RMgauss(var=SD_rho_v_male[i]^2, scale=0.1)
    #   
    #   Linf_kmeans_male = L_inf_v_male[i] + RFsimulate(model = spde_model_G_Linf_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   K_kmeans_male = K_v_male[i] + RFsimulate(model = spde_model_G_K_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   L_inf_sd_kmeans_male = L_inf_sd_v_male[i] + RFsimulate(model = spde_model_G_Linfsd_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   K_sd_kmeans_male = K_sd_v_male[i] + RFsimulate(model = spde_model_G_Ksd_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   rho_kmeans_male = rho_v_male[i] + RFsimulate(model = spde_model_G_rho_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   
    #   growth_par[,1,i,1] = Linf_kmeans_male
    #   growth_par[,2,i,1] = L_inf_sd_kmeans_male
    #   growth_par[,3,i,1] = K_kmeans_male
    #   growth_par[,4,i,1] = K_sd_kmeans_male
    #   growth_par[,5,i,1] = rho_kmeans_male
    #   
    #   spde_model_G_Linf_female   <- RMgauss(var=SD_Linf_v_female[i]^2, scale=0.1)
    #   spde_model_G_K_female      <- RMgauss(var=SD_K_v_female[i]^2, scale=0.1)
    #   spde_model_G_Linfsd_female <- RMgauss(var=SD_Linf_sd_v_female[i]^2, scale=0.1)
    #   spde_model_G_Ksd_female    <- RMgauss(var=SD_K_sd_v_female[i]^2, scale=0.1)
    #   spde_model_G_rho_female    <- RMgauss(var=SD_rho_v_female[i]^2, scale=0.1)
    #   
    #   Linf_kmeans_female = L_inf_v_female[i] + RFsimulate(model = spde_model_G_Linf_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   K_kmeans_female = K_v_female[i] + RFsimulate(model = spde_model_G_K_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   L_inf_sd_kmeans_female = L_inf_sd_v_female[i] + RFsimulate(model = spde_model_G_Linfsd_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   K_sd_kmeans_female = K_sd_v_female[i] + RFsimulate(model = spde_model_G_Ksd_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   rho_kmeans_female = rho_v_female[i] + RFsimulate(model = spde_model_G_rho_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    #   
    #   growth_par[,1,i,2] = Linf_kmeans_female
    #   growth_par[,2,i,2] = L_inf_sd_kmeans_female
    #   growth_par[,3,i,2] = K_kmeans_female
    #   growth_par[,4,i,2] = K_sd_kmeans_female
    #   growth_par[,5,i,2] = rho_kmeans_female
    # }
    
    ###################################################################################
    # equilibrium state with no fishing mortality; M happens before growth
          # n_g:  1-immature; 2-mature; 3-new shell; 4-old shell
    
    Plot_equil = FALSE
    n_t_equil = 10 # make sure n_t = n_t_equil otherwise growth and M don't match with the dimension
    N_at_size_grid_female = array (data = 0, dim = c(n_s,n_p,n_t_equil,n_g))
    N_at_size_grid_male   = array (data = 0, dim = c(n_s,n_p,n_t_equil,n_g))
    N_grid_total = array (data = 0, dim = c(n_s,n_t_equil,n_sex)) # 1-male; 2-female
    N_at_size_grid_female[,,1,1] = R_female_at_size_s[,,1]
    N_at_size_grid_male[,,1,1]   = R_male_at_size_s[,,1]
    N_grid_total[,1,1] = apply(N_at_size_grid_male[,,1,1],1,sum)
    N_grid_total[,1,2] = apply(N_at_size_grid_female[,,1,1],1,sum)
      
    for (i in 2:n_t_equil)
      {
        for (s in 1:n_s)
          {
              #growth_matrix_female = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,2][1],growth_par[s,,i-1,2][2],growth_par[s,,i-1,2][3],growth_par[s,,i-1,2][4],growth_par[s,,i-1,2][5])
              growth_matrix_female = growth_trans(growth_par[s,,i-1,2][1],growth_par[s,,i-1,2][2],growth_par[s,,i-1,2][3],low_s,up_s,bin_width)
              N_at_size_grid_female[s,,i,1] = R_female_at_size_s[s,,i] + ((N_at_size_grid_female[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]))%*%growth_matrix_female)*(1-mat_at_size_female)
              N_at_size_grid_female[s,,i,2] = ((N_at_size_grid_female[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]))%*%growth_matrix_female)*mat_at_size_female+N_at_size_grid_female[s,,i-1,2]*exp(-M_at_size_s[s,,i-1])
              N_at_size_grid_female[s,,i,3] = (N_at_size_grid_female[s,,i,1] + N_at_size_grid_female[s,,i,2])* (1-prob_old_female)
              N_at_size_grid_female[s,,i,4] = (N_at_size_grid_female[s,,i,1] + N_at_size_grid_female[s,,i,2])* prob_old_female
              N_grid_total[s,i,2] = sum((N_at_size_grid_female[s,,i,1] + N_at_size_grid_female[s,,i,2]))
              
              #growth_matrix_male = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],growth_par[s,,i-1,1][4],growth_par[s,,i-1,1][5])
              growth_matrix_male = growth_trans(growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],low_s,up_s,bin_width)
              N_at_size_grid_male[s,,i,1] = R_male_at_size_s[s,,i] + ((N_at_size_grid_male[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]))%*%growth_matrix_male)*(1-mat_at_size_male)
              N_at_size_grid_male[s,,i,2] = ((N_at_size_grid_male[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]))%*%growth_matrix_male)*mat_at_size_male+N_at_size_grid_male[s,,i-1,2]*exp(-M_at_size_s[s,,i-1])
              N_at_size_grid_male[s,,i,3] = (N_at_size_grid_male[s,,i,1] + N_at_size_grid_male[s,,i,2])* (1-prob_old_male)
              N_at_size_grid_male[s,,i,4] = (N_at_size_grid_male[s,,i,1] + N_at_size_grid_male[s,,i,2])* prob_old_male
              N_grid_total[s,i,1] = sum((N_at_size_grid_male[s,,i,1] + N_at_size_grid_male[s,,i,2]))
          }
     }
    
    N_at_size_equil_female = N_at_size_grid_female[,,n_t_equil,1]+N_at_size_grid_female[,,n_t_equil,2]-R_female_at_size_s[,,n_t_equil]
    N_at_size_equil_male   = N_at_size_grid_male[,,n_t_equil,1]+N_at_size_grid_male[,,n_t_equil,2]-R_male_at_size_s[,,n_t_equil]
    
    if(Plot_equil)
      {
        setwd(SimFile)
        col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
        min=min(N_grid_total);max=max(N_grid_total)
          for (i in 1:n_t_equil)
              {
                  P_male = P1; P_female = P1
                  P_male$male    = N_grid_total[,i,1]
                  P_female$femal = N_grid_total[,i,2]
                      
                  png(paste('Year',i,'.png',sep=''), height = 6, width = 12, units = 'in', res=900)
                  par(mfrow=c(1,2))
                  plot(ocean,col="dark gray",axes=T,xlim=c(-180,-165),ylim=c(55,60),main='Male')
                  image(SpatialPixelsDataFrame(P_male,tolerance=0.3,P_male@data),col=col(10),
                        axes=TRUE,zlim = c(min,max),add=T,
                        xlim=c(-180,-150),ylim=c(55,65))
                  mtext(paste('Year',i,sep=''),3,-2,adj=0.2)
                  
                  plot(ocean,col="dark gray",axes=T,xlim=c(-180,-165),ylim=c(55,60),main='Female')
                  image(SpatialPixelsDataFrame(P_female,tolerance=0.3,P_female@data),col=col(10),
                        axes=TRUE,zlim = c(min,max),add=T,
                        xlim=c(-180,-150),ylim=c(55,65))
                  mtext(paste('Year',i,sep=''),3,-2,adj=0.2)
                  
                  dev.off()
                  
                  img_temp = image_read(paste('Year',i,'.png',sep=''))
                  
                  if (i==1) img = img_temp
                  img = c(img,img_temp)
              }
          
          img <- image_scale(img, "1200")
          animation = image_animate(img, fps = 2, dispose = "previous")
          image_write(animation, "abundance.gif")
      }
    
    ############################################################
  
    # initial condition
        # n_g: 1-immature; 2-mature; 3-new shell; 4-old shell
    
    N_at_size_ns_female = array (data = 0, dim = c(n_s,n_p,n_t,n_g))
    N_at_size_ns_male   = array (data = 0, dim = c(n_s,n_p,n_t,n_g))
    N_ns_total          = array (data = 0, dim = c(n_s,n_t,n_sex))
    C_at_size_ns_male   = array (data = 0, dim = c(n_s,n_p,n_t,n_f)) # total catch; retained catch; bycatch
    C_at_size_ns_female = array (data = 0, dim = c(n_s,n_p,n_t)) # total catch - bycatch from trawl
    
    N_at_size_ns_female[,,1,]   = N_at_size_grid_female[,,n_t_equil,] 
    N_at_size_ns_female[,,1,1]  = N_at_size_ns_female[,,1,1] - R_female_at_size_s[,,n_t_equil] + R_female_at_size_s[,,1]
    N_at_size_ns_male[,,1,]    = N_at_size_grid_male[,,n_t_equil,] 
    N_at_size_ns_male[,,1,1]  = N_at_size_ns_male[,,1,1] - R_male_at_size_s[,,n_t_equil] + R_male_at_size_s[,,1]

    # population dynamics - July 1st, survey -> fishery -> growth
    
    alpha             = 0.62 # fishery time point  F_at_size_s
    F_bycatch_at_size = 0.00 * sels_bycatch
    
    for (i in 2:n_t)
      {   
          for (s in 1:n_s)
            {
                # female
                #growth_matrix_female = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,2][1],growth_par[s,,i-1,2][2],growth_par[s,,i-1,2][3],growth_par[s,,i-1,2][4],growth_par[s,,i-1,2][5])
                growth_matrix_female = growth_trans(growth_par[s,,i-1,2][1],growth_par[s,,i-1,2][2],growth_par[s,,i-1,2][3],low_s,up_s,bin_width)
                N_female_temp = (N_at_size_ns_female[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,2]))%*%growth_matrix_female
                N_at_size_ns_female[s,,i,1] = R_female_at_size_s[s,,i] + N_female_temp*(1-mat_at_size_female)
                N_at_size_ns_female[s,,i,2] = N_female_temp*mat_at_size_female + N_at_size_ns_female[s,,i-1,2]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,2])
                N_at_size_ns_female[s,,i,3] = (N_at_size_ns_female[s,,i,1] + N_at_size_ns_female[s,,i,2])* (1-prob_old_female)
                N_at_size_ns_female[s,,i,4] = (N_at_size_ns_female[s,,i,1] + N_at_size_ns_female[s,,i,2])* prob_old_female
                N_grid_total[s,i,2] = sum((N_at_size_ns_female[s,,i,1] + N_at_size_ns_female[s,,i,2]))
                
                # male
                #growth_matrix_male = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],growth_par[s,,i-1,1][4],growth_par[s,,i-1,1][5])
                growth_matrix_male = growth_trans(growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],low_s,up_s,bin_width)
                N_male_temp = (N_at_size_ns_male[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,1]))%*%growth_matrix_male
                N_at_size_ns_male[s,,i,1] = R_male_at_size_s[s,,i] + N_male_temp*(1-mat_at_size_male)
                N_at_size_ns_male[s,,i,2] = N_male_temp*mat_at_size_male+N_at_size_ns_male[s,,i-1,2]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,1])
                N_at_size_ns_male[s,,i,3] = (N_at_size_ns_male[s,,i,1] + N_at_size_ns_male[s,,i,2])* (1-prob_old_male)
                N_at_size_ns_male[s,,i,4] = (N_at_size_ns_male[s,,i,1] + N_at_size_ns_male[s,,i,2])* prob_old_male
                N_grid_total[s,i,1] = sum((N_at_size_ns_male[s,,i,1] + N_at_size_ns_male[s,,i,2]))
                
                # catch
                C_at_size_ns_female[s,,i-1] = (1-exp(-F_at_size_s[s,,i-1,2]-F_bycatch_at_size))*(N_at_size_ns_female[s,,i-1,1]+N_at_size_ns_female[s,,i-1,2])*exp(-alpha*M_at_size_s[s,,i-1])
                C_at_size_ns_male[s,,i-1,1] = (F_at_size_s[s,,i-1,1]/(F_at_size_s[s,,i-1,1]+F_bycatch_at_size))*(1-exp(-F_at_size_s[s,,i-1,1]-F_bycatch_at_size))*(N_at_size_ns_male[s,,i-1,1]+N_at_size_ns_male[s,,i-1,2])*exp(-alpha*M_at_size_s[s,,i-1])
                C_at_size_ns_male[s,,i-1,2] = C_at_size_ns_male[s,,i-1,1]*sels_male_re
                C_at_size_ns_male[s,,i-1,3] = (F_bycatch_at_size/(F_at_size_s[s,,i-1,1]+F_bycatch_at_size))*(1-exp(-F_at_size_s[s,,i-1,1]-F_bycatch_at_size))*(N_at_size_ns_male[s,,i-1,1]+N_at_size_ns_male[s,,i-1,2])*exp(-alpha*M_at_size_s[s,,i-1])
                
                if (i == n_t)
                  {
                    C_at_size_ns_female[s,,i] = (1-exp(-F_at_size_s[s,,i,2]-F_bycatch_at_size))*(N_at_size_ns_female[s,,i,1]+N_at_size_ns_female[s,,i,2])*exp(-alpha*M_at_size_s[s,,i])
                    C_at_size_ns_male[s,,i,1] = (F_at_size_s[s,,i,1]/(F_at_size_s[s,,i,1]+F_bycatch_at_size))*(1-exp(-F_at_size_s[s,,i,1]-F_bycatch_at_size))*(N_at_size_ns_male[s,,i,1]+N_at_size_ns_male[s,,i,2])*exp(-alpha*M_at_size_s[s,,i])
                    C_at_size_ns_male[s,,i,2] = C_at_size_ns_male[s,,i,1]*sels_male_re
                    C_at_size_ns_male[s,,i,3] = (F_bycatch_at_size/(F_at_size_s[s,,i,1]+F_bycatch_at_size))*(1-exp(-F_at_size_s[s,,i,1]-F_bycatch_at_size))*(N_at_size_ns_male[s,,i,1]+N_at_size_ns_male[s,,i,2])*exp(-alpha*M_at_size_s[s,,i])
                  }
            }
      }
    
  # save results
    
    Sim = list('N_at_size_ns_male'=N_at_size_ns_male,'C_at_size_ns_male'=C_at_size_ns_male,'F_at_size_s'=F_at_size_s)
    save(Sim, file=paste0(SimFile,'/',"Sim.RData"))
    
# produce images
    
    Dir = SimFile
    P= P1
    data = N_at_size_ns_male[,,,1] + N_at_size_ns_male[,,,2]
    
   # data = C_at_size_ns_male[,,,1]
    name = 'den'
   # Plot_annimation(P,data,n_t,Dir,name)  
   
    Plot_annimation = function(P,data,n_t,Dir,name,breaks)
    {
      setwd(paste(Dir))
      for (i in 1:n_t){
        P$data = apply(data[,,i],1,sum)
        #min=min(P$data);max=max(P$data)
        col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
        png(paste(name,'-Year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=600)
        rast.temp <- rasterize(P, rast, P$data, fun = mean)
        plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
        image(rast.temp,col=col(12),axes=TRUE,breaks=breaks,add=T,xlim=c(-180,-150),ylim=c(55,65))
        mtext(paste('Density of all size classes','','-','','Year',i,sep=''),1,-2,adj=0.2)
        dev.off()
        img_temp = image_read(paste(name,'-Year',i,'.png',sep=''))
        if (i==1) img = img_temp
        img = c(img,img_temp)
      }
      img <- image_scale(img, "600")
      animation = image_animate(img, fps = 0.5, dispose = "previous")
      image_write(animation, "abundance.gif")
    }
    
    # rast <- raster(ncol=200,nrow=50)
    # extent(rast) <- extent(coords)
    # rast2 <- rasterize(P, rast, P$data, fun = mean)
    # 
    # par(mfrow=c(4,3))
    # for (i in 1:n_t){
    #   P$data = data[,,i]
    #   min=min(P$data);max=max(P$data)
    #   col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
    #   #png(paste(name,'-Year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=600)
    #   plot(ocean,col="dark gray",axes=T,xlim=c(-180,-155),ylim=c(55,60))
    #   image(SpatialPixelsDataFrame(P,tolerance=0.51,P@data),col=col(10),
    #         axes=TRUE,zlim = c(min,max),add=T,pch=16,
    #         xlim=c(-180,-150),ylim=c(55,65))
    #   mtext(paste('Year',i,sep=''),3,-2,adj=0.2)
    #   #dev.off()
    # }
    
############################################################################################
    library(raster)
# density plots
    area = 13.71962
    Dir = SimFile
    name = 'den_sim' 
    P= P1
    den.data = N_at_size_ns_male[,,,1] + N_at_size_ns_male[,,,2]
    rast <- raster(ncol=200,nrow=100)
    extent(rast) <- extent(coords)
    col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
    breakpoints <- c(0,0.3,0.5,0.7,0.9,1.1,1.3,1.5,2,2.5,3,3.5,4,5)
    
    
    break_tot <- c(1,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,15)
    Plot_annimation(P,den.data/area,n_t,Dir,name,breaks=break_tot)
    
    setwd(Dir)
    png(paste(name,'-',Sys.Date(),'.png',sep=''), height = 9, width = 18, units = 'in', res=600)
    par(mfrow=c(n_p,n_t))
    par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    for(p in 1:n_p){
      for(t in 1:n_t){
        P$data = den.data[,p,t]/area
        rast.temp <- rasterize(P, rast, P$data, fun = mean)
        plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
        image(rast.temp,col=col(13),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
        mtext(paste('Size class',p,',','Year',t,sep=''),1,-2,adj=1,cex=0.5)
      }
    }
    dev.off()
    
    for(p in 1:n_p){
      png(paste(name,'_size-',p,'.png',sep=''), height = 8, width = 6, units = 'in', res=600)
      par(mfrow=c(4,3))
      par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
      par(tcl = -0.25)
      par(mgp = c(2, 0.6, 0))
    
        for(t in 1:n_t){
          P$data = den.data[,p,t]/area
          rast.temp <- rasterize(P, rast, P$data, fun = mean)
          plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
          image(rast.temp,col=col(13),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
          mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
        }
      dev.off()
    }
    
    # recruitment
    png(paste('recruitment.png'), height = 8, width = 6, units = 'in', res=600)
    par(mfrow=c(4,3))
    par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0)) 
    
    for(t in 1:n_t){
      P$data = R_male_at_size_s[,1,t]/area
      rast.temp <- rasterize(P, rast, P$data, fun = mean)
      plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
      image(rast.temp,col=col(10),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
      mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
    }
    dev.off()
    
    # fishing mortality
    breakpoints_f=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.5,2)
    
    png(paste('F_sim','.png',sep=''), height = 4, width = 8, units = 'in', res=600)
    par(mfrow=c(2,5))
    par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    
      for(t in 1:n_t){
        P$data = F_kmeans[,t]
        rast.temp <- rasterize(P, rast, P$data, fun = mean)
        plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
        image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints_f,add=T,xlim=c(-180,-150),ylim=c(55,65))
        mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
      }
    
    dev.off()
    
    # animation 
    for (i in 1:n_t){
      P$data = F_kmeans[,i]
      col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
      png(paste('F_sim','_year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=900)
      par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
      par(tcl = -0.25)
      par(mgp = c(2, 0.6, 0))
      rast.temp <- rasterize(P, rast, P$data, fun = mean)
      plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
      image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints_f,add=T,xlim=c(-180,-150),ylim=c(55,65))
      mtext(paste('Year',i,sep=''),1,-2,adj=0.95,cex=0.8)
      dev.off()
      img_temp = image_read(paste('F_sim','_year',i,'.png',sep=''))
      if (i==1) img = img_temp
      img = c(img,img_temp)
    }
    img <- image_scale(img, "900")
    animation = image_animate(img, fps = 0.5, dispose = "previous")
    image_write(animation, "F_sim.gif")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
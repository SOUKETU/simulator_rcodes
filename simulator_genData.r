
# generate survey data
library(Rlab)
n_loc_sample = n_loc_sample
sampling_error=FALSE
ObsModel_sim = c("PL", "Lognormal")[1]
obs = switch(ObsModel_sim,"PL"=1, "Lognormal"=2)
#logsigma_p = c(0.1,0.1,0.1,0.1,0.1)
logsigma_p = c(0.25,0.25,0.25,0.25,0.25)

data_male = c()
data_female = c()
catch_data =c()
area = 13.71962

for (i in 1:n_t){
  sample_id = sample(1:n_s, n_loc_sample, replace=FALSE)
  for (p in 1:n_p){
    temp_male = (N_at_size_ns_male[,,i,1] + N_at_size_ns_male[,,i,2])[,p]
    temp_female = (N_at_size_ns_female[,,i,1] + N_at_size_ns_female[,,i,2])[,p]
    for (k in 1:n_loc_sample){
      temp_n = temp_male[sample_id[k]]/area
      if(sampling_error){
        if(obs==2){
          temp_n=rlnorm(1,log(temp_n)-logsigma_p[p]^2,logsigma_p[p])
          #print(temp_n)
        }
        if(obs==1){
          encout = 1-exp(-temp_n)
          bp = rbern(1,encout)
          #print(bp)
            if(bp==0){
              temp_n = 0
            }else{
              temp_n = rlnorm(1,log(temp_n/encout)-logsigma_p[p]^2,logsigma_p[p])
            }
        }
      }
      temp1 = c(p,i,temp_n,coords[sample_id[k],],area)
      #temp1 = c(p,i,round(temp_male[sample_id[k]]/area,0),coords[sample_id[k],],area)
      data_male = rbind(data_male,temp1)
      

      temp_f = temp_female[sample_id[k]]/area
      if(sampling_error){
        if(obs==2){
          temp_f = rlnorm(1,log(temp_f),logsigma_p[p])
        }
        if(obs==1){
          bp = rbern(1,1-exp(-temp_n))
          if(bp==0){
            temp_f = 0
          }else{
            temp_f = rlnorm(1,log(temp_f),logsigma_p[p])
          }
        }
      }
      temp2 = c(n_p+p,i,temp_f,coords[sample_id[k],],area)
      #temp2 = c(n_p+p,i,round(temp_female[sample_id[k]]/area,0),coords[sample_id[k],],area)
      data_female = rbind(data_female,temp2)
    }
  }
}

index = matrix(0,nrow=n_t,ncol=n_p)
for(t in 1:n_t){
  for(p in 1:n_p){
    index[t,p] = sum(N_at_size_ns_male[,p,t,1] + N_at_size_ns_male[,p,t,2])
  }
}

survey_data = rbind(data_male,data_female)
if(R_sex_ratio == 1) survey_data = data_male
colnames(survey_data)=c('size_class','year','count','lon','lat','area_swept')

write.csv(survey_data,file = '/Users/jiecao/Desktop/SSST/simulation/survey_data.csv',row.names = FALSE)

catch_data_temp = C_at_size_ns_male[,,,1]

for (i in 1:n_t){
  catch_data = rbind(catch_data,catch_data_temp[,,i])
}
year = rep(seq(1,n_t,1),each=n_s)
lat = rep(coords[,2],n_t)
lon = rep(coords[,1],n_t)
catch_data = cbind(lat,lon,year,catch_data)

write.csv(catch_data,file = '/Users/jiecao/Desktop/SSST/simulation/catch_data.csv')

















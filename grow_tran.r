
# growth transition matrix 

# parameters int, slope, and beta (gamma - shape and scale)

growth_trans <- function(int,slope,beta,low_s,up_s,bin_width){
  n_p = (up_s-low_s)/bin_width # number of bins
  s_mid_points = seq(low_s+bin_width/2, up_s-bin_width/2, bin_width)
  growmat <- matrix(0, nrow = n_p, ncol = n_p)
  for (i in 1:n_p){
    mean = int + slope*s_mid_points[i] 
    alpha = mean/beta
    for (k in i:n_p){
      growmat[i,k] = pgamma((s_mid_points[k]+bin_width/2),alpha,scale=beta)-pgamma((s_mid_points[k]-bin_width/2),alpha,scale=beta)
    }
  }
  growmat <- growmat/rowSums(growmat, na.rm = T)
  return(growmat)
}
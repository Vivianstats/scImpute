### root-finding equation
fn = function(alpha, target){
  log(alpha) - digamma(alpha) - target
}

### update parameters in gamma distribution
update_gmm_pars = function(x, wt){
  tp_s = sum(wt)
  tp_t = sum(wt * x)
  tp_u = sum(wt * log(x))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v, 
                      extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}

### estimate parameters in the mixture distribution
get_mix = function(xdata, point){
  inits = rep(0, 5)
  inits[1] = sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {inits[1] = 0.01}
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > point]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {inits[5] = 0}
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0
  
  while(eps > 0.5) {
    wt = calculate_weight(xdata, paramt)
    paramt[1] = sum(wt[, 1])/nrow(wt)
    paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
    paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))
    paramt[2:3] = update_gmm_pars(x=xdata, wt=wt[,1])
    
    loglik = sum(log10(dmix(xdata, paramt)))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100) 
      break
  }
  return(paramt)
}

get_mix_parameters <-
function (count, point = log10(1.01), path, ncores = 8) 
{
    count = as.matrix(count)
    null_genes = which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
    parslist = mclapply(1:nrow(count), function(ii) {
      if (ii%%500 == 0) {
        gc()
        print(ii)
      }
      if (ii %in% null_genes) {
        return(rep(NA, 5))
      }
      xdata = count[ii, ]
      paramt = try(get_mix(xdata, point), silent = TRUE)
      if (class(paramt) == "try-error"){
        paramt = rep(NA, 5)
      }
      return(paramt)
    }, mc.cores = ncores)
    save(parslist, file = path)
    parslist = Reduce(rbind, parslist)
    colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
    saveRDS(parslist, file = path)
    return(0)
}


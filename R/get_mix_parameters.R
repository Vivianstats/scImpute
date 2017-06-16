get_mix_parameters <-
function (count, point = log10(1.01), path, ncores = 8) 
{
    count = as.matrix(count)
    null_genes = which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
    parslist = mclapply(1:nrow(count), function(ii) {
        # print(ii)
        if (ii%%500 == 0) {
            gc()
            print(ii)
        }
        if (ii %in% null_genes) {
            return(rep(NA, 5))
        }
        xdata = count[ii, ]
        inits = rep(0, 5)
        inits[1] = sum(xdata == point)/length(xdata)
        if (inits[1] == 0) {inits[1] = 0.01}
        inits[2:3] = c(0.5, 1)
        xdata_rm = xdata[xdata > point]
        inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
        if (is.na(inits[5])) {inits[5] = 0}
        paramt = inits
        A = diag(1, 2)
        B = matrix(0, ncol = 1, nrow = 2)
        eps = 1
        iter = 0
        loglik = 0
        loglik_old = 0
        while (eps > 0.5) {
            wt = calculate_weight(xdata, paramt)
            paramt_old = paramt
            paramt[1] = sum(wt[, 1])/nrow(wt)
            paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
            paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 
                2]))
            mstep_gamma = maxLik(fn_gamma, grad = fn_gamma_dev, 
                start = paramt[2:3], constraints = list(ineqA = A, 
                  ineqB = B), x = xdata, wt = wt[, 1])
            paramt[2:3] = mstep_gamma$estimate
            loglik = sum(log10(dmix(xdata, paramt)))
            eps = (loglik - loglik_old)^2
            iter = iter + 1
            if (iter > 100) 
                break
        }
        return(paramt)
    }, mc.cores = ncores)
    save(parslist, file = path)
    parslist = Reduce(rbind, parslist)
    colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
    saveRDS(parslist, file = path)
    return(0)
}

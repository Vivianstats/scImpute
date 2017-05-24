fn_gamma_dev <-
function (pars, x, wt) 
{
    alp = pars[1]
    beta = pars[2]
    dev_alp = sum(wt * (log(beta) - digamma(alp) + log(x)))
    dev_beta = sum(wt * (alp/beta - x))
    return(c(dev_alp, dev_beta))
}

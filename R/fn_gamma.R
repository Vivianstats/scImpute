fn_gamma <-
function (pars, x, wt) 
{
    alp = pars[1]
    beta = pars[2]
    tp = wt * (alp * log(beta) - lgamma(alp) + (alp - 1) * log(x) - 
        beta * x)
    return(sum(tp))
}

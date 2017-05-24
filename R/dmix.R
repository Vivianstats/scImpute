dmix <-
function (x, pars) 
{
    pars[1] * dgamma(x, shape = pars[2], rate = pars[3]) + (1 - 
        pars[1]) * dnorm(x, mean = pars[4], sd = pars[5])
}

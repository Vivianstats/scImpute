calculate_weight <-
function (x, paramt) 
{
    pz1 = paramt[1] * dgamma(x, shape = paramt[2], rate = paramt[3])
    pz2 = (1 - paramt[1]) * dnorm(x, mean = paramt[4], sd = paramt[5])
    pz = pz1/(pz1 + pz2)
    pz[pz1 == 0] = 0
    return(cbind(pz, 1 - pz))
}

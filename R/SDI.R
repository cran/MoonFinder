SDI <-
function(x){
    n = length(x)
    y = rep(1/n,n)
    z1 = 1/sum(x^2)
    s = 1/sum(y^2)
    z = round(z1/s,4)
    return(z)
}

moonlightingCoefficient <-
function(rna2mod,simMat){
    RNAs = rownames(rna2mod)
    rnaStat = matrix(0,length(RNAs),3)
    for (j in 1:length(RNAs)){
        index = which(rna2mod[j,]==1)
        tmpSimMat = simMat[index,index]
        EigenCor <- eigen(tmpSimMat)
        EV <- EigenCor$values
        EVratio = EV/sum(EV)
        # Simpson's diversity index
        sdiScore = SDI(EVratio) # 1/sum(x^2)
        rnaStat[j,] = c(RNAs[j],length(index),sdiScore)
     }
    colnames(rnaStat) = c("ncRNA","Module Number","SDI")
    return(rnaStat)
}

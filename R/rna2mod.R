rna2mod <-
function(rna2prot,modlist,pCutoff = 0.01,bgProtNum){
    moduleNum = length(modlist)
    names(modlist) = 1:moduleNum
    uniRNA = names(which(sort(table(rna2prot[,1]),decreasing=TRUE)>1))
    rnaNum = length(uniRNA)
    RNA2moduleMat = matrix(0,rnaNum,moduleNum)
    colnames(RNA2moduleMat) = 1:moduleNum
    for (i in 1:rnaNum){
        tmpRNA = uniRNA[i]
        targets = rna2prot[rna2prot[,1]%in%tmpRNA,2]
        targetNum = length(targets)
        for (j in 1:moduleNum){ 
            moduleProt = modlist[[j]]
            moduleSize = length(moduleProt)
            hitNum = length(intersect(moduleProt,targets))
            RNA2moduleMat[i,j]=1-phyper(hitNum-1,moduleSize,
            bgProtNum-moduleSize,targetNum)
        }
    }

    # return a mapping matrix with rows represent ncRNAs 
    # while columns represent modules
    RNA2moduleMatrix = matrix(0,rnaNum,moduleNum)
    colnames(RNA2moduleMatrix)=1:moduleNum
    rownames(RNA2moduleMatrix)=uniRNA
    RNA2moduleMatrix[RNA2moduleMat < pCutoff]=1
    modNum0 = apply(RNA2moduleMatrix,1,sum)
    rnaNum0 = apply(RNA2moduleMatrix,2,sum)
    rna2module = RNA2moduleMatrix[modNum0>1,rnaNum0>0]
    colnames(rna2module)=1:ncol(rna2module)

    # return a list of modules targeted by at least one ncRNA
    index1 = as.numeric(colnames(rna2module))
    moduleList = modlist[index1]
    names(moduleList) = paste0("mapped_",1:length(moduleList))
    return(list(rna2module,moduleList))
}

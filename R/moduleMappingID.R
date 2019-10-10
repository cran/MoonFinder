moduleMappingID <-
function(inputList,ID1="UNIPROT",ID2="SYMBOL",thdMS=5){
    tmpList = list()
    for(i in 1:length(inputList)){
        unip_symb = bitr(inputList[[i]], fromType=ID1, toType=ID2,
        OrgDb="org.Hs.eg.db")
        tmpList[[i]] = unip_symb[,2]
    }
    mLen = sapply(tmpList,length)
    index1 = which(mLen>=thdMS)
    outputList = unique(tmpList[index1])
    return(outputList)
}

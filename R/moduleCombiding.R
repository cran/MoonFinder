moduleCombiding <-
function(inputList,OScutoff = 0.5){
    moduleLen = length(inputList)
    matchMatrix = matrix(0,moduleLen,moduleLen)
    for (i in 1:moduleLen){
        a = length(inputList[[i]])
        for(j in 1:moduleLen){
            b = length(inputList[[j]])
            olpLen = length(intersect(inputList[[i]],inputList[[j]]))
            matchMatrix[i,j] = olpLen/sqrt(a*b)
        }
    }

    gg  <- graph.adjacency(matchMatrix >= OScutoff,diag=FALSE,
	    mode="undirected")
    gg.cluster = clusters(gg)
    clusterNum = length(gg.cluster$csize)
    clusterList = list()
    for(i in 1:clusterNum){
        index = which(gg.cluster$membership==i)
        clusterList[[i]] = unique(unlist(inputList[index]))
    }
    moduleSizes = sapply(clusterList,length)
    index1 = order(moduleSizes,decreasing=TRUE)
    outputList = clusterList[index1]
    return(outputList)
}

moduleHeatmap <-
function(rna,rna2mod,modList,net){
    color1 = brewer.pal(n=11,name="RdYlBu")
    modules = modList[which(rna2mod[rownames(rna2mod)==rna,]==1)]
    mNum = length(modules)
    names(modules) = paste0("M",1:mNum)
    moduleProt = unlist(modules)
    protNum = length(moduleProt)
    moduleProt1 = make.names(moduleProt, unique=TRUE)
    uniProt = names(V(net))
    gm <- induced.subgraph(net,intersect(moduleProt,uniProt))
    links = get.edgelist(gm)
    adjMat = matrix(0,protNum,protNum)
    rownames(adjMat) = moduleProt
    colnames(adjMat) = moduleProt
    simMat = adjMat
    x1 = match(links[,1],moduleProt)
    x2 = match(links[,2],moduleProt)
    adjMat[cbind(x1,x2)] = 1
    adjMat[cbind(x2,x1)] = 1
    diag(adjMat) = 1

    title = paste0(rna," -- ",mNum," modules")
    anno_row <- data.frame(Module = rep(paste0("M",1:mNum), 
        sapply(modules,length)))
    row.names(anno_row) <- moduleProt1
    # dev.new(width=8, height=8)
    color2 = c("white",color1[4])
    pheatmap(adjMat,annotation_row=anno_row,cluster_rows=FALSE,
        cluster_cols=FALSE,show_colnames=FALSE,
        color=colorRampPalette(color2)(5),main=title,legend=FALSE)
}

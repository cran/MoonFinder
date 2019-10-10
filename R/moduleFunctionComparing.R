moduleFunctionComparing <-
function(rna,rna2mod,modList){
    modules = modList[which(rna2mod[rownames(rna2mod)==rna,]==1)]
    names(modules) = paste0("M",1:length(modules))
    cgo <- compareCluster(geneClusters = modules,keyType="SYMBOL",
        OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff=0.01)
    # head(as.data.frame(cgo))
    dotplot(cgo)
}

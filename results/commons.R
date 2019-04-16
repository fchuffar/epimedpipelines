
plot_expr_hm = function(data, nb_grp_row=4, nb_grp_col=4, rsc, csc, CLUST_COL=TRUE, main="", ZSCORE=FALSE) {
  # Remove lines with no variation
  data = data[apply(data, 1, function(l) {length(unique(l))})>1, ]
  
  # colnames(data) = s$exp_grp[colnames(data),]$tissue_group_level1
  if (ZSCORE) {
    data = data - apply(data, 1, mean)
    data = data / apply(data, 1, sd)    
  }

  # hc_col = hclust(dist(1 - cor(data, method="pe")), method="complete")
  # foo = s$exp_grp[colnames(data)[hc_col$order],]
  # foo$tmp = foo$tissue_group_level1
  # foo$tmp = factor(foo$tmp, levels=unique(foo$tmp))
  # data = data[,rownames(foo)[order(foo$tmp)]]


  if (CLUST_COL) {
    # # clustering base on eucl. dist for tissues
    tmp_d = t(data)
    tmp_d = tmp_d[,!apply(is.na(tmp_d), 2, any)]
    hc_col = hclust(dist(tmp_d), method="complete")
    Colv = as.dendrogram(hc_col)
    dendrogram="col"


    # # clustering base on correlation for tissues
    # tmp_d = data
    # tmp_d = tmp_d[!apply(is.na(tmp_d), 1, any), ]
    # tmp_d = t(tmp_d) - apply(tmp_d, 2, mean)
    # tmp_d = t(tmp_d)
    # tmp_d = cor(tmp_d, method="pe")
    # # dim(tmp_d)
    # hc_col = hclust(dist(1 - tmp_d), method="complete")
    # Colv = as.dendrogram(hc_col)
    # dendrogram="col"

    cs = rep(2:10, each=100)
    res = sapply(cs, function(c) {
      k = kmeans(tmp_d, centers=c);  
      ret = k$betweenss / k$totss
      return(ret)
    })

    best_score = NULL
    for (i in 1:10) {
      k = kmeans(tmp_d, centers=nb_grp_col)
      score = k$betweenss / k$totss
      print(score)
      if (is.null(best_score)) {
        best_k = k
        best_score = score
      } else if (score > best_score) {
        best_k = k
        best_score = score
      }
    }

    ColSideColors = palette(RColorBrewer::brewer.pal(n=8, "Dark2"))[best_k$cluster[colnames(data)]]    
    names(ColSideColors) = colnames(data)    

    if (!missing(csc)) {
      ColSideColors = csc
    }

    # PCA on tissues
    pca = prcomp(tmp_d, scale=FALSE)
    v = pca$sdev * pca$sdev
    p = v / sum(v) * 100
    layout(matrix(1:6,2, byrow=FALSE), respect=TRUE)
    barplot(p)
    i=3
    j=2
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    i=1
    j=3
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    i=1
    j=2
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    i=4
    j=5
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])

    plot(jitter(cs), jitter(res))
    points(nb_grp_col, best_score, col=2)
    
  } else {
    Colv = NULL
    ColSideColors = rep("white", ncol(data))
    hc_col=NULL
  }

  # stop("EFN")

  # # clustering base on eucl. dist. for feats
  # tmp_d = data
  # tmp_d = tmp_d[,!apply(is.na(tmp_d), 2, any)]
  # d = dist(tmp_d)
  # hc_row = hclust(d, method="complete")
  # Rowv = as.dendrogram(hc_row)
  # dendrogram="row"

  ## features
  # clustering base on eucl. dist. for feats
  tmp_d = data
  tmp_d = tmp_d[!apply(is.na(tmp_d), 1, any), ]
  tmp_d = cor(t(tmp_d), method="pe")
  # dim(tmp_d)
  hc_row = hclust(dist(1 - tmp_d), method="complete")
  Rowv = as.dendrogram(hc_row)
  if (!is.null(Colv)) {
    dendrogram="both"
  } else {
    dendrogram="row"  
  }

  if (missing(rsc)) {
    grps = list()
    ct = cutree(hc_row, nb_grp_row)
    for (i in 1:nb_grp_row) {
      grps[[palette()[i]]] = names(ct)[ct==i]
    }
    # print(grps)
    RowSideColors = palette()[ct[rownames(data)]]    
    names(RowSideColors) = rownames(data)    
  } else {
    RowSideColors = rep("white", nrow(data))
    names(RowSideColors) = rownames(data)
    idx = intersect(rownames(data), names(rsc))
    RowSideColors[idx] = rsc[idx]
  }



  # col
  # colors = c("green", "black", "red")
  colors = c("cyan", "black", "red")
  # colors = c("blue", "yellow", "red")
  # colors = rev(RColorBrewer::brewer.pal(n=11, "RdYlBu"))
  cols = colorRampPalette(colors)(20)
  foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, trace="none", col=cols, main=paste0(main, " (", nrow(data), " features x ", ncol(data), " tissues)"), mar=c(10,5), useRaster=TRUE, RowSideColors=RowSideColors, ColSideColors=ColSideColors, cex.axis=0.5)  
  return(list(rsc=RowSideColors, csc=ColSideColors, hc_row=hc_row, hc_col=hc_col))  
}

reduce_rows = function(tmp_meth_data, map,  indicator_func2=mean, ...) {
  dim(tmp_meth_data)
  meth_by_tissues_by_feature = epimedtools::monitored_apply(mod=10, t(t(names(map))), 1, function(f) {
    # print(f)
    # f = features[1,]
    # f = features["HACD4",]
    probe_idx = intersect(rownames(tmp_meth_data), map[[f]])
    if (length(probe_idx) == 0) {    
      # print(f)
      tmp_meth_by_tissues_by_feature = rep(NA,ncol(tmp_meth_data))
      names(tmp_meth_by_tissues_by_feature) = colnames(tmp_meth_data)
    } else if (length(probe_idx) > 1) {
      tmp_meth_by_tissues_by_feature = apply(tmp_meth_data[probe_idx,], 2, indicator_func2, ...)    
    } else {
      tmp_meth_by_tissues_by_feature  = sapply(tmp_meth_data[probe_idx,], indicator_func2, ...) 
    }
    return(tmp_meth_by_tissues_by_feature)
  })
  colnames(meth_by_tissues_by_feature) = names(map)
  meth_by_tissues_by_feature = t(meth_by_tissues_by_feature)
  return(meth_by_tissues_by_feature)
}

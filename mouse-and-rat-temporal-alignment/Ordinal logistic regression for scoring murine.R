library(psupertime)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(cowplot)
library(data.table)
library(scales)
library(xfun)

#psudotime in stemcell in mouse
data_integrated@active.assay <- "integrated"
OLR <- as.SingleCellExperiment(subset(subset(data_integrated, subset = species %in% "mouse"), subset = Celltype %in% c("Neuroepithelium cell (NEC)","Radial glia cell (RGC)","Intermediate progenitor (IP)")))
OLR_v1 <- as.SingleCellExperiment(subset(subset(data_integrated, subset = species %in% "rat"), subset = Celltype %in% c("Neuroepithelium cell (NEC)","Radial glia cell (RGC)","Intermediate progenitor (IP)")))

scoremodel <- psupertime(OLR, str_remove(OLR$day, "E"), n_folds = 10)
images <- psupertime_plot_all(scoremodel, label_name = "Developmental stage (days)", ext = "pdf")
plot_labels_over_psupertime_homemade(scoremodel, palette = "Blues")
plot_identified_genes_over_psupertime_homemade (scoremodel, n_to_plot=6)
goanalysis <- psupertime_go_analysis_homemade(scoremodel, 'org.Mm.eg.db', k = 3, sig_cutoff = 5)
plot_profiles_of_gene_clusters_homemade(goanalysis)
plot_heatmap_of_gene_clusters_homemade(goanalysis)
plot_go_results_homemade(goanalysis, sig_cutoff = 2, p_cutoff = 0.23)
plot_specified_genes_over_psupertime(scoremodel, "Ccnd1", palette = "Blues")

view(scoremodel$beta_dt)
view(goanalysis$clusters_dt)
view(scoremodel$x_data)
write.csv(scoremodel$beta_dt, "gene candidates in all stem cell in mouse.csv")

#psudotime in stemcell in rat
scoremodel_v1 <- psupertime(OLR_v1, str_remove(OLR_v1$day, "E"), n_folds = 10)
images <- psupertime_plot_all(scoremodel_v1, label_name = "Developmental stage (days)", ext = "pdf")
plot_labels_over_psupertime_homemade(scoremodel_v1, palette = "Reds")
plot_identified_genes_over_psupertime_homemade (scoremodel_v1, n_to_plot=6)
goanalysis_v1 <- psupertime_go_analysis_homemade(scoremodel_v1, 'org.Mm.eg.db', k = 3, sig_cutoff = 5)
plot_profiles_of_gene_clusters_homemade(goanalysis_v1)
plot_heatmap_of_gene_clusters_homemade(goanalysis_v1)
plot_go_results_homemade(goanalysis_v1, sig_cutoff = 2, p_cutoff = 0.23)
plot_specified_genes_over_psupertime(scoremodel_v1, "Ccnd1", palette = "Blues")

view(scoremodel_v1$beta_dt)
view(goanalysis_v1$clusters_dt)
view(scoremodel_v1$x_data)
write.csv(scoremodel_v1$beta_dt, "gene candidates in all stem cell in rat.csv")


#psudotime double projection
Geneshared <- intersect(scoremodel$beta_dt$symbol,scoremodel_v1$beta_dt$symbol)
write(Geneshared, "shared clock genes.csv")
#Overlapped with E14-E16 DEG
intersect(rownames(head(mono.de, n = 50)), Geneshared)

mutualprojection <- psupertime::double_psupertime(scoremodel,scoremodel_v1, run_names= c("mouse", "rat"))
plot_double_psupertime_histogram(double_obj = mutualprojection)
plot_double_psupertime_xy(double_obj = mutualprojection)



#identify color function
.make_col_vals <- function(y_labels, palette='Blues') {
  n_labels 	= length(levels(y_labels))
  max_col 	= 11
  if (n_labels==1) {
    col_vals 	= brewer.pal(4, palette)
    col_vals 	= col_vals[1]
  } else if (n_labels==2) {
    col_vals 	= brewer.pal(4, palette)
    col_vals 	= col_vals[-2]
  } else if (n_labels<=max_col) {
    col_vals 	= brewer.pal(n_labels, palette)
  } else {
    col_pal 	= brewer.pal(max_col, palette)
    col_vals 	= rev(colorRampPalette(col_pal)(n_labels))
  }
  
  return(col_vals)
}
#homemade plot specific gene
plot_specified_genes_over_psupertime_homemade <- function(psuper_obj, extra_genes, label_name='Developmental stages', palette='Blues 2', plot_ratio=1.25, ylimit = c(-1,3)) {
  # unpack
  proj_dt     = psuper_obj$proj_dt
  beta_dt     = psuper_obj$beta_dt
  x_data      = psuper_obj$x_data
  ps_params 	= psuper_obj$ps_params
  
  # restrict to specified genes
  extra_genes = intersect(extra_genes, colnames(x_data))
  if (length(extra_genes)==0) {
    warning('genes not found; did not plot')
    return()
  }
  
  # set up data
  plot_wide   = cbind(proj_dt, data.table(x_data[, extra_genes, drop=FALSE]))
  plot_dt     = melt.data.table(
    plot_wide, 
    id 				= c("psuper", "label_input", "label_psuper"), 
    measure 		= extra_genes, 
    variable.name 	= "symbol"
  )
  plot_dt[, `:=`(symbol, factor(symbol, levels = extra_genes))]
  
  # set up plot
  col_vals 	= .make_col_vals(plot_dt$label_input, palette)
  n_genes 	= length(extra_genes)
  ncol 		= ceiling(sqrt(n_genes*plot_ratio))
  nrow 		= ceiling(n_genes/ncol)
  
  # plot
  g =	ggplot(plot_dt) +
    aes( x=psuper, y=value ) +
    geom_point( size=1, aes(colour=label_input) ) +
    geom_smooth(se=FALSE, colour='black') +
    scale_colour_manual( values=col_vals ) +
    scale_x_continuous( ) +
    scale_y_continuous( limits = ylimit) +
    facet_wrap( ~ symbol, scales='free_y', nrow=nrow, ncol=ncol ) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(
      axis.text.x = element_blank()
    ) +
    labs(
      x 		= 'Pseudotime'
      ,y 		= 'z-scored log2 expression'
      ,colour = label_name
    )
  return(g)
}
#homemade psupertime label over time
plot_labels_over_psupertime_homemade <- function(psuper_obj, label_name='Embryonic stage (day)', palette='RdBu') {
  # unpack
  proj_dt 		= psuper_obj$proj_dt
  cuts_dt 		= psuper_obj$cuts_dt
  
  # make nice colours
  col_vals 		= .make_col_vals(proj_dt$label_input, palette)
  
  # plot
  g = ggplot(proj_dt) +
    aes( x=psuper, fill=label_input, colour=label_input ) +
    geom_density( alpha=0.5 ) +
    scale_fill_manual( values=col_vals ) +
    #   geom_vline( data=cuts_dt, aes(xintercept=psuper, colour=label_input) ) +
    scale_colour_manual( values=col_vals ) +
    guides(
      fill 	= guide_legend(override.aes = list(alpha=1))
      ,colour = FALSE
    ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    labs(
      x 		= 'Psudotime score'
      ,y 		= 'Density'
      ,fill 	= label_name
    ) +
    theme_void()
  
  return(g)
}
#homemade identified genes
plot_identified_genes_over_psupertime_homemade <- function(psuper_obj, label_name='Embryonic stage (day)', n_to_plot=30, palette='Blues', plot_ratio=1) {
  # unpack
  proj_dt 	= psuper_obj$proj_dt
  #set shared gene beta_dt 	= subset(psuper_obj$beta_dt, symbol %in% Geneshared) 
  beta_dt 	= psuper_obj$beta_dt 
  x_data 		= psuper_obj$x_data
  ps_params 	= psuper_obj$ps_params
  # aset
  beta_nzero 	= beta_dt[ abs_beta > 0 ]
  n_nzero 	= nrow(beta_nzero)
  top_genes 	= as.character(beta_nzero[1:min(n_to_plot, nrow(beta_nzero))]$symbol)
  # set up data for plotting
  plot_wide 	= cbind(proj_dt, data.table(x_data[, top_genes, drop=FALSE]))
  plot_dt 	= melt.data.table(plot_wide, id=c('cell_id', 'psuper', 'label_input', 'label_psuper'), measure=top_genes, variable.name='symbol')
  plot_dt[, symbol := factor(symbol, levels=top_genes)]
  # get colours
  col_vals 	= .make_col_vals(plot_dt$label_input, palette)
  n_genes 	= length(top_genes)
  ncol 		= ceiling(sqrt(n_genes*plot_ratio))
  nrow 		= ceiling(n_genes/ncol)
  # plot
  g =	ggplot(plot_dt) +
    aes( x=psuper, y=value) +
    geom_point( size=1, aes(colour=label_input) ) +
    geom_smooth(se=FALSE, colour='black') +
    scale_colour_manual( values=col_vals ) +
    scale_shape_manual( values=c(1, 16) ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    facet_wrap( ~ symbol, scales='free_y', nrow=nrow, ncol=ncol) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(), axis.title=element_text(size=12), legend.text = element_text(size=13), legend.title=element_text(size=13)
    ) +
    labs(
      x 		= 'Psudotime score'
      ,y 		= 'z-scored log2 expression'
      ,colour = label_name
    ) 
  return(g)
}
#homemade GO in sum
psupertime_go_analysis_homemade <- function(psuper_obj, org_mapping, k=5, sig_cutoff=5) {
  if ( !requireNamespace("topGO", quietly=TRUE) ) {
    message('topGO not installed; not doing GO analysis')
    return()
  }
  if ( !requireNamespace("fastcluster", quietly=TRUE) ) {
    message('fastcluster not installed; not doing GO analysis')
    return()
  }
  
  # unpack
  glmnet_best 	= psuper_obj$glmnet_best
  best_lambdas 	= psuper_obj$best_lambdas
  proj_dt 		= copy(psuper_obj$proj_dt)
  x_data 			= copy(psuper_obj$x_data)
  beta_dt 		= psuper_obj$beta_dt
  cuts_dt 		= psuper_obj$cuts_dt
  
  # put cells in nice order, label projections
  rownames(x_data) 		= sprintf('cell_%04d', 1:nrow(x_data))
  set(proj_dt, i=NULL, 'cell_id', rownames(x_data))
  setorder(proj_dt, psuper)
  
  # do clustering on symbols
  message('clustering genes')
  hclust_obj 		= fastcluster::hclust(dist(t(x_data)), method='ward.D')
  
  # extract clusters from them
  clusters_dt 	= .calc_clusters_dt(hclust_obj, x_data, proj_dt, k)
  go_results 		= .do_topgo_for_cluster(clusters_dt, sig_cutoff, org_mapping)
  
  # make plot_dt
  plot_dt 		= .make_plot_dt(x_data, hclust_obj, proj_dt, clusters_dt)
  
  # assemble outputs
  go_list = list(
    clusters_dt 	= clusters_dt,
    go_results 		= go_results,
    plot_dt 		= plot_dt,
    cuts_dt 		= copy(psuper_obj$cuts_dt)
  )
  
  return(go_list)
}
#homemade GO supplementary
.calc_clusters_dt <- function(hclust_obj, x_data, proj_dt, k=5) {
  # make thing
  clusters_dt 	= data.table( h_clust=cutree(hclust_obj, k=k), symbol=colnames(x_data))
  # add clustering
  clusters_dt[, N:=.N, by=h_clust ]
  
  # order by correlation with psupertime
  temp_dt 			= data.table(melt(x_data, varnames=c('cell_id', 'symbol')))
  temp_dt 			= clusters_dt[ temp_dt, on='symbol' ]
  means_dt 			= temp_dt[, list(mean=mean(value)), by=list(cell_id, h_clust) ]
  means_dt 			= proj_dt[ means_dt, on='cell_id' ]
  corrs_dt 			= means_dt[, list( cor=cor(mean, psuper) ), by=h_clust]
  setorder(corrs_dt, cor)
  corrs_dt[, clust := 1:.N ]
  corrs_dt[, clust := factor(clust)]
  
  # add clusters ordered by size back in
  clusters_dt			= corrs_dt[ clusters_dt, on='h_clust' ]
  clusters_dt[, clust_label := factor(sprintf('%02d (%d genes)', clust, N)) ]
  clusters_dt[, h_clust := NULL ]
  setorder(clusters_dt, clust, symbol)
  
  return(clusters_dt)
}
.do_topgo_for_cluster <- function(clusters_dt, sig_cutoff, org_mapping) {
  # set up
  all_clusters 	= unique(clusters_dt[N>=sig_cutoff]$clust)
  go_results 		= data.table()
  
  # loop through clusters
  message(sprintf('calculating GO enrichments for %d clusters:', length(all_clusters)))
  for (c in all_clusters) {
    message('.', appendLF=FALSE)
    gene_list 			= factor( as.integer(clusters_dt$clust == c) )
    names(gene_list) 	= clusters_dt$symbol
    
    # make topGO object
    suppressMessages({
      topGO_data 	= new("topGOdata", 
                        description 		= c, 
                        allGenes 			= gene_list, 
                        # geneSelectionFun 	= function(x) {x==TRUE},
                        annot 				= annFUN.org, 
                        mapping 			= org_mapping, 
                        ontology 			= 'BP',
                        ID 					= 'symbol'
      )
    })
    # run enrichment tests on these, extract results
    suppressMessages({go_weight 	= runTest(topGO_data, algorithm = "weight", statistic = "fisher")})
    n_terms 		= length(go_weight@score)
    temp_results 	= data.table(GenTable(topGO_data, 
                                        p_go 	= go_weight, 
                                        orderBy 		= 'p_go', 
                                        ranksOf 		= 'p_go', 
                                        topNodes 		= n_terms
    ))
    temp_results[ , cluster := c ]
    
    # store
    go_results 	= rbind(go_results, temp_results)
  }
  message('')
  
  # tidy up
  setnames(go_results, 'p_go', 'tmp')
  go_results[, p_go := as.numeric(tmp) ]
  go_results[ tmp == '< 1e-30', p_go := 9e-31 ]
  go_results[ , tmp := NULL ]
  go_results[ , cluster := factor(cluster, levels=all_clusters) ]
  
  return(go_results)
}
.make_plot_dt <- function(x_data, hclust_obj, proj_dt, clusters_dt) {
  # plot
  plot_dt 			= data.table(melt(x_data, varnames=c('cell_id', 'symbol')))
  
  # nice ordering
  symbol_order 		= colnames(x_data)[hclust_obj$order]
  plot_dt[, symbol 	:= factor(symbol, levels=symbol_order)]
  plot_dt[, cell_id 	:= factor(cell_id, levels=proj_dt$cell_id)]
  
  # put this into plotting 
  plot_dt 			= clusters_dt[ plot_dt, on='symbol' ]
  plot_dt 			= proj_dt[ plot_dt, on='cell_id' ]
  
  return(plot_dt)
}
plot_go_results_homemade <- function(go_list, sig_cutoff=5, p_cutoff=0.1) {
  # unpack
  go_results 	= go_list$go_results
  
  # set up
  plot_dt 	= go_results[ Significant>=sig_cutoff & p_go<p_cutoff ]
  setorder(plot_dt, cluster, -p_go)
  plot_dt[, N := .N, by=Term ]
  plot_dt[	  , term_n := Term ]
  plot_dt[ N > 1, term_n := paste0(Term, '_', 1:.N), by=Term ]
  plot_dt[, term_n := factor(term_n, levels=plot_dt$term_n) ]
  
  # plot
  g = ggplot(plot_dt) +
    aes( x=term_n, y=-log10(p_go) ) +
    geom_col() +
    scale_y_continuous( breaks=pretty_breaks() ) +
    facet_grid( cluster ~ ., scales='free_y', space='free_y') +
    coord_flip() +
    labs(
      x 	= NULL
      ,y 	= '-log10( p-value )'
    ) 
  return(g)
}
plot_heatmap_of_gene_clusters_homemade <- function(go_list) {
  # unpack
  plot_dt 	= go_list$plot_dt
  
  # plot
  g = ggplot(plot_dt) +
    aes( x=cell_id, y=symbol, fill=value ) +
    geom_tile() +
    scale_fill_distiller( palette='RdBu', limits=c(-3, 3) ) +
    facet_grid( clust_label ~ ., scale='free_y', space='free_y' ) +
    theme_bw() +
    theme(
      axis.text 	 = element_blank()
      ,axis.ticks  = element_blank()
    ) +
    labs(
      x 		= 'Cell'
      ,y 		= 'Symbol'
      ,fill 	= 'z-scored gene\nexpression'
    )
  return(g)
}
plot_profiles_of_gene_clusters_homemade <- function(go_list, label_name='Embryonic stage (day)', palette='RdBu') {
  # unpack
  plot_dt 	= go_list$plot_dt
  cuts_dt 	= go_list$cuts_dt
  
  # set up what to plot
  means_dt 	= plot_dt[, list(value=mean(value)), by=list(psuper, clust_label)]
  
  # make nice colours
  col_vals 	= .make_col_vals(cuts_dt$label_input, palette)
  
  # plot
  g = ggplot(means_dt) +
    geom_vline(data=cuts_dt, aes(xintercept=psuper, colour=label_input)) +
    scale_colour_manual( values=col_vals ) +
    geom_smooth( colour='black', span=0.2, method='loess', aes( x=psuper, y=value ) )
  n_cells 	= length(unique(plot_dt$cell_id))
  if ( n_cells<=2000 ) {
    rug_dt 		= unique(plot_dt[, list(psuper, cell_id)])
    g = g + geom_rug(data=rug_dt, sides='b', alpha=0.1, aes(x=psuper) )
  }
  g = g + facet_grid( clust_label ~ ., scales='free_y' ) +
    theme_bw() +
    theme(
      axis.text 	= element_blank()
    ) +
    labs(
      x 		= 'Cell (psudotime score)'
      ,y 		= 'z-scored gene expression'
      ,colour = label_name
    )
  return(g)
}
#homemade double plot function
plot_double_psupertime_histogram <- function(double_obj=NULL, psuper_1=NULL, psuper_2=NULL, run_names=NULL, process=FALSE) {
  # check inputs
  if (is.null(double_obj)) {
    if ( is.null(psuper_1) | is.null(psuper_2) ) {
      stop('either a double_obj must be given, or psuper_1 and psuper_2 must both be given')
    } else {
      double_obj 		= double_psupertime(psuper_1, psuper_2, run_names, process)
    }
  }
  # unpack
  run_names 		= double_obj$run_names
  label_x 		= run_names[[1]]
  label_y 		= run_names[[2]]
  doubles_wide 	= double_obj$doubles_wide
  # make colours
  col_vals 		= .make_col_vals(doubles_wide$label_input)
  # add facet labels
  plot_dt 		= copy(doubles_wide)
  plot_dt[, input_label := paste0('Input data: ', input) ]
  # do some plotting
  g = ggplot(plot_dt) +
    aes_string(
      x 			= paste0('psuper_', label_x)
      ,fill   = paste0('label_psuper_', label_x) 
      ,colour	= paste0('label_psuper_', label_x)
    ) +
    geom_density( alpha=0.5 ) +
    scale_fill_manual(  values = brewer.pal(8, "Blues") ) +
    scale_colour_manual( values = rev(brewer.pal(8, "Blues")) ) +
    facet_grid(input_label ~ .) +
    theme_bw() +
    guides(color=guide_legend("Predicted identity")) + 
    theme_bw() +
    labs(
        x 			= paste0('Psudotime (trained on mouse cells)')
        ,y 			= paste0('Density')
      )+
    theme(legend.position="none")
  return(g)
}
plot_double_psupertime_xy <- function(double_obj=NULL, psuper_1="mouse", psuper_2="rat", run_names=NULL, process=FALSE) {
  # check inputs
  if (is.null(double_obj)) {
    if ( is.null(psuper_1) | is.null(psuper_2) ) {
      stop('either a double_obj must be given, or psuper_1 and psuper_2 must both be given')
    } else {
      double_obj 		= double_psupertime(psuper_1, psuper_2, run_names, process)
    }
  }
  # unpack
  run_names 		= double_obj$run_names
  label_x 		= run_names[[1]]
  label_y 		= run_names[[2]]
  doubles_wide 	= double_obj$doubles_wide
  # make colours
  col_vals 		= .make_col_vals(doubles_wide$label_input)
  # add facet labels
  plot_dt 		= copy(doubles_wide)
  plot_dt[, input_label := paste0('Input data: ', input) ]
  # do some plotting
  g = ggplot(plot_dt) +
    aes_string(
      x 			= paste0('psuper_', label_x)
      ,y 			= paste0('psuper_', label_y)
      ,colour 	= paste0('label_input')
    ) +
    geom_point(alpha = 0.25) +
    scale_colour_manual( values=col_vals ) +
    facet_grid( . ~ input_label) +
    theme_dark() +
    labs(
      x 			= paste0('Psudotime (trained on mouse cells)')
      ,y 			= paste0('Psudotime (trained on rat cells)')
    )+
    theme(legend.position="none")
  return(g)
}
#Plot contour
plot_double_psupertime_contour <- function(double_obj=NULL, psuper_1=NULL, psuper_2=NULL, run_names=NULL) {
  # check run_names
  if ( is.null(run_names) ) {
    run_names 	= c('1','2')
    message('using default values for run_names:', paste(run_names, sep=', '))
  } else {
    if ( !is.character(run_names) | length(unique(run_names))!=2 ) {
      stop('run_names must be character vector of length two with no repeated values')
    }
  }
  # check inputs
  if (is.null(double_obj)) {
    if ( is.null(psuper_1) | is.null(psuper_2) ) {
      stop('either a double_obj must be given, or psuper_1 and psuper_2 must both be given')
    } else {
      double_obj 		= double_psupertime(psuper_1, psuper_2, run_names)
    }
  }
  
  # unpack
  run_names 		= double_obj$run_names
  label_x 		= run_names[[1]]
  label_y 		= run_names[[2]]
  doubles_wide 	= double_obj$doubles_wide
  
  # do some plotting
  g = ggplot(doubles_wide) +
    aes_string(
      x 			= paste0('psuper_', label_x)
      ,y 			= paste0('psuper_', label_y)
      ,colour 	= 'input'
    ) +
    geom_density2d() +
    scale_colour_brewer( palette='Set1' ) +
    theme_bw() +
    labs(
      x 			= paste0('Scoring model trained on mouse cells')
      ,y 			= paste0('Scoring model trained on rat cells')
      ,colour 	= 'Input\ndata'
    )
  
  return(g)
}


psupertime::plot_double_psupertime_confusion(mutualprojection)
psupertime::plot_double_psupertime_contour(mutualprojection)
plot_double_psupertime_contour(mutualprojection)
psupertime::plot_double_psupertime_genes(scoremodel_v2, scoremodel_v3)
psupertime::project_onto_psupertime(scoremodel)

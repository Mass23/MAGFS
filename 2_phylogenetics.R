library(mgcv)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(matrixStats)
library(reshape2)
library(gratia)
library(ggpubr)
library(ape)
library(poolr)
library(foreach)
library(doMC)

set.seed(23)
setwd('~/Documents/PhD/MAGFS')

GetPhyloRandomMeans <- function(all_data, tree, resp_var, covariate){
    colnames(all_data)[colnames(all_data) == resp_var] = 'resp_var'
    colnames(all_data)[colnames(all_data) == covariate] = 'covariate'
    all_data = all_data %>% select(MAG, sample, rel_abundance, 
                                   resp_var, covariate, latitude, longitude, 
                                   Category, p, c, o, f, g, glacier)
    
    coph_dists = cophenetic.phylo(tree)
    hclust = hclust(as.dist(coph_dists), method='average')
    
    phylo_means = data.frame() 
    for (h_step in round(seq(0.00, 1.00, by = 0.05), digits = 3)){
      print(h_step)
      clusters = cutree(hclust, h=h_step*max(hclust$height))
      all_data$cluster = map_chr(all_data$MAG, function(x) as.character(clusters[paste0(c(x, 'fa'), collapse = '.')]))
      data_glom = all_data %>% filter(!is.na(cluster)) %>% 
        select(sample, glacier, latitude, longitude, covariate, rel_abundance, resp_var, cluster)
      
      phylo_i_means = data.frame()
      for (i in 1:20){
        set.seed(i)
        data_glom = data_glom %>% group_by(cluster, sample) %>% mutate(rel_abundance = sample(rel_abundance)) %>% ungroup()
        data_i = data_glom %>% group_by(sample, glacier, latitude, longitude, covariate) %>% 
          summarise(mean_resp = weighted.mean(resp_var, rel_abundance))

        i_model = bam(data = data_i, formula = mean_resp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + covariate, 
                         method = 'REML', family = gaussian())
        coef = summary(i_model)$p.coeff['covariate']
        p_val =  summary(i_model)$p.pv['covariate']
        
        i_df = data.frame(h=h_step, i=i, coef=coef, pval=p_val, covariate=covariate, resp_var=resp_var, n_clu=length(unique(clusters)))
        phylo_i_means = rbind(phylo_i_means, i_df)}
      phylo_means = rbind(phylo_means, phylo_i_means)}
  return(phylo_means)}
    
SummarisePhyloRandomMeans <- function(phylo_means){
    phylo_summ = phylo_means %>% group_by(h) %>% summarise(mean_coef = mean(coef),
                                                           sd_coef = sd(coef),
                                                           p = stouffer(pval)$p)
    phylo_summ$padj = p.adjust(phylo_summ$p, method = 'holm')
  return(phylo_summ)}  
  
TaxonomyPhyloDepth <- function(all_data, tree){
  coph_dists = cophenetic.phylo(tree)
  # Within genera
  within_genus = c()
  for (genus in unique(all_data$g)){
    if (genus != 'g__'){
    genus_mags = all_data %>% filter(g == genus) %>% pull(MAG)  %>% unique()
    if (length(genus_mags) > 4){
    genus_mags = map_chr(genus_mags, function(x) paste0(c(x, '.fa'), collapse = ''))
    dists = melt(coph_dists[rownames(coph_dists) %in% genus_mags, colnames(coph_dists) %in% genus_mags])
    dist = dists[dists$Var1 != dists$Var2,] %>% pull(value) %>% median(na.rm=T)
    within_genus = c(within_genus, dist)}}}
  
  # Within families
  within_family = c()
  for (family in unique(all_data$f)){
    if (family != 'f__'){
    family_mags = all_data %>% filter(f == family) %>% pull(MAG)  %>% unique()
    if (length(family_mags) > 4){
    family_mags = map_chr(family_mags, function(x) paste0(c(x, '.fa'), collapse = ''))
    dists = melt(coph_dists[rownames(coph_dists) %in% family_mags, colnames(coph_dists) %in% family_mags])
    dist = dists[dists$Var1 != dists$Var2,] %>% pull(value) %>% median()
    within_family = c(within_family, dist)}}}
  
  # Within orders
  within_order = c()
  for (order in unique(all_data$o)){
    if (order != 'o__'){
    order_mags = all_data %>% filter(o == order) %>% pull(MAG)  %>% unique()
    if (length(order_mags) > 4){
    order_mags = map_chr(order_mags, function(x) paste0(c(x, '.fa'), collapse = ''))
    dists = melt(coph_dists[rownames(coph_dists) %in% order_mags, colnames(coph_dists) %in% order_mags])
    dist = dists[dists$Var1 != dists$Var2,] %>% pull(value) %>% median()
    within_order = c(within_order, dist)}}}
  
  # Within classes
  within_class = c()
  for (class in unique(all_data$c)){
    if (order != 'o__'){
      class_mags = all_data %>% filter(c == class) %>% pull(MAG)  %>% unique()
      if (length(class_mags) > 4){
        class_mags = map_chr(class_mags, function(x) paste0(c(x, '.fa'), collapse = ''))
        dists = melt(coph_dists[rownames(coph_dists) %in% class_mags, colnames(coph_dists) %in% class_mags])
        dist = dists[dists$Var1 != dists$Var2,] %>% pull(value) %>% median()
        within_class = c(within_class, dist)}}}
  return(list(genus=within_genus, family=within_family, order=within_order, class=within_class))}

PhyloLOCO <- function(all_data, tree, h_value, resp_var, covariate){
    phylo_loco = data.frame() 
    coph_dists = cophenetic.phylo(tree)
    hclust = hclust(as.dist(coph_dists), method='average')
    
    clusters = cutree(hclust, h=h_value*max(hclust$height))
    all_data$cluster = map_chr(all_data$MAG, function(x) as.character(clusters[paste0(c(x, 'fa'), collapse = '.')]))
    colnames(all_data)[colnames(all_data) == resp_var] = 'resp_var'
    colnames(all_data)[colnames(all_data) == covariate] = 'covariate'
    
    for (clu in unique(clusters)){
      print(clu)
      phylo_clu_means = data.frame()
      loco_data = all_data %>% filter(cluster != clu)
      
      for (i in 1:20){
        set.seed(i)
        data_glom = loco_data %>% select(sample, glacier, latitude, longitude, 
                                         resp_var, rel_abundance, covariate, cluster) %>% 
          group_by(cluster, sample) %>% mutate(rel_abundance = sample(rel_abundance)) %>% ungroup()
        
        data_i = data_glom %>% group_by(sample, glacier, latitude, longitude, covariate) %>% 
          summarise(mean_resp = weighted.mean(resp_var, rel_abundance))

        i_model = bam(data = data_i, formula = mean_resp ~ s(latitude, longitude, bs='sos', k=-1) + covariate, 
                        method = 'REML', family = gaussian())
        coef_cov = as.numeric(summary(i_model)$p.coeff['covariate'])

        phylo_clu_means = rbind(phylo_clu_means, data.frame('h'=h_value, 'i'=i, 'coef'=coef_cov, 'loco_cluster'=clu, 'resp_var'=resp_var, 'covariate'=covariate))}
      phylo_loco = rbind(phylo_loco, phylo_clu_means)}
    return(phylo_loco)}

PermutationClusterTest <- function(phylo_loco, phylo_null_means, h_value, tree){
  # check sign clusters
  significant_clusters = data.frame()
  for (clu in unique(phylo_loco$loco_cluster)){
    loco_values = as.vector(phylo_loco$coef[phylo_loco$loco_cluster == clu])
    null_values = as.vector(phylo_null_means$coef[phylo_null_means$h == h_value])

    wt = wilcox.test(loco_values, null_values)
    median_effect = (median(null_values) - median(loco_values)) / median(null_values)

    significant_clusters = rbind(significant_clusters, data.frame('cluster' = clu, 'p' = wt$p.value, 'median_effect' = median_effect))}
  significant_clusters$padj = p.adjust(significant_clusters$p, method = 'bonferroni')
  return(significant_clusters)}

################################################################################################################################################################################################

all_data = read.csv('Data/all_data.csv')
tree = read.tree('Data/treeBacteria.tree')

# Redundancy index
phylo_means_red_gldist = GetPhyloRandomMeans(all_data, tree, 'redundancy_index', 'gl_dist')
summ_phylo_red_gldist = SummarisePhyloRandomMeans(phylo_means_red_gldist)
summ_phylo_red_gldist %>% filter(padj < 0.05) %>% arrange(-h)
summ_phylo_red_gldist$resp_var = 'Redundancy index'
summ_phylo_red_gldist$covariate = 'Distance to the glacier'

phylo_means_red_glindex = GetPhyloRandomMeans(all_data, tree, 'redundancy_index', 'gl_index')
summ_phylo_red_glindex = SummarisePhyloRandomMeans(phylo_means_red_glindex)
summ_phylo_red_glindex %>% filter(padj < 0.05) %>% arrange(-h)
summ_phylo_red_glindex$resp_var = 'Redundancy index'
summ_phylo_red_glindex$covariate = 'Glacier index'

phylo_means_red_watemp = GetPhyloRandomMeans(all_data, tree, 'redundancy_index', 'water_temp')
summ_phylo_red_watemp = SummarisePhyloRandomMeans(phylo_means_red_watemp)
summ_phylo_red_watemp %>% filter(padj < 0.05) %>% arrange(-h)
summ_phylo_red_watemp$resp_var = 'Redundancy index'
summ_phylo_red_watemp$covariate = 'Water temperature'

# Chlorophyll-a
phylo_means_trnas_chla = GetPhyloRandomMeans(all_data, tree, 'norm_tRNAs', 'chla')
summ_phylo_trnas_chla = SummarisePhyloRandomMeans(phylo_means_trnas_chla)
summ_phylo_trnas_chla %>% filter(padj < 0.05) %>% arrange(-h)
summ_phylo_trnas_chla$resp_var = 'tRNA number'
summ_phylo_trnas_chla$covariate = 'Chlorophyll-a'

phylo_means_size_chla = GetPhyloRandomMeans(all_data, tree, 'norm_size', 'chla')
summ_phylo_size_chla = SummarisePhyloRandomMeans(phylo_means_size_chla)
summ_phylo_size_chla %>% filter(padj < 0.05) %>% arrange(-h)
summ_phylo_size_chla$resp_var = 'Genome size'
summ_phylo_size_chla$covariate = 'Chlorophyll-a'

phylo_means_gene_chla = GetPhyloRandomMeans(all_data, tree, 'norm_gene_number', 'chla')
summ_phylo_gene_chla = SummarisePhyloRandomMeans(phylo_means_gene_chla)
summ_phylo_gene_chla %>% filter(padj < 0.05) %>% arrange(-h)
summ_phylo_gene_chla$resp_var = 'Gene number'
summ_phylo_gene_chla$covariate = 'Chlorophyll-a'

taxa_depth = TaxonomyPhyloDepth(all_data, tree)

chla_phylo_summ = rbind(summ_phylo_gene_chla, summ_phylo_trnas_chla, summ_phylo_size_chla)
chla_phylo_summ = chla_phylo_summ %>% group_by(resp_var, covariate) %>% mutate(max_coef = max(abs(mean_coef))) %>% ungroup() %>% mutate(mean_effect = abs(mean_coef)/max_coef,
                                                                                                                                        mean_effect_sd = abs(sd_coef)/max_coef)
chla_phylo_summ$formula = paste0(chla_phylo_summ$resp_var, ' ~ ', chla_phylo_summ$covariate)

chla_phylo_summ$formula <- factor(chla_phylo_summ$formula, levels=c('tRNA number ~ Chlorophyll-a',
                                                                    'Gene number ~ Chlorophyll-a',
                                                                    'Genome size ~ Chlorophyll-a'))

p1 = ggplot(chla_phylo_summ) + 
  geom_hline(yintercept = 0, colour='dimgrey', linetype='dashed') +
  geom_vline(xintercept = median(taxa_depth$genus), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Genus', x=0.25, y=-0.3), colour='black', size=4) +
  geom_vline(xintercept = median(taxa_depth$family), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Family', x=0.46, y=-0.3), colour='black', size=4) +
  geom_vline(xintercept = median(taxa_depth$order), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Order', x=0.58, y=-0.3), colour='black', size=4) +
  geom_vline(xintercept = median(taxa_depth$class), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Class', x=0.76, y=-0.3), colour='black', size=4) +
  geom_ribbon(mapping = aes(x=h, ymin=mean_effect-mean_effect_sd, ymax=mean_effect+mean_effect_sd, fill=formula), alpha=0.3) +
  geom_line(mapping = aes(x=h, y=mean_effect, colour=formula), size=2) + theme_linedraw() + xlim(0,1) + xlab('') + ylab('Signal') +
  theme(panel.grid = element_line(colour='grey'), legend.position = c(0.8, 0.85), legend.key.size = unit(0.5, 'cm'), legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0)) + 
  guides(fill=guide_legend(title=''), colour=guide_legend(title='')) + 
  scale_colour_manual(values=c('#007BE0','#004F8F','#00223D'))  + 
  scale_fill_manual(values=c('#007BE0','#004F8F','#00223D'))


red_phylo_summ = rbind(summ_phylo_red_gldist, summ_phylo_red_watemp, summ_phylo_red_glindex)
red_phylo_summ = red_phylo_summ %>% group_by(resp_var, covariate) %>% mutate(max_coef = max(abs(mean_coef))) %>% ungroup() %>% mutate(mean_effect = abs(mean_coef)/max_coef,
                                                                                                                                      mean_effect_sd = abs(sd_coef)/max_coef)
                              
red_phylo_summ$formula = paste0(red_phylo_summ$resp_var, ' ~ ', red_phylo_summ$covariate)

red_phylo_summ$formula <- factor(red_phylo_summ$formula, levels=c('Redundancy index ~ Distance to the glacier',
                                                                  'Redundancy index ~ Glacier index',
                                                                  'Redundancy index ~ Water temperature'))

p2 = ggplot(red_phylo_summ) + 
  geom_hline(yintercept = 0, colour='dimgrey', linetype='dashed') +
  geom_vline(xintercept = median(taxa_depth$genus), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Genus', x=0.25, y=0.9), colour='black', size=4) +
  geom_vline(xintercept = median(taxa_depth$family), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Family', x=0.46, y=0.9), colour='black', size=4) +
  geom_vline(xintercept = median(taxa_depth$order), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Order', x=0.58, y=0.9), colour='black', size=4) +
  geom_vline(xintercept = median(taxa_depth$class), colour='dimgrey', linewidth=0.75) + 
  geom_text(aes(label='Class', x=0.76, y=0.9), colour='black', size=4) +
  geom_ribbon(mapping = aes(x=h, ymin=mean_effect-mean_effect_sd, ymax=mean_effect+mean_effect_sd, fill=formula), alpha=0.3) +
  geom_line(mapping = aes(x=h, y=mean_effect, colour=formula), size=2) + theme_linedraw() + xlim(0,1) + xlab('Phylogenetic height') + ylab('Signal') +
  theme(panel.grid = element_line(colour='grey'), legend.position = c(0.75, 0.7), legend.key.size = unit(0.5, 'cm'), legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0)) + 
  scale_colour_manual(values=c('#D66976', '#69B5AC', '#D6A756')) + scale_fill_manual(values=c('#D66976', '#69B5AC', '#D6A756')) + 
  guides(fill=guide_legend(title=''), colour=guide_legend(title=''))


p = ggarrange(p2 , p1, nrow = 2, align = 'v', labels = c('A', 'B'))
ggsave('Figure_2.pdf', p, width = 7, height = 7)

#########################
coph_dists = cophenetic.phylo(tree)
hclust = hclust(as.dist(coph_dists), method='average')
clusters = cutree(hclust, h=0.6*max(hclust$height))
all_data$cluster_pc2 = map_chr(all_data$MAG, function(x) as.character(clusters[paste0(c(x, 'fa'), collapse = '.')]))
all_data$cluster = all_data$cluster_pc2

ploco_size_chla = PhyloLOCO(all_data, tree, 0.600, 'norm_size', 'chla')
test_size_chla = PermutationClusterTest(ploco_size_chla, phylo_means_size_chla, 0.600, tree)
test_size_chla %>% filter(padj < 0.05)
#   cluster            p median_effect         padj
#       1   1.341202e-07     0.4926574 4.425966e-06  c__Gammaproteobacteria
#      30   7.937835e-05     0.2880802 2.619486e-03  "p__Acidobacteriota" "p__Desulfobacterota" "p__Desulfobacterota_B" 
                                                    #"p__Desulfobacterota_E" "p__Desulfobacterota_I" "p__Myxococcota" 
                                                    #"p__Myxococcota_A" "p__Nitrospirota"  

ploco_gene_chla = PhyloLOCO(all_data, tree, 0.600, 'norm_gene_number', 'chla')
test_gene_chla = PermutationClusterTest(ploco_gene_chla, phylo_means_gene_chla, 0.600, tree)
test_gene_chla %>% filter(padj < 0.05)
#   cluster            p median_effect         padj
#         1 6.531902e-08     0.7048170 2.155528e-06  c__Gammaproteobacteria
#        30 5.832835e-06     0.4693787 1.924835e-04  "p__Acidobacteriota" "p__Desulfobacterota" "p__Desulfobacterota_B" 
                                                    #"p__Desulfobacterota_E" "p__Desulfobacterota_I" "p__Myxococcota" 
                                                    #"p__Myxococcota_A" "p__Nitrospirota"  

ploco_trna_chla = PhyloLOCO(all_data, tree, 0.600, 'norm_tRNAs', 'chla')
test_trna_chla = PermutationClusterTest(ploco_trna_chla, phylo_means_trnas_chla, 0.600, tree)
test_trna_chla %>% filter(padj < 0.05)
#   cluster            p median_effect         padj
#         1 5.411816e-09     0.7069186 1.785899e-07 c__Gammaproteobacteria
#         5 2.579738e-04    -0.4238117 8.513136e-03 "p__Bacteroidota"   "p__Fibrobacterota"
#        28 3.364318e-05     0.2401113 1.110225e-03 "p__Cyanobacteria"

######################################################################
clusters = cutree(hclust, h=0.200*max(hclust$height))
all_data$cluster_pc1 = map_chr(all_data$MAG, function(x) as.character(clusters[paste0(c(x, 'fa'), collapse = '.')]))
all_data$cluster = all_data$cluster_pc1
ploco_red_gldist = PhyloLOCO(all_data, tree, 0.2, 'redundancy_index', 'gl_dist')
test_red_gldist = PermutationClusterTest(ploco_red_gldist, phylo_means_red_gldist, 0.200, tree)
test_red_gldist %>% filter(padj < 0.05)
#   cluster            p median_effect         padj. consensus         
#         2 1.192103e-04     0.6561007 4.696886e-02  o__Burkholderiales
#         4 7.937835e-05     0.8690717 3.127507e-02  o__Burkholderiales
#        15 1.192103e-04     0.9337818 4.696886e-02  o__Rhizobiales
#        18 1.135908e-05     0.8691132 4.475478e-03  o__Rhizobiales
#        20 9.104785e-05     0.7587717 3.587285e-02  o__Rhizobiales o__Rhizobiales_A
#        22 7.937835e-05     0.8684029 3.127507e-02  o__Micropepsales
#        45 2.897784e-05     0.7680832 1.141727e-02  o__Palsa-1295
#        49 1.192103e-04     0.6676426 4.696886e-02  o__UBA10030
#        60 1.830634e-05     0.7712549 7.212700e-03  o__Cytophagales
#        62 6.909214e-05     0.6659444 2.722230e-02  o__Cytophagales
#       126 1.564634e-05     0.8650131 6.164659e-03  o__UBA1135
#       133 1.654449e-06     0.9652294 6.518528e-04  o__Chthoniobacterales
#       145 6.003771e-05     0.9942862 2.365486e-02  o__Omnitrophales
#       147 5.208212e-05     1.0221004 2.052036e-02  o__UBA10015
#       153 1.042636e-04     0.9435636 4.107987e-02  o__Treponematales
#       156 1.564634e-05     0.7171809 6.164659e-03  o__Deinococcales
#       190 7.937835e-05     0.8075204 3.127507e-02  o__UBA9983_A
#       232 6.003771e-05     0.7407982 2.365486e-02  o__Saccharimonadales
#       245 4.510181e-05     0.8697875 1.777011e-02  o__Anaerolineales
#       258 3.898887e-05     0.8147188 1.536161e-02  o__IMCC26256
#       264 3.364318e-05     0.7674151 1.325541e-02  o__Acidimicrobiales
#       282 2.402846e-06     0.7458309 9.467214e-04  o__Baltobacterales
#       289 6.909214e-05     0.5865411 2.722230e-02  o__Armatimonadales
#       295 3.898887e-05     0.5445854 1.536161e-02  o__Mor1
#       297 1.564634e-05     0.7305707 6.164659e-03  o__Pyrinomonadales
#       314 9.924080e-09     1.6662588 3.910088e-06  o__Nitrospirales
#       318 1.135908e-05     0.9092941 4.475478e-03  o__Desulfobulbales
#       326 4.909097e-06     0.8926769 1.934184e-03  o__UBA2466
#       329 6.003771e-05     0.6387693 2.365486e-02  o__Bdellovibrionales
#       331 2.137650e-05     0.8894864 8.422340e-03  o__Bdellovibrionales
#       332 1.192103e-04     0.7997320 4.696886e-02  o__Bdellovibrionales
#       348 4.510181e-05     0.6721270 1.777011e-02  o__Oligoflexales
#       354 9.104785e-05     0.7154731 3.587285e-02  o__Polyangiales
#       357 1.042636e-04     0.8949298 4.107987e-02  o__Polyangiales
#       370 4.510181e-05     0.8501970 1.777011e-02  o__Xanthomonadales
#       393 2.491213e-05     0.7596084 9.815377e-03  o__Burkholderiales

ploco_red_watemp = PhyloLOCO(all_data, tree, 0.200, 'redundancy_index', 'water_temp')
test_red_watemp = PermutationClusterTest(ploco_red_watemp, phylo_means_red_watemp, 0.200, tree)
test_red_watemp %>% filter(padj < 0.05)
#   cluster            p median_effect         padj.  consensus
#       314 1.367071e-06     0.5498099 0.000538626    o__Nitrospirales

all_data$cluster = NULL
write.csv(all_data, 'Data/all_data_phylo_cluster_added.csv', quote = F, row.names = F)

test_size_chla$test = 'norm_size ~ chla'
test_gene_chla$test = 'norm_gene ~ chla'
test_trna_chla$test = 'norm_trna ~ chla'
test_red_gldist$test = 'KO redundancy ~ gl_dist'
test_red_watemp$test = 'KO redundancy ~ watemp'
all_tests_chla = rbind(test_size_chla, test_gene_chla, test_trna_chla)
all_tests_red = rbind(test_red_gldist, test_red_watemp)
write.csv(all_tests_chla, file = 'Statistics/LOCO_statistics_chla.csv', quote = F, row.names = F)
write.csv(all_tests_red, file = 'Statistics/LOCO_statistics_red.csv', quote = F, row.names = F)









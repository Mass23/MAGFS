library(mgcv)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(matrixStats)
library(reshape2)
library(gratia)
library(ggpubr)
library(performance)
library(factoextra)
library(ggExtra)

setwd('~/Documents/PhD/MAGFS')

############################
# Community averages
GetCommunityAverages <- function(all_data){
   sample_data = all_data %>% group_by(sample, glacier, mountain_range, latitude, longitude,
                                      water_temp, chla, gl_area, gl_dist, gl_cov, gl_index) %>% 
    summarise(mean_norm_size = weighted.mean(norm_size, rel_abundance),
              mean_GC = weighted.mean(GC, rel_abundance),
              mean_norm_gene_number = weighted.mean(norm_gene_number, rel_abundance),
              mean_norm_tRNAs = weighted.mean(norm_tRNAs, rel_abundance),
              mean_redundancy_index = weighted.mean(redundancy_index, rel_abundance),
              mean_coding_density = weighted.mean(coding_density, rel_abundance)) %>% ungroup()
    return(sample_data)}

GlacierInfluencePCA <- function(comm_avgs){
    pca_data = comm_avgs %>% select(mountain_range, water_temp, chla, gl_area, gl_dist, gl_cov, gl_index)
    pca_res = prcomp(pca_data %>% select(-mountain_range), center = T, scale. = T)
    eigens = get_eigenvalue(pca_res)
    pca_df = as.data.frame(pca_res$x)
    pca_df$mountain_range = pca_data$mountain_range
    pca_arr = as.data.frame(pca_res$rotation)
    rownames(pca_arr) = c('Water temperature', 'Chlorophyll-a', 'Glacier area', 
                          'Distance to the glacier', 'Glacier coverage', 'Glacier index')
    p = ggplot() + geom_point(pca_df, mapping = aes(y=PC1, x=PC2, shape=mountain_range)) +
      geom_segment(pca_arr, mapping = aes(y=0,yend=PC1*4,x=0,xend=PC2*4)) +
      geom_text(pca_arr, mapping = aes(x=PC2*3.2,y=PC1*5,label = rownames(pca_arr))) + 
      theme_linedraw() + 
      xlab(paste0('PC1 (', as.character(round(eigens$variance.percent[1],1)), '%)')) + 
      ylab(paste0('PC2 (', as.character(round(eigens$variance.percent[2],1)), '%)')) + 
      labs(shape='Mountain range') + 
      scale_colour_brewer(palette='Paired') + scale_shape_manual(values=c(19, 0, 17, 3, 8, 7, 4, 15, 6)) + 
      theme(panel.grid.major = element_line(colour='darkgrey'),
            panel.grid.minor = element_line(colour='darkgrey'))
    
    return(p)}

Main <- function(){
  all_data = read.csv('Data/all_data.csv')
  comm_avgs = GetCommunityAverages(all_data)
  p1 = GlacierInfluencePCA(comm_avgs)
  
  p2a = ggplot(comm_avgs, aes(x=mean_norm_size/1000000)) + geom_boxplot() + theme_minimal() + xlab('Norm. genome size (mbp)') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p2b = ggplot(comm_avgs, aes(x=mean_norm_gene_number)) + geom_boxplot() + theme_minimal() + xlab('Norm. gene number (#)') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p2c = ggplot(comm_avgs, aes(x=mean_norm_tRNAs)) + geom_boxplot() + theme_minimal() + xlab('Norm. tRNA number (#)') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p2d = ggplot(comm_avgs, aes(x=mean_GC)) + geom_boxplot() + theme_minimal() + xlab('GC content (%)') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p2e = ggplot(comm_avgs, aes(x=mean_redundancy_index)) + geom_boxplot() + theme_minimal() + xlab('Redundancy index') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p2f = ggplot(comm_avgs, aes(x=mean_coding_density)) + geom_boxplot() + theme_minimal() + xlab('Coding density (%)') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p2 = ggarrange(p2a, p2b, p2c, p2d, p2e, p2f, ncol=1)
  
  models_out = data.frame()
  models = list()
  for (covariate in c('water_temp','chla','gl_area','gl_dist','gl_cov','gl_index')){
    for (response_variable in c('mean_norm_size','mean_GC','mean_norm_gene_number','mean_norm_tRNAs',
                                'mean_redundancy_index','mean_coding_density')){
      mod_data = as.data.frame(comm_avgs)
      form = paste('scale(', response_variable, ') ~', "s(latitude, longitude, bs='sos', k=-1, m=1)", '+', covariate)
      model = bam(data = mod_data, formula = reformulate(form), method = 'REML', family = gaussian())
      models = append(models, list(model))
      #stats
      summ = summary(model)
      coef = summ$p.coeff[2]
      coef_se = summ$se[2]
      t = summ$p.t[2]
      p_val = summ$p.pv[2]*2
      spatial_f = summ$s.table[3]
      spatial_p = summ$s.table[4]*2
      models_out = rbind(models_out, data.frame(covariate=covariate, resp_var=response_variable, t=t,
                                                coef=coef, coef_se=coef_se, p=p_val, spatial_f=spatial_f, spatial_p=spatial_p))}}
  models_out$padj = p.adjust(models_out$p, method = 'holm')
  models_out %>% filter(padj < 0.05)
  write.csv(models_out, file = 'Statistics/PC_models.csv', quote = F, row.names = F)
  
  models_out$covariate[models_out$covariate == 'chla'] = "Chlorophyll-a"
  models_out$covariate[models_out$covariate == 'gl_dist'] = 'Distance to the glacier'
  models_out$covariate[models_out$covariate == 'gl_area'] = 'Glacier area'
  models_out$covariate[models_out$covariate == 'water_temp'] = 'Water temperature'
  models_out$covariate[models_out$covariate == 'gl_index'] = 'Glacier index'
  
  models_out$resp_var[models_out$resp_var == 'mean_redundancy_index'] = 'Redundancy index'
  models_out$resp_var[models_out$resp_var == 'mean_norm_size'] = 'Genome size'
  models_out$resp_var[models_out$resp_var == 'mean_norm_gene_number'] = 'Gene number'
  models_out$resp_var[models_out$resp_var == 'mean_norm_tRNAs'] = 'tRNA number'
  
  models_out$formula = paste0(models_out$resp_var, ' ~ ', models_out$covariate)
  
  models_out$formula <- factor(models_out$formula, levels=c('Redundancy index ~ Distance to the glacier',
                                                            'Redundancy index ~ Glacier index',
                                                            'Redundancy index ~ Water temperature',
                                                            'tRNA number ~ Chlorophyll-a',
                                                            'Gene number ~ Chlorophyll-a',
                                                            'Genome size ~ Chlorophyll-a'))

  p3 = ggplot(models_out %>% filter(padj < 0.05), aes(y=formula, x=coef, colour=formula)) + 
    geom_vline(xintercept = 0, colour='dimgrey', linetype='dashed') + 
    geom_point(position=position_dodge(width=0.9), size=3) + 
    geom_errorbarh(aes(xmin=coef-coef_se, xmax=coef+coef_se), position="dodge", linewidth=1) + 
    theme_minimal() + ylab('') + xlab('Regression coefficient') + 
    theme(legend.position = 'none') + scale_colour_manual(values = c('#D66976', '#69B5AC', '#D6A756','#007BE0','#004F8F','#00223D'))
  sp = ggarrange(p2, p3, ncol=2, nrow=1, labels = c('B', 'C'), widths = c(0.33, 0.66))
  
  p = ggarrange(p1, sp, nrow=2, labels=c('A',''), heights = c(0.55,0.45))
  ggsave('figure_1.pdf', p, width = 8, height = 8.5)
}

Main()




library(mgcv)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(matrixStats)
library(reshape2)
library(gratia)
library(ggpubr)
library(compositions)

setwd('~/Documents/PhD/MAGFS')
all_data = read.csv('Data/all_data_phylo_cluster_added.csv')

strain_data = all_data %>% select(MAG, completeness, contamination, Category,
                                  coding_density, norm_size, norm_gene_number, norm_tRNAs, GC,
                                  ko_per_mbp, ko_per_gene, redundancy_index, rel_abundance, cluster_pc2) %>%
  group_by(MAG, completeness, contamination, Category,
           coding_density, norm_size, norm_gene_number, norm_tRNAs, GC,
           redundancy_index, cluster_pc2) %>% summarise(mean_abundance = mean(rel_abundance), prevalence = mean(rel_abundance > 0))

strain_data$norm_size = strain_data$norm_size / 1000000
strain_data$Clade = 'Others'
strain_data$Clade[strain_data$cluster_pc2 == 1] = 'Gammaproteobacteria'

wtest_df = data.frame()
for (var in c('norm_gene_number', 'norm_size', 'norm_tRNAs', 'coding_density', 'redundancy_index', 'GC')){
  vals_gamma = strain_data %>% filter(Clade == 'Gammaproteobacteria') %>% pull(var)
  vals_other = strain_data %>% filter(Clade != 'Gammaproteobacteria') %>% pull(var)
  wtest = wilcox.test(vals_gamma, vals_other)
  med_gamma = median(vals_gamma)
  med_others = median(vals_other)
  med_diff = med_gamma - med_others
  wtest_df = rbind(wtest_df, data.frame(Variable=var, med_gamma=med_gamma, med_others=med_others, med_diff=med_diff, p=wtest$p.value))}
wtest_df$padj = p.adjust(wtest_df$p, method = 'holm')
wtest_df %>% filter(padj < 0.05)

p1 = ggplot(strain_data, aes(y=Clade, x=norm_tRNAs, colour=Clade)) + geom_boxplot() + xlab('tRNA number [N] ***') + ylab('') + theme_linedraw() + scale_x_log10() +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none')
p2 = ggplot(strain_data, aes(y=Clade, x=coding_density, colour=Clade)) + geom_boxplot() + xlab('Coding density [%] ***') + ylab('') + theme_linedraw() +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none')
p3 = ggplot(strain_data, aes(y=Clade, x=norm_size, colour=Clade)) + geom_boxplot() + xlab('Genome size [mbp]') + ylab('') + theme_linedraw() +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none')
p4 = ggplot(strain_data, aes(y=Clade, x=norm_gene_number, colour=Clade)) + geom_boxplot() + xlab('Gene number [N]') + ylab('') + theme_linedraw() +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none')
p5 = ggplot(strain_data, aes(y=Clade, x=redundancy_index, colour=Clade)) + geom_boxplot() + xlab('Redundancy index ***') + ylab('') + theme_linedraw() +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none')
p6 = ggplot(strain_data, aes(y=Clade, x=GC, colour=Clade)) + geom_boxplot() + xlab('GC content [%] ***') + ylab('') + theme_linedraw() +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none')

sp2 = ggarrange(p3, p4, p1, p2, p5, p6, nrow=6) 

####################################################################################
PredictGamMedian <- function(model, lat, lon, covariate, val){
  new_data = data.frame(latitude=lat, longitude = lon)
  new_data[,covariate] = val
  preds = predict.gam(model, newdata = new_data, se.fit = T)
  return(list(fit=median(preds$fit),se=median(preds$se.fit)))}

ComputeAverages <- function(cluster_data, taxonomy){
  averages = cluster_data %>% group_by(sample, chla, gl_area, latitude, longitude) %>% 
    summarise_at(vars(norm_size, norm_tRNAs, norm_gene_number, redundancy_index, coding_density, GC), funs(weighted.mean(. , w=rel_abundance)))
  averages$sum_rel_ab = map_dbl(averages$sample, function(x) sum(cluster_data$rel_abundance[cluster_data$sample == x]))   
  averages$Cluster = taxonomy
  return(averages)}

ModelAverages <- function(cluster_data, new_data, resp_var, covariate){
  form = paste0(resp_var, " ~ s(latitude, longitude, bs='sos', k=-1) + ", covariate, collapse = ' ')
  model = bam(cluster_data, formula = as.formula(form), method='REML', family=gaussian())
  new_data$pred = map_dbl(new_data[,colnames(new_data) == covariate], function(val) PredictGamMedian(model, cluster_data$latitude, cluster_data$longitude, covariate, val)$fit)
  new_data$se = map_dbl(new_data[,colnames(new_data) == covariate], function(val) PredictGamMedian(model, cluster_data$latitude, cluster_data$longitude, covariate, val)$se)
  new_data$Cluster = unique(cluster_data$Cluster)
  return(list(data=new_data,model=model))}

####################################################################################
table(all_data$p[all_data$cluster_pc2 == 1]) # Proteobacteria only
table(all_data$c[all_data$cluster_pc2 == 1]) # Gammaproteobacteria only (all of them)
table(all_data$o[all_data$cluster_pc2 == 1]) # multiples
table(all_data$c[all_data$cluster_pc2 != 1]) # all except Gammaproteobacteria

#######################################################################################
gamma_data = all_data %>% filter(cluster_pc2 == 1, rel_abundance > 0) %>% select(sample, norm_size, norm_tRNAs, norm_gene_number, coding_density, GC, ko_per_mbp, ko_per_gene, redundancy_index, rel_abundance, chla, gl_area, latitude, longitude) 
other_data = all_data %>% filter(cluster_pc2 != 1, rel_abundance > 0) %>% select(sample, norm_size, norm_tRNAs, norm_gene_number, coding_density, GC, ko_per_mbp, ko_per_gene, redundancy_index, rel_abundance, chla, gl_area, latitude, longitude) 

gamma_avgs = ComputeAverages(gamma_data, 'Gammaproteobacteria') %>% mutate(sum_rel_ab = log10(sum_rel_ab))
other_avgs = ComputeAverages(other_data, 'Others') %>% mutate(sum_rel_ab = log10(sum_rel_ab))

#################################
# Plot relative abundance
new_data = data.frame(chla = seq(-10, 0.5, 0.1))
gamma_plot = ModelAverages(gamma_avgs, new_data, 'sum_rel_ab', 'chla')
other_plot = ModelAverages(other_avgs, new_data, 'sum_rel_ab', 'chla')
plot_data = rbind(gamma_plot$data, other_plot$data)
clu_data = rbind(gamma_avgs, other_avgs)

sp1 = ggplot() + 
  geom_ribbon(plot_data, mapping=aes(x=chla, ymin=pred-se, ymax=pred+se, fill=Cluster), alpha=0.5) + 
  geom_point(clu_data, mapping=aes(x=chla, y=sum_rel_ab, colour=Cluster), alpha = 0.75) + 
  geom_line(plot_data, mapping=aes(x=chla, y=pred, colour=Cluster)) + 
  facet_wrap(~Cluster, nrow = 1, scales = 'fixed') + theme_linedraw() + xlab('ln chlorophyll-a') + ylab('log10 rel. abundance') + 
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), legend.position = 'none') + 
  scale_fill_manual(values=c('#987AB2','#DF844F')) + 
  scale_colour_manual(values=c('#987AB2','#DF844F'))

#################################
# Plot effects
reg_coefs = data.frame()

gamma_plot = ModelAverages(gamma_avgs, new_data, 'norm_tRNAs', 'chla')
other_plot = ModelAverages(other_avgs, new_data, 'norm_tRNAs', 'chla')
reg_coefs = rbind(reg_coefs, data.frame(Clade='Gammaproteobacteria', `Response variable`='tRNA number [N]', coef=summary(gamma_plot$model)$p.coeff['chla'], p=summary(gamma_plot$model)$p.pv['chla'], se=summary(gamma_plot$model)$se['chla']))
reg_coefs = rbind(reg_coefs, data.frame(Clade='Others',              `Response variable`='tRNA number [N]', coef=summary(other_plot$model)$p.coeff['chla'], p=summary(other_plot$model)$p.pv['chla'], se=summary(other_plot$model)$se['chla']))

gamma_plot = ModelAverages(gamma_avgs, new_data, 'norm_gene_number', 'chla')
other_plot = ModelAverages(other_avgs, new_data, 'norm_gene_number', 'chla')
reg_coefs = rbind(reg_coefs, data.frame(Clade='Gammaproteobacteria', `Response variable`='Gene number [N]', coef=summary(gamma_plot$model)$p.coeff['chla'], p=summary(gamma_plot$model)$p.pv['chla'], se=summary(gamma_plot$model)$se['chla']))
reg_coefs = rbind(reg_coefs, data.frame(Clade='Others',              `Response variable`='Gene number [N]', coef=summary(other_plot$model)$p.coeff['chla'], p=summary(other_plot$model)$p.pv['chla'], se=summary(other_plot$model)$se['chla']))

gamma_plot = ModelAverages(gamma_avgs, new_data, 'norm_size', 'chla')
other_plot = ModelAverages(other_avgs, new_data, 'norm_size', 'chla')
reg_coefs = rbind(reg_coefs, data.frame(Clade='Gammaproteobacteria', `Response variable`='Genome size [mbp]', coef=summary(gamma_plot$model)$p.coeff['chla'], p=summary(gamma_plot$model)$p.pv['chla'], se=summary(gamma_plot$model)$se['chla']))
reg_coefs = rbind(reg_coefs, data.frame(Clade='Others',              `Response variable`='Genome size [mbp]', coef=summary(other_plot$model)$p.coeff['chla'], p=summary(other_plot$model)$p.pv['chla'], se=summary(other_plot$model)$se['chla']))

sp3 = ggplot(reg_coefs, aes(x=coef, y=Response.variable, fill=Clade, colour=Clade))  + geom_vline(xintercept = 0, colour='dimgrey', linetype='dashed') + 
  geom_point(position=position_dodge(width=0.9), size=3) + geom_errorbarh(aes(xmin=coef-se, xmax=coef+se), position="dodge", linewidth=1) + facet_wrap(Response.variable~., nrow = 6, ncol = 1, scales = 'free') +
  theme_minimal() + theme(axis.text.y = element_blank()) + ylab('') + xlab('Regression coefficient') + scale_colour_manual(values = c('#987AB2', '#DF844F'))

sp = ggarrange(sp2, sp3, ncol = 2, legend = 'bottom', labels = c('B', 'C'))
p = ggarrange(sp1, sp, nrow=2, labels = c('A',''), heights = c(0.4, 0.6))
ggsave('figure_3.pdf', p, width=6.5, height = 8)


#################################
dir.create('Data/Pangenome')
alpha_mags = all_data %>% filter(c == 'c__Alphaproteobacteria', completeness >= 0.9, contamination <= 0.1) %>% mutate(MAG=paste0('Data/MAGs/',MAG,'.faa')) %>% pull(MAG) %>% unique()
gamma_mags = all_data %>% filter(cluster_pc2 == 1, completeness >= 0.9, contamination <= 0.1) %>% mutate(MAG=paste0('Data/MAGs/',MAG,'.faa')) %>% pull(MAG) %>% unique()

comp_motupan = all_data %>% select(MAG, completeness, contamination) %>% distinct() %>% 
  mutate(`Bin Id`= MAG, Completeness=completeness*100, Contamination=contamination*100) %>% select(`Bin Id`, Completeness, Contamination)

dir.create('Data/Pangenome/all_pangenome')
file.copy(c(alpha_mags, gamma_mags), 'Data/Pangenome/all_pangenome')
write.table(comp_motupan, file='Data/Pangenome/all_pangenome/completeness.tsv', sep='\t', quote = F, row.names = F)

dir.create('Data/Pangenome/gamma_pangenome')
file.copy(gamma_mags, 'Data/Pangenome/gamma_pangenome')
write.table(comp_motupan, file='Data/Pangenome/gamma_pangenome/completeness.tsv', sep='\t', quote = F, row.names = F)

# Run in command line:
# mOTUpan.py --output ALL --faas all_pangenome/*.faa --checkm all_pangenome/completeness.tsv --seed 90 --name ALL --threads 6
# mOTUpan.py --output GAMMA --faas gamma_pangenome/*.faa --checkm gamma_pangenome/completeness.tsv --seed 90 --name GAMMA --threads 6






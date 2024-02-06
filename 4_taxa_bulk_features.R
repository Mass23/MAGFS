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
library(performance)

setwd('~/Documents/PhD/MAGFS')
all_data = read.csv('Data/all_data_phylo_cluster_added.csv')

strain_data = all_data %>% select(MAG, completeness, contamination, Category,
                                  coding_density, norm_size, norm_gene_number, norm_tRNAs, GC,
                                  ko_per_mbp, ko_per_gene, redundancy_index, rel_abundance, cluster_pc2) %>%
  group_by(MAG, completeness, contamination, Category,
           coding_density, norm_size, norm_gene_number, norm_tRNAs, GC,
           ko_per_mbp, ko_per_gene, redundancy_index, cluster_pc2) %>% summarise(mean_abundance = mean(rel_abundance), prevalence = mean(rel_abundance > 0))

strain_data$norm_size = strain_data$norm_size / 1000000
strain_data$Clade = 'Others'
strain_data$Clade[strain_data$cluster_pc2 == 1] = 'Gammaproteobacteria'

#### Eco success & bulk features
hist(strain_data$mean_abundance)
hist(strain_data$prevalence)
strain_data = strain_data %>% filter(mean_abundance > 0)
strain_data$mean_abundance = log10(strain_data$mean_abundance + (min(strain_data$mean_abundance)/2))
strain_data$prevalence = log10(strain_data$prevalence)
hist(strain_data$mean_abundance)
hist(strain_data$prevalence)

wtest_df = data.frame()
for (var in c('mean_abundance', 'prevalence')){
  vals_gamma = strain_data %>% filter(Clade == 'Gammaproteobacteria') %>% pull(var)
  vals_other = strain_data %>% filter(Clade != 'Gammaproteobacteria') %>% pull(var)
  wtest = wilcox.test(vals_gamma, vals_other)
  med_gamma = median(vals_gamma)
  med_others = median(vals_other)
  med_diff = med_gamma - med_others
  wtest_df = rbind(wtest_df, data.frame(Variable=var, med_gamma=med_gamma, med_others=med_others, med_diff=med_diff, p=wtest$p.value))}
wtest_df$padj = p.adjust(wtest_df$p, method = 'holm')
wtest_df %>% filter(padj < 0.05)


p1 = ggplot(strain_data, aes(y=Clade, x=mean_abundance, colour=Clade)) + geom_boxplot() + xlab('') + ylab('') + theme_linedraw() + xlab(bquote(log[10]~mean~relative~abundance)) +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none') + xlim(-6.75, -1.5)
p2 = ggplot(strain_data, aes(y=Clade, x=prevalence, colour=Clade)) + geom_boxplot() + xlab('') + ylab('') + theme_linedraw() + xlab(bquote(log[10]~prevalence)) +
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + theme(legend.position = 'none') + xlim(-2.25, 0.25)


####################################################################
m0_size_ab = bam(data = strain_data, formula = norm_size ~ s(mean_abundance, k=3, bs='cs'))
m1_size_ab = bam(data = strain_data, formula = norm_size ~ Clade + s(mean_abundance, k=3, bs='cs', by=as.factor(Clade)))
test_performance(m0_size_ab, m1_size_ab)

m0_size_pr = bam(data = strain_data, formula = norm_size ~ s(prevalence, k=5, bs='ts'))
m1_size_pr = bam(data = strain_data, formula = norm_size ~ Clade + s(prevalence, k=5, bs='ts', by = as.factor(Clade)))
test_performance(m0_size_ab, m1_size_ab)

m0_codden_ab = bam(data = strain_data, formula = coding_density ~ s(mean_abundance, k=3, bs='cs'))
m1_codden_ab = bam(data = strain_data, formula = coding_density ~ Clade + s(mean_abundance, k=3, bs='cs', by = as.factor(Clade)))
test_performance(m0_codden_ab, m1_codden_ab)

m0_codden_pr = bam(data = strain_data, formula = coding_density ~ s(prevalence, k=5, bs='ts'))
m1_codden_pr = bam(data = strain_data, formula = coding_density ~ Clade + s(prevalence, k=5, bs='ts', by = as.factor(Clade)))
test_performance(m0_codden_pr, m1_codden_pr)

# Size ~ abundance
new_data = data.frame(mean_abundance = c(seq(-6.5,-1.75, 0.05),seq(-6.5,-1.75, 0.05)),
                      Clade = c(rep('Gammaproteobacteria', 96), rep('Others', 96)))
new_data$norm_size = predict.gam(m1_size_ab, newdata = new_data, se.fit = T)$fit
new_data$norm_size_se = predict.gam(m1_size_ab, newdata = new_data, se.fit = T)$se.fit
p3 = ggplot() + geom_point(strain_data, mapping=aes(x=mean_abundance, y=norm_size, colour=Clade), alpha=0.05) +
  geom_ribbon(new_data, mapping = aes(x=mean_abundance, ymin=norm_size-norm_size_se, ymax=norm_size+norm_size_se, fill=Clade), alpha=0.5) +
  geom_line(new_data, mapping = aes(x=mean_abundance, y=norm_size, colour=Clade), size=1.2) + theme_linedraw() + 
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  xlab('') + ylab('Genome size [mbp]') + guides(colour = guide_legend(override.aes = list(alpha=1))) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + xlim(-6.75, -1.5)

# Size ~ prevalence
new_data = data.frame(prevalence = c(seq(-2.2, 0, 0.02),seq(-2.2, 0, 0.02)),
                      Clade = c(rep('Gammaproteobacteria', 111), rep('Others', 111)))
new_data$norm_size = predict.gam(m1_size_pr, newdata = new_data, se.fit = T)$fit
new_data$norm_size_se = predict.gam(m1_size_pr, newdata = new_data, se.fit = T)$se.fit
p4 = ggplot() + geom_point(strain_data, mapping=aes(x=prevalence, y=norm_size, colour=Clade), alpha=0.05) +
  geom_ribbon(new_data, mapping = aes(x=prevalence, ymin=norm_size-norm_size_se, ymax=norm_size+norm_size_se, fill=Clade), alpha=0.5) +
  geom_line(new_data, mapping = aes(x=prevalence, y=norm_size, colour=Clade), size=1.2) + theme_linedraw() + 
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  xlab('') + ylab('') + guides(colour = guide_legend(override.aes = list(alpha=1))) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F'))  + xlim(-2.25, 0.25)

# kombp ~ abundance
new_data = data.frame(mean_abundance = c(seq(-6.5,-1.75, 0.05),seq(-6.5,-1.75, 0.05)),
                      Clade = c(rep('Gammaproteobacteria', 96), rep('Others', 96)))
new_data$norm_size = predict.gam(m1_codden_ab, newdata = new_data, se.fit = T)$fit
new_data$norm_size_se = predict.gam(m1_codden_ab, newdata = new_data, se.fit = T)$se.fit
p5 = ggplot() + geom_point(strain_data, mapping=aes(x=mean_abundance, y=coding_density, colour=Clade), alpha=0.05) +
  geom_ribbon(new_data, mapping = aes(x=mean_abundance, ymin=norm_size-norm_size_se, ymax=norm_size+norm_size_se, fill=Clade), alpha=0.5) +
  geom_line(new_data, mapping = aes(x=mean_abundance, y=norm_size, colour=Clade), size=1.2) + theme_linedraw() + 
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey')) + 
  xlab(bquote(log[10]~mean~relative~abundance)) + ylab('Coding density [%]') + guides(colour = guide_legend(override.aes = list(alpha=1))) +
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + xlim(-6.75, -1.5)

# kombp ~ prevalence
new_data = data.frame(prevalence = c(seq(-2.2, 0, 0.02),seq(-2.2, 0, 0.02)),
                      Clade = c(rep('Gammaproteobacteria', 111), rep('Others', 111)))
new_data$norm_size = predict.gam(m1_codden_pr, newdata = new_data, se.fit = T)$fit
new_data$norm_size_se = predict.gam(m1_codden_pr, newdata = new_data, se.fit = T)$se.fit
p6 = ggplot() + geom_point(strain_data, mapping=aes(x=prevalence, y=coding_density, colour=Clade), alpha=0.05) +
  geom_ribbon(new_data, mapping = aes(x=prevalence, ymin=norm_size-norm_size_se, ymax=norm_size+norm_size_se, fill=Clade), alpha=0.5) +
  geom_line(new_data, mapping = aes(x=prevalence, y=norm_size, colour=Clade), size=1.2) + theme_linedraw() + 
  theme(panel.grid.major = element_line(colour='grey'), panel.grid.minor = element_line(colour='grey')) + 
  xlab(bquote(log[10]~prevalence)) + ylab('') + guides(colour = guide_legend(override.aes = list(alpha=1))) + 
  scale_fill_manual(values=c('#987AB2', '#DF844F')) + 
  scale_colour_manual(values=c('#987AB2', '#DF844F')) + xlim(-2.25, 0.25)


p = ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, common.legend = T, legend.grob = get_legend(p1), legend = 'bottom', align='v',
                heights = c(0.3,0.35,0.35), labels = c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave('figure_4.pdf', p, width = 6, height = 7.3)



















##################################################################
BinsAnalysis <- function(str_data, bins_var, clade_test, resp_var){
  quantiles = quantile(str_data %>% pull(bins_var), probs=seq(0,1,0.25))
  out_df = data.frame()
  for (bin in 2:5){
    min_val = quantiles[bin-1]
    max_val = quantiles[bin]
    bin_data_with = str_data %>% filter((!!sym(bins_var)) > min_val, (!!sym(bins_var)) <= max_val) %>% pull(resp_var) 
    bin_data_without = str_data %>% filter((!!sym(bins_var)) > min_val, (!!sym(bins_var)) <= max_val, Clade != clade_test) %>% pull(resp_var) 
    wtest = wilcox.test(bin_data_with, bin_data_without)
    median_difference = median(bin_data_with) - median(bin_data_without)
    out_df = rbind(out_df, data.frame(bins_var=bins_var, bin=bin-1, resp_var, p=wtest$p.value, median_diff=median_difference))}
  return(out_df)}

BinsAnalysis(strain_data, 'mean_abundance', 'Gammaproteobacteria', 'norm_size')
BinsAnalysis(strain_data, 'prevalence', 'Gammaproteobacteria', 'norm_size')
BinsAnalysis(strain_data, 'mean_abundance', 'Gammaproteobacteria', 'ko_per_mbp')
BinsAnalysis(strain_data, 'prevalence', 'Gammaproteobacteria', 'ko_per_mbp')


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

LoadStrainData <- function(){
  tax_data = read.csv('Data/NOMIS_MAGS_tax.tsv', sep='\t') %>% filter(d == 'd__Bacteria') %>% mutate(MAG = MAGs) %>% select(-MAGs)
  
  # Extract genome size completeness contamination information
  comp_data = read.csv('Data/checkm2_all.tsv', sep='\t') %>% mutate(MAG = gsub(pattern = ".fasta", "", Name)) %>% filter(MAG %in% tax_data$MAG)
  size_data = read.table('Data/NOMIS_MAGs_statsTable.txt', sep='\t', header = T) %>% mutate(MAG = MAGs) %>% filter(MAG %in% tax_data$MAG)
  tax_data$completeness = map_dbl(tax_data$MAG, function(x) comp_data$Completeness[comp_data$MAG == x]/100)
  tax_data$contamination = map_dbl(tax_data$MAG, function(x) comp_data$Contamination[comp_data$MAG == x]/100)
  tax_data$size = map_dbl(tax_data$MAG, function(x) size_data$Length[size_data$MAG == x])
  tax_data$GC = map_dbl(tax_data$MAG, function(x) size_data$GC[size_data$MAG == x])
  tax_data$gene_number = map_dbl(tax_data$MAG, function(x) size_data$CDSs[size_data$MAG == x])
  tax_data$coding_density = map_dbl(tax_data$MAG, function(x) size_data$coding.density[size_data$MAG == x])
  tax_data$tRNAs = map_dbl(tax_data$MAG, function(x) size_data$tRNAs[size_data$MAG == x])
  
  tax_data$norm_size = tax_data$size * (1/tax_data$completeness) * (1-tax_data$contamination)
  tax_data$norm_gene_number = tax_data$gene_number * (1/tax_data$completeness) * (1-tax_data$contamination)
  tax_data$norm_tRNAs = tax_data$tRNAs * (1/tax_data$completeness) * (1-tax_data$contamination)
  
  # Specialists/generalist
  spec_data = read.csv('Data/spec_gen.tsv', sep='\t')
  tax_data$Category = map_chr(tax_data$MAG, function(x) ifelse(x %in% spec_data$MAGs, spec_data$sign[spec_data$MAGs == x], 'Non sign.'))
  tax_data$Category[tax_data$Category == 'SPECIALIST'] = 'Specialist'
  tax_data$Category[tax_data$Category == 'GENERALIST'] = 'Generalist'
  tax_data$Category[tax_data$Category == 'NON SIGNIFICANT'] = 'Non sign.'
  
  kegg_data = read.csv('Data/NOMIS_MAGs_Contigs_Genes_KO.txt', sep='\t') %>% mutate(MAG = MAGs) %>% select(-MAGs)
  kegg_data = kegg_data %>% group_by(MAG) %>% summarise(KoN = n(), KoUnique = length(unique(KEGG_ko))) %>% mutate(KoRed = KoN/KoUnique)
  tax_data$ko_number = map_int(tax_data$MAG, function(x) as.integer(kegg_data$KoN[kegg_data$MAG == x]))
  tax_data$ko_unique_number = map_int(tax_data$MAG, function(x) as.integer(kegg_data$KoUnique[kegg_data$MAG == x]))
  tax_data$redundancy_index = tax_data$ko_number/tax_data$ko_unique_number
  return(tax_data)}

log_const <- function(x){return(log(x + (min(x[which(x > 0)])/2)))}

LoadStreamData <- function(){
  stream_data = read.csv('Data/nomis-20230320-0855-db.csv', sep=',') %>% mutate(sample = map_chr(patch, function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_')))
  stream_data = stream_data %>% rename(mountain_range=mountain_range, water_temp=water_temp..C., gl_area=gl_sa..km2., gl_dist=sn_sp_dist..m., 
                                       gl_cov=gl_cov...., chla=chla..ug.g.1., turb=turb..NTU., latitude=lat_sp..DD., longitude=lon_sp..DD.) %>% 
                                select(sample, mountain_range, water_temp, gl_area, gl_dist, gl_cov, chla, latitude, longitude) %>%
                                mutate(gl_index =  sqrt(gl_area) / (gl_dist + sqrt(gl_area))) %>% 
                                group_by(sample, mountain_range) %>% 
                                summarise_all(mean)
 
  # Transform skewed variables
  stream_data$gl_index = log_const(stream_data$gl_index)
  stream_data$gl_area = log_const(stream_data$gl_area)
  stream_data$gl_dist = log_const(stream_data$gl_dist)
  
  stream_data$chla = log_const(stream_data$chla)
  stream_data$water_temp = log_const(stream_data$water_temp)
  
  stream_data$glacier = map_chr(stream_data$sample, function(x) strsplit(as.character(x), '_')[[1]][1])
  return(stream_data)}
  
LoadAbundanceData <- function(x){
  cov_data = read.csv('Data/MAGs_cov_norm.txt', sep='\t')
  rownames(cov_data) = cov_data$MAGs
  cov_data$MAGs = NULL
  cov_data[cov_data < 100] = 0
  cov_data = sweep(cov_data, 2, colSums(cov_data),'/')
  
  colnames(cov_data) = gsub('UpB', 'UP', colnames(cov_data))
  colnames(cov_data) = gsub('Up', 'UP', colnames(cov_data))
  colnames(cov_data) = gsub('DownB', 'DN', colnames(cov_data))
  colnames(cov_data) = gsub('Down', 'DN', colnames(cov_data))
  
  colnames(cov_data)[colnames(cov_data) == 'GL140_1'] = 'GL140_UP'
  colnames(cov_data)[colnames(cov_data) == 'GL140_2'] = 'GL140_UP'
  colnames(cov_data)[colnames(cov_data) == 'GL140_3'] = 'GL140_UP'
  colnames(cov_data)[colnames(cov_data) == 'GL140_4'] = 'GL140_DN'
  colnames(cov_data)[colnames(cov_data) == 'GL140_5'] = 'GL140_DN'
  colnames(cov_data)[colnames(cov_data) == 'GL140_6'] = 'GL140_DN'
  
  cov_data = as.data.frame(t(rowsum(t(cov_data), names(cov_data))/c(table(names(cov_data)))))
  
  all_data = melt(as.matrix(cov_data), varnames = c('MAG', 'sample'))
  all_data = all_data %>% filter(!(sample %in% c("GLR10_GL11_UP_3","GLR12_GL15_UP_2","GLR17_GL16_UP_1","GLR18_GL16_UP_2",
                                                 "GLR19_GL16_UP_3","GLR1_GL5_UP_1","GLR20_GL16_DN_1","GLR21_GL16_DN_2","GL15_Down",
                                                 "GLR22_GL16_DN_3","GLR2_GL5_UP_2","GLR3_GL5_UP_3","GLR56_GL46_UP_1","GLR58_GL46_UP_3",
                                                 "GLR59_GL47_UP_1","GLR68_GL53_UP_2","GLR80_GL56_UP_1","GLR81_GL56_UP_2",
                                                 "GLR8_GL11_UP_1","GLR9_GL11_UP_2","X336R","X337R","X338R","X339R","X340R","X341R")))
  
  
  all_data$rel_abundance = map_dbl(1:nrow(all_data), function(i) ifelse(all_data$sample[i] %in% colnames(cov_data), 
                                                                    cov_data[rownames(cov_data) == all_data$MAG[i], all_data$sample[i]],
                                                                    -999))
  return(all_data %>% select(MAG, sample, rel_abundance))}

Main <- function(){
  strain_data = LoadStrainData()
  stream_data = LoadStreamData()
  abund_data = LoadAbundanceData()
  all_data = purrr::reduce(list(strain_data,abund_data), dplyr::left_join, by = 'MAG')
  all_data = purrr::reduce(list(all_data,stream_data), dplyr::left_join, by = 'sample')
  write.csv(all_data, 'Data/all_data.csv', quote = F, row.names = F)}

Main()

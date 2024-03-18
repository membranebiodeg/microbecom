library


pathabun_unstrat<-read.delim(file="/home/dlaw16/metagenomics/processing/humann/pathabun/allsamples_pathabundance-cpm_unstratified.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
filter_pathabun_unstrat<-pathabun_unstrat[-(1:2),]
colnames(filter_pathabun_unstrat)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_pathabun_unstrat))
metadata_maaslin <- read.delim(file="/home/dlaw16/metagenomics/processing/profiling/functional/maaslin/metadata_maaslin.txt", sep="\t",header = TRUE, stringsAsFactors = FALSE, row.names = 1)
pathway_fit = Maaslin2(input_data     = filter_pathabun_unstrat, 
                       input_metadata = metadata_maaslin, 
                       output         = "pathway_output", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label"),
                       random_effects = c("site"),
                       reference      = c("label,C0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)

#uniref50
filter_pathabun_unstrat50<-pathabun50[-(1:2),]
colnames(filter_pathabun_unstrat50)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_pathabun_unstrat50))
pathway_fit50 = Maaslin2(input_data     = filter_pathabun_unstrat50, 
                       input_metadata = metadata_maaslin, 
                       output         = "pathway_output50", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label"),
                       random_effects = c("site"),
                       reference      = c("label,c0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)

#KO
gene_ko_unstrat<-read.delim(file="/home/dlaw16/metagenomics/processing/humann/gene/allsamples_ko50-cpm_unstratified.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
filter_gene_ko_unstrat<-gene_ko_unstrat[-(1:2),]
colnames(filter_gene_ko_unstrat)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_gene_ko_unstrat))
gene_ko_fit50 = Maaslin2(input_data   = filter_gene_ko_unstrat, 
                       input_metadata = metadata_maaslin, 
                       output         = "ko_output50", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label"),
                       random_effects = c("site"),
                       reference      = c("label,c0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)

#EC
gene_ec_unstrat<-read.delim(file="/home/dlaw16/metagenomics/processing/humann/gene/allsamples_ec50-cpm_unstratified.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
filter_gene_ec_unstrat<-gene_ec_unstrat[-(1:2),]
colnames(filter_gene_ec_unstrat)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_gene_ec_unstrat))
gene_ec_fit50 = Maaslin2(input_data   = filter_gene_ec_unstrat, 
                       input_metadata = metadata_maaslin, 
                       output         = "ec_output50", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label"),
                       random_effects = c("site"),
                       reference      = c("label,c0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)

#GO
gene_go_unstrat<-read.delim(file="/home/dlaw16/metagenomics/processing/humann/gene/allsamples_go50-cpm_unstratified.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
filter_gene_go_unstrat<-gene_go_unstrat[-(1:2),]
colnames(filter_gene_go_unstrat)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_gene_go_unstrat))
gene_go_fit50 = Maaslin2(input_data   = filter_gene_go_unstrat, 
                       input_metadata = metadata_maaslin, 
                       output         = "go_output50", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label"),
                       random_effects = c("site"),
                       reference      = c("label,c0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)

#eggnog
eggnog_unstrat<-read.delim(file="/home/dlaw16/metagenomics/processing/humann/gene/allsamples_eggnog50-cpm_unstratified.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

filter_eggnog_unstrat<-eggnog_unstrat[-(1:2),]

colnames(filter_eggnog_unstrat)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_eggnog_unstrat))

eggnog_fit50 = Maaslin2(input_data   = filter_eggnog_unstrat, 
                       input_metadata = metadata_maaslin, 
                       output         = "eggnog_output50", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label"),
                       random_effects = c("site"),
                       reference      = c("label,c0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)


# Gene families
genefam_unstrat<-read.delim(file="/home/dlaw16/metagenomics/processing/humann/gene/allsamples_rxn-cpm_named_unstratified.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
filter_genefam_unstrat<-genefam_unstrat[-(1:2),]
colnames(filter_genefam_unstrat)=gsub("_merged_trimmed_Abundance.CPM","", colnames(filter_genefam_unstrat))
genefam_fit = Maaslin2(input_data     = filter_genefam_unstrat, 
                       input_metadata = metadata_maaslin, 
                       output         = "gene_output", 
                       normalization  = "NONE", 
                       fixed_effects  = c("label","treated"),
                       random_effects = c("site"),
                       reference      = c("label,c0"),
                       min_abundance  = 0.1,
                       min_prevalence = 0.1)

# Filtering metabolic pathways
sigf_pathway <- read.delim(file="/home/dlaw16/pathway_output50/significant_results.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE);

#to arrange from smallest q values 
top_pathway <- sigf_pathway %>% group_by(value) %>% slice_max(order = -qval, n = 100)
top_pathway <- top_pathway[, !(names(top_pathway) %in% c('metadata','stderr', 'N','N.not.0', 'pval'))]

# reshape your data frame from a long format to a wide format (tidyr package)
wide_top_pathway <- top_pathway %>% pivot_wider(names_from = c(value), values_from = c(value, coef, qval))
wide_top_pathway <- wide_top_pathway[, !(names(wide_top_pathway) %in% c('value_cf', 'value_t1', 'value_t2', 'value_t3','value_t4', 'value_t5'))]

names(wide_top_pathway) <- c('pathway','cf_coef','t1_coef','t2_coef','t3_coef','t4_coef','t5_coef','cf','t1','t2','t3','t4','t5')
wide_top_pathway <- replace(wide_top_pathway, is.na(wide_top_pathway), 0)
wide_top_pathway <- as.data.frame(wide_top_pathway)
row.names(wide_top_pathway)<-wide_top_pathway$pathway

wide_top_pathway <- wide_top_pathway[, !(names(wide_top_pathway) %in% c('pathway'))]

write.table(wide_top_pathway, file = "top_pathway.tsv", sep = "\t", row.names = TRUE)

# plotting in heatmap (top 80)
hmap_pathway <- read.delim(file="/home/dlaw16/pathway_output50/top_pathway_heatmap.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE);

row.names(hmap_pathway) <- hmap_pathway$X
hmap_pathway <- hmap_pathway[, !(names(hmap_pathway) %in% c('X'))]

pheatmap(hmap_pathway, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 25,
         cellheight = 10,
         fontsize_number = 6, 
         number_color = "black",
         #col = brewer.pal(9, "PiYG"),
         #display_numbers = T,
         border_color = "grey60")

# plotting heatmap for chosen enzymes
hmap_ec <- read.delim(file="/home/dlaw16/ec_output50/ec_heatmap.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
row.names(hmap_ec) <- hmap_ec$Enzyme
hmap_ec_edit <- hmap_ec[, !(names(hmap_ec) %in% c('Enzyme','Label'))]

row_annotate <- data.frame(Activity = c("Chitosan Degrading", "Chitosan Degrading", "Glycolysis", "Glycolysis","Glycolysis","Glycolysis","Nitrification","Nitrification","Nitrification"))
row.names(row_annotate) = hmap_ec$Enzyme

pheatmap(hmap_ec_edit, 
         annotation_row = row_annotate,
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 25,
         cellheight = 10,
         fontsize_number = 6, 
         number_color = "black",
         #col = brewer.pal(9, "PiYG"),
         #display_numbers = T,
         border_color = "grey60")

# plotting heatmap for chosen genes
hmap_ko <- read.delim(file="/home/dlaw16/ko_output50/ko_heatmap.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
row.names(hmap_ko) <- hmap_ko$Gene
hmap_ko_edit <- hmap_ko[, !(names(hmap_ko) %in% c('Gene','Activity'))]

row_annotate_ko <- data.frame(Activity = c("Assimilatory nitrate reduction", "Assimilatory nitrate reduction", "Denitrification"))
row.names(row_annotate_ko) = hmap_ko$Gene

pheatmap(hmap_ko_edit, 
         annotation_row = row_annotate_ko,
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 25,
         cellheight = 10,
         fontsize_number = 6, 
         number_color = "black",
         #col = brewer.pal(9, "PiYG"),
         #display_numbers = T,
         border_color = "grey60")

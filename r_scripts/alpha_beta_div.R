library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)
library(phyloseq)
library(tidyverse)
library(vegan)
library(dplyr)

# ALPHA DIVERSITY
merged_df<-read.table(file ="/home/dlaw16/metagenomics/processing/profiling/metaphlan/metaphlan_merged_table.txt", sep="\t", header = TRUE)
colnames(merged_df)<-gsub("d_metagenome","",colnames(merged_df))

# Remove first column
rownames(merged_df) <- merged_df[,1]
merged_df <- merged_df[,-1]

# Remove samples with all unknowns
removed <- which(colSums(merged_df) == 0)

# Data transformation
merged_df <- merged_df / 100 

richness_val<- richness(merged_df)
dominance_val <- dominance(merged_df, index = "all")
rarity_val <- rarity(merged_df, index = "all")
even_val <- evenness(merged_df, "all")

shannon<-microbiome::alpha(merged_df, index = c("diversity_shannon"))
gini_simpson<-microbiome::alpha(merged_df, index = c("diversity_gini_simpson"))
div_invsimp<-microbiome::alpha(merged_df, index = c("diversity_inverse_simpson"))


metadata_div<-read.delim(file="/home/dlaw16/metagenomics/processing/profiling/metaphlan/metadata.tsv", sep="\t", header = TRUE)

metadata_div$observed<-richness_val$observed
metadata_div$chao<-richness_val$chao1
metadata_div$shannon<-shannon$diversity_shannon
metadata_div$inv_simpson<-div_invsimp$diversity_inverse_simpson
metadata_div$lowabun<-rarity_val$low_abundance
metadata_div$rareabun<-rarity_val$rare_abundance

p_obv <- ggplot(metadata_div, aes(x=type,y=observed,color=type, fill=type))+geom_boxplot(alpha=0.3)+geom_point()+theme_light()+labs(y= "Observed", x = "Samples")+theme(legend.position = "none")
p_chao <- ggplot(metadata_div, aes(x=type,y=chao,color=type, fill=type))+geom_boxplot(alpha=0.3)+geom_point()+theme_light()+labs(y= "Chao", x = "Samples")+theme(legend.position = "none")
p_shannon <- ggplot(metadata_div, aes(x=type,y=shannon,color=type, fill=type))+geom_boxplot(alpha=0.3)+geom_point()+theme_light()+labs(y= "Shannon", x = "Samples")+theme(legend.position = "none")
p_simpson <- ggplot(metadata_div, aes(x=type,y=inv_simpson,color=type, fill=type))+geom_boxplot(alpha=0.3)+geom_point()+theme_light()+labs(y= "Simpson", x = "Samples")+theme(legend.position = "none")
p_lowabun <- ggplot(metadata_div, aes(x=type,y=lowabun,color=type, fill=type))+geom_boxplot(alpha=0.3)+geom_point()+theme_light()+labs(y= "Low abundance", x = "Samples")+theme(legend.position = "none")
p_rareabun <- ggplot(metadata_div, aes(x=type,y=rareabun,color=type, fill=type))+geom_boxplot(alpha=0.3)+geom_point()+theme_light()+labs(y= "Rare abundance", x = "Samples")+theme(legend.position = "none")

alphadiv_combi <- gather(metadata_div, key="measure", value="value", c("chao","shannon")) 
alphadiv_combi$ordered_measure <- factor(alphadiv_combi$measure, levels = c("chao","shannon"))
ggplot(alphadiv_combi, aes(x=type, y=value))+ geom_boxplot()+facet_wrap(~ordered_measure,scales="free_y",labeller = as_labeller(c(chao='Chao1',shannon='Shannon'))) +labs(y= "Alpha diversity", x = "Samples")+theme(strip.text.x = element_text(size = 20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))


# Ordination: PCoA with Bray-Curtis distance (beta diversity)
srs_df<-read.table(file ="/home/dlaw16/metagenomics/processing/profiling/metaphlan/merged_abundance_table_species.txt", sep="\t", header = TRUE)
srs_mat<-srs_df |> column_to_rownames("clade_name") |> as.matrix() |> t()
dist_mat<-vegdist(srs_mat, method="bray")
cmd_res<-cmdscale(dist_mat,k=(nrow(srs_mat) -1), eig = TRUE)
str(cmd_res)
pcoa_df<- tibble(PC1 = cmd_res$points[,1],PC2 = cmd_res$points[,2])
metadata<-read.delim(file="/home/dlaw16/metagenomics/processing/profiling/metaphlan/metadata.tsv", sep="\t", header = TRUE)
pcoa_meta=bind_cols(pcoa_df,metadata)
rownames(pcoa_meta)<-pcoa_meta$id;
pcoa_meta$Category <- c("control","control","control","treated", "treated","treated","treated","treated","treated","treated","treated","treated","treated","treated",
                        "treated","treated","treated","treated","control","control","control");
pcoa_var_explained <- round(cmd_res$eig / sum(cmd_res$eig) * 100, 2);

pal_edit<- brewer.pal(11, "PuOr");
pal_edit<-pal_edit[-c(3,4,5,6)];

ggplot(pcoa_meta, aes(x = PC1, y = PC2, color=type, shape = Category)) +geom_point(size =8)+ guides(color=guide_legend(title="Type", ncol =7),
      shape = guide_legend(position = "top"))+ scale_color_manual(values = pal_edit)+theme_light() +theme(legend.position='bottom',legend.text = element_text(size=20),
      legend.key.size = unit(1, 'cm'),legend.title = element_text(size=30), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30),axis.title = element_text(size = 30)) +
      labs(x = paste0("PCoA1 (", pcoa_var_explained[1], "%)"), y = paste0("PCoA2 (", pcoa_var_explained[2], "%)"));

# PERMANOVA (beta diversity statistics)
srs_mat_copy<-srs_mat
label<-rep(c("C0","S1","S2","S3","S4","S5","Cfinal"),each=3)
permstat<-adonis2(srs_mat_copy ~ label,method="bray",perm=999)

#within sample (replicates)
c0_data<-srs_mat_copy[1:3,]
label_rep<-c("1","2","3")
stat_c0<-adonis2(c0_data ~ label_rep,method="bray",perm=999)

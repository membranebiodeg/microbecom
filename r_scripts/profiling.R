library(tidyverse)
library(fs)
library(wesanderson)

# colors for the taxnomical composition figure
name_wespal<-names(wes_palettes)
#name_wespal<-name_wespal[-c(1,3,7,9,12,13,14,19)]
listpal=list()

for (x in 1:length(name_wespal)) {
  i = name_wespal[x]
  valuepal=wes_palette(i)
  listpal=c(listpal,valuepal)
}
wespal<-unlist(listpal)


path <-"/home/dlaw16/metagenomics/processing/profiling/metaphlan"
file_dir<- dir_ls(path = path, recurse = TRUE)
#list all the path in that directory for every files
file_list <- file_dir[grepl("metagenome.txt$",file_dir)]
#in this case file_dir and file_list is the same because there is only txt files in the dir

# Composition plot at the phylum level
results_metaphyl<-data.frame() #empty dataframe
results_metaphyl10<-data.frame()
for (x_metaphyl in 1:length(file_list)) {
  i_metaphyl <- file_list[x_metaphyl] #storing the path in variable i
  name_metaphyl <- gsub(".*/(.*)_profiled_metagenome.txt","\\1",i_metaphyl)
  #\\ reads 1 literally (environment variable)
  
  meta_out<-read.table(file =i_metaphyl, sep="\t")
  colnames(meta_out)<-c("clade_name","NCBI_tax_id","relative_abundance","additional_species")
  
  metaphyl_extract<-meta_out[grepl("p__[A-Za-z0-9]+$", meta_out$clade_name),] 
  metaphyl_extract$Sample<-name_metaphyl #add a column named sample
  metaphyl_extract$label=sapply(strsplit(metaphyl_extract$clade_name,'p__'),"[",2)
  #strsplit(metaphlan_family$clade_name,'f__') will split it when it found f__. 
  #Since there is only one f__, it will split it into 2 columns. 
  #sapply(fileNameSplit, `[`, 1) can be rewritten as sapply(fileNameSplit, function(x) x[1])
  #?sapply or ?'[' for more info
  
  results_metaphyl<-rbind(results_metaphyl,metaphyl_extract)
  
}

results_metaphyl$Sample<-str_replace(results_metaphyl$Sample, "t0_1", "c0_1")
results_metaphyl$Sample<-str_replace(results_metaphyl$Sample, "t0_2", "c0_2")
results_metaphyl$Sample<-str_replace(results_metaphyl$Sample, "t0_3", "c0_3")
results_metaphyl$Sample<-str_replace(results_metaphyl$Sample, "t6_1", "cf_1")
results_metaphyl$Sample<-str_replace(results_metaphyl$Sample, "t6_2", "cf_2")
results_metaphyl$Sample<-str_replace(results_metaphyl$Sample, "t6_3", "cf_3")

#control comparison
controlphyl<-results_metaphyl[results_metaphyl$Sample == 'c0_1'|results_metaphyl$Sample == 'c0_2'|results_metaphyl$Sample == 'c0_3'|results_metaphyl$Sample=='cf_1'|results_metaphyl$Sample=='cf_2'|results_metaphyl$Sample=='cf_3',]
n_phyl1<-length(unique(controlphyl$label))
ggplot(controlphyl,aes(y=relative_abundance,x=Sample))+geom_bar(aes(fill = label),stat = "identity")+labs(y= "Relative Abundance", x = "Samples",fill = "Bacteria phylum")+theme(legend.position='bottom',legend.text = element_text(size=15),
legend.key.size = unit(1.5, 'cm'),legend.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))+scale_fill_manual(values=sample(wespal,n_phyl1))

#control_t0 vs samples
metaphyl_controlinit_sample<-results_metaphyl[results_metaphyl$Sample !='cf_1'& results_metaphyl$Sample !='cf_2'& results_metaphyl$Sample != 'cf_3',]
n_phyl2<-length(unique(metaphyl_controlinit_sample$label))
ggplot(metaphyl_controlinit_sample,aes(y=relative_abundance,x=Sample))+geom_bar(aes(fill = label),stat = "identity")+labs(y= "Relative Abundance", x = "Samples",fill = "Bacteria phylum")+theme(legend.position='bottom',
legend.text = element_text(size=15),legend.key.size = unit(1.5, 'cm'),legend.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))+scale_fill_manual(values=sample(wespal,n_phyl2))

#control_final vs samples
metaphyl_controlfinal_sample<-results_metaphyl[results_metaphyl$Sample !='c0_1'& results_metaphyl$Sample !='c0_2'& results_metaphyl$Sample != 'c0_3',]
n_phyl3<-length(unique(metaphyl_controlfinal_sample$label))
ggplot(metaphyl_controlfinal_sample,aes(y=relative_abundance,x=Sample))+geom_bar(aes(fill = label),stat = "identity")+labs(y= "Relative Abundance", x = "Samples",fill = "Bacteria phylum")+theme(legend.position='bottom',
legend.text = element_text(size=15),legend.key.size = unit(1.5, 'cm'),legend.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))+scale_fill_manual(values=sample(wespal,n_phyl3))

# Composition plot at the genus level
metagen_new=data.frame()
for (x in 1:length(names)){
metagen_temp<-results_metagen[results_metagen$finallabel != "Others" & results_metagen$Sample == names[x],]
sum_metagen<-sum(metagen_temp$relative_abundance)
metagen_others<-100-sum_metagen
othersgen_row <- data.frame(clade_name="Others",
                         NCBI_tax_id="",
                         relative_abundance=metagen_others,
                         additional_species="",
                         Sample= names[x],
                         label="Others",
                         finallabel="",
                         labelfam="")
metagen_new<-rbind(metagen_new,metagen_temp,othersgen_row)

}
metagen_new$Sample<-str_replace(metagen_new$Sample, "t0_1", "c0_1")
metagen_new$Sample<-str_replace(metagen_new$Sample, "t0_2", "c0_2")
metagen_new$Sample<-str_replace(metagen_new$Sample, "t0_3", "c0_3")
metagen_new$Sample<-str_replace(metagen_new$Sample, "t6_1", "cf_1")
metagen_new$Sample<-str_replace(metagen_new$Sample, "t6_2", "cf_2")
metagen_new$Sample<-str_replace(metagen_new$Sample, "t6_3", "cf_3")

#control comparison
controlgen<-metagen_new[metagen_new$Sample == 'c0_1'|metagen_new$Sample == 'c0_2'|metagen_new$Sample == 'c0_3'|metagen_new$Sample=='cf_1'|metagen_new$Sample=='cf_2'|metagen_new$Sample=='cf_3',]
orderedlabelgen<-c("Others",sort(unique(controlgen$label)[!unique(controlgen$label) %in% "Others"]))
controlgen$orderedlabelgen<-factor(controlgen$label,levels = orderedlabelgen) 
controlgen$orderedlabelgen <- gsub("GGB63470","Nitrospiraceae_unclassified", controlgen$orderedlabelgen);
controlgen$orderedlabelgen <- gsub("GGB66002","Proteobacteria_unclassified", controlgen$orderedlabelgen);
n_gen1<-length(unique(controlgen$orderedlabelgen))
ggplot(controlgen,aes(y=relative_abundance,x=Sample))+geom_bar(aes(fill = orderedlabelgen),stat = "identity")+labs(y= "Relative Abundance", x = "Samples",fill = "Bacteria genus")+theme(legend.position='bottom',
legend.text = element_text(size=12),legend.key.size = unit(1, 'cm'),legend.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))+scale_fill_manual(values=sample(wespal,n_gen1))

#control_t0 vs samples
metagen_controlinit_sample<-metagen_new[metagen_new$Sample !='cf_1'& metagen_new$Sample !='cf_2'& metagen_new$Sample != 'cf_3',]
orderedlabelgen<-c("Others",sort(unique(metagen_controlinit_sample$label)[!unique(metagen_controlinit_sample$label) %in% "Others"]))
metagen_controlinit_sample$orderedlabelgen<-factor(metagen_controlinit_sample$label,levels = orderedlabelgen) 
n_gen2<-length(unique(metagen_controlinit_sample$orderedlabelgen))
ggplot(metagen_controlinit_sample,aes(y=relative_abundance,x=Sample))+geom_bar(aes(fill = orderedlabelgen),stat = "identity")+labs(y= "Relative Abundance", x = "Samples",fill = "Bacteria genus")+theme(legend.position='bottom',
legend.text = element_text(size=10),legend.key.size = unit(0.5, 'cm'),legend.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))+scale_fill_manual(values=sample(wespal,n_gen2))

#control_final vs samples
metagen_controlfinal_sample<-metagen_new[metagen_new$Sample !='c0_1'& metagen_new$Sample !='c0_2'& metagen_new$Sample != 'c0_3',]
orderedlabelgen<-c("Others",sort(unique(metagen_controlfinal_sample$label)[!unique(metagen_controlfinal_sample$label) %in% "Others"]))
metagen_controlfinal_sample$orderedlabelgen<-factor(metagen_controlfinal_sample$label,levels = orderedlabelgen) 
n_gen3<-length(unique(metagen_controlfinal_sample$orderedlabelgen))
ggplot(metagen_controlfinal_sample,aes(y=relative_abundance,x=Sample))+geom_bar(aes(fill = orderedlabelgen),stat = "identity")+labs(y= "Relative Abundance", x = "Samples",fill = "Bacteria genus")+theme(legend.position='bottom',
legend.text = element_text(size=10),legend.key.size = unit(0.5, 'cm'),legend.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),axis.title = element_text(size = 20))+scale_fill_manual(values=sample(wespal,n_gen3))

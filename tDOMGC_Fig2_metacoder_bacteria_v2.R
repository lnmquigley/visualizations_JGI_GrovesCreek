require(metacoder)
install.packages("ggrepel")
require(ggrepel)
bacteria <- read.csv("bacteria_order_75_mc.csv")
bacteria <- subset(bacteria, sum >= 1600)
head(bacteria)
metadata <- read.csv("july_metadata.csv")

bac <- parse_tax_data(bacteria, datasets = list(meta=metadata), mappings=c("{{name}}"="sample_id"),
                            class_cols = "lineage", class_sep = ";",
                            class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.+)$")

bac
bac$data$tax_abund <- calc_taxon_abund(bac, "tax_data",
                                             cols=metadata$sample_id)
bac$data$tax_occ <- calc_n_samples(bac, "tax_abund", groups=metadata$tide)
bac$data$diff_table <- compare_groups(bac, dataset = "tax_abund",
                                            cols=metadata$sample_id,
                                            groups=metadata$tide)
bact <- heat_tree(bac, node_label=tax_name, node_size=n_obs, 
                  make_node_legend=FALSE)
print(bact)
####metagenome#####
bacteria_mg <- subset(bacteria, select=c(MG_01, MG_09B, MG_13B, MG_02, MG_04, MG_05A, MG_03,
                                   MG_05C, MG_06, MG_07, MG_09A, MG_10, MG_11, MG_12A,
                                   MG_12B, MG_13A, lineage))
head(bacteria_mg)

meta_mg <- subset(metadata, type=="metagenome")
head(meta_mg)

bac_metag <- parse_tax_data(bacteria_mg, datasets = list(meta=meta_mg), mappings=c("{{name}}"="sample_id"),
                        class_cols = "lineage", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
bac_metag
bac_metag$data$tax_abund <- calc_taxon_abund(bac_metag, "tax_data",
                                         cols=meta_mg$sample_id)
bac_metag$data$tax_occ <- calc_n_samples(bac_metag, "tax_abund", groups=meta_mg$tide)
bac_metag$data$diff_table <- compare_groups(bac_metag, dataset = "tax_abund",
                                        cols=meta_mg$sample_id,
                                        groups=meta_mg$tide)
print(bac_metag$data$diff_table, n=126)
print(bac_metag$data$tax_occ, n=100)
print(bac_metag$data$tax_abund)
bac_metag
?heat_tree
bac_tree <- heat_tree(bac_metag, 
                      node_label=ifelse(wilcox_p_value <= 0.05, tax_name, NA),
                      node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                      node_color_interval=c(-2,2),
                      node_color_range=c("darkviolet", "gray", "darkorange"), 
                      overlap_avoidance=1, 
                      make_node_legend=FALSE)
print(bac_tree)


#####metatranscriptome#####
bacteria_mt <- subset(bacteria, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07,
                                         MT_08, MT_09A, MT_10, MT_11, MT_12A, MT_12B, MT_13A,
                                         MT_13C, lineage))
head(bacteria_mt)

metadata <- read.csv("july_metadata.csv")
meta_mt <- subset(metadata, type=="metatranscriptome")
head(meta_mt)
meta_mt <- subset(meta_mt, sample_id != "MT_09B")

bac_metat <- parse_tax_data(bacteria_mt, datasets = list(meta=meta_mt), mappings=c("{{name}}"="sample_id"),
                            class_cols = "lineage", class_sep = ";",
                            class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.+)$")
bac_metat
bac_metat$data$tax_abund <- calc_taxon_abund(bac_metat, "tax_data",
                                             cols=meta_mt$sample_id)
bac_metat$data$tax_occ <- calc_n_samples(bac_metat, "tax_abund", groups=meta_mt$tide)
bac_metat$data$diff_table <- compare_groups(bac_metat, dataset = "tax_abund",
                                            cols=meta_mt$sample_id,
                                            groups=meta_mt$tide)
print(bac_metat$data$diff_table, n=100)
print(bac_metat$data$tax_occ, n=100)
print(bac_metat$data$tax_abund)
bac_metat
?heat_tree
bac_tree_mt <- heat_tree(bac_metat, 
                      node_label=ifelse(wilcox_p_value <=0.05, tax_name, NA),
                      node_color=log2_median_ratio,
                      node_color_interval=c(-2,2),
                      node_color_range=c("darkviolet", "gray", "darkorange"), 
                      make_node_legend=FALSE)
print(bac_tree_mt)

######16S######

r16S <- subset(bacteria, select=c(r16_01, r16_02, r16_03, r16_04, r16_05A, r16_05C, r16_06, r16_07,
                                         r16_08, r16_09A, r16_09B, r16_10, r16_11, r16_12A, r16_12B, r16_13A,
                                         r16_13B, lineage))
head(r16S)

metadata <- read.csv("july_metadata.csv")
meta_16s <- subset(metadata, type=="16S")
head(meta_16s)

bac_16S <- parse_tax_data(r16S, datasets = list(meta=meta_16s), mappings=c("{{name}}"="sample_id"),
                            class_cols = "lineage", class_sep = ";",
                            class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.+)$")
bac_16S
bac_16S$data$tax_abund <- calc_taxon_abund(bac_16S, "tax_data",
                                             cols=meta_16s$sample_id)
bac_16S$data$tax_occ <- calc_n_samples(bac_16S, "tax_abund", groups=meta_16s$tide)
bac_16S$data$diff_table <- compare_groups(bac_16S, dataset = "tax_abund",
                                            cols=meta_16s$sample_id,
                                            groups=meta_16s$tide)
print(bac_16S$data$diff_table, n=100)
print(bac_16S$data$tax_occ, n=100)
print(bac_16S$data$tax_abund)
bac_16S
?heat_tree
bac_tree_16S <- heat_tree(bac_16S, 
                         node_label=ifelse(wilcox_p_value <=0.05, tax_name, NA),
                         node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                         node_color_interval=c(-2,2),
                         node_color_range=c("darkviolet", "gray", "darkorange"), 
                         make_node_legend=FALSE)
print(bac_tree_16S)




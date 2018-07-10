require(metacoder)
install.packages("ggrepel")
require(ggrepel)
archaea <- read.csv("archaea_order_75_mc.csv")
archaea <- subset(archaea, sum >= 4)
head(archaea)
metadata <- read.csv("july_metadata.csv")

arc <- parse_tax_data(archaea, datasets = list(meta=metadata), mappings=c("{{name}}"="sample_id"),
                            class_cols = "lineage", class_sep = ";",
                            class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.+)$")

arc
arc$data$tax_abund <- calc_taxon_abund(arc, "tax_data",
                                             cols=metadata$sample_id)
arc$data$tax_occ <- calc_n_samples(arc, "tax_abund", groups=metadata$tide)
arc$data$diff_table <- compare_groups(arc, dataset = "tax_abund",
                                            cols=metadata$sample_id,
                                            groups=metadata$tide)
arct <- heat_tree(arc, node_label=tax_name, node_size=n_obs, 
                  make_node_legend=FALSE)
print(arct)
####metagenome#####
arc_mg <- subset(archaea, select=c(MG_01, MG_09B, MG_13B, MG_02, MG_04, MG_05A, MG_03,
                                   MG_05C, MG_06, MG_07, MG_09A, MG_10, MG_11, MG_12A,
                                   MG_12B, MG_13A, lineage))
head(arc_mg)

meta_mg <- subset(metadata, type=="metagenome")
head(meta_mg)

arc_metag <- parse_tax_data(arc_mg, datasets = list(meta=meta_mg), mappings=c("{{name}}"="sample_id"),
                        class_cols = "lineage", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
arc_metag
arc_metag$data$tax_abund <- calc_taxon_abund(arc_metag, "tax_data",
                                         cols=meta_mg$sample_id)
arc_metag$data$tax_occ <- calc_n_samples(arc_metag, "tax_abund", groups=meta_mg$tide)
arc_metag$data$diff_table <- compare_groups(arc_metag, dataset = "tax_abund",
                                        cols=meta_mg$sample_id,
                                        groups=meta_mg$tide)
print(arc_metag$data$diff_table)
print(arc_metag$data$tax_occ, n=100)
print(arc_metag$data$tax_abund)
arc_metag
arc_tree <- heat_tree(arc_metag, 
                      node_label=ifelse(wilcox_p_value <= 0.05, tax_name, NA),
                      node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                      node_color_interval=c(-2,2),
                      node_color_range=c("darkviolet", "gray", "darkorange"), 
                      overlap_avoidance=1,
                      node_color_axis_label = "Log-Fold Change across Tidal Cycle")
print(arc_tree)


#####metatranscriptome#####
arc_mt <- subset(archaea, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07,
                                         MT_08, MT_09A, MT_10, MT_11, MT_12A, MT_12B, MT_13A,
                                         MT_13C, lineage))
head(arc_mt)

metadata <- read.csv("july_metadata.csv")
meta_mt <- subset(metadata, type=="metatranscriptome")
head(meta_mt)
meta_mt <- subset(meta_mt, sample_id != "MT_09B")

arc_metat <- parse_tax_data(arc_mt, datasets = list(meta=meta_mt), mappings=c("{{name}}"="sample_id"),
                            class_cols = "lineage", class_sep = ";",
                            class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.+)$")
arc_metat
arc_metat$data$tax_abund <- calc_taxon_abund(arc_metat, "tax_data",
                                             cols=meta_mt$sample_id)
arc_metat$data$tax_occ <- calc_n_samples(arc_metat, "tax_abund", groups=meta_mt$tide)
arc_metat$data$diff_table <- compare_groups(arc_metat, dataset = "tax_abund",
                                            cols=meta_mt$sample_id,
                                            groups=meta_mt$tide)
print(arc_metat$data$diff_table)
print(arc_metat$data$tax_occ, n=100)
print(arc_metat$data$tax_abund)
arc_metat
arc_tree_mt <- heat_tree(arc_metat, 
                      node_label=ifelse(wilcox_p_value <=0.05, tax_name, NA),
                      node_color=log2_median_ratio,
                      node_color_interval=c(-2,2),
                      node_color_range=c("darkviolet", "gray", "darkorange"),
                      repel_force=2,
                      overlap_avoidance=5, 
                      make_node_legend=FALSE)
print(arc_tree_mt)

#####16s#####
arc_16s <- subset(archaea, select=c(r16_01, r16_02, r16_03, r16_04, r16_05A, r16_05C, r16_06, r16_07,
                                   r16_08, r16_09A, r16_09B, r16_10, r16_11, r16_12A, r16_12B, r16_13A,
                                   r16_13B, lineage))
head(arc_16s)

metadata <- read.csv("july_metadata.csv")
meta_16s <- subset(metadata, type=="16S")
head(meta_16s)
#meta_mt <- subset(meta_mt, sample_id != "MT_09B")

arc_16s <- parse_tax_data(arc_16s, datasets = list(meta=meta_mt), mappings=c("{{name}}"="sample_id"),
                            class_cols = "lineage", class_sep = ";",
                            class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.+)$")
arc_16s
arc_16s$data$tax_abund <- calc_taxon_abund(arc_16s, "tax_data",
                                             cols=meta_mt$sample_id)
arc_16s$data$tax_occ <- calc_n_samples(arc_16s, "tax_abund", groups=meta_mt$tide)
arc_16s$data$diff_table <- compare_groups(arc_16s, dataset = "tax_abund",
                                            cols=meta_mt$sample_id,
                                            groups=meta_mt$tide)
print(arc_16s$data$diff_table)
print(arc_16s$data$tax_occ, n=100)
print(arc_16s$data$tax_abund)
arc_16s
arc_tree_16s <- heat_tree(arc_16s, 
                         node_label=ifelse(wilcox_p_value <=0.05, tax_name, NA),
                         node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                         node_color_interval=c(-2,2),
                         node_color_range=c("darkviolet", "gray", "darkorange"),
                         repel_force=2,
                         overlap_avoidance=5, 
                         make_node_legend=FALSE)
print(arc_tree_16s)




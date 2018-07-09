require(metacoder)
bacteria <- read.csv("bacteria_order.csv", row.names=11)
bac_july <- subset(bacteria, select=c(MT_05BTJ, MT_05CTJ, MT_12ATJ, MT_12BTJ))
bac_july$sum <- rowSums(bac_july)
head(bac_july)
bac_july <- subset(bac_july, select=-c(occur, sum))
bac_july <- subset(bac_july, sum >0)
write.csv(bac_july, "bacteria_order_july.csv")

bac_april <- subset(bacteria, select=-c(MT_05BTJ, MT_05CTJ, MT_12ATJ, MT_12BTJ))
bac_april$sum <- rowSums(bac_april)
head(bac_april)
bac_april <- subset(bac_april, sum > 0)
write.csv(bac_april, "bacteria_order_april.csv")

bac_july <- read.csv("bacteria_order_july.csv")
head(bac_july)
bac_july <- subset(bac_july, sum >=800)

bac_july_5 <- subset(bac_july, select=c(X, MT_05BTJ, MT_05CTJ, sum))
bac_july_5 <- subset(bac_july_5, sum >= 800)
bac_july_12 <- subset(bac_july, select=c(X, MT_12ATJ, MT_12BTJ, sum))
bac_july_12 <- subset(bac_july_12, sum >= 800)

metadata <- read.csv("april_meta.csv")
meta5 <- subset(metadata, timepoint=="J-1")
meta12 <- subset(metadata, timepoint=="J-2")
meta <- subset(metadata, month=="July")

july <- parse_tax_data(bac_july, datasets = list(meta=meta), mappings=c("{{name}}"="Sample_ID"),
                        class_cols = "X", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
july
july$data$tax_abund <- calc_taxon_abund(july, "tax_data",
                                         cols=meta$Sample_ID)

july_tree <- heat_tree(july, node_label=tax_name,
                        overlap_avoidance=1, 
                        layout="fruchterman-reingold",
                       initial_layout="da",
                        make_node_legend=FALSE)
print(july_tree)



####July 5####
july5 <- parse_tax_data(bac_july_5, datasets = list(meta=meta5), mappings=c("{{name}}"="Sample_ID"),
                        class_cols = "X", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
july5
july5$data$tax_abund <- calc_taxon_abund(july5, "tax_data",
                                         cols=meta5$Sample_ID)
july5$data$diff_table <- compare_groups(july5, dataset = "tax_abund",
                                        cols=meta5$Sample_ID,
                                        groups=meta5$replicate)
print(july5$data$diff_table, n=101)
print(july5$data$tax_abund)

july5_tree <- heat_tree(july5, node_label=ifelse(log2_median_ratio > 1 | log2_median_ratio < -1, tax_name, NA),
                      node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                      node_color_range=c("#00746F", "#F0EDE3", "#ABC178"),
                      node_color_interval=c(-2, 2),
                      overlap_avoidance=1, 
                      make_node_legend=TRUE,
                      node_color_axis_label="Log-fold Change")
print(july5_tree)

####July12####
july12 <- parse_tax_data(bac_july_12, datasets = list(meta=meta12), mappings=c("{{name}}"="Sample_ID"),
                        class_cols = "X", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
july12
july12$data$tax_abund <- calc_taxon_abund(july12, "tax_data",
                                         cols=meta12$Sample_ID)
july12$data$diff_table <- compare_groups(july12, dataset = "tax_abund",
                                        cols=meta12$Sample_ID,
                                        groups=meta12$replicate)
print(july12$data$diff_table, n=101)
print(july12$data$tax_abund)
july12_tree <- heat_tree(july12, node_label=ifelse(log2_median_ratio > 1 | log2_median_ratio < -1, tax_name, NA),
                        node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                        node_color_range=c("#00746F", "#F0EDE3", "#ABC178"), 
                        overlap_avoidance=1, 
                        node_color_interval=c(-2,2),
                        node_color_axis_label="Log-fold Change",
                        make_node_legend=TRUE)
print(july12_tree)














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


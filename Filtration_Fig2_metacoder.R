require(metacoder)
bacteria <- read.csv("bacteria_order.csv", row.names=11)
bac_july <- subset(bacteria, select=c(MT_05BTJ, MT_05CTJ, MT_12ATJ, MT_12BTJ))
bac_july$sum <- rowSums(bac_july)
head(bac_july)
bac_july <- subset(bac_july, sum >0)
write.csv(bac_july, "bacteria_order_july.csv")

bac_april <- subset(bacteria, select=-c(MT_05BTJ, MT_05CTJ, MT_12ATJ, MT_12BTJ))
bac_april$sum <- rowSums(bac_april)
head(bac_april)
bac_april <- subset(bac_april, sum > 0)
write.csv(bac_april, "bacteria_order_april.csv")

bac_april <- read.csv("bacteria_order_april.csv")
head(bac_april)

bac_april <- subset(bac_april, sum >=800)

bac_april_6 <- subset(bac_april, select=c(X, MT_06BSA, MT_06CTA, sum))
bac_april_6 <- subset(bac_april_6, sum >= 800)
bac_april_12 <- subset(bac_april, select=c(X, MT_12BSA, MT_12ATA, sum))
bac_april_12 <- subset(bac_april_12, sum >= 800)
bac_april_13 <- subset(bac_april, select=c(X, MT_13BSA, MT_13YTA, sum))
bac_april_13 <- subset(bac_april_13, sum >= 800)

metadata <- read.csv("april_meta.csv")
meta6 <- subset(metadata, timepoint=="A-1")
meta12 <- subset(metadata, timepoint=="A-2")
meta13 <- subset(metadata, timepoint=="A-3")

####guide tree####
april <- parse_tax_data(bac_april, datasets = list(meta=meta), mappings=c("{{name}}"="Sample_ID"),
                         class_cols = "X", class_sep = ";",
                         class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                         class_regex = "^(.+)__(.+)$")
april
april$data$tax_abund <- calc_taxon_abund(april, "tax_data",
                                          cols=meta$Sample_ID)

april_tree <- heat_tree(april, node_label=tax_name,
                        overlap_avoidance=1, 
                        layout="fruchterman-reingold",
                        make_node_legend=FALSE)
print(april_tree)

####april 6####
april6 <- parse_tax_data(bac_april_6, datasets = list(meta=meta6), mappings=c("{{name}}"="Sample_ID"),
                        class_cols = "X", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
april6$tax_data
april6$data$tax_abund <- calc_taxon_abund(april6, "tax_data",
                                         cols=meta6$Sample_ID)
april6$data$diff_table <- compare_groups(april6, dataset = "tax_abund",
                                        cols=meta6$Sample_ID,
                                        groups=meta6$time)
print(april6$data$diff_table, n=101)
print(april6$data$tax_abund, n=20)
april6_tree <- heat_tree(april6, node_label=ifelse(log2_median_ratio > 1 | log2_median_ratio < -1, tax_name, NA),
                      node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                      node_color_range=c("#FF8200", "#F0EDE3", "#58595B"), 
                      overlap_avoidance=1,
                      node_color_interval=c(-5,5),
                      node_color_axis_label="Log-fold Change",
                      make_node_legend=FALSE)
print(april6_tree)

####april 12####
april12 <- parse_tax_data(bac_april_12, datasets = list(meta=meta12), mappings=c("{{name}}"="Sample_ID"),
                        class_cols = "X", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")
april12
april12$data$tax_abund <- calc_taxon_abund(april12, "tax_data",
                                         cols=meta12$Sample_ID)
april12$data$diff_table <- compare_groups(april12, dataset = "tax_abund",
                                        cols=meta12$Sample_ID,
                                        groups=meta12$time)
print(april12$data$diff_table, n=101)
print(april12$data$tax_abund)
april12_tree <- heat_tree(april12, node_label=ifelse(log2_median_ratio > 1 | log2_median_ratio < -1, tax_name, NA),
                        node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                        node_color_range=c("#FF8200", "#F0EDE3", "#58595B"), 
                        overlap_avoidance=1,
                        node_color_interval=c(-5,5),
                        node_color_axis_label= "Log-fold Change")
print(april12_tree)

####april 13####
april13 <- parse_tax_data(bac_april_13, datasets = list(meta=meta12), mappings=c("{{name}}"="Sample_ID"),
                          class_cols = "X", class_sep = ";",
                          class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                          class_regex = "^(.+)__(.+)$")
april13
april13$data$tax_abund <- calc_taxon_abund(april13, "tax_data",
                                           cols=meta13$Sample_ID)
april13$data$diff_table <- compare_groups(april13, dataset = "tax_abund",
                                          cols=meta13$Sample_ID,
                                          groups=meta13$time)
print(april13$data$diff_table, n=101)
print(april13$data$tax_abund)
april13_tree <- heat_tree(april13, node_label=ifelse(log2_median_ratio > 1 | log2_median_ratio < -1, tax_name, NA),
                          node_color=ifelse(mean_diff==0, "#99FFFF", log2_median_ratio),
                          node_color_range=c("#FF8200", "#F0EDE3", "#58595B"), 
                          overlap_avoidance=1, 
                          node_color_interval=c(-5,5),
                          make_node_legend=TRUE,
                          node_color_axis_label="Log-fold Change")
print(april13_tree)



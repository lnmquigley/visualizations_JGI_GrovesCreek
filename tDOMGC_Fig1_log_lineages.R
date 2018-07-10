phylo_ex <- read.csv("phylo_family_tide.csv")
head(phylo_ex)
phylo_ex_low <- subset(phylo_ex, select=c(lineage, Ex_Low, Ave_Low))
head(phylo_ex_low)
phylo_ex_low <- subset(phylo_ex_low, Ave_Low >=50)
phylo_ex_low <- subset(phylo_ex_low, Ex_Low !="#DIV/0!")

write.csv(phylo_ex_low, "phylo_expression_lowtide.csv", row.names = FALSE)

phylo_ex_high <- subset(phylo_ex, select=c(lineage, Ex_High, Ave_High))
head(phylo_ex_high)
phylo_ex_high <- subset(phylo_ex_high, Ave_High >=50)
phylo_ex_high <- subset(phylo_ex_high, Ex_High !="#DIV/0!")

write.csv(phylo_ex_high, "phylo_expression_hightide.csv", row.names = FALSE)

high <- read.csv("phylo_expression_hightide.csv")
head(high)
require(ggplot2)
log_high_p <- ggplot(high, aes(x=lineage, y=Log_Ex)) + 
  geom_point(aes(size=Abundance), stat="identity") + 
  coord_flip() + 
  xlab("Lineage") +
  scale_y_continuous(limits=c(-1.3,2)) +
  ylab("Log Fold Expression (TPM/GPM)") +
  geom_hline(yintercept=0) +
  theme_classic()
print(log_high_p)

low <- read.csv("phylo_expression_lowtide.csv")
log_april_p <- ggplot(low, aes(x=lineage, y=Log_Ex)) + 
  geom_point(aes(size=Abundance), stat="identity") + 
  coord_flip() + 
  xlab("Lineage") +
  geom_hline(yintercept=0) +
  ylab("Log Fold Expression (TPM/GPM)") +
  theme_classic()
print(log_april_p)

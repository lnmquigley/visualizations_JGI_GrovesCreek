require(ggplot2)
require(tidyr)

require(dplyr)
early <- read.csv("bacteria_order_april_early.csv")
head(early)

early <- subset(early, select=-c(X, total))
head(early)

early_filt <- subset(early, rel_abund >=1)
head(early_filt)

early$Taxonomy<-ifelse(early$rel_abund <= 1, "other", early$Taxonomy)
early <- subset(early, select=-c(Taxonomy))
early <- merge(early, early_filt, all=TRUE)

write.csv(early, "early_april_filt.csv", row.names=FALSE)

early$Taxonomy[is.na(early$Taxonomy)] <- "other"

early_filt <- separate(early_filt, X, c("domain", "phylum", "class", "order"), sep=";", remove=TRUE)
head(early_filt)
early_filt <- subset(early_filt, select=-c(domain, phylum))
early_filt <- subset(early_filt, class!="c__Clostridia")

early_filt$class <- sub("c__", "", early_filt$class)
early_filt$order <- sub("o__", "", early_filt$order)

early_filt$Taxonomy <- paste(early_filt$class, early_filt$order, sep=";")

early_filt <- subset(early_filt, select=-c(order, class))
write.csv(early_filt, "bacteria_april_early_filt.csv", row.names=FALSE)

early_filt <- read.csv("bacteria_april_early_filt.csv")
head(early_filt)
p <- ggplot(early_filt, aes(Sample, RelAbund, fill=Taxonomy)) +
              geom_bar(stat="identity") +
              theme_classic()
print(p)
?geom_bar

april <- read.csv("early_april_filt.csv")
head(april)
april <- separate(april, Taxonomy, c("domain", "phylum", "class", "order"), sep=";", remove=TRUE)
april <- subset(april, select=-c(domain, phylum))
april$class <- sub("c__", "", april$class)
april$order <- sub("o__", "", april$order)
april$Taxonomy <- paste(april$class, april$order, sep=";")
head(april)
april <- subset(april, select=-c(class, order))



p <- ggplot(april, aes(Sample, RelAbund, fill=Taxonomy)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Time Point") +
  scale_fill_manual(values=c("pink","#EE3E80", "red", "#FF8200","#FF9966",
                             "gold2", "yellow", "palegreen", "lawngreen",
                             "springgreen3", "darkgreen", "skyblue4", "navyblue",
                             "royalblue2", "darkturquoise", "violet", "darkmagenta", "deeppink4",
                             "rosybrown4")) +
  theme_classic() +
  theme(text = element_text(size=20))
print(p)

levels(april$Taxonomy)
str(april)

april$Taxonomy[april$Taxonomy=="Under 1% of ;the community"] <- "Under 1% of the community"

boxB <- subset(aromatic_sep, type=="metatranscriptome" & KO_num=="KO:K15512")

boxB_active <- subset(boxB, GPM >= 2)

boxB <- read.csv("boxB_active.csv")

boxB <- aggregate(boxB$GPM, by=list(boxB$timepoint, boxB$habitat), FUN=sum)
head(boxB)

require(dplyr)
boxB <- rename(boxB, TimePoint = Group.1, Habitat=Group.2, TPM = x)
write.csv(boxB, "boxB_habitat.csv", row.names=FALSE)
boxB <- read.csv("boxB_habitat.csv")
boxB <- subset(boxB, Habitat!="wastewater")
require(ggplot2)

habitat_p <- ggplot(boxB, aes(TimePoint, TPM, fill=Habitat)) +
  geom_area(position="stack") +
  xlab("Time Point") +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("#B9E1E2", "#006C93", "#705550")) +
  theme_classic() +
  theme(text = element_text(size=20)) +
  theme(legend.position="bottom")
print(habitat_p)

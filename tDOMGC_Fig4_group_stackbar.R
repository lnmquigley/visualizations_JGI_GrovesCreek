require(ggplot2)
require(dplyr)
require(scales)

group <- read.csv("group1_stackbar_raw.csv")
head(group)
group <- subset(group, padj <=0.05)

group <- aggregate(group$baseMean, by=list(group$group_1, group$tide), FUN=sum)
head(group)

group <- rename(group, KO_Group1=Group.1, Tide=Group.2, Mean_Reads=x)
head(group)

group <- read.csv("group1relabund.csv")

group$RelativeAbundance <- group$RelativeAbundance/100

group_plot <- ggplot(group, aes(x=Tide, y=RelativeAbundance)) +
  geom_bar(stat="identity", aes(fill=KOGroup1)) +
  scale_fill_manual(values=c("pink", "red", "#FF8200", "lawngreen",
                             "royalblue2", "lightskyblue", "darkmagenta")) +
  scale_y_continuous(labels = scales::percent) +
  ylab(expression(atop("Relative Abundance of Reads", paste("within KO Groups")))) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  guides(fill=guide_legend(ncol=2)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="bottom")

print(group_plot)


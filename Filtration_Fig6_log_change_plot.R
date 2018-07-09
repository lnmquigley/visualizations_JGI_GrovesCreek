require(ggplot2)
log <- read.csv("log_fold_change_ko_long.csv")
head(log)
log <- subset(log, group_1!="Unclassified")
log <- subset(log, group_2!="Signaling molecules and interactions" & 
                group_2!="Metabolism of terpenoids and polyketides" &
                group_2!="Metabolism of other amino acids" &
                group_2!="Cellular community - eukaryotes")


log_april <- subset(log, month=="April")
log_april$timepoint <- factor(log_april$timepoint)
head(log_april)
log_april_p <- ggplot(log_april, aes(x=group_2, y=log)) + 
  geom_bar(aes(fill=timepoint), position="dodge", stat="identity") + 
  coord_flip() + 
  xlab("KO Group") +
  ylab("log-fold change") +
  theme_classic()
print(log_april_p)

log_july <- subset(log, month=="July")
log_july$timepoint <- factor(log_july$timepoint)
head(log_july)
log_july_p <- ggplot(log_july, aes(x=group_2, y=log)) + 
  geom_bar(aes(fill=timepoint), position="dodge", stat="identity") + 
  coord_flip() + 
  xlab("KO Group") +
  scale_y_continuous(limits=c(-1, 1)) +
  ylab("log-fold change") +
  theme_classic()
print(log_july_p)

log_p <- ggplot(log, aes(x=group_2, y=log)) + 
  geom_bar(aes(fill=Time_Point), position="dodge", stat="identity") + 
  coord_flip() + 
  facet_wrap(~ month) +
  xlab("KO Group") +
  #scale_y_continuous(limits=c(-1, 1)) +
  ylab("Log-fold Change") +
  theme_classic() +
  theme(text = element_text(size=20))
print(log_p)


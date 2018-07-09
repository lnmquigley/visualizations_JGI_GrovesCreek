require(ggplot2)

ko <- read.csv("log_fold_change_long_early.csv")
ko <- subset(ko, ProcessingTime!="Late")
head(ko)
ko$log <- log10(ko$norm)
ko <- subset(ko, group_2!="Cellular community - eukaryotes")
ko <- subset(ko, group_2!="Cellular processes and signaling")
ko <- subset(ko, group_2!="Genetic information processing")
ko <- subset(ko, group_2!="Metabolism")
ko <- subset(ko, group_2!="Metabolism of other amino acids")
ko <- subset(ko, group_2!="Metabolism of terpenoids and polyketides")
ko <- subset(ko, group_2!="Poorly characterized")
ko <- subset(ko, group_2!="Viral protein family")

kop <- ggplot(ko, aes(Month, norm)) +
  geom_jitter(aes(colour=TimePoint), size=6) +
  facet_wrap(~group_2, scales="free", nrow=7) +
  ylab("Average TPM per KO number within each Group") +
  theme_classic() +
  theme(text = element_text(size=20))

print(kop)

?geom_jitter

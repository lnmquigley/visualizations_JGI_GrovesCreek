require(dplyr)
require(splitstackshape)
require(tidyr)
require(ggplot2)
CAZy <- read.csv("CAZy_GPM.csv")
CAZy_anno <- read.delim("CAZyDB-ec-info.txt", sep="\t", header=FALSE)
CAZy <- cSplit(CAZy, "CAZy", "|")
head(CAZy)
CAZy <- subset(CAZy, select=-c(CAZy_3, CAZy_4))
CAZy <- rename(CAZy, genbank=CAZy_1, family=CAZy_2)
head(CAZy)

CAZy_anno <- rename(CAZy_anno, genbank=V1, family=V2, EC=V3, enzyme_name=V4)
head(CAZy_anno)

CAZy <- merge(CAZy, CAZy_anno)

head(CAZy)
CAZy_by_class <- aggregate(CAZy$GPM, by=list(CAZy$type, CAZy$time, CAZy$class), FUN=sum)
head(CAZy_by_class)
write.csv(CAZy_by_class, "cazy_by_class.csv", row.names=FALSE)
CAZy_by_class <- read.csv("cazy_by_class.csv")
require(ggplot2)

cazy_class_p <- ggplot(CAZy_by_class, aes(Time, Count, fill=Class)) +
  geom_area() +
  facet_grid(~Type) +
  theme_classic() +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("pink", "coral", "springgreen3", "turquoise", "violet", "rosybrown4")) +
  theme(text = element_text(size=20)) +
  ylab("Genes or Transcripts per Million") +
  xlab("Time Point") +
  theme(legend.position="bottom")
print(cazy_class_p)


CAZy_sep <- separate(CAZy, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                            sep=";", remove=TRUE)
CAZy_sep <- subset(CAZy_sep, select=-c(phylum, class, order, family, genus, species, strain))
head(CAZy_sep)
CAZy_domain <- aggregate(CAZy_sep$GPM, by=list(CAZy_sep$type, CAZy_sep$time, CAZy_sep$domain), FUN=sum)
head(CAZy_domain)
tail(CAZy_domain)
CAZy_domain <- rename(CAZy_domain, Type=Group.1, Time=Group.2, Domain=Group.3, Count=x)
CAZy_domain_family <- aggregate(CAZy_sep$GPM, by=list(CAZy_sep$type, CAZy$time, CAZy_sep$domain, CAZy$family), FUN=sum)
head(CAZy_domain_family)
CAZy_domain_family <- rename(CAZy_domain_family, Type=Group.1, Time=Group.2, Domain=Group.3, CAZy_Family=Group.4, Count=x)
CAZy_domain_family <- subset(CAZy_domain_family, CAZy_Family!="PL0" & CAZy_Family!="PL11" & CAZy_Family!="PL1" & CAZy_Family!="AA3" & CAZy_Family!="AA4"
                             & CAZy_Family!="GT0" & CAZy_Family!="CE0" & CAZy_Family!="GH0")


cazy_domain_p <- ggplot(CAZy_domain, aes(Time, Count, fill=Domain)) +
  geom_area() +
  facet_wrap(~Type) +
  theme_classic() +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("coral", "springgreen3", "turquoise", "violet")) +
  theme(text = element_text(size=20)) +
  ylab("Genes or Transcripts per Million") +
  xlab("Time Point") +
  theme(legend.position="bottom")
print(cazy_domain_p)


Aux <- subset(CAZy, class=="Auxiliary Activities")
head(Aux)
tail(Aux)
Aux_by_acitivity <- aggregate(Aux$GPM, by=list(Aux$type, Aux$time, Aux$activity), FUN=sum)
head(Aux_by_acitivity)
Aux_by_acitivity <- rename(Aux_by_acitivity, Type=Group.1, time=Group.2, activity=Group.3, GPM=x)


oxo_p <- ggplot(Aux_by_acitivity, aes(time, GPM)) +
  geom_point(aes(shape=Type, colour=Type), size=8) +
  facet_wrap(~ activity, scales="free_y") +
  theme_classic() +
  scale_colour_manual(values=c("springgreen", "darkgreen")) +
  scale_x_continuous(breaks=c(4,8,12)) +
  theme(text = element_text(size=24)) +
  xlab("Time Point") +
  ylab("Genes or Transcripts per Million") +
  theme(legend.position="bottom")
print(oxo_p)

Auxcommcomp <- subset(Aux, family=="AA10" | family=="AA2" | family=="AA6")
Auxcommcomp_sep <- separate(Auxcommcomp, taxonomy, c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                         sep=";", remove=TRUE)
Auxcommcomp_sep <- subset(Auxcommcomp_sep, select=-c(genus, species, strain))
Auxcommcomp_sep$taxonomy <- paste(Auxcommcomp_sep$domain, Auxcommcomp_sep$phylum, Auxcommcomp_sep$class, Auxcommcomp_sep$order,
                               Auxcommcomp_sep$family, sep=";")
Auxcommcomp_set <- subset(Auxcommcomp_sep, select=-c(domain, phylum, class, order, family))
head(Auxcommcomp_set)
Aux_by_family <- aggregate(Auxcommcomp_sep$GPM, by=list(Auxcommcomp_sep$taxonomy, Auxcommcomp_sep$activity, Auxcommcomp_sep$time, Auxcommcomp_sep$type), FUN=sum)
head(Aux_by_family)
Aux_by_family <- rename(Aux_by_family, Taxonomy=Group.1, Activity=Group.2, Time=Group.3, Type=Group.4, GPM=x)

Aux_by_family_pad <- expand(Aux_by_family, nesting(Taxonomy, Activity, Time, Type), Time)
head(Aux_by_family_pad)
head(Aux_by_family)
Aux_by_family <- rename(Aux_by_family, Time1=Time)
Aux_by_family_pad <- subset(Aux_by_family_pad, select=-c(Time))

Aux_by_family_padded <- merge(Aux_by_family_pad, Aux_by_family, by=c("Taxonomy", "Activity", "Type", "Time1"), all=TRUE)
head(Aux_by_family_padded)
tail(Aux_by_family_padded)
Aux_by_family_padded[is.na(Aux_by_family_padded)] <- 0
Aux_by_family_padded <- unique(Aux_by_family_padded)

write.csv(Aux_by_family_padded, "dbCAN_oxidases_commcomp.csv", row.names = FALSE)

AA <- read.csv("dbCAN_oxidases_commcomp_under1.csv")
head(AA)
LPMO <- subset(AA, Activity=="lytic polysaccharide monooxygenases (LPMO)")
per <- subset(AA, Activity=="class II lignin-modifying peroxidases")
bq <- subset(AA, Activity=="1,4-benzoquinone reductases")

LPMO <- aggregate(LPMO$GPM, by=list(LPMO$Taxonomy, LPMO$Time, LPMO$Type), FUN=sum)
head(LPMO)
LPMO <- rename(LPMO, Taxonomy=Group.1, Time=Group.2, Type=Group.3, Count=x)
LPMO <- separate(LPMO, Taxonomy, c("domain", "phylum", "class", "order", "family"),
                            sep=";", remove=TRUE)
head(LPMO)
LPMO <- subset(LPMO, select=-c(class, order))
LPMO$Taxonomy <- paste(LPMO$domain, LPMO$phylum, LPMO$family, sep=";")

LPMO_p <- ggplot(LPMO, aes(Time, Count, fill=Taxonomy)) +
  geom_area() +
  facet_grid(~Type) +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("pink","#EE3E80", "red", "#FF8200", "tomato", 
                             "#FF9966", "gold2", "khaki", "yellow", "lawngreen", 
                             "springgreen3", "seagreen", "darkgreen", "navyblue", "royalblue2",
                             "turquoise", "powderblue", "thistle", "violet", "purple",
                             "darkmagenta", "violetred4", "deeppink", "palevioletred", "rosybrown4",
                             "mistyrose4", "snow3", "lightskyblue", "palegreen", "wheat",
                             "lightsalmon")) +
  theme_classic() +
  guides(fill=guide_legend(ncol=4)) +
  theme(text = element_text(size=17.5)) +
  xlab("Time Point") +
  ylab("Genes or Transcripts per Million") +
  theme(legend.position="bottom")
print(LPMO_p)

per <- aggregate(per$GPM, by=list(per$Taxonomy, per$Time, per$Type), FUN=sum)
head(per)
per <- rename(per, Taxonomy=Group.1, Time=Group.2, Type=Group.3, Count=x)
per <- separate(per, Taxonomy, c("domain", "phylum", "class", "order", "family"),
                 sep=";", remove=TRUE)
head(per)
per <- subset(per, select=-c(class, order))
per$Taxonomy <- paste(per$domain, per$phylum, per$family, sep=";")

per_p <- ggplot(per, aes(Time, Count, fill=Taxonomy)) +
  geom_area() +
  facet_grid(~Type) +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("pink","#EE3E80", "red", "#FF8200", "tomato", 
                             "#FF9966", "gold2", "khaki", "yellow", "lawngreen", 
                             "springgreen3", "seagreen", "darkgreen", "navyblue", "royalblue2",
                             "turquoise", "powderblue", "thistle", "violet", "purple",
                             "darkmagenta", "violetred4", "deeppink", "palevioletred", "rosybrown4",
                             "mistyrose4")) +
  theme_classic() +
  xlab("Time Point") +
  ylab("Genes or Transcripts per Million") +
  guides(fill=guide_legend(ncol=3)) +
  theme(text = element_text(size=20)) +
  theme(legend.position="bottom")
print(per_p)

bq <- aggregate(bq$GPM, by=list(bq$Taxonomy, bq$Time, bq$Type), FUN=sum)
head(bq)
bq <- rename(bq, Taxonomy=Group.1, Time=Group.2, Type=Group.3, Count=x)
bq <- separate(bq, Taxonomy, c("domain", "phylum", "class", "order", "family"),
                 sep=";", remove=TRUE)
head(bq)
bq <- subset(bq, select=-c(class, order))
bq$Taxonomy <- paste(bq$domain, bq$phylum, bq$family, sep=";")

bq_p <- ggplot(bq, aes(Time, Count, fill=Taxonomy)) +
  geom_area() +
  facet_grid(~Type) +
  theme_classic() +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("pink","#EE3E80", "red","#FF9966",
                             "gold2", "yellow", "palegreen", "lawngreen",
                             "springgreen3", "darkgreen", "skyblue4", "navyblue",
                             "royalblue2", "darkturquoise", "violet", "darkmagenta", "deeppink4",
                             "rosybrown4")) +
  xlab("Time Point") +
  ylab("Genes or Transcripts per Million") +
  guides(fill=guide_legend(ncol=3)) +
  theme(legend.position="bottom")
  theme(text = element_text(size=20))
print(bq_p)

plantderived <- subset(CAZy, family=="GH19" | family=="GH23" | family=="GH29" | family=="GH2" | family=="GH51" |
                         family=="GH43" | family=="GT2" | family=="GT4" | family=="CE4" | family=="CBM48")
head(plantderived)
plantderived <- aggregate(plantderived$GPM, by=list(plantderived$type, plantderived$time, plantderived$family), FUN=sum)
plantderived <- rename(plantderived, Type=Group.1, time=Group.2, family=Group.3, GPM=x)

pd_p <- ggplot(plantderived, aes(time, GPM)) +
  geom_point(aes(shape=Type, colour=Type), size=8) +
  facet_wrap(~family, scales="free_y", ncol=2) +
  theme_classic() +
  scale_colour_manual(values=c("springgreen", "darkgreen")) +
  scale_x_continuous(breaks=c(4,8,12)) +
  theme(text = element_text(size=24)) +
  xlab("Time Point") +
  ylab("Genes or Transcripts per Million") +
  theme(legend.position="bottom")
print(pd_p)

chitinase <- subset(CAZy, family=="GH19")
chitinase <- aggregate(chitinase$GPM, by=list(chitinase$taxonomy, chitinase$time, chitinase$type), FUN=sum)
head(chitinase)
tail(chitinase)
chitinase <- rename(chitinase, Taxonomy=Group.1, Time=Group.2, Type=Group.3, Count=x)
chitinase <- separate(chitinase, Taxonomy, c("domain", "phylum", "class", "order", "family"),
               sep=";", remove=TRUE)
head(chitinase)
chitinase <- subset(chitinase, select=-c(class, order))
chitinase$Taxonomy <- paste(chitinase$domain, chitinase$phylum, chitinase$family, sep=";")
chitinase <- subset(chitinase, select=-c(domain, phylum, family))

chitinase_by_family_pad <- expand(chitinase, nesting(Taxonomy, Time, Type), Time)
head(chitinase_by_family_pad)
chitinase <- rename(chitinase, Time1=Time)
chitinase_by_family_pad <- subset(chitinase_by_family_pad, select=-c(Time))

chitinase_by_family_padded <- merge(chitinase_by_family_pad, chitinase, by=c("Taxonomy", "Type", "Time1"), all=TRUE)
head(chitinase_by_family_padded)
tail(chitinase_by_family_padded)
chitinase_by_family_padded[is.na(chitinase_by_family_padded)] <- 0
chitinase_by_family_padded <- unique(chitinase_by_family_padded)
write.csv(chitinase_by_family_padded, "chitinase.csv", row.names=FALSE)


catA_mg <- read.csv("catA_mg_commcomp_mod.csv")
catA_mg$name <- "catA"
catA_mt <- read.csv("catA_mt_commcomp_mod.csv")
catA_mt$name <- "catA"
catE_mg <- read.csv("catE_mg_commcomp_mod.csv")
catE_mg$name <- "catE"
catE_mt <- read.csv("catE_mt_commcomp_mod.csv")
catE_mt$name <- "catE"
boxA_mg <- read.csv("boxA_mg_commcomp_mod.csv")
boxA_mg$name <- "boxA"
boxA_mt <- read.csv("boxA_mt_commcomp_mod.csv")
boxA_mt$name <- "boxA"
boxB_mg <- read.csv("boxB_mg_commcomp_mod.csv")
boxB_mg$name <- "boxB"
boxB_mt <- read.csv("boxB_mt_commcomp_mod.csv")
boxB_mt$name <- "boxB"
gdoA_mg <- read.csv("gentisate_mg_commcomp_mod.csv")
gdoA_mg$name <- "gdoA"
gdoA_mt <- read.csv("gentisate_mt_commcomp_mod.csv")
gdoA_mt$name <- "gdoA"
pcaG_mg <- read.csv("pcaG_mg_commcomp_mod.csv")
pcaG_mg$name <- "pcaG"
pcaG_mt <- read.csv("pcaG_mt_commcomp_mod.csv")
pcaG_mt$name <- "pcaG"
pcaH_mg <- read.csv("pcaH_mg_commcomp_mod.csv")
pcaH_mg$name <- "pcaH"
pcaH_mt <- read.csv("pcaH_mt_commcomp_mod.csv")
pcaH_mt$name <- "pcaH"

require(ggplot2)
require(dplyr)
aromatic_ko_commcomp <- rbind(boxA_mg, boxB_mg, catA_mg, catE_mg, gdoA_mg, pcaG_mg, pcaH_mg, boxA_mt, boxB_mt, catA_mt, catE_mt, gdoA_mt, pcaG_mt, pcaH_mt)
write.csv(aromatic_ko_commcomp, "aromatic_ko_commcomp_051418.csv", row.names=FALSE)

aromatic_ko_commcomp <- read.csv("aromatic_ko_commcomp_051418.csv")
aromatic_ko_commcomp <- subset(aromatic_ko_commcomp, name !="boxA" & name!="catE" & name!="pcaG")
write.csv(aromatic_ko_commcomp, "aromatic_ko_commcomp_reducded_051418.csv")
aromatic_ko_commcomp<-read.csv("aromatic_ko_commcomp_reducded_051418.csv")
head(aromatic_ko_commcomp)

aromatic_ko_commcomp <- aggregate(aromatic_ko_commcomp$GPM, by=list(aromatic_ko_commcomp$taxonomy, aromatic_ko_commcomp$timepoint, aromatic_ko_commcomp$name, aromatic_ko_commcomp$type), FUN=sum)
write.csv(aromatic_ko_commcomp, "ko_aromatic_workinprogress.csv", row.names = FALSE)

aromatic_ko_commcomp <- rename(aromatic_ko_commcomp, Taxonomy=Group.1, Timepoint=Group.2, Gene=Group.3, Type=Group.4, GPM=x)
aromatic_ko_commcomp <- read.csv("ko_aromatic_workinprogress.csv")

p <- ggplot(aromatic_ko_commcomp, aes(Timepoint, GPM, fill=Taxonomy)) + 
            geom_area(position="stack") + 
            scale_fill_manual(values=c("red","#EE3E80", "pink","#FF9966", "#FF8200",
                             "yellow", "lawngreen", "springgreen3", "darkgreen", "palegreen", "skyblue2", "navyblue",
                             "royalblue2", "darkturquoise", "violet", "deeppink4", "darkmagenta",
                             "seashell4")) +
            scale_x_continuous(breaks=c(4,8,12)) +
            facet_wrap(Type~Gene, scales="free", nrow=2) + 
            theme_classic() +
            guides(fill=guide_legend(ncol=3)) +
            theme(text=element_text(size=24)) +
            theme(legend.position = "bottom")
print(p)

box <- subset(aromatic_ko_commcomp, Gene=="boxB")
box_p <- ggplot(box, aes(Timepoint, GPM, fill=Taxonomy)) +
  geom_area(position="stack") +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("pink", "magenta", "#FF9966", "yellow", 
                             "lawngreen", "springgreen4", "royalblue2", "lightskyblue",
                            "thistle2", "purple", "seashell4")) +
  facet_wrap(~Type) +
  theme_classic() +
  xlab("Time Point") +
  ylab("GPM or TPM") +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "bottom") +
  theme(text=element_text(size=24))
print(box_p)

proto <- subset(aromatic_ko_commcomp, Gene=="pcaH")
proto_p <- ggplot(proto, aes(Timepoint, GPM, fill=Taxonomy)) +
  geom_area(position="stack") +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("pink", "hotpink1", "orange1", "yellow", 
                             "darkorchid4", "seashell4")) +
  facet_wrap(~Type) +
  theme_classic() +
  xlab("Time Point") +
  ylab("GPM or TPM") +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "bottom") +
  theme(text=element_text(size=24))
print(proto_p)

gen <- subset(aromatic_ko_commcomp, Gene=="gdoA")
gen_p <- ggplot(gen, aes(Timepoint, GPM, fill=Taxonomy)) +
  geom_area(position="stack") +
  scale_x_continuous(breaks=c(4,8,12)) +
  scale_fill_manual(values=c("red", "pink", "coral", "yellow", 
                            "lawngreen", "purple", "seashell4")) +
  facet_wrap(~Type) +
  xlab("Time Point") +
  ylab("GPM or TPM") +
  theme_classic() +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "bottom") +
  theme(text=element_text(size=24))
print(gen_p)

cat <- subset(aromatic_ko_commcomp, Gene=="catA")
cat_p <- ggplot(cat, aes(Timepoint, GPM, fill=Taxonomy)) +
  geom_area(position="stack") +
  scale_x_continuous(breaks=c(4,8,12)) +
 scale_fill_manual(values=c("violetred4", "pink", "coral", "darkturquoise", "seashell4")) +
  facet_wrap(~Type) +
  xlab("Time Point") +
  ylab("GPM or TPM") +
  theme_classic() +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "bottom") +
  theme(text=element_text(size=24))
print(cat_p)

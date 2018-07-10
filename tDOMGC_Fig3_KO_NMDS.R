require(vegan)
require(ggplot2)
#require(phyloseq)
options(max.print=999999)
abund_table_ko <- read.csv("MGMT_KOnum_GPM.csv", row.names = 1, check.names = F)
head(abund_table_ko)
meta_data <- read.csv("july_metadata.csv", row.names = 1, check.names=F)

abund_table_ko$rowsum <- rowSums(abund_table_ko)
head(abund_table_ko)
abund_table_ko <- subset(abund_table_ko, abund_table_ko$rowsum > 0)
abund_table_ko <- subset(abund_table_ko, select= -c(rowsum))
head(abund_table_ko)
abund_table_mg_ko <- subset(abund_table_ko, select=c(MG_01, MG_09B, MG_13B, MG_02, MG_04, MG_05A, MG_03, MG_05C, MG_06, MG_07, MG_09A, MG_10, MG_11, MG_12A, MG_12B, MG_13A))
head(abund_table_mg_ko)
abund_table_mg_ko$rowsum <- rowSums(abund_table_mg_ko !=0)
mg_single_ko <- subset(abund_table_mg_ko, rowsum==1)
head(abund_table_mg_ko)
abund_table_mg_ko <- subset(abund_table_mg_ko, abund_table_mg_ko$rowsum > 1)
abund_table_mg_ko <- subset(abund_table_mg_ko, select= -c(rowsum))


abund_table_mt_ko <- subset(abund_table_ko, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07, MT_08, MT_09A, MT_10, MT_11, MT_12A, MT_12B, MT_13A, MT_13C))
head(abund_table_mt_ko)
abund_table_mt_ko$rowsum <- rowSums(abund_table_mt_ko !=0)
head(abund_table_mt_ko)
mt_single_ko <- subset(abund_table_mt_ko, rowsum==1)
abund_table_mt_ko <- subset(abund_table_mt_ko, abund_table_mt_ko$rowsum > 1)
abund_table_mt_ko <- subset(abund_table_mt_ko, select= -c(rowsum))


abund_table_ko <- t(abund_table_ko)
abund_table_mg_ko <- t(abund_table_mg_ko)
abund_table_mt_ko <- t(abund_table_mt_ko)

july_dist <- vegdist(abund_table_ko, method ="bray", na.rm=T)
july_dist_mg <- vegdist(abund_table_mg_ko, method="bray", na.rm=T)
july_dist_mt <- vegdist(abund_table_mt_ko, method="bray", na.rm=T)
? metaMDS
july_MDS <- metaMDS(july_dist, distance="bray", k=2, trymax=100, trace=T, plot=T)
july_MDS_mg <- metaMDS(july_dist_mg, distance="bray", k=2, trymax=100)
july_MDS_mt <- metaMDS(july_dist_mt, distance="bray", k=2, trymax=100)
stressplot(july_MDS_mg)
stressplot(july_MDS_mt)
head(july_MDS)
plot(july_MDS)
attach(meta_data)
head(meta_data)
tail(meta_data)
meta_data <- t(meta_data)

meta_data_mg <- subset(meta_data, select=c(MG_01, MG_09B, MG_13B, MG_02, MG_04, MG_05A, MG_03, MG_05C, MG_06, MG_07, MG_09A, MG_10, MG_11, MG_12A, MG_12B, MG_13A))
meta_data_mg <- t(meta_data_mg)

meta_data_mt <- subset(meta_data, select=c(MT_01, MT_02, MT_03, MT_05B, MT_05C, MT_06, MT_07, MT_08, MT_09A, MT_10, MT_11, MT_12A, MT_12B, MT_13A, MT_13C))
meta_data_mt <- t(meta_data_mt)

NMDS <- data.frame(x=july_MDS$point[,1], y=july_MDS$point[,2], timepoint=as.factor(meta_data[,1]),
                   depth=as.factor(meta_data[,2]), temperature=as.factor(meta_data[,3]),
                   salinity=as.factor(meta_data[,4]), a254=as.factor(meta_data[,5]),
                   type=as.factor(meta_data[,6]), tide=as.factor(meta_data[,7]))

NMDS_mg <- data.frame(x=july_MDS_mg$point[,1], y=july_MDS_mg$point[,2], timepoint=as.factor(meta_data_mg[,1]),
                   depth=as.factor(meta_data_mg[,2]), temperature=as.factor(meta_data_mg[,3]),
                   salinity=as.factor(meta_data_mg[,4]), a254=as.factor(meta_data_mg[,5]),
                   type=as.factor(meta_data_mg[,6]), tide=as.factor(meta_data_mg[,7]))

NMDS_mt <- data.frame(x=july_MDS_mt$point[,1], y=july_MDS_mt$point[,2], timepoint=as.factor(meta_data_mt[,1]),
                      depth=as.factor(meta_data_mt[,2]), temperature=as.factor(meta_data_mt[,3]),
                      salinity=as.factor(meta_data_mt[,4]), a254=as.factor(meta_data_mt[,5]),
                      type=as.factor(meta_data_mt[,6]), tide=as.factor(meta_data_mt[,7]))

#ord<-ordiellipse(july_MDS, as.factor(NMDS[,9]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
#?ordiellipse
#dev.off()

#veganCovEllipse <-function(cov, center = c(0,0), scale = 1, npoints = 100) {
    ## Basically taken from the 'car' package: The Cirlce
  #  theta <- (0:npoints) * 2 * pi/npoints
   # Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    #Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    #o <- attr(Q, "pivot")
    #t(center + scale * t(Circle %*% Q[,o])) }

#df_ell <- data.frame()
#for(g in levels(NMDS$tide)){
  #if(g!="" && (g %in% names(ord))){
   # df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$tide==g,],
                                                     #veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                 # ,tide=g))}}
  
#head(df_ell)
#NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$tide),mean)
#head(NMDS)
#tail(NMDS)
#str(NMDS)

NMDS$depth=as.numeric(levels(NMDS$depth))[NMDS$depth]
NMDS$timepoint=as.numeric(levels(NMDS$timepoint))[NMDS$timepoint]
NMDS$a254=as.numeric(levels(NMDS$a254))[NMDS$a254]
NMDS$salinity=as.numeric(levels(NMDS$salinity))[NMDS$salinity]

NMDS_mg$depth=as.numeric(levels(NMDS_mg$depth))[NMDS_mg$depth]
NMDS_mg$timepoint=as.numeric(levels(NMDS_mg$timepoint))[NMDS_mg$timepoint]

NMDS_mt$depth=as.numeric(levels(NMDS_mt$depth))[NMDS_mt$depth]
NMDS_mt$timepoint=as.numeric(levels(NMDS_mt$timepoint))[NMDS_mt$timepoint]

p<-ggplot(NMDS, aes(x = x, y = y)) +
  geom_point(aes(colour = depth, shape = type), size=8) +
  geom_text(aes(x = x, y = y, label = as.character(timepoint))) + 
  ###geom_count(aes(colour=depth, shape=type, size=timepoint)) +
  scale_color_distiller(palette="Spectral") +
  theme_bw()
print(p)

p_mg <-ggplot(NMDS_mg, aes(x = x, y = y)) +
  geom_point(aes(colour = depth), size=8) +
  geom_text(aes(x = x, y = y, label = as.character(timepoint))) + 
  ###geom_count(aes(colour=depth, shape=type, size=timepoint)) +
  scale_color_distiller(palette="Spectral") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(text=element_text(size=20))
print(p_mg)

p_mt <-ggplot(NMDS_mt, aes(x = x, y = y)) +
  geom_point(aes(colour = depth), shape=17, size=8) +
  geom_text(aes(x = x, y = y, label = as.character(timepoint))) + 
  ###geom_count(aes(colour=depth, shape=type, size=timepoint)) +
  scale_color_distiller(palette="Spectral") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(text=element_text(size=20))
print(p_mt)

head(NMDS)
tail(NMDS)
p_time <-ggplot(NMDS, aes(timepoint, depth)) +
  geom_point(aes(colour=depth), shape=15, size=8) +
  geom_line() +
  xlab("Timepoint") +
  ylab("Depth (m)") +
  scale_x_continuous(breaks=c(4, 8, 12)) +
  scale_color_distiller(palette="Spectral") +
  theme_classic() +
  theme(text=element_text(size=20))
print(p_time)

head(meta_data)
head(meta_data_mt)
head(meta_data_mg)
meta_data_mg <- as.data.frame(meta_data_mg)
meta_data_mt <- as.data.frame(meta_data_mt)


meta_data_mg <- subset(meta_data, type="metagenome")
adonis(formula=july_dist_mg ~ Lignin, data=meta_data_mg)

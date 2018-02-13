require(vegan)
require(ggplot2)
require(phyloseq)
options(max.print=999999)
abund_table_ko <- read.csv("MGMT_KOnum_GPMTPM.csv", row.names = 1, check.names = F)
head(abund_table_ko)
meta_data <- read.csv("july_metadata.csv", row.names = 1, check.names=F)

abund_table_ko$rowsum <- rowSums(abund_table_ko)
head(abund_table_ko)

abund_table_ko <- subset(abund_table_ko, abund_table_ko$rowsum > 0)
abund_table_ko <- subset(abund_table_ko, select= -c(rowsum))
abund_table_ko <- t(abund_table_ko)

write.csv(abund_table_ko, "ko_mt_no0.csv")
dim(abund_table_ko)
tail(abund_table_ko)

july_dist <- vegdist(abund_table_ko, method ="bray", na.rm=T)
? metaMDS
july_MDS <- metaMDS(july_dist, distance="bray", k=2, trymax=100, trace=T, plot=T)
stressplot(july_MDS)
head(july_MDS)
plot(july_MDS)
attach(meta_data)
head(meta_data)

NMDS <- data.frame(x=july_MDS$point[,1], y=july_MDS$point[,2], timepoint=as.factor(meta_data[,1]),
                   depth=as.factor(meta_data[,2]), temperature=as.factor(meta_data[,3]),
                   salinity=as.factor(meta_data[,4]), a254=as.factor(meta_data[,5]),
                   type=as.factor(meta_data[,6]), tide=as.factor(meta_data[,7]))
head(NMDS)

ord<-ordiellipse(july_MDS, as.factor(NMDS[,9]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
?ordiellipse
dev.off()

veganCovEllipse <-function(cov, center = c(0,0), scale = 1, npoints = 100) {
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    o <- attr(Q, "pivot")
    t(center + scale * t(Circle %*% Q[,o])) }

df_ell <- data.frame()
for(g in levels(NMDS$tide)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$tide==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,tide=g))}}
  
head(df_ell)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$tide),mean)
head(NMDS)
tail(NMDS)
str(NMDS)

NMDS$depth=as.numeric(levels(NMDS$depth))[NMDS$depth]
NMDS$timepoint=as.numeric(levels(NMDS$timepoint))[NMDS$timepoint]

p<-ggplot(NMDS, aes(x, y)) +
  geom_point(aes(colour=depth, shape=type)) +
  geom_count(aes(colour=depth, shape=type, size=timepoint)) +
  scale_color_distiller(palette="Spectral") +
  theme_bw()
print(p)

geom_count(aes(colour=depth, shape=type, size=timepoint)) +

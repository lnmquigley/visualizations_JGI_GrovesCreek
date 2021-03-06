require(vegan)
require(ggplot2)
#require(phyloseq)
options(max.print=999999)
abund_table_ko <- read.csv("KOnum_TPM.csv", row.names = 1, check.names = F)
head(abund_table_ko)
meta_data <- read.csv("april_meta.csv", row.names = 1, check.names=F)

abund_table_ko$rowsum <- rowSums(abund_table_ko)
head(abund_table_ko)
abund_table_ko <- subset(abund_table_ko, abund_table_ko$rowsum > 0)
abund_table_ko <- subset(abund_table_ko, select= -c(rowsum))
head(abund_table_ko)


abund_table_ko <- t(abund_table_ko)
abund_table_ko <- subset(abund_table_ko, select=c(-2))

filt_dist <- vegdist(abund_table_ko, method ="bray", na.rm=T)

filt_MDS <- metaMDS(filt_dist, distance="bray", k=2, trymax=100, trace=T, plot=T)

stressplot(filt_MDS)

attach(meta_data)
head(meta_data)
tail(meta_data)
meta_data <- t(meta_data)

?adonis
adonis(filt_dist~time*month, data=meta_data)

NMDS <- data.frame(x=filt_MDS$point[,1], y=filt_MDS$point[,2], timepoint=as.factor(meta_data[1,]),
                   time=as.factor(meta_data[2,]), month=as.factor(meta_data[3,]))

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
#NMDS$timepoint=as.numeric(levels(NMDS$timepoint))[NMDS$timepoint]
NMDS$a254=as.numeric(levels(NMDS$a254))[NMDS$a254]
NMDS$salinity=as.numeric(levels(NMDS$salinity))[NMDS$salinity]

p<-ggplot(NMDS, aes(x = x, y = y)) +
  geom_point(aes(colour = time, shape=month), size=12) +
  scale_color_manual(values=c("#FF8200", "#58595B"))+
  geom_text(aes(x = x, y = y, label = as.character(timepoint)), colour="white") + 
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic() +
  theme(text = element_text(size=20))
print(p)


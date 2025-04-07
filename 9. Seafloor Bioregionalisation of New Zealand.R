##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 9 - DEVELOPMENT OF A SEAFLOOR BIOREGIONALISATION FOR NEW ZEALAND

##------ Authors: Fabrice Stephenson
##------ Start date : 01/07/2019
##------ End date : 01/07/2023

# The code to develop a Seafloor Bioregionaliation of New Zealand was funded by 
# the New Zealand Department of Conservation and Fisheries New Zealand 
# (Ministry for Primary Industries) described in:
# Stephenson et al. (2023). A seafloor bioregionalisation for New Zealand. 
#     Ocean & Coastal Management 242, 106688.

##============================================================================##

# DESCRIPTION: Classification of compositional turnover into a low number of 
# groups representing the lowest number of groups needed to describe geographic 
# regions that were relatively homogeneous and distinct in terms of their 
# environmental conditions and in turn biological contents. The discrimination 
# across classification levels was assessed using the biological data included 
# in the GF models in an analysis of similarities test (ANOSIM)

# 1.  Load files and packages
# 2.  Heirarchical classification at different levels
# 3.  Export Classification as rasters
# 4.  Taxa and environmental values within groups
# 5.  Taxa anosim analysis to explore group descrimination 

####==========    1. LOAD FILES AND PACKAGES   =============================####
require(raster); require(cluster); require(devEMF); require(tidyverse)
require(vegan); require(ecodist); require(parallel)

# Load MPI projection
MPIproj<- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# load master predictor stack 
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

imp.vars <- c("Bathy","BedDist", "BotOxy", "BotNi", "BotPhos","BotSal",
              "BotSil", "BotTemp","BPI_broad", "BPI_fine", "ChlAGrad",
              "DET", "PB555nm","SeasTDiff", "Slope", "SSTGrad", "sed.class", 
              "TC", "POCFlux","Ebed")

setwd(paste0(dir, "/Test data"))
mask <- raster("Template_1km.tif")
mask[values(mask)>0] <- 1
plot(mask)

# load 1km predictors for whole area; extract for imp.var
load("Pred_1km.CMB.source")
head(Pred_1km.CMB)
Pred_1km <- Pred_1km.CMB; rm(Pred_1km.CMB)

# generate a depth mask if needed
Pred_1km$Bathy[Pred_1km$Bathy < 0] <- 0

# bring in the turnover object for multi-level classification
load("Pred_EEZ_CMB.source")
Turnover_CMB <- Pred_EEZ_CMB; rm(Pred_EEZ_CMB)

# Load biological data for taxa of interest
load("DF.source") # DF bio 
load("RF.source") # RF bio
load("MA.source") # MA bio
load("BI.source") # BI bio

####==========    2. HEIRARCHICAL CLASSIFICATION AT DIFFERENT LEVELS    ====####
# loop through Transformed environmental space and extract each class 
# 4 - 20 at 1 km

# setwd("O:/DOC19208/Working/R/COMBINING_MODELS/classification")
load("EEZClaraClassification.source")

# now create a 500 row dataset summarising the transformed envrionmental attributes for the groups from clara
EEZMedoidMeans <- matrix(0, nrow = 500, ncol = length(imp.vars))

dimnames(EEZMedoidMeans)[[1]] <- paste('Grp_',c(1:500),sep='')
dimnames(EEZMedoidMeans)[[2]] <- imp.vars

for (i in c(1:length(imp.vars))) EEZMedoidMeans[,i] <- tapply(Turnover_CMB[,imp.vars[[i]]],
                                                              EEZClaraClassification[[4]],mean)

summary(EEZMedoidMeans)
# and apply a hierarchical classification to it using agnes
EEZMedoidAgnesClassification <- agnes(EEZMedoidMeans, 
                                      metric = 'manhattan',
                                      method = 'gaverage',
                                      par.method = -0.1)

# plot(EEZMedoidAgnesClassification, which.plots = 2)
# rect.hclust(EEZMedoidAgnesClassification, k=25, border="red")
# rect.hclust(EEZMedoidAgnesClassification, k=75, border="blue")

# now reduce to a smaller number of groups based on the clara results
ClaraGroupExpansion <- cutree(EEZMedoidAgnesClassification,1:500)

# head(ClaraGroupExpansion)
i <- match(EEZClaraClassification$clustering,ClaraGroupExpansion[,500])
# summary(i)

EEZClaraClassification <- EEZClaraClassification[-c(11:length(EEZClaraClassification))]

Group.Num <- seq(3,20,1)
Group.Name <- paste("Grp_",Group.Num, sep = "")

for (j in 1:length(Group.Num)){
  index <- j + 10
  EEZClaraClassification[[index]] <- ClaraGroupExpansion[i,Group.Num[j]]
  names(EEZClaraClassification)[index] <- Group.Name[j]
}

# setwd("O:/DOC19208/Working/The Seafloor Community Classification (SCC)/Classification_assessment")
save(EEZClaraClassification,file = 'EEZClaraClassification_3_20.source')

####==========    3.EXPORT CLASSIFICATION AS RASTERS    ====================####
load('EEZClaraClassification_3_20.source')

Group.Num <- seq(3,20,1)
Group.Name <- paste("Grp_",Group.Num, sep = "")

setwd(paste0(dir, "/Outputs"))
# k = 5
for (k in 1:length(Group.Num)){
  Means <- matrix(0, nrow = Group.Num[k], ncol = length(imp.vars))
  dimnames(Means)[[1]] <- paste('Grp_', c(1:Group.Num[k]),sep='')
  dimnames(Means)[[2]] <- imp.vars
  
  for (i in c(1:length(imp.vars))) Means[,i] <- tapply(Turnover_CMB[,imp.vars[[i]]],
                                                               EEZClaraClassification[[Group.Name[k]]],mean)
  EEZClusterPCA <- prcomp(Means)
  
  # set up colours using the same PCA space
  a1 <- EEZClusterPCA$x[, 1]
  a2 <- EEZClusterPCA$x[, 2]
  a3 <- EEZClusterPCA$x[, 3]
  r <- a1 + a2
  g <- -a2
  b <- a3 + a2 - a1
  r <- (r - min(r))/(max(r) - min(r)) * 255
  g <- (g - min(g))/(max(g) - min(g)) * 255
  b <- (b - min(b))/(max(b) - min(b)) * 255
  
  #####   3.1 Plot and save PCA   ----------------------------------------------
  nvs <- dim(EEZClusterPCA$rotation)[1]
  vec <- imp.vars[c(1,3,6,8,9,19)]
  # vec <- imp.vars[c(1:20)]
  lv <- length(vec)
  vind <- rownames(EEZClusterPCA$rotation) %in% vec
  
  scal <- 15
  xrng <- range(EEZClusterPCA$x[, 1], EEZClusterPCA$rotation[, 1]/scal) * + 1.1
  yrng <- range(EEZClusterPCA$x[, 2], EEZClusterPCA$rotation[, 2]/scal) * + 1.1
  
  variableCEX <- (as.numeric(table(EEZClaraClassification[[Group.Name[k]]]))^0.110) * 2
  
  # save as EMF
  dir.create(paste(Group.Name[k],sep =""))
  setwd(paste(dir, "/Outputs/", Group.Name[k],sep =""))
  emf(file = paste(Group.Name[k], "_EEZ_PCA.emf"), emfPlus = T)
  plot((EEZClusterPCA$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = variableCEX, 
       col = rgb(r, g, b, max = 255), asp = 1)
  arrows(rep(0, lv), rep(0, lv), EEZClusterPCA$rotation[vec,1]/scal, EEZClusterPCA$rotation[vec, 2]/scal, length = 0.06210)
  jit <- 0.00110
  
  text(EEZClusterPCA$rotation[vec, 1]/scal + jit * sign(EEZClusterPCA$rotation[vec, 1]), 
       EEZClusterPCA$rotation[vec, 2]/scal + jit * sign(EEZClusterPCA$rotation[vec, 2]), labels = vec)
  text(EEZClusterPCA$x[,1], EEZClusterPCA$x[,2], seq(1,Group.Num[k],1),
       cex = variableCEX/7, #0.8,
       pos = 4,
       adj = c(0.2,-0.1))
  dev.off() # finishes the plot and saves
  
  #####   3.2 Spatial prediction   --------------------------------------------- 
  # create a dataframe with spatial information for th classification 
  CMB_DF.1km <- cbind(Turnover_CMB[,c(1:3)],EEZClaraClassification[[Group.Name[k]]])
  
  # export as raster
  CMB_DF.1km.R <- rasterFromXYZ(data.frame(x = CMB_DF.1km[,1],
                                           y = CMB_DF.1km[,2],
                                           z = CMB_DF.1km[,4]),
                                       crs = MPIproj)
  plot(CMB_DF.1km.R)

  ##### PLOT MAP 
  breakpoints <- seq(1,Group.Num[k],1)
  rgb <- (cbind(r,g,b)/255)
  colors <- rgb(rgb)
 
  jpeg(filename = paste(Group.Name[k], "_EEZ_Spatial.jpeg", sep = ""), width = 11000, height = 11000, quality = 100, bg = "white", res = 1000)
  par(mar=c(2,2,2,2))
  plot(CMB_DF.1km.R,col=colors, legend = F)
  title(paste(Group.Name[k]))
  dev.off()
  
  writeRaster(CMB_DF.1km.R, filename= paste(Group.Name[k], ".tif", sep = ""), 
              format = "GTiff", 
              overwrite = TRUE)
  
  #####   3.3 Tag biological data    -------------------------------------------
  # DF
  DF.class <- as.data.frame(raster::extract(CMB_DF.1km.R,DF[,c("X","Y")]))
  colnames(DF.class) <- Group.Name[k]
  DF <- cbind(DF, DF.class)
  # RF
  RF.class <- as.data.frame(raster::extract(CMB_DF.1km.R,RF[,c("X","Y")]))
  colnames(RF.class) <- Group.Name[k]
  RF <- cbind(RF, RF.class)
  # MA
  MA.class <- as.data.frame(raster::extract(CMB_DF.1km.R,MA[,c("X","Y")]))
  colnames(MA.class) <- Group.Name[k]
  MA <- cbind(MA, MA.class)
  # DF
  BI.class <- as.data.frame(raster::extract(CMB_DF.1km.R,BI[,c("X","Y")]))
  colnames(BI.class) <- Group.Name[k]
  BI <- cbind(BI, BI.class)
  
  print(paste("FINISHED: GROUP ",Group.Num[k], sep ="" ))
}

setwd(paste(dir, "/Outputs/", sep =""))
save(DF, file = "DF.class.source")
save(RF, file = "RF.class.source")
save(MA, file = "MA.class.source")
save(BI, file = "BI.class.source")

####==========    4. TAXA AND ENVIRONMENTAL VALUE SUMMARY METRICS    =======####
# ENVIRONMENTAL SUMMARY
setwd(paste(dir, "/Outputs/", "Grp_7",sep =""))
BR <- raster("Grp_7.tif")
plot(BR)
BR_env <- as.data.frame(raster::extract(BR,Pred_1km[,c(1,2)]))
BR_env <- round(BR_env,0)
BR_env <- as.factor(BR_env[,1])
BR_pred_1km <- cbind(BR_env, Pred_1km)

# boxplot of variables by group
BR_pred_1km_plt <- BR_pred_1km[,c(1,4,6,11,22)]
devEMF::emf(file = "Bioregion/Env_conditions.emf", emfPlus = T)
par(mfrow=c(1,3))
par(mar=c(2,4.2,1,1))
boxplot(Bathy ~ BR_env, data=BR_pred_1km_plt,outline=FALSE,
        xlab="Group number", ylab="Bathy (m)",par(cex.lab=1.0),
        ylim = c(0,1500))
boxplot(BotOxy ~ BR_env, data=BR_pred_1km_plt,outline=FALSE,
        xlab="Group number", ylab= expression('Dissolved oxygen at depth (ml l ' ^-1* ")"),par(cex.lab=1.1))
boxplot(BotTemp ~ BR_env, data=BR_pred_1km_plt,outline=FALSE,
        xlab="Group number", ylab=expression("Temperature at depth ("*~degree*C* " km"^-1*")"),par(cex.lab=1.0))
dev.off()

# Median ± 5 - 95% quantile of environmental values
BR_pred_1km_MEAN <- aggregate(BR_pred_1km[,imp.vars], list(BR_pred_1km$BR_env), median)
BR_pred_1km_MEAN[,2:21] <- round(BR_pred_1km_MEAN[,2:21], 3)
BR_pred_1km_quant5 <- aggregate(BR_pred_1km[,imp.vars], list(BR_pred_1km$BR_env), FUN = 'quantile', probs=c(0.25))
BR_pred_1km_quant95 <- aggregate(BR_pred_1km[,imp.vars], list(BR_pred_1km$BR_env), FUN = 'quantile', probs=c(0.75))
colnames(BR_pred_1km_MEAN)[1] <- "GF_Class"

setwd(paste(dir, "/Outputs/", sep =""))
write.csv(BR_pred_1km_MEAN, file = "BR_ENV_MedianF.csv")
write.csv(BR_pred_1km_quant5, file = "Bioregion/BR_ENV5.csv")
write.csv(BR_pred_1km_quant95, file = "Bioregion/BR_ENV95.csv")

# SPECIES SUMMARY
species.DF <- colnames(DF[,23:342])
BR_fac <- as.data.frame(raster::extract(BR, DF[,c(1,2)]))
BR_fac <- as.factor(BR_fac[,1])
BR_samp <- cbind(BR_fac, DF)

# loop through and create means by factor classifcation for summary per group
BR_spe <- BR_samp[,species.DF]
indx <- sapply(BR_spe, is.factor)
BR_spe[indx] <- lapply(BR_spe[indx], function(x) as.numeric(as.character(x)))
BR_preds <- BR_samp[,imp.vars]
BR_spe <- cbind(BR_fac, BR_spe)

# calculat the mean species occurrence
mean.Pres <- BR_spe %>%   #
  group_by(BR_fac) %>%    #Should only leave GF
  summarise_all(funs(mean), na.rm=T)

setwd(paste(dir, "/Outputs/", "Grp_7",sep =""))
write.csv(mean.Pres, file = "Outputs/BR_MEAN.Pres.csv")

# count number of samples per group
FishFullGF.cnt <- as.data.frame(BR_samp[,1])
names(FishFullGF.cnt) <- "GF_grp"
FishFullGF.cnt$n.samples <- 1

mean.cnt1 <- FishFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

# count number of cells per GF group (extent)
FishFullGF.cnt <- as.data.frame(BR_pred_1km[,1])
names(FishFullGF.cnt) <- "GF_grp"
FishFullGF.cnt$extent <- 1

mean.cnt2 <- FishFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

mean.cnt_F<- merge(mean.cnt1,mean.cnt2, all = TRUE)

setwd(paste(dir, "/Outputs/", "Grp_7",sep =""))
write.csv(mean.cnt_F, file = "BR.mean.cnt_F.csv")

# count number of raster cells in each group  (extent)
BR_FINAL <- cbind(BR_preds, BR_spe)
BR_MEAN <- aggregate(BR_FINAL[,1:ncol(BR_FINAL)], list(BR_samp$BR_fac), mean)
colnames(BR_MEAN)[1] <- "GF_Class"

setwd(paste(dir, "/Outputs/", "Grp_7",sep =""))
write.csv(BR_MEAN, file = "BR_MEAN.csv")


####==========    5. TAXA ANOSIM FOR EACH GROUP     ========================####
source("./Pairwise_adonis.R")

# biological groups with class number tags
setwd(paste(dir, "/Outputs/", sep =""))
load("DF.class.source") # DF bio 
load("RF.class.source") # RF bio
load("MA.class.source") # MA bio
load("BI.class.source") # BI bio

Group.Num <- seq(3,20,1)
Group.Name <- paste("Grp_",Group.Num, sep = "")

####    5.1 Reef fish   --------------------------------------------------------
RF <- na.omit(RF)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)

for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(RF[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  RF.cut <- RF[RF[,Group.Name[j]] %in% grps.sum,]
  
  RF.spe <- RF.cut[,c(24:(ncol(RF.cut)-(length(Group.Num)+1)))]
  
  for (i in 1:ncol(RF.spe)){RF.spe[,i] <- as.numeric(RF.spe[,i])}
  RF.spe[RF.spe==1] <- 0
  RF.spe[RF.spe==2] <- 1
  RF.spe <- RF.spe[,colSums(RF.spe[,1:length(RF.spe)]) > 0]
  RF.spe <- RF.spe[rowSums(RF.spe[1:nrow(RF.spe),]) > 0,]
  
  RF.class <- RF.cut[,c((ncol(RF.cut)-length(Group.Num)):ncol(RF.cut))]
  RF.class <- RF.class[rownames(RF.class) %in% rownames(RF.spe), ]
  
  RF.dist <- vegdist(RF.spe, distance="jaccard", binary = T)
  
  # ANOSIM
  perm <- anosim(RF.dist, RF.class[,Group.Name[j]], permutations = 100)
  # PERMANOVA
  # perm <- adonis(RF.dist ~ RF.class[,Group.Name[j]], data = RF.class, permutations = 100)
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  # if (any(j == c(1,5,10,15,25))){
  #   Perm.pair <- pairwise.adonis(RF.dist,
  #                                RF.class[,Group.Name[j]], perm = 25, parallel = 5)
  #   gc()
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # } else {
  #   Summ.table[j,3] <- 0
  # }
  # Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  
  Summ.table[j,4] <- perm$statistic
  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  Summ.table[j,5] <- perm$signif
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}

emf(file = "RF_R2.emf", emfPlus = T)
plot(Summ.table$Num.Classes[2:18], Summ.table$R2[2:18], xlab="Group number ", ylab="R2", pch=19)
dev.off()
write.csv(Summ.table, file = "RF_scores.csv")

# dispersion if you need to look at this
disp.grp <- betadisper(RF.dist, RF.class$Grp_150, type = "median")
pairwise <- permutest(disp.grp, pairwise=TRUE, permutations=100)
pairwise$pairwise[1]
# plot(disp.grp, ellipse = F, hull = T, conf = 0.40) # 90% data ellips
# boxplot(disp.grp)

####    5.2 Macroalgae    ------------------------------------------------------
MA <- na.omit(MA)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)

for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(MA[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  MA.cut <- MA[MA[,Group.Name[j]] %in% grps.sum,]
  
  MA.spe <- MA.cut[,c(24:(ncol(MA.cut)-(length(Group.Num)+1)))]
  
  for (i in 1:ncol(MA.spe)){MA.spe[,i] <- as.numeric(MA.spe[,i])}
  MA.spe[MA.spe==1] <- 0
  MA.spe[MA.spe==2] <- 1
  MA.spe <- MA.spe[,colSums(MA.spe[,1:length(MA.spe)]) > 0]
  MA.spe <- MA.spe[rowSums(MA.spe[1:nrow(MA.spe),]) > 0,]
  
  MA.class <- MA.cut[,c((ncol(MA.cut)-length(Group.Num)):ncol(MA.cut))]
  MA.class <- MA.class[rownames(MA.class) %in% rownames(MA.spe), ]
  
  MA.dist <- bcdist(MA.spe)
  perm <- anosim(MA.dist, MA.class[,Group.Name[j]], permutations = 100, parallel = 5)
  # perm <- adonis(MA.dist ~ MA.class[,Group.Name[j]], data = MA.class, permutations = 100, parallel = 10)
 
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  # if (any(j == c(1,5,10,15,25))){
  #   Perm.pair <- pairwise.adonis(MA.dist,
  #                                MA.class[,Group.Name[j]], perm = 25, parallel = 5)
  #   gc()
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # } else {
  #   Summ.table[j,3] <- 0
  # }
  # Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)

  Summ.table[j,4] <- perm$statistic
  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  Summ.table[j,5] <- perm$signif
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}
# setwd("O:/DOC19208/Working/The Seafloor Community Classification (SCC)/Classification_assessment")
emf(file = "MA_R2.emf", emfPlus = T)
plot(Summ.table$Num.Classes[2:18], Summ.table$R2[2:18], xlab="Group number ", ylab="R2", pch=19)
dev.off()
write.csv(Summ.table, file = "MA_scores.csv")

####    5.3 Demersal fish   ----------------------------------------------------
DF <- na.omit(DF)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)
start <- Sys.time() # recording the time

# for (j in 1:length(Group.Num)){
#   grps.sum <- as.data.frame(table(DF[,Group.Name[j]]))
#   grps.sum <- grps.sum[grps.sum[,2] > 4,1]
#   
#   DF.cut <- DF[DF[,Group.Name[j]] %in% grps.sum,]
#   
#   DF.spe <- DF.cut[,c(24:(ncol(DF.cut)-(length(Group.Num)+1)))]
#   for (i in 1:ncol(DF.spe)){DF.spe[,i] <- as.numeric(DF.spe[,i])}
#   DF.spe[DF.spe==1] <- 0
#   DF.spe[DF.spe==2] <- 1
#   DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]
#   DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,]
#   
#   DF.class <- DF.cut[,c((ncol(DF.cut)-length(Group.Num)):ncol(DF.cut))]
#   DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe), ]
#   
#   DF.dist <- bcdist(DF.spe, rmzero = T)
#   gc()
#   perm <- anosim(DF.dist, DF.class[,Group.Name[j]], permutations = 100, parallel = 5)
#   # perm <- adonis(DF.dist ~ DF.class[,Group.Name[j]], permutations = 100, parallel = 5)
#   gc()
#   
#   Summ.table[j,1] <- Group.Num[j]
#   Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
#   
#   # if (any(j == c(1,5,10,15,25))){
#   #   Perm.pair <- pairwise.adonis(DF.dist,
#   #                                DF.class[,Group.Name[j]], perm = 25, parallel = 5)
#   #   gc()
#   #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
#   # } else {
#   #   Summ.table[j,3] <- 0
#   # }
# 
#   # Summ.table[j,4] <- perm$aov.tab$R2[1]
#   Summ.table[j,4] <- perm$statistic
#   Summ.table[j,5] <- perm$signif
#   # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
#   
#   print(paste("Iteration finished for ", Group.Num[j], sep =""))
# }
# end <- Sys.time()
# end - start

# SUBSAMPLING NEEDED

start <- Sys.time()
for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(DF[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]

  DF.cut <- DF[DF[,Group.Name[j]] %in% grps.sum,]

  for (k in 1:length(grps.sum)){
    grp <- grps.sum[k]
    count <- as.data.frame(table(DF[,Group.Name[j]]))
    samp <- DF.cut[DF.cut[,Group.Name[j]] == grp,]
    if(nrow(samp) > 1500){
      sub.samp <- samp[sample(nrow(samp), 1500, replace = F), ]
    } else {
      sub.samp <-samp
    }
    if (k == 1){DF.samp <-sub.samp} else {DF.samp <- rbind(DF.samp, sub.samp)}
    }

  DF.cut <- DF.samp

  DF.spe <- DF.cut[,c(24:(ncol(DF.cut)-(length(Group.Num)+1)))]
  for (i in 1:ncol(DF.spe)){DF.spe[,i] <- as.numeric(DF.spe[,i])}
  DF.spe[DF.spe==1] <- 0
  DF.spe[DF.spe==2] <- 1
  DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]
  DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,]

  DF.class <- DF.cut[,c((ncol(DF.cut)-length(Group.Num)):ncol(DF.cut))]
  DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe), ]

  DF.dist <- bcdist(DF.spe, rmzero = T)
  gc()
  perm <- anosim(DF.dist, DF.class[,Group.Name[j]], permutations = 100, parallel = 5)
  # perm <- adonis(DF.dist ~ DF.class[,Group.Name[j]], permutations = 100, parallel = 5)
  gc()
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  
  # if (any(j == c(1,5,10,15,25))){
  #   Perm.pair <- pairwise.adonis(DF.dist,
  #                                DF.class[,Group.Name[j]], perm = 25, parallel = 5)
  #   gc()
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # } else {
  #   Summ.table[j,3] <- 0
  # }
  
  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  Summ.table[j,4] <- perm$statistic
  Summ.table[j,5] <- perm$signif
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}
end <- Sys.time()
end - start

emf(file = "DF_R2.emf", emfPlus = T)
plot(Summ.table$Num.Classes[2:18], Summ.table$R2[2:18], xlab="Group number ", ylab="R2", pch=19)
dev.off()
write.csv(Summ.table, file = "DF_scores.csv")

####    5.4 Benthic invertebrates   --------------------------------------------
BI <- na.omit(BI)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)
start <- Sys.time() # recording the time

# for (j in 1:length(Group.Num)){
#   grps.sum <- as.data.frame(table(BI[,Group.Name[j]]))
#   grps.sum <- grps.sum[grps.sum[,2] > 4,1]
#   
#   BI.cut <- BI[BI[,Group.Name[j]] %in% grps.sum,]
#   
#   BI.spe <- BI.cut[,c(3:(ncol(BI.cut)-(length(Group.Num)+1)))]
#   for (i in 1:ncol(BI.spe)){BI.spe[,i] <- as.numeric(BI.spe[,i])}
#   # BI.spe[BI.spe==1] <- 0
#   # BI.spe[BI.spe==2] <- 1
#   BI.spe <- BI.spe[,colSums(BI.spe[,1:length(BI.spe)]) > 0]
#   BI.spe <- BI.spe[rowSums(BI.spe[1:nrow(BI.spe),]) > 0,]
#   
#   BI.class <- BI.cut[,c((ncol(BI.cut)-29):ncol(BI.cut))]
#   BI.class <- BI.class[rownames(BI.class) %in% rownames(BI.spe), ]
#   
#   BI.dist <- bcdist(BI.spe, rmzero = T)
#   perm <- anosim(BI.dist, BI.class[,Group.Name[j]], permutations = 100, parallel = 5)
#   # perm <- adonis(BI.dist ~ BI.class[,Group.Name[j]], permutations = 100, parallel = 5)
#   
#   Summ.table[j,1] <- Group.Num[j]
#   Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
#   Summ.table[j,4] <- perm$statistic
#   Summ.table[j,5] <- perm$signif
#   
#   # Summ.table[j,4] <- perm$aov.tab$R2[1]
#   # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
#   # 
#   # if (any(j == c(1,2,3,4,5,10,15,20,30))){
#   #   for (k in 1:length(grps.sum)){
#   #     grp <- grps.sum[k]
#   #     count <- as.data.frame(table(BI.cut[,Group.Name[j]]))
#   #     samp <- BI.cut[BI.cut[,Group.Name[j]] == grp,]
#   #     if(nrow(samp) > 100){
#   #       sub.samp <- samp[sample(nrow(samp), 100, replace = F), ]
#   #     } else {
#   #       sub.samp <-samp
#   #     }
#   #     if (k == 1){BI.samp <-sub.samp} else {BI.samp <- rbind(BI.samp, sub.samp)}
#   #   }
#     
#     # BI.cut <- BI.samp
#     
#   #   BI.spe <- BI.cut[,c(3:(ncol(BI.cut)-32))]
#   #   for (i in 1:ncol(BI.spe)){BI.spe[,i] <- as.numeric(BI.spe[,i])}
#   #   # BI.spe[BI.spe==1] <- 0
#   #   # BI.spe[BI.spe==2] <- 1
#   #   BI.spe <- BI.spe[,colSums(BI.spe[,1:length(BI.spe)]) > 0]
#   #   BI.spe <- BI.spe[rowSums(BI.spe[1:nrow(BI.spe),]) > 0,]
#   #   
#   #   BI.class <- BI.cut[,c((ncol(BI.cut)-29):ncol(BI.cut))]
#   #   BI.class <- BI.class[rownames(BI.class) %in% rownames(BI.spe), ]
#   #   
#   #   BI.dist <- bcdist(BI.spe, rmzero = T)
#   #   Perm.pair <- pairwise.adonis(BI.dist,BI.class[,Group.Name[j]], perm = 100)
#   #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
#   #   } else {
#   #   Summ.table[j,3] <- 0
#   # }
#   
#   print(paste("Iteration finished for ", Group.Num[j], sep =""))
# }
# end <- Sys.time()
# end - start

# SUBSAMPLING NEEDED

for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(BI[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  BI.cut <- BI[BI[,Group.Name[j]] %in% grps.sum,]
 
  for (k in 1:length(grps.sum)){
       grp <- grps.sum[k]
       count <- as.data.frame(table(BI.cut[,Group.Name[j]]))
       samp <- BI.cut[BI.cut[,Group.Name[j]] == grp,]
       if(nrow(samp) > 1500){
         sub.samp <- samp[sample(nrow(samp), 1500, replace = F), ]
       } else {
         sub.samp <-samp
       }
       if (k == 1){BI.samp <-sub.samp} else {BI.samp <- rbind(BI.samp, sub.samp)}
     }

   BI.cut <- BI.samp
   

     BI.spe <- BI.cut[,c(3:(ncol(BI.cut)-(length(Group.Num)+1)))]
     BI.spe <- BI.spe[,colSums(BI.spe[,1:length(BI.spe)]) > 0]
     BI.spe <- BI.spe[rowSums(BI.spe[1:nrow(BI.spe),]) > 0,]

     BI.class <- BI.cut[,c((ncol(BI.cut)-length(Group.Num)):ncol(BI.cut))]
     BI.class <- BI.class[rownames(BI.class) %in% rownames(BI.spe), ]
  
     BI.dist <- bcdist(BI.spe, rmzero = T)
  
     perm <- anosim(BI.dist, BI.class[,Group.Name[j]], permutations = 100, parallel = 5)
     # perm <- adonis(BI.dist ~ BI.class[,Group.Name[j]], permutations = 100, parallel = 5)
     
     Summ.table[j,1] <- Group.Num[j]
     Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
     Summ.table[j,4] <- perm$statistic
     Summ.table[j,5] <- perm$signif
     
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}
end <- Sys.time()
end - start

emf(file = "BI_R2.emf", emfPlus = T)
plot(Summ.table$Num.Classes[2:18], Summ.table$R2[2:18], xlab="Group number ", ylab="R2", pch=19)
dev.off()
write.csv(Summ.table, file = "BI_scores.csv")

library(ecospat)

#load data
dat<-read.csv("NicheOverlapData.csv",h=T,sep=";")

#subset wanted season (calving, summer or winter)
dat.sub<-subset(dat,season=="calving")

#subset reindeer in given season
rein.dat<-subset(dat.sub,cat.col=="reindeer")  
#subset other wanted species for niche analyses
spec.dat<-subset(dat.sub,cat.col=="moose")

    
#set names of columns to keep in the analyses
covars<-names(rein.dat)[c(3:8)]
    
#Run PCA
pca.rm <- dudi.pca(rbind(rein.dat,spec.dat)[,covars],scannf=F,nf=2)

        #redefine function with new title
        ecospat.plot.contrib.mod<-function (contrib, eigen) 
        {
          if (ncol(contrib) == 1) {
            h <- c(unlist(contrib))
            n <- row.names(contrib)
            barplot(h, space = 0, names.arg = n)
            title(main = "variable contribution")
          }
          if (ncol(contrib) == 2) {
            s.corcircle(contrib[, 1:2]/max(abs(contrib[, 1:2])), 
                        grid = FALSE)
            title(main = "Reindeer and moose, calving", sub = paste("axis1 = ", 
                  round(eigen[1]/sum(eigen) * 100, 2), "%", "axis2 = ", 
                  round(eigen[2]/sum(eigen) * 100, 2), "%"))
          }
        }
        
#Plot variable contribution
ecospat.plot.contrib.mod(contrib=pca.rm$co, eigen=pca.rm$eig)

        
#PCA scores for the whole study area
scores.glob <- pca.rm$li

#PCA scores for reindeer calving area occurence data
scores.sp.rein.dat <- suprow(pca.rm,rein.dat[which(rein.dat[,9]==1),covars])$li
#PCA scores for moose calving area occurence data
scores.sp.spec.dat <- suprow(pca.rm,spec.dat[which(spec.dat[,9]==1),covars])$li

#PCA scores for the whole reindeer calving environment
scores.glob.rein.dat <- suprow(pca.rm,rein.dat[,covars])$li
#PCA scores for the whole moose calving environment
scores.glob.spec.dat <- suprow(pca.rm,spec.dat[,covars])$li
        
# gridding reindeer seasonal niche
grid.rein.dat <- ecospat.grid.clim.dyn(glob=scores.glob,
                                        glob1=scores.glob.rein.dat,
                                        sp=scores.sp.rein.dat, R=100,
                                        th.sp=0)

#plot
par(mfrow=c(1,2))
plot(grid.rein.dat$z.uncor,main="Reindeer, calving",sub="density of occurrence")
plot(grid.rein.dat$z.cor,sub="occupancy of the environment")

        
# gridding moose seasonal niche
grid.spec.dat <- ecospat.grid.clim.dyn(glob=scores.glob,
                                        glob1=scores.glob.spec.dat,
                                        sp=scores.sp.spec.dat, R=100,
                                        th.sp=0)

par(mfrow=c(1,2))
plot(grid.spec.dat$z.uncor,main="Moose, calving",sub="density of occurrence")
plot(grid.spec.dat$z.cor,sub="occupancy of the environment")

  
#niche overlap 
ecospat.niche.overlap(grid.rein.dat, grid.spec.dat, cor = T)

#visualizing the niche dynamics between species
par(mfrow=c(1,1))
ecospat.plot.niche.dyn(grid.rein.dat, grid.spec.dat, quant=c(0.1), 
                        transparency=60,
                        title= paste("Niche overlap reindeer and moose, calving"), 
                        name.axis1="Environmental gradient 1",
                        name.axis2="Environmental gradient 2",
                        col.unf = "orange",
                        col.exp = "olivedrab",
                        col.stab = "sienna",
                        colZ1 = "orange",
                        colZ2 = "olivedrab")

        
#niche dynamics along one environmental gradient

#change "vari" to the other variables in "covars" for more plots
vari<-"elevation"
        
# gridding reindeer niche
grid.uni.rein.dat <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(rein.dat,spec.dat)[,vari]),
                                            glob1=as.data.frame(rein.dat[,vari]),
                                            sp=as.data.frame(rein.dat[which(rein.dat[,9]==1),vari]),
                                            R=1000, th.sp=0)

# gridding moose niche
grid.uni.spec.dat <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(rein.dat,spec.dat)[,vari]),
                                            glob1=as.data.frame(spec.dat[,vari]),
                                            sp=as.data.frame(spec.dat[which(spec.dat[,9]==1),vari]),
                                            R=1000, th.sp=0)

ecospat.plot.niche.dyn(grid.uni.rein.dat, grid.uni.spec.dat, quant=0,
                        transparency=30,
                        title= "Niche Overlap",
                        name.axis1="Elevation",
                       col.unf = "sienna3",
                       col.exp = "olivedrab",
                       col.stab = "sienna",
                       colZ1 = "tan2",
                       colZ2 = "olivedrab")
        

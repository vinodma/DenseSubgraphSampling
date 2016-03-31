setwd("/Users/vmangipudi/graphsampling/data")
dat <- read.csv("rmat_pg.csv")
dat<-dat[-1,]
dat$Samplingrate <- round(dat$Samplingrate,2)

#dat$Samplingrate <- dat$Samplingrate/100.0 #uncomment for large graph dataset

dat$samplerate <- dat$Samplingrate
dat$Samplingrate <- as.factor(dat$Samplingrate)
dat$totalnodes <- dat$numNodes

dat$numNodes <- as.factor(dat$numNodes)
#Plot of sampling rate vs sampled edges
g <- ggplot(dat,aes(x=Samplingrate,y=sampedges,color=numNodes)) + geom_boxplot() + scale_y_continuous("# of Edges Sampled") +theme(legend.position = "bottom")
ggsave(filename = "/Users/vmangipudi/largeedgesvssamplerate.png",g)  
  
g <- ggplot(dat,aes(x=Samplingrate,y=sampnodes,color=numNodes)) + geom_boxplot() + scale_y_continuous("# of Nodes remaning in sampled graph") +theme(legend.position = "bottom")
ggsave(filename = "/Users/vmangipudi/largenodesvssamplerate.png",g)  

#plot of density vs Sampling rate
dat$factoreddensity <- dat$SampledDensity/dat$samplerate

g<-ggplot(dat,(aes(x=Samplingrate,y=factoreddensity,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=ActualDensity,color=numNodes),shape=5,size=3 )  +theme(legend.position = "bottom")+scale_y_continuous("Sampled density/Sampling rate")
ggsave(filename = "/Users/vmangipudi/largefactoreddensityvssamplerate.png",g) 

#plot of sampling rate vs #of communities

g<-ggplot(dat,(aes(x=Samplingrate,y=SampledCommunitysize,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=ActualCommunitysize,color=numNodes),shape=5,size=3 ) +theme(legend.position = "bottom") + scale_y_continuous("Number of communities in sampled Graph")
ggsave(filename = "/Users/vmangipudi/numcommunitiesvssamplerate.png",g) 
#plot of sampling rate vs modularity
g<-ggplot(dat,(aes(x=Samplingrate,y=sampledmodularity,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=actualmodularity,color=numNodes),shape=5,size=3 ) +theme(legend.position = "bottom") + scale_y_continuous("Modularity in sampled Graph")
ggsave(filename = "/Users/vmangipudi/modularityvssamplerate.png",g) 

#plot of sampling rate vs 1st largest community

g<-ggplot(dat,(aes(x=Samplingrate,y=Lcommsize_samp1,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=Lcommsize_act1,color=numNodes),shape=5,size=3 ) +theme(legend.position = "bottom") +scale_y_continuous("Size of Largest community")
ggsave(filename = "/Users/vmangipudi/1commsizevssamplerate.png",g) 
#plot of sampling rate vs 2nd largest community

g<-ggplot(dat,(aes(x=Samplingrate,y=Lcommsize_samp2,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=Lcommsize_act2,color=numNodes),shape=5,size=3 ) +theme(legend.position = "bottom") +scale_y_continuous("Size of 2nd Largest community")
ggsave(filename = "/Users/vmangipudi/2commsizevssamplerate.png",g) 
#plot of sampling rate vs 3nd largest community

g<-ggplot(dat,(aes(x=Samplingrate,y=Lcommsize_samp3,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=Lcommsize_act3,color=numNodes),shape=5,size=3 ) +theme(legend.position = "bottom")+scale_y_continuous("Size of 3rd Largest community")
ggsave(filename = "/Users/vmangipudi/3commsizevssamplerate.png",g) 

#plot of sampling rate vs clustering coefficient 
g<-ggplot(dat,(aes(x=Samplingrate,y=Clustcoeff_samp,fill=numNodes))) + geom_boxplot() + geom_point(aes(y=Clustcoeff_act,color=numNodes),shape=5,size=3 ) +theme(legend.position = "bottom") + scale_y_continuous("Clustering coefficient in sampled Graph")
ggsave(filename = "/Users/vmangipudi/ccvssamplerate.png",g) 
#plot of sampling rate vs clustering coeficient distribution (D-Statistic)
g<-ggplot(dat,(aes(x=Samplingrate,y=ks_cc,fill=numNodes))) + geom_boxplot() +theme(legend.position = "bottom") + scale_y_continuous("Clustering coefficient (D-statistic)")
ggsave(filename = "/Users/vmangipudi/ks_ccvssamplerate.png",g) 

#plot of sampling rate vs page rank distribution (D-statitistic)
g<-ggplot(dat,(aes(x=Samplingrate,y=ks_pg,fill=numNodes))) + geom_boxplot() +theme(legend.position = "bottom") + scale_y_continuous("Page rank (D-statistic)")
ggsave(filename = "/Users/vmangipudi/ks_pgvssamplerate.png",g) 
#plot of sampling rate vs degree distribution (D-statistic)
g<-ggplot(dat,(aes(x=Samplingrate,y=ks_deg,fill=numNodes))) + geom_boxplot() +theme(legend.position = "bottom") + scale_y_continuous("Degree distribution (D-statistic)")
ggsave(filename = "/Users/vmangipudi/ks_degvssamplerate.png",g) 


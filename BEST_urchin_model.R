
setwd("C:/Users/kcresswe/OneDrive - University of Tasmania/WORK/IMAS/Sea_urchins/Bioeconomic_model/working model/Paper 1/Best/FINALFINAL")
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)
library(patchwork)
library(devtools)
library(StanHeaders)
library(doBy)

getwd()
        
zmax <- 9
fcmax <- 51
endtime<- 140                           #140 or 163, 2100 or 2123
fitsum <- array(0,dim=c(zmax, 4,8))
biomean <- array(0,dim=c(zmax,fcmax,6,endtime))
exploitmean <- array(0,dim=c(zmax,fcmax,6,endtime))
nevermean <- array(0,dim=c(zmax,6,endtime))
fitsummean <- array(0,dim=c(zmax,4,8))
catchmean <- array(0,dim=c(zmax,fcmax,6,endtime))
fishedvalsmean <- array(0,dim=c(zmax,5000))
neverfishedvalsmean <- array(0,dim=c(zmax,5000))
impactfishedvalsmean <- array(0,dim=c(zmax,5000))
fishedvals <- array(0,dim=c(zmax,5000))
neverfishedvals <- array(0,dim=c(zmax,5000))
impactfishedvals <- array(0,dim=c(zmax,5000))

# parameters survey dat
allsurveydat <- read.csv("survey_data.csv", header=T)
surveydat <- allsurveydat[,c(3,4,5)]  
surveydat <- as.matrix(surveydat)
reefareadat <- allsurveydat[,2]*1000*1000


historicalcatch <- read.csv("historical_catch.csv", header=T)
str(historicalcatch)
 
ncatchyr <- 15 ################## 
zeromatrix <- rep(0,(ncatchyr*9)) # from 2009 to 2023
zeromatrix <- as.data.frame(zeromatrix)
zeromatrix$Year <- rep(2009:2023,each=9)
zeromatrix$Region <- rep(1:9, by=ncatchyr)
zeromatrix <- zeromatrix[,-1]

historicalcatchdat <- left_join(zeromatrix,historicalcatch)
historicalcatchdat[is.na(historicalcatchdat)] <- 0

histcatch <- array(0,dim=c(zmax,ncatchyr))
histeffort <- array(0,dim=c(zmax,ncatchyr))
for (z in 1:9){
  for (t in 1:ncatchyr){
    row <- (t-1)*9+z
      histcatch[z,t] <- historicalcatchdat[row,3]
      histeffort[z,t] <- historicalcatchdat[row,4]
  }}

#test diameter maximum in mm
TDmax <- 150
TDmin <- 25
TDkg <- c(1:50)
sizebin <- (TDmax-TDmin)/50
tdmid <- c(1:50)
# convert test diameter in mm to wet weight in kg from Keane unpublished data 2020
for (i in 1:50){
  tdmid[i] <- TDmin + (sizebin)*i - sizebin/2 
  TDkg[i] <- 0.0035*(tdmid[i]^2.5257)/1000 # Keane factory data
  if (TDkg[i]<0) TDkg[i]=0
}

#growth transition matrix read in
G <- read.csv("size_transition_matrix.csv", header=T)
G <- G[,-1]
G <- as.matrix(G)

#stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(beepr)
beep()


pooledthing <- stan_model(file = "best_stan_urchin_model.stan", 
                          model_name = "GrowthCurve")
beep()

for (z in c(9,2,4,1,3,5,6,7,8)) {
  datanew = list(
  "fcmax" = fcmax,
  "histcatch" = histcatch[z,],
  "weight" = TDkg[1:50],
  "reefarea" = reefareadat[z],
  "G" = G,
  "N" = 5,                                        ########### 
  "x" =
    c(1, 20,  # years 1961 and 1980
      42, # 2002 
      57, #2017
      61), # 2021 
  "Y" =
    c(0.00000001, # made up near zero in 1961 (kg urchin/m2)
      0.00001, # made up - very low in 1980 when first urchin discovered
      surveydat[z,1],
      surveydat[z,2],
      surveydat[z,3]),
  "newt" = as.integer(endtime),
  "newx" = as.integer(c(1:endtime))
)


fit4 <- sampling(pooledthing, iter=3000,warmup=500,
                 thin=2,cores=4,chains=4, 
                 data=datanew)
print("z = ")
print(z)
print(fit4, "m")
print(fit4, c("apram", "bpram", "cpram", "sigma"))

#############Fished
Yffish <- rstan::extract(fit4, "Y_meanfished")
Yffishpred <- rstan::extract(fit4, "Y_predfished")
for (fc in 1:fcmax){
  for (t in 1:endtime) {
    vec1 <- which(Yffish$Y_meanfished[,fc,t]>=0)
    vec2 <- Yffish$Y_meanfished[vec1,fc,t]
    biomean[z,fc,c(1),t] <- mean(vec2)
    if (fc==1&t==63) fishedvals[z,1:length(vec2)] <- vec2
    biomean[z,fc,c(3,5),t] <- quantile(vec2,c(0.05,0.95))
    vec2 <- Yffishpred$Y_predfished[vec1,fc,t]
    biomean[z,fc,c(2),t] <-  max(mean(vec2),0)
    biomean[z,fc,c(4),t] <- max(quantile(vec2,c(0.05)),0)
    biomean[z,fc,c(6),t] <- max(quantile(vec2,c(0.95)),0)
  }  
}

# posterior prediction using best fit 
# for mean urchin density if no fishing had ever occurred,
#############NEVERFished - extract the biomass
Ynever <- rstan::extract(fit4, "Y_neverfished")
  for (t in 1:endtime) {
    vec1 <- which(Ynever$Y_neverfished[,t]>=0)
    vec2 <- Ynever$Y_neverfished[vec1,t]
    nevermean[z,c(1),t] <- mean(vec2)
    nevermean[z,c(3,5),t] <- quantile(vec2,c(0.05,0.95))
    if (t==63) {
      neverfishedvals[z,1:length(vec2)] <- vec2
      impactfishedvals[z,1:length(vec2)] <- (neverfishedvals[z,1:length(vec2)]-fishedvals[z,1:length(vec2)])/neverfishedvals[z,1:length(vec2)]
    }
  }  

# posterior for the exploitable biomass
#############Fished - extract the exploitable biomass
Yexploit <- rstan::extract(fit4, "Exploitbio")
for (fc in 1:fcmax){
  for (t in 1:endtime) {
    vec1 <- which(Yexploit$Exploitbio[,fc,t]>=0)
    vec2 <- Yexploit$Exploitbio[vec1,fc,t]/reefareadat[z]
    exploitmean[z,fc,c(1),t] <- mean(vec2)
    exploitmean[z,fc,c(3,5),t] <- quantile(vec2,c(0.05,0.95))
  }  }

# posterior of the actual catch removed over time depending on fishing levels fc
#############Fished - extract the catch
Ycatch <- rstan::extract(fit4, "Timecatch")
for (fc in 1:fcmax){
  for (t in 1:endtime) {
    vec1 <- which(Ycatch$Timecatch[,fc,t]>=0)
    vec2 <- Ycatch$Timecatch[vec1,fc,t]
    catchmean[z,fc,c(1),t] <- mean(vec2)
    catchmean[z,fc,c(3,5),t] <- quantile(vec2,c(0.05,0.95))
  }  
  }

fitsummary <- as.array(summary(fit4, pars = c("apram", "bpram","cpram","sigma"), probs = c(0.05, 0.5,0.95))$summary)
fitsum[z,,] <- fitsummary

}


beep()
mtext(side=1, line=2, "Years since 1970", outer=T, font=2,cex=1.4)
mtext(side=2, line=1, "Urchin biomass density (kg/m2)",  outer=T, font=2, cex=1.4)

#write - run this if writing to file
# fname <- "best_subsite_final_impact_4ch"
# for (z in c(1:9)) {
#   name <- paste0("results/",fname,"/biomean_z",z,".csv")
#   write.csv(biomean[z,,,],name)
#   name <- paste0("results/",fname,"/exploitmean_z",z,".csv")
#   write.csv(exploitmean[z,,,],name)
#   name <- paste0("results/",fname,"/nevermean_z",z,".csv")
#   write.csv(nevermean[z,,],name)
#   name <- paste0("results/",fname,"/catchmean_z",z,".csv")
#   write.csv(catchmean[z,,,],name)
#   name <- paste0("results/",fname,"/fitsum_",z,".csv")
#   write.csv(fitsum[z,,],name)
#   
#   name <- paste0("results/",fname,"/neverfishedvals_",z,".csv")
#   write.csv(neverfishedvals[z,],name)
#   name <- paste0("results/",fname,"/fishedvals_",z,".csv")
#   write.csv(fishedvals[z,],name)
#   name <- paste0("results/",fname,"/impactfishedvals_",z,".csv")
#   write.csv(impactfishedvals[z,],name)
#   
# }

#read
fname2 <- "best_subsite_final_impact_4ch"
for (z in c(1:9)) {
  name <- paste0("results/",fname2,"/biomean_z",z,".csv")
  bioread <- as.data.frame(read.csv(name))
  bioread <- bioread[-1]
  name <- paste0("results/",fname2,"/exploitmean_z",z,".csv")
  expread <- as.data.frame(read.csv(name))
  expread <- expread[-1]
  name <- paste0("results/",fname2,"/catchmean_z",z,".csv")
  catchread <- as.data.frame(read.csv(name))
  catchread <- catchread[-1]
  name <- paste0("results/",fname2,"/nevermean_z",z,".csv")
  neverread <- as.data.frame(read.csv(name))
  neverread <- neverread[-1]
  
  name <- paste0("results/",fname2,"/fishedvals_",z,".csv")
  fishedvalsR <- as.data.frame(read.csv(name))
  fishedvalsR <- fishedvalsR[-1]
  for (n in 1:5000){
    fishedvalsmean[z,n] <- fishedvalsR[n,1]
  }                                 
  name <- paste0("results/",fname2,"/neverfishedvals_",z,".csv")
  neverfishedvalsR <- as.data.frame(read.csv(name))
  neverfishedvalsR <- neverfishedvalsR[-1]
  for (n in 1:5000){
    neverfishedvalsmean[z,n] <- neverfishedvalsR[n,1]
  }                
  name <- paste0("results/",fname2,"/impactfishedvals_",z,".csv")
  impactfishedvalsR <- as.data.frame(read.csv(name))
  impactfishedvalsR <- impactfishedvalsR[-1]
  for (n in 1:5000){
    impactfishedvalsmean[z,n] <- impactfishedvalsR[n,1]
  }                
  
  name <- paste0("results/",fname2,"/fitsum_",z,".csv")
  fitsumread <- as.data.frame(read.csv(name))
  fitsumread <- fitsumread[-1]
  
  for (f in 1:fcmax){
    for (p in 1:6){
      for (n in 1:endtime){
        col <- (n-1)*6+p
        biomean[z,f,p,n] <- bioread[f,col]
        exploitmean[z,f,p,n] <- expread[f,col]#/reefareadat[z]
        catchmean[z,f,p,n] <- catchread[f,col]
        if (f==1) {
          nevermean[z,p,n] <- neverread[p,n]
          if (p<=4&n<=8) fitsummean[z,p,n] <- fitsumread[p,n]
        }
      }}}
  }# z


{
nabiomean <- biomean
nabiomean[is.na(nabiomean)] <- -100
nanevermean <- nevermean
nanevermean[is.na(nanevermean)] <- -100
naexpmean <- exploitmean
naexpmean[is.na(naexpmean)] <- -100
nacatchmean <- catchmean
nacatchmean[is.na(nacatchmean)] <- -100
}


####################### PLOT of chosen zones, with different catch levels ##################
####################### showing urchin DENSITY, CPUE & CATCH ##################
{

q[] <- .0006 # best estimate based on catch data and model data
timestart<-49   #2013 = 53    49 is 2009 model starts in 1 1961
timeend <-68   #2034 = 74   
fcmaxmax <- rep(0,9)
timetotal <- timeend-timestart+1
plotnum <-0
df <-as.data.frame(rep(c(timestart:timeend),fcmax))
names(df)[1] <- "x"
breaks <- seq(4+timestart,timeend,by=5)
labels <- (2010+seq(3,(timeend-timestart),by=5))
df$Y<-0
df$Ycatch<-0
df$fc <- 0
df$ExpB<- 0
df$Effort <- 0
df$cpue <- 0
df$truecpue <- 0
bottomzone <- 2
selectfc<-.5
bestfit <- 1
y1 <- array(0,dim=c(9,15))
y2 <- array(0,dim=c(9,15))
for(z in c(1:9)){
  df$Y<-0
  df$Ycatch<-0
  df$fc <- 0
  df$ExpB<- 0
  df$Effort <- 0
  df$cpue <- 0
  df$truecpue <- 0
  for (fc in 1:(fcmax)){ # for different levels of fishing from 1 to fcmax
    if (nabiomean[z,fc,bestfit,timeend]>=0) {
      #fc <- (f*2)-1
      fcmaxmax[z] <- fc
      df$Y[c(1:(timetotal))+(fc-1)*(timetotal)] <- nabiomean[z,fc,bestfit,timestart:timeend]
      df$ExpB[c(1:(timetotal))+(fc-1)*(timetotal)] <- naexpmean[z,fc,bestfit,timestart:timeend]
      df$Ycatch[c(1:(ncatchyr))+(fc-1)*(timetotal)] <- histcatch[z,]/1000
      df$Ycatch[c((ncatchyr+1):(timetotal))+(fc-1)*(timetotal)] <- nacatchmean[z,fc,bestfit,64:timeend]/1000
      df$Effort[c(1:(ncatchyr))+(fc-1)*(timetotal)] <- histeffort[z,]
      df$Effort[c((ncatchyr+1):(timetotal))+(fc-1)*(timetotal)] <- nacatchmean[z,fc,bestfit,64:timeend]/(naexpmean[z,fc,bestfit,64:timeend]*reefareadat[z]*q[z])
      for (t in 1:ncatchyr) if (histeffort[z,t]>0) df$truecpue[c(1:(ncatchyr))+(fc-1)*(timetotal)] <- histcatch[z,]/histeffort[z,]
      #df$cpue[c(15:(timetotal))+(fc-1)*(timetotal)] <- naexpmean[z,fc,1,63:timeend]*reefareadat[z]*q[z]
      df$cpue[c(1:(timetotal))+(fc-1)*(timetotal)] <- naexpmean[z,fc,bestfit,timestart:timeend]*reefareadat[z]*q[z]
      df$fc[c(1:(timetotal))+(fc-1)*(timetotal)] <- rep((fc-1)*0.02,(timetotal))
    } else break }
  
  df$effort2 <- df$Ycatch*1000/(df$ExpB*reefareadat[z]*q[z])
  df$fcfact <- as.factor(df$fc)
  df[is.na(df)] <- 0
  y1[z,] <- df[(df$fc==0)&(df$x<=63),2]
  
  #density
  pa <- ggplot() + labs(x=NULL, y=expression("Urchin density (kg/m"^2*")"),fill="F")+ 
    geom_hline(yintercept=c(0.1,0.2),colour="grey94")+
    geom_line(data=df[(df$Y>0)&(df$fc<=1),],aes(x=x,y=Y,color=fc,group=fcfact),linewidth=.8)+
    geom_line(show.legend = TRUE) + ylim(0,.24)+ geom_hline(yintercept = .07,col="green3",linetype="longdash",linewidth=1)+
    #geom_point(aes(x=c(49:63),y=y1),colour="red")+
     #geom_point(data=surveypoints,aes(x=Year,y=dat)) +
    theme(axis.text.x = element_blank()) + theme(axis.text.y = element_text(size=18))
  if (z==bottomzone)  pa <- pa + theme(axis.text.x = element_text(size=18,angle=45,hjust=1))+ labs(x="Year")
  assign(paste0("p",z,"d"),pa)
  
  #cpue
  pb <- ggplot()+labs(x=NULL, y=expression("Catch rate (kg/hr)"))+
    geom_line(data=df[(df$Y>0)&(df$fc<=1),],aes(x=x,y=cpue,color=fc,group=fcfact),show.legend = FALSE) + ylim(0,850) + 
    theme(axis.text.x = element_blank())+ theme(axis.text.y = element_text(size=18))
  if (z==bottomzone)  pb <- pb + theme(axis.text.x = element_text(size=18,angle=45,hjust=1))
  assign(paste0("p",z,"cp"),pb)
  
  #catch
  pc <- ggplot() + geom_line(data=df[(df$Y>0)&(df$fc<=1),],aes(x=x,y=(Ycatch),color=fc,group=fcfact),linewidth=.5)+ ylim(0,420)+labs(x=NULL, y="Annual catch (tonnes)")+
    geom_line(show.legend = FALSE) + #geom_point(aes(x=c(49:63),y=y2),colour="red")+
    theme(axis.text.x = element_blank())+ theme(axis.text.y = element_text(size=18))
  if (z==bottomzone)  pc <- pc + theme(axis.text.x = element_text(size=18,angle=45,hjust=1))
  assign(paste0("p",z,"ca"),pc)
  assign(paste0("p",z,"catch"),y2)
  
  #effort
  y2[z,] <- df[(df$fc==1)&(df$x<=63),6]
  pe <- ggplot() + geom_hline(yintercept=c(500,1000,1500,2000),colour="grey94")+
    geom_line(data=df[(df$Y>0)&(df$fc<=1),],aes(x=x,y=Effort,color=fc,group=fcfact),linewidth=1)+labs(x=NULL, y=expression("Annual effort (hr)"))+
    geom_line(show.legend = FALSE) + ylim(0,2000) + #geom_hline(yintercept = histeffort[z,ncatchyr],col="black",linetype="dotted",linewidth=1)+
    theme(axis.text.x = element_blank()) + theme(axis.text.y = element_text(size=18))
  if (z==bottomzone)  pe <- pe + theme(axis.text.x = element_text(size=18,angle=45,hjust=1)) + labs(x="Year")
  assign(paste0("p",z,"e"),pe)
}


p2d <- p2d + geom_line(aes(x=c(49:63),y=y1[2,]),col="grey",linewidth = 1) + geom_point(aes(x=c(49:63),y=y1[2,])) # + geom_hline(yintercept=c(0.1,0.2),colour='grey')
p2e <- p2e  + geom_line(aes(x=c(49:63),y=y2[2,]),col="red",linewidth = 1)+ geom_point(aes(x=c(49:63),y=y2[2,]))
p4d <- p4d + geom_line(aes(x=c(49:63),y=y1[4,]),col="grey",linewidth = 1) + geom_point(aes(x=c(49:63),y=y1[4,]))
p4e <- p4e + geom_line(aes(x=c(49:63),y=y2[4,]),col="red",linewidth = 1) + geom_point(aes(x=c(49:63),y=y2[4,]))
p9d <- p9d + geom_line(aes(x=c(49:63),y=y1[9,]),col="grey",linewidth = 1) + geom_point(aes(x=c(49:63),y=y1[9,]))
p9e <- p9e + geom_line(aes(x=c(49:63),y=y2[9,]),col="red",linewidth = 1) + geom_point(aes(x=c(49:63),y=y2[9,]))

combined <- (p4d)  /  (p9d  ) /(p2d ) &
  theme(legend.position = "right",text = element_text(size=18),
        panel.background = element_rect(fill = "white", colour = "grey50"))&
  theme( axis.title.y = element_text(size=18, angle=90))&
  theme( axis.title.x = element_text(size=18))&
  scale_x_continuous(breaks=breaks, labels=labels)&
  labs(main="b")&
  scale_color_gradient(low = "dark blue", high = "orange",name="H")
  pdf("final_figs/Figure_3_new.pdf", width = 12, height = 12)  # Adjust width and height as needed
    
  combined + plot_layout(guides="collect")
  dev.off()

}


####################################################################################################
regionmean <- summaryBy(kg.m2~c(time,site_no), data=subsitedensdat, FUN = mean)

#### best fit graph

pdf("final_figs/best_fit.pdf", width = 10, height = 8)  # Adjust width and height as needed


par(mfrow=c(3,3), mar=c(1,1,0.1,0), oma=c(5,5,1,1),cex=1.1)
timestart <- 0
# all regions fit model to data
timeend<-73 #73 is 2033 # 140 is 2100
names <- c("(a)",'(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)')
names <- c("(1)",'(2)','(3)','(4)','(5)','(6)','(7)','(8)','(9)')
for(z in c(1:9)){
  z
  
  plot(nabiomean[z,1,1,timestart:timeend],type="l", col="blue",ylim=c(0,.45),xaxt="n",yaxt = "n",font.main = 1)
  lines(nabiomean[z,1,3,timestart:timeend],col="red",lty="dashed")
  lines(nabiomean[z,1,5,timestart:timeend],col="red",lty="dashed")
  
  abline(h=c(0.1,0.2,0.3,0.4),col="grey94")
  points(regionmean[c(z,z+9,z+18),3]~c(43-timestart,58-timestart,62-timestart),col="aquamarine4",pch=16)
  if (z %in% c(7,8,9)) axis(1, at=seq(1,(timeend-timestart+1),by=5),las=2, labels=(seq(1,(timeend-timestart+1),by=5)+(1960+timestart-1)))
  if (z %in% c(1,4,7)) axis(2, at=c(0,0.1,0.2,0.3,0.4),las=2, labels=c(0,0.1,0.2,0.3,0.4))
  if (z %in% c(2,5,8,3,6,9)) axis(2, at=c(0,0.1,0.2,0.3,0.4),las=2, labels=c("","","","",""))
  
  if (z==1) {legend("left", c('Best fit',"90% prediction interval", "Surveyed urchin density"), 
                    lty=c(1,2,0),pch=c(0,0,16),col=c("blue","red","aquamarine4"), cex=.8,bty='n',lwd=c(1,1,0))}
  rect(xleft = 49, xright = 63, ybottom = par("usr")[3], ytop = par("usr")[4], 
       border = NA, col = adjustcolor("light grey", alpha = 0.3))
  title(names[z],line=-1.5,adj=0.1,font.main=1)
} 
mtext(side=1, line=2.5, "Year",  font=1,cex=1.6,outer=T)
mtext(side=2, line=1.5, expression("Urchin density (kg/m"^2*")"),  font=1, cex=1.6,outer=T)
mean(subsitedat[2,,3])*3

dev.off()



####################################################################################################

###########   NO FISHING   #############################################################
#IMPACT of fishing then

yrs <- 5
impact<- as.data.frame(c(1:(yrs*18)))
names(impact)[1] <- "density"
for (t in 1:yrs){
  timeend <- 43+(t-1)*5
  rowstart <- (t-1)*18
  impact$density[(rowstart+1):(rowstart+9)]  <- nabiomean[,1,1,timeend]
  impact$density[(rowstart+10):(rowstart+18)] <- nanevermean[,1,timeend]-nabiomean[,1,1,timeend]
}

impact$density[impact$density<0] <- 0
impact$Region <- rep(c(1:9), by=(yrs*2))
impact$Year <- 2003
impact$fished <- "Density with historical fishing"
for (t in 1:yrs){
  timeend <- 2003+(t-1)*5
  rowstart = (t-1)*18
  impact$fished[(rowstart+1):(rowstart+9)] <- "Density if no fishing had occurred"
  impact$Year[(rowstart+1):(rowstart+18)] <- timeend 
}


dat02 <- filter(impact,Year==2003)
dat07 <- filter(impact,Year==2008)
dat12 <- filter(impact,Year==2013)
dat17 <- filter(impact,Year==2018)
dat22 <- filter(impact,Year==2023)

dev.off()
barwidth <- 0.15
gapwidth <- 0.02
ggplot() +
  geom_col(data=dat02[order(dat02$fished),],
           mapping = aes(x=Region ,y=density,fill=fished),
           # stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat07[order(dat07$fished),],
           mapping = aes(x=Region + barwidth+gapwidth,y=density,fill=fished),
           # stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat12[order(dat12$fished),],
           mapping = aes(x=Region + barwidth*2+gapwidth*2,y=density,fill=fished),
           #stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat17[order(dat17$fished),],
           mapping = aes(x=Region + barwidth*3+gapwidth*3,y=density,fill=fished),
           #stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat22[order(dat22$fished),],
           mapping = aes(x=Region + barwidth*4+gapwidth*4,y=density,fill=fished),
           #stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  labs(x="", y="Urchin density (kg/m2)") + ylim(0,0.31) + 
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=16, angle=90, face="bold"), 
        axis.title.y=element_text(size=20,face="bold",margin=margin(r=10))) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_fill_manual(values = c("lightseagreen", "black"))

write.csv(impact,"impact_fishing.csv")

(impact[83,1]+impact[74,1])/impact[74,1]
(1.705407e-01-1.258862e-01)/1.258862e-01

(impact[c(82:90),1]+impact[c(73:81),1])/impact[c(73:81),1]*100

# reduction in density due to fishing
impact[c(82:90),1]/(impact[c(82:90),1]+impact[c(73:81),1])*100


###########################################################################
impactupper <- as.data.frame(c(1:(yrs*18)))
names(impactupper)[1] <- "density"
for (t in 1:yrs){
  timeend <- 43+(t-1)*5
  rowstart <- (t-1)*18
  impactupper$density[(rowstart+1):(rowstart+9)]  <- nanevermean[,3,timeend]
  impactupper$density[(rowstart+10):(rowstart+18)] <- nabiomean[,1,3,timeend]
}

impactupper$density[impactupper$density<0] <- 0
impactupper$Region <- rep(c(1:9), by=(yrs*2))
impactupper$Year <- 2003
impactupper$fished <- "Density with historical fishing"
for (t in 1:yrs){
  timeend <- 2003+(t-1)*5
  rowstart = (t-1)*18
  impactupper$fished[(rowstart+1):(rowstart+9)] <- "Density if no fishing had occurred"
  impactupper$Year[(rowstart+1):(rowstart+18)] <- timeend 
}

# reduction in density due to fishing
impactupper
(impactupper[c(73:81),1]-impactupper[c(82:90),1])/impactupper[c(73:81),1]*100

############################################################################################

#3 is quantile 0.05, 5 is quantile 0.95

impactlower <- as.data.frame(c(1:(yrs*18)))
names(impactlower)[1] <- "density"
for (t in 1:yrs){
  timeend <- 43+(t-1)*5
  rowstart <- (t-1)*18
  impactlower$density[(rowstart+1):(rowstart+9)]  <- nanevermean[,5,timeend]
  impactlower$density[(rowstart+10):(rowstart+18)] <- nabiomean[,1,5,timeend]
}

impactlower$density[impactupper$density<0] <- 0
impactlower$Region <- rep(c(1:9), by=(yrs*2))
impactlower$Year <- 2003
impactlower$fished <- "Density with historical fishing"
for (t in 1:yrs){
  timeend <- 2003+(t-1)*5
  rowstart = (t-1)*18
  impactlower$fished[(rowstart+1):(rowstart+9)] <- "Density if no fishing had occurred"
  impactlower$Year[(rowstart+1):(rowstart+18)] <- timeend 
}

# reduction in density due to fishing
(impactlower[c(73:81),1]-impactlower[c(82:90),1])/impactlower[c(73:81),1]*100

####################################################################################################

#plot mean impact and error bars
# Example data frame
data <- data.frame(
  region = factor(c(1:9)),  # 9 regions
  mean_density = impact[c(82:90),1]/(impact[c(73:81),1]+impact[c(82:90),1])*100,  # Example mean densities
  lower_error = (impactlower[c(73:81),1]-impactlower[c(82:90),1])/impactlower[c(73:81),1]*100,  # Example lower error (sd or se)
  upper_error = (impactupper[c(73:81),1]-impactupper[c(82:90),1])/impactupper[c(73:81),1]*100   # Example upper error
)

# Calculate ymin and ymax for error bars
data$ymin <- data$lower_error
data$ymax <- data$upper_error

impactvals <- fishedvalsmean
impactvals[1,] <- neverfishedvals[1,]-fishedvals[1,]

ndat <- 5000
impactdatframe <- as.data.frame(c(1:(ndat*9)))
impactdatframe[,1] <- rep(c(1:9),each=ndat)
names(impactdatframe)[1] <- "region"
impactdatframe$modelpoints <- 0
for (z in 1:9){
  row <- (z-1)*ndat+1
  impactdatframe$modelpoints[c(row:(row+ndat-1))] <- impactfishedvalsmean[z,c(1:ndat)]*100
}
length(is.na(impactfishedvalsmean)==TRUE)
impactdatframe$region <- as.factor(impactdatframe$region)

pdf("final_figs/impact_reduction_mean.pdf", width = 10, height = 8)  # Adjust width and height as needed
# Your plot code here

library(dplyr)

# Remove the top 500 data points per region
filtered_data <- impactdatframe %>%
  group_by(region) %>%
  arrange(desc(modelpoints)) %>%  # Sort by descending order of modelpoints
  slice((51:n())) %>%  # Retain rows from the 501st to the last
  ungroup()


pdf("final_figs/impact_reduction_mean_2.pdf", width = 8, height = 6)  # Adjust width and height as needed

ggplot(filtered_data, aes(x = region, y = modelpoints)) + 
  geom_jitter(width = 0.2, alpha = 0.2, color = "coral2") +  # Points first
  geom_violin(trim = TRUE, fill = "cadetblue1", color = "black", alpha = 0.8) +  # Violin next
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.05),
        ymax = quantile(x, 0.95)
      )
    },
    geom = "errorbar", width = 0.2, color = "black", size = .3
  ) +  # Custom boxplot with 5th and 95th percentiles
  stat_summary(
    fun = median, geom = "point", size = .3, color = "black"
  ) +  # Add median as a point
  theme_classic(base_size = 14) +
  labs(x = "Region", y = "Reduction in density from fishing (%)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )+
  ylim(-1, 78) 

dev.off()


########## to inform the values in the harvest rate table ##################################################################


#what was the harvest rate in 2023
histcatch[,15]
Htarget <- as.data.frame(c(1:9))
fiveyeartargetH0.06 <- rep(0,9)
fiveyeartargetH0.1 <- rep(0,9)
maintaintarget <- reduce20target <- deplete30UF <- rep(0,9)
names(Htarget)[1] <- "Region"
Htarget$harvest23 <- histcatch[,15]/(naexpmean[,2,1,63]*reefareadat[])
#what harvest rate needed to reach ecological target by 5 years max
for (z in 1:9) fiveyeartargetH0.06[z] <- ((which.min(nabiomean[z,,1,68]>0.07)))
for (z in 1:9) fiveyeartargetH0.1[z] <- ((which.min(nabiomean[z,,1,68]>0.15)))
for (z in 1:9) maintaintarget[z] <- ((which.min(nabiomean[z,,1,68]>nabiomean[z,,1,63])))
for (z in 1:9) reduce20target[z] <- ((which.min(nabiomean[z,,1,68]>nabiomean[z,,1,63]*.8)))
for (z in 1:9) deplete30UF[z] <- ((which.min(nabiomean[z,,1,68]/nanevermean[z,1,68]>.3)))
(deplete30UF-1)*1/50
z<-7
nabiomean[z,deplete30UF[z],1,68]
nanevermean[z,1,68]
nabiomean[z,,1,68]/nanevermean[z,1,68]
Htarget$fiveyear0.07 <- (fiveyeartargetH0.06-1)*1/50
Htarget$fiveyear0.15 <- (fiveyeartargetH0.1-1)*1/50
Htarget$fiveyearmaintain <- (maintaintarget-1)*1/50


plot(nabiomean[1,50,1,],ylim=c(0,.4))
lines(nanevermean[1,1,])
lines(nabiomean[1,50,1,]/nanevermean[1,1,])

#effort
effort5year0.07 <- effort5yrmaintain <- effort5yrdeplete30 <- array(0,dim=c(zmax,5))
Htarget$effort2023 <- histeffort[,15]
for (z in 1:9) effort5year0.07[z,] <- nacatchmean[z,fiveyeartargetH0.06[z],1,64:68]/(naexpmean[z,fiveyeartargetH0.06[z],1,64:68]*reefareadat[z]*q[z])
Htarget$fiveyr0.07effort <- effort5year0.07[,5]
for (z in 1:9) effort5yrmaintain[z,] <- nacatchmean[z,maintaintarget[z],1,64:68]/(naexpmean[z,maintaintarget[z],1,64:68]*reefareadat[z]*q[z])
Htarget$fiveyr0.07maintain <- effort5yrmaintain[,5]
for (z in 1:9) effort5yrdeplete30[z,] <- nacatchmean[z,deplete30UF[z],1,64:68]/(naexpmean[z,deplete30UF[z],1,64:68]*reefareadat[z]*q[z])
Htarget$fiveyr0.07maintain <- effort5yrdeplete30[,5]
effort5year0.07[,1]/Htarget$effort2023
effort5yrdeplete30[z,1]/Htarget$effort2023
effort5yrmaintain[z,1]/Htarget$effort2023

Htarget$fiveyr0.07effort/Htarget$effort2023

#target depleted
Htargetdepleted <- as.data.frame(c(1:9))
fiveyear_20pc <- rep(0,9)
names(Htargetdepleted)[1] <- "Region"
Htargetdepleted$harvest23 <- histcatch[,15]/(naexpmean[,2,1,63]*reefareadat[])
for (z in 1:9) fiveyear_20pc[z] <- ((which.min(nabiomean[z,,1,68]/nanevermean[z,1,68]>0.5)))
z2deplete <- rep(0,140)
for (t in 60:140) z2deplete[t] <- which.min(nabiomean[2,,1,t]/nanevermean[2,1,t]>0.5)
(z2deplete-1)*1/50
plot(nabiomean[2,22,1,1:73])
Htargetdepleted$fiveyr50pc <- (fiveyear_20pc-1)*1/50
-log(1-Htargetdepleted$fiveyr50pc)
-log(1-seq(0.0,1,by=0.02))
-log(1-Htargetdepleted$harvest23)*100




###########   NO FISHING   #############################################################
#IMPACT of fishing then

yrs <- 5
impact<- as.data.frame(c(1:(yrs*18)))
names(impact)[1] <- "density"
for (t in 1:yrs){
  timeend <- 43+(t-1)*5
  rowstart <- (t-1)*18
  impact$density[(rowstart+1):(rowstart+9)]  <- nabiomean[,1,1,timeend]
  impact$density[(rowstart+10):(rowstart+18)] <- nanevermean[,1,timeend]-nabiomean[,1,1,timeend]
}

impact$density[impact$density<0] <- 0
impact$Region <- rep(c(1:9), by=(yrs*2))
impact$Year <- 2003
impact$fished <- "Density with historical fishing"
for (t in 1:yrs){
  timeend <- 2003+(t-1)*5
  rowstart = (t-1)*18
  impact$fished[(rowstart+1):(rowstart+9)] <- "Density if no fishing had occurred"
  impact$Year[(rowstart+1):(rowstart+18)] <- timeend 
}


dat02 <- filter(impact,Year==2003)
dat07 <- filter(impact,Year==2008)
dat12 <- filter(impact,Year==2013)
dat17 <- filter(impact,Year==2018)
dat22 <- filter(impact,Year==2023)

dev.off()
barwidth <- 0.15
gapwidth <- 0.02
ggplot() +
  geom_col(data=dat02[order(dat02$fished),],
           mapping = aes(x=Region ,y=density,fill=fished),
          # stat="identity",
          position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat07[order(dat07$fished),],
           mapping = aes(x=Region + barwidth+gapwidth,y=density,fill=fished),
          # stat="identity",
          position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat12[order(dat12$fished),],
         mapping = aes(x=Region + barwidth*2+gapwidth*2,y=density,fill=fished),
         #stat="identity",
         position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat17[order(dat17$fished),],
           mapping = aes(x=Region + barwidth*3+gapwidth*3,y=density,fill=fished),
           #stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  geom_col(data=dat22[order(dat22$fished),],
           mapping = aes(x=Region + barwidth*4+gapwidth*4,y=density,fill=fished),
           #stat="identity",
           position=position_stack(reverse = TRUE),width = barwidth)+
  labs(x="", y="Urchin density (kg/m2)") + ylim(0,0.31) + 
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=16, angle=90, face="bold"), 
        axis.title.y=element_text(size=20,face="bold",margin=margin(r=10))) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_fill_manual(values = c("lightseagreen", "black"))

write.csv(impact,"impact_fishing.csv")

(impact[83,1]+impact[74,1])/impact[74,1]
(1.705407e-01-1.258862e-01)/1.258862e-01

# actual values for impact here
(impact[c(82:90),1]+impact[c(73:81),1])/impact[c(73:81),1]*100



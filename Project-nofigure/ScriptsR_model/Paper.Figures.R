##Introduction: checking the value for 
library(stringr)
#install.packages("deming")
library(deming)
library(plot3D)

par(mfrow=c(1,1),pin=c(8,7),mai=c(0.5,0.75,0.5,0.25))
#Analysis of Bar-Even11 data
##I Importing the global data
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Data")
enzyme_data<-read.csv(file="BarEven11.txt",sep="\t",header=TRUE)#,header=TRUE)
enzyme_data$kcat<-as.numeric(enzyme_data$kcat)
head(enzyme_data)

#1.Testing for a trade-off (non log) between kcat and KM
toff01<-cor.test(enzyme_data$kcat,enzyme_data$KM)
toff01#Conclusion: No linear relationship R^2~0, although significant light increase
toff02<-lm(enzyme_data$KM~enzyme_data$kcat)
summary(toff02)
plot(enzyme_data$KM,enzyme_data$kcat,col="red",pch=19)

enzyme_data$log10kcat<-log10(enzyme_data$kcat)
enzyme_data$log10KM<-log10(enzyme_data$KM)
enzyme_data$log10kcat_KM<-log10(enzyme_data$kcat/enzyme_data$KM*10^6)#*10^6 to visualize in MicroMolar

##2.kcat vs KM trade-off on log-scale
###y-linear regression
toff1<-lm(enzyme_data$log10KM~enzyme_data$log10kcat)
toff2<-lm(enzyme_data$log10kcat_KM~enzyme_data$log10kcat)
summary(toff1)#Conclusion: Light apprent trade-off with linear relationship, R^2=0.09
summary(toff2)#Conclusion: no trade-off at all, R^2=0.4
par(mfrow=c(1,2))
plot(enzyme_data$log10kcat,enzyme_data$log10KM,col="red",pch=19)
plot(enzyme_data$log10kcat_KM,enzyme_data$log10kcat,col="red",pch=19)

##3.log analysis of the relationship between kcat and kcat/KM
#0.Opening the multiple data files #Need to be put in the same file #setwd(dir="Desktop/ScriptsR_modeles/Data")
module_data<-read.csv(file="ID_modules.csv",sep=";",header=TRUE)
reaction_data<-read.csv(file="ID_reaction.csv",sep=";",header=TRUE)
enzyme_data<-read.csv(file="KM_kcat.csv",sep=";",header=TRUE)
organism_data<-read.csv(file="Organism_ID.csv",sep=";",header=TRUE)
#Renaming the columns
colnames(enzyme_data)<-c("Reaction_ID","Organism_ID","KM(muM)","kcat(1/s)")
colnames(organism_data)<-c("Organism_ID","Organism_name")

#a.Attributing module type to reaction
module_Type<-c()
for (r in 1:length(reaction_data$Reaction_ID)){
  module_Type<-c(module_Type,
                 toString(module_data$Type[module_data$ID_MOD==reaction_data$Module_ID[r]]))
}
reaction_data<-cbind(reaction_data,module_Type)

#b.Attributing module type to enzyme
module_Type<-c()
for (r in 1:length(enzyme_data$Reaction_ID)){
  module_Type<-c(module_Type,
                 toString(reaction_data$module_Type[reaction_data$Reaction_ID==enzyme_data$Reaction_ID[r]]))
}

enzyme_data<-cbind(enzyme_data,module_Type)

#c.Removing data if module type unknown
enzyme_data<-enzyme_data[enzyme_data$module_Type!="",]

##c.Attributing organism name to enzyme 
organism_ID<-c()
organism_Name<-c()

for (s in 1:length(enzyme_data$Organism_ID)){
  organism_ID<-c(organism_ID,enzyme_data$Organism_ID[s])
  organism_Name<-c(organism_Name,
                   toString(organism_data[organism_data$Organism_ID==toString(enzyme_data$Organism_ID[s]),]$Organism_name[1]))
}

enzyme_data<-cbind(enzyme_data,organism_Name)

##c.For enzymes involved in many modules, picking for the putatively most demanding (1.CEM, 2.AAFAN, 3.INTER,4.SEC) 
for (s in 1:length(enzyme_data$Organism_ID)){
  if(str_detect(enzyme_data$module_Type[s],"CEM")){
    enzyme_data$module_Type[s]="CEM"
  }
  else if(str_detect(enzyme_data$module_Type[s],"AA,FA,N")){
    enzyme_data$module_Type[s]="AA,FA,N"
  }
  else if(str_detect(enzyme_data$module_Type[s],"INTER")){
    enzyme_data$module_Type[s]="INTER"
  }
  else if(str_detect(enzyme_data$module_Type[s],"SEC")){
    enzyme_data$module_Type[s]="SEC"
  }
}


#1.kcat vs KM trade-off
##a.Analysing data with regards to categories
###With regards to Metabolic modules

#log transforming the data
log10kcat<-c()
log10KM<-c()
log10kcat_KM<-c()
for (s in 1:length(enzyme_data$Organism_ID)){
  log10kcat<-c(log10kcat,log(as.numeric(as.character(enzyme_data$`kcat(1/s)`[s])),10))
  log10KM<-c(log10KM,log(as.numeric(as.character(enzyme_data$`KM(muM)`[s])),10))
  log10kcat_KM<-c(log10kcat_KM,6+log(as.numeric(as.character(enzyme_data$`kcat(1/s)`[s]))/as.numeric(as.character(enzyme_data$`KM(muM)`[s])),10))
}
enzyme_data<-cbind(enzyme_data,log10kcat,log10KM,log10kcat_KM)

##Main findings (of our own)
setwd(dir="~")
setwd(dir="enzyme-evolution/ScriptsR_model")
source("AP1.Plotting_multiple_images.R")
source("AP2.RaphsonNewton_Long.R") #To find equilibrium for 2 variables problems
require(Rmpfr)
library(RColorBrewer)
setwd(dir="Data")

options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
ncontour=5
#palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)


##0.Estimates of cell minimum sizes

#Considering only basal metabolism
#For 1M
SE<-1
G_ATP_yield<-30
CM<-3.9*10^8/((60*60)*G_ATP_yield)#(V^0.88)in Glucose
SE_molecpermicrometer<-SE*6.02*10^23/(10^15)#conversion to molecule/micrometer^3
P<-10^-6
Max_size<-P*SE_molecpermicrometer/CM*4*pi/(4/3*pi)^1.88#Max size corresponds to R^3.64
Max_R<-Max_size^(1/3.64)
Max_V<-4/3*pi*Max_R^3

##I Resuts for the model case FD
###All results are presented for concentrations of 1M and KT=0.1M, except when expressly mentioned

#1a.Analysis with mechanistic values (not completely mechanistic however) and Vm variations
N_reso=100 #Resolution to build graphs
#Creation of two sets of parameters
KT=10^-4.5##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-4.5) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=c(10*KT)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(10^3) #An intermediate value seems reasonable to show general results
E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list() #Equilibrium product flux
log10_P_Phi_eq<-list() #log10 of equilibrium product flux

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      print(paste(i))
      kr=kr_set[l]
      for (j in 1:N_reso){
        a<- mpfr((Vm[s]*kf[i]+kf[i]*E_tot_conc*kcat[j]*(1+Se/KT)),120)
        b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf[i])+kf[i]*kcat[j]*E_tot_conc*(KT+Se)),120)
        c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),120)
        delta<-mpfr(b^2-4*a*c,120)
        S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
        ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
        tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10)
  }
}

#1.Plotting the results in 3D
plot.new()
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

par(mfrow=c(1,1),pin=c(8,7),mai=c(0.5,0,0.5,0.25))
x_log10kf<-seq(0,10,length=N_reso)
y_log10kcat<-seq(-4,6,length=N_reso)
fluxfacet <- (P_Phi_eq[[1]][-1, -1] + P_Phi_eq[[1]][-1, -N_reso] + P_Phi_eq[[1]][-N_reso, -1] + P_Phi_eq[[1]][-N_reso, -N_reso])/4
P_cut<-cut(fluxfacet,ncol)

##B&W figure
#persp(x_log10kf,y_log10kcat,xlab="log10 (kf)",ylab="log10 (kcat)",zlab="Equilibrium flux (M/s)",P_Phi_eq[[1]],theta = -30, phi = 30, expand = 0.6,
#      d=0.5,shade=0.25,border="black",axes = T, ticktype="detailed",nticks = 5,
#      box=T,cex.lab=0.75,cex.axis=0.4)

##Coloured figure
persp3D(x_log10kf,y_log10kcat,xlab="log10 kf",ylab="log10 kcat",zlab="Equilibrium flux (M/s)",zlim=c(0,10^-4.5),P_Phi_eq[[1]],theta = -20, phi = 20, expand = 1
      ,clab = c("Flux(M/s)"),shade=0.1,border="black",axes = T, ticktype="simple",ticks = FALSE,col=palet,#palet[P_cut]
      box=TRUE,cex.lab=1,cex.axis=0.8,r=5,lwd=0.35,lty=1,legend.shrink=0.4,colkey = list(plot = FALSE),contour = list(side = c("z"),levels=c(0.99*max(P_Phi_eq[[1]]),
                                                                                                                               0.999*max(P_Phi_eq[[1]]),
                                                                                                                               0.9999*max(P_Phi_eq[[1]]),
                                                                                                                               0.99999*max(P_Phi_eq[[1]])),colvar=P_Phi_eq[[1]],col="black",lwd=1,labcex=1))
ncut_leg=5
at_vec<-c()
lab_at<-c()
for (i in 1:ncut_leg){
  at_vec[i]=i/ncut_leg*max(P_Phi_eq[[1]])
  lab_at[i]<-c(round(i/ncut_leg,2))
}
persp3D(x_log10kf,y_log10kcat,xlab="log10 (kf)",ylab="log10 (kcat)",zlab="Equilibrium flux (M/s)",P_Phi_eq[[1]],theta = -20, zlim=c(0,1.2*10^-6),phi = 18, expand = 0.8
        ,shade=0.1,border="black",axes = T, ticktype="detailed",nticks = 5,col=palet,#palet[P_cut]
        box=FALSE,cex.lab=1,cex.axis=0.8,r=5,lwd=0.15,lty=1,legend.shrink=0.5,colkey = FALSE,add=TRUE,image=FALSE,contour = list(side = c("z"),levels=c(0.9*max(P_Phi_eq[[1]])),colvar=P_Phi_eq[[1]],col="white",lwd=1,lty=1))

colkey(clim=c(0,max(P_Phi_eq[[1]])),col=c("white"),side=4,length=0.8,width=0.7,cex.axis=0.75,side.clab=2,cex.clab=0.8,dist=-0.105,add=T,tck=2,hadj=2.5)
colkey(at=at_vec,labels=as.character(lab_at),clim=c(0,max(P_Phi_eq[[1]])),col=palet,side=4,length=0.8,cex.axis=0.75,side.clab=2,cex.clab=0.8,dist=-0.11,add=T)


#text(-0.153,-0.04,"0.2",srt=-60,font=1,col=1,cex=0.8,family="serif")
#text(-0.153,0.003,"0.4",srt=-40,font=1,col=1,cex=0.8,family="serif")
#text(-0.141,0.042,"0.6",srt=-40,font=1,col=1,cex=0.8,family="serif")
#text(-0.13,0.08,"0.8",srt=-40,font=1,col=1,cex=0.8,family="serif")
text(-0.125,0.105,"0.9*",srt=-45,font=1,col=1,cex=1.25,family="serif")
text(-0.091,0.121,"0.99",srt=-25,font=1,col=1,cex=0.7,family="serif")
text(-0.066,0.126,"0.999",srt=-25,font=1,col=1,cex=0.7,family="serif")
text(-0.042,0.128,"0.9999",srt=-25,font=1,col=1,cex=0.7,family="serif")
text(-0.019,0.130,"0.99999",srt=-25,font=1,col=1,cex=0.7,family="serif")

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_generic<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    w_kf_FD_generic[[s+length(Vm)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=20)
  }
}
##2D plot
source("AP11.Plotting_multiple_images.R")
multiplePlot("","","","","","",ncol,log10kf_var,log10kcat_var,w_kf_FD_generic,
             abs="log10 kf",ord="log10 kcat",palette=pal,image=TRUE,scale="AUTO",labcex=1,legcex=1,axcex=1,colorkey=TRUE,globcex=0.5) ##The unit and the value of Vm need to be checked and explained
contour(w_kf_FD_generic[[1]],lev=c(0.9),add=TRUE,lty=1,labcex=0.6,method='edge',col="black",labels=c("0.9"))
contour(w_kf_FD_generic[[1]],lev=c(0.05),add=TRUE,lty=4,labcex=1,method='edge',col="green",drawlabels=FALSE)
contour(w_kf_FD_generic[[1]],lev=c(0.9999),add=TRUE,lty=4,labcex=1,method='edge',col="red",drawlabels=FALSE)
setwd(dir="~")
setwd(dir="Desktop/")
dev.print(device = jpeg, file = "FullSensitivityAnalysis.jpeg", width = 560*3,height=475*3,res=100*3,type="cairo")


#2.Sensitivity study

##Specific values
N_reso=100
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5,5*10^-3)#
tab_P_Phi_eq_sens <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_spec<-list() #Equilibrium product flux
log10_P_Phi_eq_sens_spec<-list() #log10 of equilibrium product flux

for (l in 1:(length(KT_set))){
    kr=kr_set[1]
    Vm=Vm_set[l]
    KT=KT_set[l]
    Se=10*KT
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      for (j in 1:N_reso){
        a<- mpfr((Vm*kf[i]+kf[i]*E_tot_conc*kcat[j]*(1+Se/KT)),120)
        b<- mpfr((Vm*(kr+kcat[j]-Se*kf[i])+kf[i]*kcat[j]*E_tot_conc*(KT+Se)),120)
        c<- mpfr(-Vm*Se*(kr+kcat[j]),120)
        delta<-mpfr(b^2-4*a*c,120)
        S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
        ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
        tab_P_Phi_eq_sens[i,j]=as.numeric(kcat[j]*ES_conc_eq)
      }
    }
    P_Phi_eq_sens_spec[[l]]<-tab_P_Phi_eq_sens
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10)
}

w_kf_FD_sens_spec<-list()
for (l in 1:(length(KT_set))){
    w_kf_FD_sens_spec[[l]]<-round(P_Phi_eq_sens_spec[[l]]/max(P_Phi_eq_sens_spec[[l]]),digits=20)
}

##Sensitivity study
N_reso=100
Vm_set=c(10^-6,10^(-4.5),10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(10^-5,10^-1)#
tab_P_Phi_eq_sens <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens<-list() #Equilibrium product flux
log10_P_Phi_eq_sens<-list() #log10 of equilibrium product flux

for (l in 1:(length(KT_set))){
  for (s in 1:(length(Vm_set))){
    kr=kr_set[1]
    Vm=Vm_set[s]
    KT=KT_set[l]
    Se=10*KT
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      for (j in 1:N_reso){
        a<- mpfr((Vm*kf[i]+kf[i]*E_tot_conc*kcat[j]*(1+Se/KT)),120)
        b<- mpfr((Vm*(kr+kcat[j]-Se*kf[i])+kf[i]*kcat[j]*E_tot_conc*(KT+Se)),120)
        c<- mpfr(-Vm*Se*(kr+kcat[j]),120)
        delta<-mpfr(b^2-4*a*c,120)
        S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
        ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
        tab_P_Phi_eq_sens[i,j]=as.numeric(kcat[j]*ES_conc_eq)
      }
    }
    P_Phi_eq_sens[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq_sens
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10)
  }
}

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Vm_set))){
    w_kf_FD_sens[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq_sens[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq_sens[[s+length(Vm_set)*(l-1)]]),digits=20)
  }
}

#Plot for the SM
#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=20)
Affinity_set<-c("High affinity","Moderate affinity","Low affinity")
Flux_set<-c("Weak flux","Moderate flux","High flux")
addtxt<-list(l=0.95,h=0.95,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
#subtitle<-list("Weak flux, high affinity",
 #              "Moderate flux, high affinity",
  #             "High flux, high affinity",
   #            "Weak flux, moderate affinity",
   #            "Moderate flux, moderate affinity",
     #          "High flux, moderate affinity",
     #          "Weak flux, low affinity",
       #        "Moderate flux, low affinity",
       #        "High flux, low affinity")
multiplePlot("","","","",Affinity_set,Flux_set,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens,
             abs="log10 kf",ord="log10 kcat",lev=c(0.9),palette=pal,image=TRUE,scale="AUTO",
             labcex=1.5,TEXT_to_Add=addtxt,cextext=1.5,meth = "edge",legcex=0.575,axcex=0.7,colorkey="COMMON",globcex=0.5,contcex=1,contourlab=TRUE)
##Plot above finally put in appendix
##Plotting the different isoclines for FD on the same plot

#KT_set_mod<-c(10^-5,10^-1)

par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75,1.1, 0.75, 0.85))
plot.new()
contour(w_kf_FD_sens[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method = "edge",xlim=c(0,1),drawlabels =FALSE)
title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
for (l in 1:(length(KT_set))){
  for (s in 1:length(Vm_set)){
      col_set<-c(2,3,1)
      contour(w_kf_FD_sens[[(l-1)*length(Vm_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=col_set[s],lty=l+0,method = "edge",add=TRUE,xlim=c(0,1),drawlabels =FALSE)
      axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
           labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
      axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
           labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
      print((l-1)*length(KT_set)+s)
  }
}
sublegend<-c("Low",
             "Moderate",
             "High",
             "High",
             "Low")

L1=legend("bottomleft",legend=sublegend,lty=c(rep(1,3),1,2),col=c(2,3,1,1,1),ncol=2,cex=0.65,inset=.04,bty='n', x.intersp=0.5)
L2=legend(x = L1$rect$left-0.05, y = L1$rect$top+0.05, legend = c("Flux level","Affinity"), col=rep(NA,2), lty=c(1), ncol=2, x.intersp = -0.2, bg = NA,cex=0.65,,bty='n')
rect(L1$rect$left-0.02,L2$rect$top-L2$rect$h-L1$rect$h+0.02, L1$rect$left+L1$rect$w-0.04, L2$rect$top,
     col = "white",border="black", lty = 1, lwd = 1)
legend("bottomleft",legend=sublegend,lty=c(rep(1,3),1,2),col=c(2,3,1,1,1),ncol=2,cex=0.65,inset=.04,bty='n', x.intersp=0.5)
legend(x = L1$rect$left-0.05, y = L1$rect$top+0.05, legend = c("Flux level","  Affinity"), col=rep(NA,2), lty=c(1), ncol=2, x.intersp = -0.2, bg = NA,cex=0.65,,bty='n')
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFit_Flux&Sat.jpeg", width = 575*3,height=500*3,res=300,type="cairo")

par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75,1.1, 0.75, 0.85))
contour(w_kf_FD_sens_spec[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=0,col=c(rgb(255,0,0, maxColorValue=255, alpha=20)),lty=1,lwd=60,method = "edge",xlim=c(0,1),drawlabels =FALSE)
contour(w_kf_FD_sens_spec[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=rgb(10,10,10, maxColorValue=255, alpha=20),lty=1,lwd=60,method = "edge",xlim=c(0,1),drawlabels =FALSE,add=TRUE)
axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
     labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
     labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
text(0.41,0.52,"AA-N",srt=0,font=1,cex=1.3,family="serif",col="tomato")
text(0.3375,0.73,"CEM",srt=0,font=1,cex=1.3,family="serif",col="dark grey")
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)

title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFit_Flux&Trans.jpeg", width = 575*3,height=500*3,res=300,type="cairo")

##3.Analysis of saturation

options(digits=20)
N_reso=100
Vm_set=c(10^-4.5,10^-4.5) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(10^-5,10^-1)#
Se_set=c(1,0.01,100)
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
tab_P_Phi_eq_sens_sat <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_sat<-list() #Equilibrium product flux
log10_P_Phi_eq_sens_sat<-list() #log10 of equilibrium product flux

for (l in 1:(length(Vm_set))){
  for (s in 1:(length(Se_set))){
    kr=kr_set[1]
    Vm=Vm_set[l]
    KT=KT_set[l]
    Se=Se_set[s]*KT
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      for (j in 1:N_reso){
        a<- mpfr((Vm*kf[i]+kf[i]*E_tot_conc*kcat[j]*(1+Se/KT)),120)
        b<- mpfr((Vm*(kr+kcat[j]-Se*kf[i])+kf[i]*kcat[j]*E_tot_conc*(KT+Se)),120)
        c<- mpfr(-Vm*Se*(kr+kcat[j]),120)
        delta<-mpfr(b^2-4*a*c,120)
        S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
        ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
        tab_P_Phi_eq_sens_sat[i,j]=as.numeric(kcat[j]*ES_conc_eq)
      }
    }
    P_Phi_eq_sens_sat[[s+length(Se_set)*(l-1)]]<-tab_P_Phi_eq_sens_sat
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10)
  }
}
dev.off()
par(mfrow=c(2,3),mai=c(1,1,1,1))
for (i in 1:6){
  image.plot(P_Phi_eq_sens_sat[[i]])
}

multiplePlot("","","","",c(""),KT_set,ncol,log10kf_var,log10kcat_var,P_Phi_eq_sens_sat,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",lev=c(0.9),palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",
             globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)


##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_sat<-list()

for (l in 1:(length(Vm_set))){
  for (s in 1:(length(Se_set))){
    w_kf_FD_sens_sat[[s+length(Se_set)*(l-1)]]<-
      round(P_Phi_eq_sens_sat[[s+length(Se_set)*(l-1)]]/max(P_Phi_eq_sens_sat[[s+length(Se_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)


ncol=64
pal<-colorRampPalette(c("white"))(ncol)
palet<-list(pal)

par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75,1.1, 0.75, 0.85))
contour(w_kf_FD_sens_sat[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method = "edge",xlim=c(0,1),drawlabels =FALSE)
title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
for (l in 1:(length(Vm_set))){
  for (s in 1:length(Se_set)){
  contour(w_kf_FD_sens_sat[[(l-1)*length(Se_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s,lty=l+0,method = "edge",add=TRUE,xlim=c(0,1),drawlabels =FALSE)
  axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
     labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
  axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
     labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
  print((l-1)*length(Se_set)+s)
  }
}

list(expression(paste("[",S[out],"]/KT=0.01")),"[Sout]/KT=1","[Sout]/KT=100","  High","  Low")
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=1)

sublegend<-c(expression(paste(S[out],"/",K[T],"=",10^-2)),
             expression(paste(S[out],"/",K[T],"=",1)),
             expression(paste(S[out],"/",K[T],"=",10^2)),
             "High",
             "Low")

L1=legend("bottomleft",legend=sublegend,lty=c(rep(1,3),1,2),col=c(2,1,3,1,1),ncol=2,cex=0.65,inset=.04,bty='n', x.intersp=0.5)
L2=legend(x = L1$rect$left-0.05, y = L1$rect$top+0.05, legend = c("Saturation level","Affinity"), col=rep(NA,2), lty=c(1), ncol=2, x.intersp = -0.2, bg = NA,cex=0.65,,bty='n')
rect(L1$rect$left-0.02,L2$rect$top-L2$rect$h-L1$rect$h+0.02, L1$rect$left+L1$rect$w-0.08, L2$rect$top+0.02,
     col = "white",border="black", lty = 1, lwd = 1)
legend("bottomleft",legend=sublegend,lty=c(rep(1,3),1,2),col=c(2,1,3,1,1),ncol=2,cex=0.65,inset=.04,bty='n', x.intersp=0.5)
legend(x = L1$rect$left-0.05, y = L1$rect$top+0.05, legend = c("Saturation level","Affinity"), col=rep(NA,2), lty=c(1), ncol=2, x.intersp = -0.2, bg = NA,cex=0.65,,bty='n')
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFit_NonSatToSat.jpeg", , width = 575*3,height=500*3,res=300,type="cairo")

##4.Best value for the first enzyme followed by a variable second enzyme with sensitivy study
##With degradation in the pathway
N_reso=100
#Defining parameters
Vm_set=rep(10^-6,3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^10
kcat1<-10^6
kr1<-10^3
#Ki<-10^-6
eta_d_set<-c(10^-6,10^-4,10^-2)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(KT_set))){
  for (s in 1:(length(eta_d_set))){
    #kr=kr_set[1]
    Vm=Vm_set[s]
    eta_d<-eta_d_set[s]
    KT=KT_set[l]
    Se=10*KT
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

  
#b.Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(eta_d_set))){
    w_kf_FD_sens_2Enz[[s+length(eta_d_set)*(l-1)]]<-round(P_Phi_eq_sens_2Enz[[s+length(eta_d_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz[[s+length(eta_d_set)*(l-1)]]),digits=20)
  }
}

#image.plot(P_Phi_eq[[1]],legend.shrink=1)
options(digits=10)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B","C"),srt = 0,font=2,col="black",cex=1)

options(digits=3)
deg=round(eta_d_set*Vm_set[1],5)
round(deg,5)
options(digits=10)
subtitle<-list(substitute(paste("Low degradation rate - ",, eta," =",eta_d,"/s, irreversible first reaction"),list(eta_d=eta_d_set[1])),
               substitute(paste("Moderate degradation rate - ", eta," =",eta_d,"/s, irreversible first reaction"),list(eta_d=eta_d_set[2])),
               substitute(paste("High degradation rate - ", eta," =",eta_d,"/s, irreversible first reaction"),list(eta_d=eta_d_set[3])))

setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
FitLandscape_etahigh_noreverse<- multiplePlot("KT=","","Vm=","M/s",KT_set[1],signif(Vm_set[1],4),ncol,log10kf_var,log10kcat_var,list(w_kf_FD_sens_2Enz[[3]]),
                                              abs="log10 kf",ord="log10 kcat",lev=c(0.9),palette=pal,image=TRUE,scale="AUTO",labcex=1,globcex=0.5,axcex=1,,contcex=1,
                                              TEXT_to_Add=addtxt2,meth="flattest",colorkey="TRUE",contourlab=FALSE)+
  contour(w_kf_FD_generic[[1]],lev=c(0.9),add=TRUE,lty=4,labcex=1,method='edge',col="white",labels=c("1st enz")) 
dev.print(device = jpeg, file = "2DFitLandscape_etahigh_noreverse.jpeg", width = 575*3,height=460*3,res=300,type="cairo")

multiplePlot("","","","",c(""),c("","",""),ncol,log10kf_var,log10kcat_var, Conc_prod_eq_sens_2Enz,
             abs="log10 kf",ord="log10 kcat",palette=pal,image=TRUE,scale="AUTO",labcex=1.5,globcex=0.5,axcex=1,contcex=1,
             TEXT_to_Add=addtxt2,meth="flattest",colorkey="TRUE",contourlab=FALSE,cextext=1.5)
dev.print(device = jpeg, file = "2DConcentrationLandscape_degrates.jpeg", width = 1400*3,height=390*3,res=300,type="cairo")

###5.Effect of background enzyme on the landscape of the second enzyme (assuming no reversibility)
N_reso=50
#Defining parameters
Vm_set=c(10^-6)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^10)#10^2,
kcat1_set<-c(10^6)#10^-2,
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^-6,10^-4,10^-2)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-15


tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(eta_d_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[l]
    KT=KT_set[1]
    Se=10*KT
    kf1<-kf1_set[s]
    kcat1<-kcat1_set[s]
    kr1<-10^3#kcat1
    for (i in 1:(N_reso)){
      init<-c(10^-20,10^-20)#c(10^-4,10^-10)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
        tab_Conc_prod_eq[i,j]=SSR$y
        #init<-c(SSR$x,SSR$y)
      }
    }
    P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Maximum of fitness
max(P_Phi_eq_sens_2Enz_background_eff[[1]])/max(P_Phi_eq_sens_2Enz_background_eff[[2]])

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff<-list()
for (l in 1:(length(eta_d_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}
eff_enz<-c("Low","High")
multiplePlot("Degradation rate=","/s","First enzyme efficiency=","",signif(eta_d_set,4),eff_enz,ncol,log10kf_var,log10kcat_var,P_Phi_eq_sens_2Enz_background_eff,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey="COMMON",globcex=0.5,contcex=1,contourlab=TRUE,axcex=1,legcex=1)

multiplePlot("Degradation rate=","/s","First enzyme efficiency=","",signif(eta_d_set,4),eff_enz,ncol,log10kf_var,log10kcat_var, w_kf_FD_sens_2Enz_background_eff,
             abs="log10 (kf)",ord="log10 (kcat)",lev=c(0.9,0.99,0.999),palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey="COMMON",globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)

multiplePlot("Degradation rate=","/s","First enzyme efficiency=","",signif(eta_d_set,4),eff_enz,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey=TRUE,globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)
sublegend<-c(expression(paste(eta,"=",10^-6,"/s")),
             expression(paste(eta,"=",10^-4,"/s")),
             expression(paste(eta,"=",10^-2,"/s")))


plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1.125, 0.75, 1.125))

contour(w_kf_FD_sens_2Enz_background_eff[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
  #contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
  #contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
  title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
  for (l in 2:(length(eta_d_set))){
    for(s in 1:length(kcat1_set)){
      contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=s,lty=l+1,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
      #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
      #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
      axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
           labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
      axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
           labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
      print((l-1)*length(kcat1_set)+s)
    }
  }
    
  
addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","  Perfect"),bty="n",lty=c(rep(1,3)),col=c(1,2),ncol=1,cex=0.65)
L2=legend(0.71,0.84,title="Degradation rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","  Perfect"),bty="n",lty=c(rep(1,2)),col=c(1,2),ncol=1,cex=0.65)
L2=legend(0.71,0.84,title="  Degradation rate",legend=sublegend,bty="n",lty=c(1,3,4),col=c(1),ncol=1,cex=0.65)

L2=legend(0.72,0.975,title="  Degradation rate",legend=sublegend,bty="y",lty=c(1),col=c(2,1,3),ncol=1,cex=0.65)
x0<-c(0.1)
y0<-c(0.1)
x1<-c(0.6)
y1<-c(0.6)
arrows(x0,y0,x1,y1,code=2,lwd=2,lty=4)
mtext(text="Effect of increasing
      toxicity",at=c(0.6),line=-6.7,outer=FALSE,font=2,cex=0.75)

dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_Influence.jpeg", width = 575*3,height=475*3,res=100*3,type="cairo")
dev.print(device = jpeg, file = "2DFit_Landscape_Deg_Poster.jpeg", width = 560*3,height=475*3,res=100*3,type="cairo")

###5ter.Plotting itness with a sigmoid toxicity
##High flux
N_reso=100
#Defining parameters
Vm_set=c(10^-3)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-3)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^5,10^10)
kcat1_set<-c(10^-2,10^1,10^6)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^-15)#degradation rate #Sensitivity study to be done relatively to Vm
T_half_set<-c(10^-4,10^-2)#Toxicity penality

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-15


tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[1]
    T_half<-T_half_set[l]
    KT=KT_set[1]
    Se=10*KT
    kf1<-kf1_set[s]
    kcat1<-kcat1_set[s]
    kr1<-10^3#kcat1
    for (i in 1:(N_reso)){
      init<-c(10^-20,10^-20)#c(10^-4,10^-10)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)#+eta_d
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)*(1-SSR$y/(SSR$y+T_half))
        tab_Conc_prod_eq[i,j]=SSR$y
        #init<-c(SSR$x,SSR$y)
      }
    }
    P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Maximum of fitness
max(P_Phi_eq_sens_2Enz_background_eff[[1]])/max(P_Phi_eq_sens_2Enz_background_eff[[2]])

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff_tox<-list()
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff_tox[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C","D","E","F"),srt = 0,font=2,col="black",cex=1)
eff_enz<-c("Low","Moderate","High")
multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,P_Phi_eq_sens_2Enz_background_eff,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)


dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_Influence&Tox.jpeg", width = 1400*3,height=800*3,res=100*3,type="cairo")

#multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var, w_kf_FD_sens_2Enz_background_eff_tox,TEXT_to_Add = addtxt2,
   #         abs="log10 (kf)",ord="log10 (kcat)",lev=c(0.9,0.99,0.999),palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9)

#multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
             #abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey=TRUE,globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)
#sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),
             #expression(paste(eta,"=",10^-1,"/s")))


plot.new()
par(mfrow=c(1,2),oma=c(0,0,0,0),mai = c(0.75, 1.125, 0.75, 1.125))

contour(w_kf_FD_sens_2Enz_background_eff_tox[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
col_set<-c(3,1,2)
for (l in 2:(length(T_half_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff_tox[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=col_set[s],lty=l,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
    print((l-1)*length(kcat1_set)+s)
  }
}
sublegend=c(expression("T="~10^{-4}~"M"),expression("T="~10^{-2}~"M"))

addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
#L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","Moderate","Perfect"),bty="n",lty=c(rep(1,3)),col=c(1,2),ncol=1,cex=0.65)
#L2=legend(0.69,0.8,title="Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","Moderate","Perfect"),bty="n",lty=c(rep(1,3)),col=col_set,ncol=1,cex=0.65)
L2=legend(0.69,0.77,title="  Toxicity rate",legend=sublegend,bty="n",lty=c(1,2),col=c(1),ncol=1,cex=0.65)

##Low flux
N_reso=100
#Defining parameters
Vm_set=c(10^-6)#,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#,10^-3)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-c(10^2,10^5,10^10)
kcat1_set<-c(10^-2,10^1,10^6)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^-15)#degradation rate #Sensitivity study to be done relatively to Vm
T_half_set<-c(10^-4,10^-2)#Toxicity penality

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-15


tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    #kr=kr_set[1]
    Vm=Vm_set[1]
    eta_d<-eta_d_set[1]
    T_half<-T_half_set[l]
    KT=KT_set[1]
    Se=10*KT
    kf1<-kf1_set[s]
    kcat1<-kcat1_set[s]
    kr1<-10^3#kcat1
    for (i in 1:(N_reso)){
      init<-c(10^-20,10^-20)#c(10^-4,10^-10)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)#+eta_d
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)*(1-SSR$y/(SSR$y+T_half))
        tab_Conc_prod_eq[i,j]=SSR$y
        #init<-c(SSR$x,SSR$y)
      }
    }
    P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Maximum of fitness
max(P_Phi_eq_sens_2Enz_background_eff[[1]])/max(P_Phi_eq_sens_2Enz_background_eff[[2]])

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff_tox_LF<-list()
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff_tox_LF[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C","D","E","F"),srt = 0,font=2,col="black",cex=1)
eff_enz<-c("Low","Moderate","High")

#multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
#abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey=TRUE,globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)
#sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),
#expression(paste(eta,"=",10^-1,"/s")))

contour(w_kf_FD_sens_2Enz_background_eff_tox_LF[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox_LF[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox_LF[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
col_set<-c(3,1,2)
for (l in 2:(length(T_half_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff_tox_LF[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=col_set[s],lty=l,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
    print((l-1)*length(kcat1_set)+s)
  }
}
sublegend=c(expression("T="~10^{-4}~"M"),expression("T="~10^{-2}~"M"))

addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
#L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","Moderate","Perfect"),bty="n",lty=c(rep(1,3)),col=c(1,2),ncol=1,cex=0.65)
#L2=legend(0.69,0.8,title="Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","Moderate","Perfect"),bty="n",lty=c(rep(1,3)),col=col_set,ncol=1,cex=0.65)
L2=legend(0.69,0.77,title="  Toxicity rate",legend=sublegend,bty="n",lty=c(1,2),col=c(1),ncol=1,cex=0.65)

dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_Influence_Tox_Allflux.jpeg", width = 1200*3,height=480*3,res=100*3,type="cairo")

###5qua. Plotting with a given flux as the starting point
##High flux
N_reso=100
#Defining parameters
Phi_set=c(10^-6,10^-4.5,10^-3)
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
eta_d_set<-c(10^-2)#degradation rate #Sensitivity study to be done relatively to Vm
T_half_set<-c(10^-3)#Toxicity penality

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-15
tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_sens_2Enz_background_eff_FluxGen<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_background_eff<-list() #Equilibrium produt concentration

MaxconcP<-0

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(Phi_set))){
    #kr=kr_set[1]
    Phi=Phi_set[s]
    eta_d<-eta_d_set[1]
    T_half<-T_half_set[l]
    KT=KT_set[1]
    Se=10*KT
    kf1<-10^10
    kcat1<-10^6
    kr1<-10^3#kcat1
    for (i in 1:(N_reso)){
      init<-c(10^-20,10^-20)#c(10^-4,10^-10)
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kinh=10^3 if we are to introduce reverse reactions
        kf2<-kf2_set[i]
        #print(paste(i,"kf2",kf2))
        kcat2<-kcat2_set[j]
        kr2=kr2_set[1]
        #print(paste(j,"kcat2",kcat2))
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)#+eta_d
        F2e<-expression(Phi-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)*(1-SSR$y/(SSR$y+T_half))
        tab_Conc_prod_eq[i,j]=SSR$y
        #init<-c(SSR$x,SSR$y)
      }
    }
    P_Phi_eq_sens_2Enz_background_eff_FluxGen[[s+length(kcat1_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}


#Analysing the results
##Maximum of fitness
max(P_Phi_eq_sens_2Enz_background_eff[[1]])/max(P_Phi_eq_sens_2Enz_background_eff[[2]])

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff_tox_Genflux<-list()
for (l in 1:(length(T_half_set))){
  for (s in 1:(length(Phi_set))){
    w_kf_FD_sens_2Enz_background_eff_tox_Genflux[[s+length(kcat1_set)*(l-1)]]<-
      (round(P_Phi_eq_sens_2Enz_background_eff_FluxGen[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff_FluxGen[[s+length(kcat1_set)*(l-1)]]),digits=20))
  }
}
addtxt2<-list(l=0.98,h=0.98,txt=c("A","B","C","D","E","F"),srt = 0,font=2,col="black",cex=1)
eff_enz<-c("Low","Moderate","High")
multiplePlot("","","","flux",c(""),eff_enz,ncol,log10kf_var,log10kcat_var,w_kf_FD_sens_2Enz_background_eff_tox_Genflux,TEXT_to_Add = addtxt2,
             abs="log10 (kf)",ord="log10 (kcat)",lev=c(0.9),palette=pal,image=TRUE,scale="AUTO",labcex=1.5,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9,legcex=1,cextext=1.5)


dev.print(device = jpeg, file = "2DFit_Landscape_2Enz_First_Enz_GenFlux&Tox.jpeg", width = 1400*3,height=390*3,res=100*3,type="cairo")

#multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var, w_kf_FD_sens_2Enz_background_eff_tox,TEXT_to_Add = addtxt2,
#         abs="log10 (kf)",ord="log10 (kcat)",lev=c(0.9,0.99,0.999),palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey="COMMON",globcex=0.6,contcex=1,contourlab=TRUE,axcex=0.9)

#multiplePlot("Toxicity half constant=","M","First enzyme efficiency=","",T_half_set,eff_enz,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
#abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,subcex=2,colorkey=TRUE,globcex=0.5,contcex=1,contourlab=TRUE,axcex=1)
#sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),
#expression(paste(eta,"=",10^-1,"/s")))


plot.new()
par(mfrow=c(1,2),oma=c(0,0,0,0),mai = c(0.75, 1.125, 0.75, 1.125))

contour(w_kf_FD_sens_2Enz_background_eff_tox[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_kf_FD_sens_2Enz_background_eff_tox[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.25,font=2,xlab="log10 kf",ylab="log10 kcat",cex.lab=1)
col_set<-c(3,1,2)
for (l in 2:(length(T_half_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff_tox[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=col_set[s],lty=l,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1)
    print((l-1)*length(kcat1_set)+s)
  }
}
sublegend=c(expression("T="~10^{-4}~"M"),expression("T="~10^{-2}~"M"))

addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)
#L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","Moderate","Perfect"),bty="n",lty=c(rep(1,3)),col=c(1,2),ncol=1,cex=0.65)
#L2=legend(0.69,0.8,title="Toxicity rate",legend=sublegend,bty="n",lty=c(1,3),col=c(1),ncol=1,cex=0.65)
rect(L2$rect$left-0.02,L1$rect$top-L2$rect$h-L1$rect$h+0.02, L2$rect$left+L1$rect$w+0.03, L1$rect$top+0.02,col = "white",border="black", lty = 1, lwd = 1)
L1=legend(0.7,0.95,title="First enzyme",legend=c("Inefficient","Moderate","Perfect"),bty="n",lty=c(rep(1,3)),col=col_set,ncol=1,cex=0.65)
L2=legend(0.69,0.77,title="  Toxicity rate",legend=sublegend,bty="n",lty=c(1,2),col=c(1),ncol=1,cex=0.65)




###7.Mapping the data in the parameter space

setwd(dir="~")
setwd(dir="enzyme-evolution/ScriptsR_model")
source("AP11.Plotting_multiple_images.R")
#setwd(dir="Data")
require("Rmpfr")

#Defining features of plots
ncol=128
ncontour=5
palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue1","yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

##0.Importing the data
setwd(dir="~")
setwd(dir="/Users/florian/enzyme-evolution/Data")
enzyme_data<-read.csv(file="BarEven11.txt",sep="\t",header=TRUE)#,header=TRUE)
enzyme_data$log10kcat<-log10(enzyme_data$kcat)
enzyme_data$log10KM<-log10(enzyme_data$KM)
enzyme_data$log10kcat_KM<-log10(enzyme_data$kcat/enzyme_data$KM*10^6)

###8.Analysis with data and evolutionary results on phenomenological values kcat/KM and KM fitness landscapes
N_reso=250
KT_set=c(5*10^-5,5*10^-3) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^2
E_tot_conc=10^-3
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab
  }
}

addtxt11<-list(l=0.75,h=0.2,txt=expression(paste("Diffusion limit ",DL[2])),srt = 45,font=3,col="black",cex=1)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=100",s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt13<-list(l=1.05,h=0.75,txt=expression(paste(DL[3])),srt = 45,font=3,col="black",cex=1)
addtxt14<-list(l=1.05,h=0.85,txt=expression(paste(DL[4])),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
KMkcat_KmData<-as.data.frame(cbind(enzyme_data$log10kcat_KM,enzyme_data$log10kcat))
Data_list<-list(KMkcat_KmData)

multiplePlot("kr=","","Vm=","M/s",10^3,Vm[1],ncol,log10kcatKM_var,log10kcat_var,list(w_kcatKM[[1]]),
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),lev=c(0.9),
             xyData=Data_list,palette=pal,TEXT_to_Add =addtxt11,image=TRUE,col_data ="black",pcex=0.5,subcex=1,labcex=1,colorkey=TRUE,legcex=0.5,contourlab = FALSE,axcex=1)

abline(a=-0.4,b=1,col="lightgray",lty=1,lwd=5)
#abline(a=-0.3,b=1,col="grey",lty=6,lwd=2)
#abline(a=-0.2,b=1,col="grey",lty=6,lwd=2)

text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=1)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)


##Dissociating biosynthesis by end product

##Reorganizing the data to categorize enzymes
module_data<-read.csv(file="ID_modules.csv",sep=";",header=TRUE)
reaction_data<-read.csv(file="ID_reaction.csv",sep=";",header=TRUE)
enzyme_data<-read.csv(file="KM_kcat.csv",sep=";",header=TRUE)
organism_data<-read.csv(file="Organism_ID.csv",sep=";",header=TRUE)
#Renaming the columns
colnames(enzyme_data)<-c("Reaction_ID","Organism_ID","KM(muM)","kcat(1/s)")
colnames(organism_data)<-c("Organism_ID","Organism_name")

#a.Attributing module type to reaction
module_Type<-c()
for (r in 1:length(reaction_data$Reaction_ID)){
  module_Type<-c(module_Type,
                 toString(module_data$Type[module_data$ID_MOD==reaction_data$Module_ID[r]]))
}
reaction_data<-cbind(reaction_data,module_Type)

#b.Attributing module type to enzyme
module_Type<-c()
for (r in 1:length(enzyme_data$Reaction_ID)){
  module_Type<-c(module_Type,
                 toString(reaction_data$module_Type[reaction_data$Reaction_ID==enzyme_data$Reaction_ID[r]]))
}

enzyme_data<-cbind(enzyme_data,module_Type)

#c.Removing data if module type unknown
enzyme_data<-enzyme_data[enzyme_data$module_Type!="",]

##c.Attributing organism name to enzyme 
organism_ID<-c()
organism_Name<-c()

for (s in 1:length(enzyme_data$Organism_ID)){
  organism_ID<-c(organism_ID,enzyme_data$Organism_ID[s])
  organism_Name<-c(organism_Name,
                   toString(organism_data[organism_data$Organism_ID==toString(enzyme_data$Organism_ID[s]),]$Organism_name[1]))
}

enzyme_data<-cbind(enzyme_data,organism_Name)

##c.For enzymes involved in many modules, picking for the putatively most demanding (1.CEM, 2.AAFAN, 3.INTER,4.SEC) 
for (s in 1:length(enzyme_data$Organism_ID)){
  if(str_detect(enzyme_data$module_Type[s],"CEM")){
    enzyme_data$module_Type[s]="CEM"
  }
  else if(str_detect(enzyme_data$module_Type[s],"AA,FA,N")){
    enzyme_data$module_Type[s]="AA,FA,N"
  }
  else if(str_detect(enzyme_data$module_Type[s],"INTER")){
    enzyme_data$module_Type[s]="INTER"
  }
  else if(str_detect(enzyme_data$module_Type[s],"SEC")){
    enzyme_data$module_Type[s]="SEC"
  }
}

#log transforming the data
log10kcat<-c()
log10KM<-c()
log10kcat_KM<-c()
for (s in 1:length(enzyme_data$Organism_ID)){
  log10kcat<-c(log10kcat,log(as.numeric(as.character(enzyme_data$`kcat(1/s)`[s])),10))
  log10KM<-c(log10KM,log(as.numeric(as.character(enzyme_data$`KM(muM)`[s])),10))
  log10kcat_KM<-c(log10kcat_KM,6+log(as.numeric(as.character(enzyme_data$`kcat(1/s)`[s]))/as.numeric(as.character(enzyme_data$`KM(muM)`[s])),10))
}
enzyme_data<-cbind(enzyme_data,log10kcat,log10KM,log10kcat_KM)
Cat_Data_list<-list()
Cat_list<-c()
#for (i in 1:length(levels(factor(enzyme_data$module_Type)))){
for (i in 1:2){
  Cat_Data_list[[i]]<-as.data.frame(cbind(enzyme_data[enzyme_data$module_Type==levels(factor(enzyme_data$module_Type))[i],]$log10kcat_KM,
                           enzyme_data[enzyme_data$module_Type==levels(factor(enzyme_data$module_Type))[i],]$log10kcat))
  Cat_list[i]<-levels(factor(enzyme_data$module_Type))[i]
}

legtext1=c("AFN","CE")

#Defining features of plots
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=100",s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)

multiplePlot("","","","","","",ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),lev=c(0.9),
             xyData=Cat_Data_list,TEXT_to_Add =addtxt11,palette=pal,image=TRUE,
             ncat=2,pcex=0.9,subcex=1,labcex=2,legtext = legtext1,legpos=c(1.07,0.5),pch_dat=c(16,17),axcex=1,cextext=2,legtitle="Enzyme classes",colorkey=FALSE,legcex=1.75,contourlab=FALSE,globcex=0.5,contcex=1)

contour(w_kcatKM[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels=FALSE)
#text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
legend(x=1.05,y=0.8,title="Fitness isoclines",legend=c("amino acids","sugars"),lty=1,col=c(1,2),xpd=NA,cex=1.5)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
abline(a=-0.4,b=1,col="lightgray",lty=1,lwd=5)

dev.print(device = jpeg, file = "2DFitContour_DataClasses.jpeg", width = 590*3,height=460*3,res=100*3,type="cairo")

###9.Analysis of Evolutionary Results - no concentration evolution
##Low flux
kcat_i<-10^-3
kf_i<-10^2
kr_i<-10^3
E_tot_conc=10^-3

KT_set=c(10^-5,10^-1)
Vm_set=c(10^-6,10^-3)

#Flux determination
fit<-function(kf,kr,kcat,E_conc){
  a<- (Vm*kf+kf*E_conc*kcat*(1+Se/KT))
  b<- (Vm*(kr+kcat-Se*kf)+kf*kcat*E_conc*(KT+Se))
  c<- -Vm*Se*(kr+kcat)
  delta<-b^2-4*a*c
  S_conc_eq=(-b+delta^(1/2))/(2*a)
  ES_conc_eq=E_conc*(kf*S_conc_eq)/(kr+kcat+kf*S_conc_eq)
  Flux=kcat*ES_conc_eq
  return(Flux)
}

Max_val<-list()
KT=KT_set[1]
Vm=Vm_set[1]
Se=10*KT
Max_val[["W"]]<-fit(10^10,10^3,10^6,10^-3)

##a.Importing simulation results
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("W","H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(0,-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
      if(b==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi00",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/3-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/3-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][3+2*(long/3-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("W","H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:12){
    n=param_set$Ne_set[s]
    b=param_set$bias[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    bias_set_data[[f]]<-c(bias_set_data[[f]],rep(b,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","bias","fitness")
for (f in flux_set){
data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],bias_set_data[[f]],fit_evo_eq_summary[[f]]))
names(data_fit_eq[[f]]) <- col_headings
}

#Plotting on fitness landscapes
N_reso=250#250
KT_set=c(10^-5,10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^(2.8)
E_tot_conc=10^-3
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab
  }
}
#Defining features of plots
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

#Low flux
addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)

Evo_data<-list()

##Plot with nomutationl bias

f="W"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.7,subcex=1,labcex=1.8,axcex=0.9,cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(4,4)),colorkey=FALSE,globcex=0.5,legcex=1)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.375,0.959,"0.99",srt=-90,font=1,col=1,cex=1.65,family="serif")
#text(0.375,0.2,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.475,0.948,"0.999",srt=-90,font=1,col=2,cex=1.65,family="serif")
text(0.58,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.65,family="serif")
text(0.78,0.39,"0.99999",srt=-0,font=1,col=4,cex=1.65,family="serif")
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_lowF_nobias.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

##Plot with low mutationl bias
f="H"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.6,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=3,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[2]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
legend(x=1.05,y=0.9,title="Mutational bias",legend=c("b = 0","b = - 0.1","b = - 0.2"),pch=c(4,19,17),col=1,xpd=NA,cex=1.8)

#High mutational bias
##a.Importing simulation results
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("L")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
      if(b==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi00",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/3-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/3-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][3+2*(long/3-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("W")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:8){
    n=param_set$Ne_set[s]
    b=param_set$bias[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    bias_set_data[[f]]<-c(bias_set_data[[f]],rep(b,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","bias","fitness")
for (f in flux_set){
  data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],bias_set_data[[f]],fit_evo_eq_summary[[f]]))
  names(data_fit_eq[[f]]) <- col_headings
}
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

##Correlation test for simulation - between kcat and kcat_KM
Evo_data_full_high_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
cor.test(Evo_data_full_high_b$V1,Evo_data_full_high_b$V2)
Evo_data_full_lowf<-rbind(Evo_data_full_low_b,Evo_data_full_high_b)
cor.test(Evo_data_full_lowf$V1,Evo_data_full_lowf$V2)

xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
  for (i in 1:length(xyData)){
    points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
           (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
           col=i,pch=17,cex=0.8)
}

text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)

dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_highF_withbias.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

##Plot with low mutationl bias
mut_bias_set<-c(-0.1,-0.2)
f="H"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.7,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.32),
             cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=2.5,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
legend(x=1.08,y=0.9,title="Mutational bias",legend=c("b = 0","b = - 0.1","b = - 0.2"),pch=c(4,19,17),col=1,xpd=NA,cex=1.6)
legend(x=1.115,y=0.6,title="Protein cost",legend=c(expression(paste(10^{-3})),expression(paste(10^{-4}))),col=1,pch=c(6,2),xpd=TRUE,cex=1.6)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
axis(side=1,labels=FALSE,tck=FALSE)
axis(side=4,labels=FALSE,tck=FALSE)
#High mutational bias
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

##Correlation test for simulation - between kcat and kcat_KM
Evo_data_full_high_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
#cor.test(Evo_data_full_high_b$V1,Evo_data_full_high_b$V2)
Evo_data_full_lowf<-rbind(Evo_data_full_low_b,Evo_data_full_high_b)
#cor.test(Evo_data_full_lowf$V1,Evo_data_full_lowf$V2)
xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
for (i in 1:length(xyData)){
  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
         col=i,pch=17,cex=0.9)
}

#text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_highF_withbias.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

###9.Analysis of Evolutionary Results - no concentration evolution
##High flux
kcat_i<-10^-3
kf_i<-10^2
kr_i<-10^3
E_tot_conc=10^-3

KT_set=c(10^-5,10^-1)
Vm_set=c(10^-6,10^-3)

#Flux determination
fit<-function(kf,kr,kcat,E_conc){
  a<- (Vm*kf+kf*E_conc*kcat*(1+Se/KT))
  b<- (Vm*(kr+kcat-Se*kf)+kf*kcat*E_conc*(KT+Se))
  c<- -Vm*Se*(kr+kcat)
  delta<-b^2-4*a*c
  S_conc_eq=(-b+delta^(1/2))/(2*a)
  ES_conc_eq=E_conc*(kf*S_conc_eq)/(kr+kcat+kf*S_conc_eq)
  Flux=kcat*ES_conc_eq
  return(Flux)
}

Max_val<-list()
KT=KT_set[2]
Vm=Vm_set[2]
Se=10*KT
Max_val[["H"]]<-fit(10^10,10^3,10^6,10^-3)

##a.Importing simulation results
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("H")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(0,-0.1,-0.2)#,-0.3)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
      if(b==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi00",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/3-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/3-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][3+2*(long/3-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("H")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:8){
    n=param_set$Ne_set[s]
    b=param_set$bias[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    bias_set_data[[f]]<-c(bias_set_data[[f]],rep(b,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","bias","fitness")
for (f in flux_set){
  data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],bias_set_data[[f]],fit_evo_eq_summary[[f]]))
  names(data_fit_eq[[f]]) <- col_headings
}

#Plotting on fitness landscapes
N_reso=250#250
KT_set=c(10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^(2.8)
E_tot_conc=10^-3
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        print(paste(i,j))
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab
  }
}
#Defining features of plots
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

#High flux
addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)

Evo_data<-list()

##Plot with low mutational bias
f="W"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.6,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),
             cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=2.75,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.325,0.958,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.425,0.947,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.525,0.936,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.625,0.925,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
axis(side=1,labels=FALSE,tck=FALSE)
axis(side=4,labels=FALSE,tck=FALSE)
legend(x=1.08,y=0.9,title="Mutational bias",legend=c("b = - 0.1","b = - 0.2"),pch=c(19,17),col=1,xpd=NA,cex=1.75)


#High mutational bias
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

##Correlation test for simulation - between kcat and kcat_KM
Evo_data_full_high_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
cor.test(Evo_data_full_high_b$V1,Evo_data_full_high_b$V2)
Evo_data_full_lowf<-rbind(Evo_data_full_low_b,Evo_data_full_high_b)
cor.test(Evo_data_full_lowf$V1,Evo_data_full_lowf$V2)

xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
for (i in 1:length(xyData)){
  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
         col=i,pch=17,cex=0.8)
}

text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_highF_withbias.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")


##9.Manifold parameters involved: the case for reversibility
##And with degradation in the pathway

###With high affinity, testing the effect of reversibility
N_reso_ref=250
#Defining parameters
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
log10kcatKM_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso_ref)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso_ref)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^8
kcat1<-10^4
K_eq_ref<-1
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

kcatKM2_set<-kcatKM
kcat2_set<-kcat
kr2_set<-c(10^2)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso_ref, ncol=N_reso_ref)
tab_Conc_prod_eq<-matrix(nrow=N_reso_ref, ncol=N_reso_ref)
P_Phi_eq_sens_2Enz_rev<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_rev<-list() #Equilibrium produt concentration

MaxconcP<-0
for (l in 1:(length(KT_set))){
  for (s in 1:(length(K_eq_ref))){
    #kr=kr_set[1]
    Vm=Vm_set[l]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    kr1<-sqrt(K_eq_ref)*kcat1
    kinh1<-sqrt(K_eq_ref)*kf1
    for (i in 1:(N_reso_ref)){
      for (j in 1:N_reso_ref){
        print(paste(l,s,i,j))
        kr2=kr2_set[1]
        kcat2<-kcat2_set[j]
        kf2=kcatKM2_set[i]*(1+kr2/kcat2)
        if(kf2/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          Vm1_pos<-kcat1*E1_conc
          KM1_pos<-(kcat1+kr1)/kf1
          K_I<-kinh1/kf1
          Vm1_neg<-kr1*E1_conc
          KM1_neg<-(kcat1+kr1)/kinh1
          Vm2<-kcat2*E2_conc
          KM2<-(kcat2+kr2)/kf2
          T1<-KT+Se
          T2<-1+Se/KT
          #init<-c(10^-15,10^-15)
          F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
          F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
          Flist<-list(F1e,F2e)
          SSR<-RaphsonNewton(Flist,init,error)
          #print(paste(SSR$x,SSR$y))
          if(SSR$y>MaxconcP){
            MaxconcP<-SSR$y
          }
          tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
          tab_Conc_prod_eq[i,j]=SSR$y
        }
      }
    }
    P_Phi_eq_sens_2Enz_rev[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_rev[[s+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

N_reso=250
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Keq_set_sens<-c(10^-4,10^-2,10^2,10^4)
tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso, ncol=N_reso)
kcatKM2_set<-kcatKM
kcat2_set<-kcat

#for (l in 1:(length(kr2_set))){
#for (d in 1:length(eta_d_set)){
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set_sens))){
    #kr=kr_set[1]
    Vm=Vm_set[l]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    kr1<-sqrt(Keq_set_sens[s])*kcat1
    kinh1<-sqrt(Keq_set_sens[s])*kf1
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        kr2=kr2_set[1]
        kcat2<-kcat2_set[j]
        kf2=kcatKM2_set[i]*(1+kr2/kcat2)
        if(kf2/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          Vm1_pos<-kcat1*E1_conc
          KM1_pos<-(kcat1+kr1)/kf1
          K_I<-kinh1/kf1
          Vm1_neg<-kr1*E1_conc
          KM1_neg<-(kcat1+kr1)/kinh1
          Vm2<-kcat2*E2_conc
          KM2<-(kcat2+kr2)/kf2
          T1<-KT+Se
          T2<-1+Se/KT
          #init<-c(10^-15,10^-15)
          F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
          F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
          Flist<-list(F1e,F2e)
          SSR<-RaphsonNewton(Flist,init,error)
          #print(paste(SSR$x,SSR$y))
          if(SSR$y>MaxconcP){
            MaxconcP<-SSR$y
          }
          tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
          tab_Conc_prod_eq[i,j]=SSR$y
        }
      }
    }
    P_Phi_eq_sens_2Enz_rev[[(s+1)+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_rev[[(s+1)+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results

#w_kf_FD_sens_2Enz_rev[[1]]<-round(P_Phi_eq_sens_2Enz_rev[[1]]/max(P_Phi_eq_sens_2Enz_rev[[1]]),digits=20)
w_kf_FD_sens_2Enz_rev_tab<-matrix(nrow=N_reso , ncol=N_reso)
w_kf_FD_sens_2Enz_rev<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set_sens)+1)){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]][i,j]>0){
          w_kf_FD_sens_2Enz_rev_tab[i,j]<-round(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]][i,j]/max(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]]),digits=20)
        }
        else{
          w_kf_FD_sens_2Enz_rev_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kf_FD_sens_2Enz_rev[[s+length(Vm)*(l-1)]]<-w_kf_FD_sens_2Enz_rev_tab
  }
}


options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)
addtxt2<-list(l=0.95,h=0.95,txt=c("A","B","C"),srt = 0,font=2,col="black")
addtxt11<-list(l=0.75,h=0.2,txt=expression(paste("Diffusion limit ",DL[2])),srt = 45,font=3,col="black",cex=1)

multiplePlot("kr=","","Vm=","M/s",KT_set,K_eq_ref,ncol=128,log10kcatKM_var,log10kcat_var,list(w_kf_FD_sens_2Enz_rev[[1]]),
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),lev=c(0.9),TEXT_to_Add =addtxt11,palette=pal,image=TRUE,
             pcex=0.75,subcex=1,axcex=1,cextext=2,labcex=2,colorkey=FALSE,legtitle="",globcex=0.5,contourlab = FALSE,legcex=1)

for (s in 2:3){
  contour(w_kf_FD_sens_2Enz_rev[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=5-s,method = "simple",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
}
for (s in 4:5){
  contour(w_kf_FD_sens_2Enz_rev[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=s-2,method = "simple",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
}

contour(w_kf_FD_sens_2Enz_rev[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method = "simple",add=TRUE,xlim=c(0,1),drawlabels = FALSE)

text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
legend(xpd=TRUE,x=1.06,y=0.6,title=c("Equilibrium constant"),legend=c(expression("Keq="~10^{-4}),expression("Keq="~10^{-2}),expression("Keq="~1),expression("Keq="~10^{2}),expression("Keq="~10^{4})),lty=c(3,2,1,2,3),col=c(2,2,1,3,3),ncol=1,cex=1.4)


addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=100",s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt13<-list(l=1.07,h=0.75,txt=expression(paste(DL[3],"*")),srt = 45,font=3,col="black",cex=1)
addtxt14<-list(l=1.06,h=0.85,txt=expression(paste(DL[4])),srt = 45,font=3,col="black",cex=1)
#addtxt2<-list(l=0.02,h=0.98,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)


abline(a=-0.295,b=1.0,col="grey",lty=6,lwd=1)
abline(a=-0.195,b=1.0,col="grey",lty=6,lwd=1)
abline(a=-0.395,b=1.0,col="light grey",lty=1,lwd=4)


axis(1, at=seq(0,1,log10kcatKM_var[3]/abs(log10kcatKM_var[2]-log10kcatKM_var[1])),labels=FALSE,tck=-0.03,line=0)
axis(1, at=seq(0,1,log10kcatKM_var[3]/abs(log10kcatKM_var[2]-log10kcatKM_var[1])),labels=FALSE,tck=FALSE,line=0,side=4)

text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
text(addtxt13$l,addtxt13$h,addtxt13$txt[1],srt=addtxt13$srt,font=addtxt13$font,col=addtxt13$col,cex=2,xpd=NA)
text(addtxt14$l,addtxt14$h,addtxt14$txt[1],srt=addtxt14$srt,font=addtxt14$font,col=addtxt14$col,cex=2,xpd=NA)

dev.print(device = jpeg, file = "2DFitLandscape_Multiple_Reverse_exp_par.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

##10.Manifold parameters involved: the case for reversibility
##And with degradation in the pathway

###With high affinity, testing the effect of reversibility
N_reso=100
#Defining parameters
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1<-10^8
kcat1<-10^4
Keq_set_sens<-c(10^-4,10^-2,1,10^2,10^4)
#Ki<-10^-6
eta_d_set<-c(10^2)#degradation rate #Sensitivity study to be done relatively to Vm

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-c(10^3)##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-10
init<-c(10^-20,10^-20)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
tab_Conc_prod_eq<-matrix(nrow=N_reso, ncol=N_reso)
P_Phi_eq_sens_2Enz_rev<-list() #Equilibrium product flux
Conc_prod_eq_sens_2Enz_rev<-list() #Equilibrium produt concentration

MaxconcP<-0

for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set_sens))){
    #kr=kr_set[1]
    Vm=Vm_set[l]
    eta_d<-eta_d_set[1]*Vm
    KT=KT_set[l]
    Se=10*KT
    kr1<-sqrt(Keq_set_sens[s])*kcat1
    kinh1<-sqrt(Keq_set_sens[s])*kf1
    for (i in 1:(N_reso)){
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        kr2=kr2_set[1]
        kcat2<-kcat2_set[j]
        kf2=kf2_set[i]
        if(kf2/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          Vm1_pos<-kcat1*E1_conc
          KM1_pos<-(kcat1+kr1)/kf1
          K_I<-kinh1/kf1
          Vm1_neg<-kr1*E1_conc
          KM1_neg<-(kcat1+kr1)/kinh1
          Vm2<-kcat2*E2_conc
          KM2<-(kcat2+kr2)/kf2
          T1<-KT+Se
          T2<-1+Se/KT
          #init<-c(10^-15,10^-15)
          F1e<-expression(Vm1_pos*x/(KM1_pos+x+K_I*y)-Vm1_neg*y/(KM1_neg+y+x/K_I)-(Vm2/(KM2+y)+eta_d)*y)
          F2e<-expression(Vm*(Se-x)/(T1+x*T2)-(Vm2/(KM2+y)+eta_d)*y)
          Flist<-list(F1e,F2e)
          SSR<-RaphsonNewton(Flist,init,error)
          #print(paste(SSR$x,SSR$y))
          if(SSR$y>MaxconcP){
            MaxconcP<-SSR$y
          }
          tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
          tab_Conc_prod_eq[i,j]=SSR$y
        }
      }
    }
    P_Phi_eq_sens_2Enz_rev[[(s)+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    Conc_prod_eq_sens_2Enz_rev[[(s)+length(Vm_set)*(l-1)]]<-tab_Conc_prod_eq
  }
}

#Analysing the results
##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_rev_app<-list()
for (l in 1:(length(KT_set))){
  for (s in 1:(length(Keq_set_sens))){
    w_kf_FD_sens_2Enz_rev_app[[s+length(Keq_set_sens)*(l-1)]]<-round(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]]/max(P_Phi_eq_sens_2Enz_rev[[s+length(Keq_set_sens)*(l-1)]]),digits=20)
  }
}

#Defining features of plots
plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1.25, 0.75, 1))
addtxt2<-list(l=0.975,h=0.975,txt=c("A","B","C"),srt = 0,font=2,col="black")

ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

contour(w_kf_FD_sens_2Enz_rev_app[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="edge",xlim=c(0,1),drawlabels = FALSE)
title(main="",outer=FALSE,line=1.75,cex.main=2,font=2,xlab="log10 kf",ylab="",cex.lab=1)
title(main="",outer=FALSE,line=2.5,cex.main=2,font=2,xlab="",ylab="log10 kcat",cex.lab=1)
for (s in 1:2){
  contour(w_kf_FD_sens_2Enz_rev_app[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=4-s,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
}
for (s in 4:5){
  contour(w_kf_FD_sens_2Enz_rev_app[[s]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=s-2,method = "edge",add=TRUE,xlim=c(0,1),drawlabels = FALSE)
}

text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=1)

#legend(xpd=TRUE,x=1.1,y=0.6,title=c("Equilibrium constant"),legend=c(expression("Keq="~10^{-4}),expression("Keq="~10^{-2}),expression("Keq="~1),expression("Keq="~10^{2}),expression("Keq="~10^{4})),lty=c(3,2,1,2,3),col=c(2,2,1,3,3),ncol=1,cex=0.7)

#coordonnées des fleches 1
x0<-c(0.3)
y0<-c(0.2)
x1<-c(0.8)
y1<-c(0.7)
arrows(x0,y0,x1,y1,code=2,lwd=2,lty=4)
mtext(text="Effect of increasing
      reversibility",at=c(0.85),line=-4.7,outer=FALSE,font=2,cex=0.75)
axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
     labels=FALSE,cex.axis=1)
mtext(text=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),line=0.7,side=1,cex=1)
axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
     labels=FALSE,cex.axis=1)
mtext(text=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),line=0.6,side=2,cex=1)

#title(main=expression("Dissociation rate "~k[r]~"="~10^3~s^{-1})
   #   ,col.main="black",font.main=1,cex.main=1)

dev.print(device = jpeg, file = "2DFitLandscape_Multiple_Reverse.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")
#11.Analyzing with kcat[E]tot for the first enzyme

#Creation of the new set of parameter
N_reso=250
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
Etot_set=c(10^-9,10^-6,10^-3) #Quite efficient catalytic property for the enzyme
Vm_set=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
KT_set=c(5*10^-5)#
kr_set=10^3
#kcat_set<-kcat

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(Etot_set))){
  kr<-kr_set[1]
  E_tot_conc<-Etot_set[l]
  for (s in 1:(length(Vm_set))){
    Vm=Vm_set[s]
    KT=KT_set[s]
    Se=10*KT
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        print(paste(l,s,i,j))
        #kf=kcatKM[i]*(kr+kcat[j])/kcat[j]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          #k1=10^10
          #}
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<- (Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT))
          b<- (Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se))
          c<- -Vm[s]*Se*(kr+kcat[j])
          delta<-b^2-4*a*c
          S_conc_eq=(-b+delta^(1/2))/(2*a)
          ES_conc_eq=E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq)
          tab_P_Phi_eq[i,j]=kcat[j]*ES_conc_eq
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_Etot_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_Etot_kcatKM<-list()
for (l in 1:(length(Etot_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_Etot_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_Etot_kcatKM_tab[i,j]<-1.05 #For values above diffusion limit
        }
      }
    }
    w_Etot_kcatKM[[s+length(Vm)*(l-1)]]<-w_Etot_kcatKM_tab
  }
}

ncol=128
ncontour=5
#palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)
addtxt11<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)

Vm_length_empty<-c("")
kcat_length_empty<-c("")
subtitle<-list(expression(k[cat]~"="~10^{-2}~s^{-1}),
               expression(k[cat]~"="~10^{1}~s^{-1}),
               expression(k[cat]~"="~10^{4}~s^{-1}))
              
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)

multiplePlot("kr=","","Vm=","M/s",kcat_length_empty,Vm_length_empty,ncol=128,log10kcatKM_var,log10kcat_var,list(w_Etot_kcatKM[[2]]),
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),lev=c(0.9),meth="edge",TEXT_to_Add =addtxt11,palette=pal,image=TRUE,
             pcex=0.75,subcex=1,axcex=1,cextext=2,labcex=1.8,colorkey=FALSE,globcex=0.5,legcex=1,contourlab=FALSE,legtitle="")
contour(w_Etot_kcatKM[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
contour(w_Etot_kcatKM[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="edge",xlim=c(0,1),add=TRUE,drawlabels = FALSE)
abline(a=-0.29,b=1,col="light grey",lty=1,lwd=10)

axis(1, at=seq(0,1,log10kcatKM_var[3]/abs(log10kcatKM_var[2]-log10kcatKM_var[1])),labels=FALSE,tck=-0.03,line=0)
axis(1, at=seq(0,1,log10kcatKM_var[3]/abs(log10kcatKM_var[2]-log10kcatKM_var[1])),labels=FALSE,tck=FALSE,line=0,side=4)
legend(xpd=TRUE,x=1.05,y=0.6,title=c("Enzyme concentration"),legend=c(expression("[E]="~10^{-9}~"M"),expression("[E]="~10^{-6}~"M"),expression("[E]="~10^{-3}~"M")),lty=1,col=c(2,1,3),ncol=1,cex=1.4)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)

dev.print(device = jpeg, file = "Plot2DFitLandscape_Enz_conc.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")


plot.new()
contour(w_Etot_kcatKM[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="flattest",xlim=c(0,1))
contour(w_Etot_kcatKM[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="flattest",xlim=c(0,1),add=TRUE)
contour(w_Etot_kcatKM[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=3,lty=1,method ="flattest",xlim=c(0,1),add=TRUE)

title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 (kcat/KM)",ylab="log10 [E]tot",cex.lab=1)
axis(1, at=seq(0,1,log10kcatKM_var[3]/abs(log10kcatKM_var[2]-log10kcatKM_var[1])),
     labels=seq(log10kcatKM_var[1],log10kcatKM_var[2],log10kcatKM_var[3]),cex.axis=1.25)
axis(2, at=seq(0,1,log10E_tot_var[3]/abs(log10E_tot_var[2]-log10E_tot_var[1])),
     labels=seq(log10E_tot_var[1],log10E_tot_var[2],log10E_tot_var[3]),cex.axis=1.25)

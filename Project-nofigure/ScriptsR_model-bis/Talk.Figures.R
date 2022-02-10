setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/")
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

##I Resuts for the model case FD
###All results are presented for concentrations of 1M and KT=0.1M, except when expressly mentioned

#1a.Analysis with mechanistic values (not completely mechanistic however) and Vm variations
N_reso=50 #Resolution to build graphs
#Creation of two sets of parameters
KT=5*10^-5##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
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
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "yellow")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

par(mfrow=c(1,1),pin=c(8,7),mai=c(0.5,0.75,0.5,0.25))

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_generic<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    w_kf_FD_generic[[s+length(Vm)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=20)
  }
}
##2D plot

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10kf_var,log10kcat_var,w_kf_FD_generic,
             abs="log10 (kf)",ord="log10 (kcat)",ncont=0,palette=pal,image=TRUE,scale="AUTO",labcex=1.5) ##The unit and the value of Vm need to be checked and explained
contour(w_kf_FD_generic[[1]],lev=c(0.9),add=TRUE,lty=1,labcex=1,method='edge',col="black",labels=c("0.9"))
contour(w_kf_FD_generic[[1]],lev=c(0.99,0.999,0.9999,0.99999),add=TRUE,lty=4,labcex=1,method='edge',col="black")


setwd(dir="/Users/florianlabourel/Desktop/enzyme-evolution-master/ScriptsR_model")

###9.Analysis of Evolutionary Results

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
setwd(dir="Desktop/ScriptsR_modeles")
flux_set=c("W","H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(0,-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
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
N_reso=250
KT_set=c(10^-5,10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^3
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

plot.new()
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",palette=pal,TEXT_to_Add =addtxt1,image=TRUE,scale=c(0,1),labcex=1.5,subcex=1) ##The unit and the value of Vm need to be checked and explained

contour(w_kcatKM[[1]],lev=c(0.9,0.99,0.999,0.9999,0.99999),add=TRUE,lty=c(1,rep(4,4)),labcex=1,method='flattest',col="black")
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=1)

#Defining features of plots
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

#Low flux
addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)

Evo_data<-list()

##Plot with no mutationl bias
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
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=1,subcex=1.25,labcex=2.5,legtext = legtext1,legpos=c(1.09,0.5)
             ,axcex=1.25,cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)))


levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
text(0.41,1.051,"0.99",srt=-90,font=1,col=1,cex=2,family="serif",xpd=TRUE)
text(0.52,1.066,"0.999",srt=-90,font=1,col=2,cex=2,family="serif",xpd=TRUE)
text(0.62,1.081,"0.9999",srt=-90,font=1,col=3,cex=2,family="serif",xpd=TRUE)
text(0.72,1.096,"0.99999",srt=-90,font=1,col=4,cex=2,family="serif",xpd=TRUE)

##Plot with low mutationl bias
f="W"
b=mut_bias_set[2]
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
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.9,subcex=1,labcex=2.5,legtext = legtext1,legpos=c(1.09,0.5)
             ,axcex=1.25,cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)))


levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
text(0.41,1.051,"0.99",srt=-90,font=1,col=1,cex=2,family="serif",xpd=TRUE)
text(0.52,1.066,"0.999",srt=-90,font=1,col=2,cex=2,family="serif",xpd=TRUE)
text(0.62,1.081,"0.9999",srt=-90,font=1,col=3,cex=2,family="serif",xpd=TRUE)
text(0.72,1.096,"0.99999",srt=-90,font=1,col=4,cex=2,family="serif",xpd=TRUE)


#High mutational bias
b=mut_bias_set[3]
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
         col=i,pch=17,cex=1)
}

legend(x=1.05,y=0.8,title="Mutational bias",legend=c("b = - 0.1","b = - 0.2"),pch=c(19,17),col=1,xpd=NA,cex=1.5)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)

##Full data

kcat_KM<-c()
kcat<-c()
mut_bias_set_red<-c(-0.1,-0.2)

flux=c("W","H")
for (f in flux){
  for(b in mut_bias_set_red){
    for (n in Ne_set){
      kcat_KM<-c(kcat_KM,log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))))
      kcat<-c(kcat,log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))))
    }
  }
}

Evo_data_full<-as.data.frame(cbind(kcat_KM,kcat))

ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)

multiplePlot("kr=","","Vm=","M/s",kr,Vm[[1]],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),lev=c(0.9),
             xyData=list(Evo_data_full),TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=1,subcex=1,labcex=2,axcex=1.25,cextext=2,col_data=1,pch_dat = c(19))

contour(w_kcatKM[[2]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lty=1,method ="edge",xlim=c(0,1),add=TRUE)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
legend(x=1.05,y=0.5,title="Fitness isoclines",legend=c("amino acids","sugars"),lty=1,col=c(1,2),xpd=NA,cex=1.5)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[2],ncol=128,log10kcatKM_var,log10kcat_var,list(w_kcatKM[[2]]),
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=list(Evo_data_full),TEXT_to_Add =addtxt1,palette=pal,image=TRUE,pcex=1.05,subcex=1,labcex=2,axcex=1.25,cextext=2,col_data=1)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[2],ncol=128,log10kcatKM_var,log10kcat_var,list(w_kcatKM[[2]]),
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),ncat=1,
             xyData=list(Evo_data_full),TEXT_to_Add =addtxt1,palette=pal,image=TRUE,pcex=1.05,subcex=1,labcex=2,legtext = legtext1,legpos=c(1.09,0.5),axcex=1.25,cextext=2,legtitle="log10(Ne)",pch_dat = c(19))


jet.colors <- colorRampPalette(c("steelblue1","yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

##0.Importing the data
library(stringr)


setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/Data")
enzyme_data<-read.csv(file="BarEven11.txt",sep="\t",header=TRUE)#,header=TRUE)
enzyme_data$log10kcat<-log10(enzyme_data$kcat)
enzyme_data$log10KM<-log10(enzyme_data$KM)
enzyme_data$log10kcat_KM<-log10(enzyme_data$kcat/enzyme_data$KM*10^6)

###8.Analysis with data and evolutionary results on phenomenological values kcat/KM and KM fitness landscapes
N_reso=250
KT_set=c(5*10^-5,5*10^-3) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^3
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
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt13<-list(l=1.05,h=0.75,txt=expression(paste(DL[3])),srt = 45,font=3,col="black",cex=1)
addtxt14<-list(l=1.05,h=0.85,txt=expression(paste(DL[4])),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)
KMkcat_KmData<-as.data.frame(cbind(enzyme_data$log10kcat_KM,enzyme_data$log10kcat))
Data_list<-list(KMkcat_KmData)

multiplePlot("kr=","","Vm=","M/s",10^3,Vm[1],ncol,log10kcatKM_var,log10kcat_var,list(w_kcatKM[[1]]),
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),lev=c(0.9),
             xyData=Data_list,palette=pal,TEXT_to_Add =addtxt11,image=TRUE,col_data ="black",pcex=0.5,subcex=1,labcex=1.25)

abline(a=-0.3,b=1,col="grey",lty=6,lwd=2)
abline(a=-0.2,b=1,col="grey",lty=6,lwd=2)

text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=1)
#text(addtxt13$l,addtxt13$h,addtxt13$txt[1],srt=addtxt13$srt,font=addtxt13$font,col=addtxt13$col,cex=1,xpd=NA)
#text(addtxt14$l,addtxt14$h,addtxt14$txt[1],srt=addtxt14$srt,font=addtxt14$font,col=addtxt14$col,cex=1,xpd=NA)
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
for (i in 1:1){
  Cat_Data_list[[i]]<-as.data.frame(cbind(enzyme_data[enzyme_data$module_Type==levels(factor(enzyme_data$module_Type))[i],]$log10kcat_KM,
                                          enzyme_data[enzyme_data$module_Type==levels(factor(enzyme_data$module_Type))[i],]$log10kcat))
  Cat_list[i]<-levels(factor(enzyme_data$module_Type))[i]
}

legtext1=c("AA,F,N")

#Defining features of plots
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)

multiplePlot("kr=","","Vm=","M/s",kr,Vm[[1]],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Cat_Data_list,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,pcex=1,subcex=1,ncat=1,labcex=2.5,axcex=1.25,cextext=2)

contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lty=1,method ="flattest",xlim=c(0,1),add=TRUE)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
legend(x=1.05,y=0.8,title="Fitness isoclines",legend=c("amino acids","sugars"),lty=1,col=c(1,2),xpd=NA,cex=1.5)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)


legtext1=c("AFN","CE","IM","SM")


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
kf1_set<-c(10^2,10^10)
kcat1_set<-c(10^-2,10^6)
#kr1_set<-c(10^5,10^6,10^3,10^4)
#kinh1_set<-c(10^9,10^10,10^7,10^8)
#Ki<-10^-6
eta_d_set<-c(10^2,10^4)#degradation rate #Sensitivity study to be done relatively to Vm

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
    eta_d<-eta_d_set[l]*Vm
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
max(P_Phi_eq_sens_2Enz_background_eff[[2]])/max(P_Phi_eq_sens_2Enz_background_eff[[4]])

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_FD_sens_2Enz_background_eff<-list()
for (l in 1:(length(eta_d_set))){
  for (s in 1:(length(kcat1_set))){
    w_kf_FD_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]<-
      round(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]/max(P_Phi_eq_sens_2Enz_background_eff[[s+length(kcat1_set)*(l-1)]]),digits=20)
  }
}

#multiplePlot("KT=","","Vm=","M/s",signif(Vm_set,4),kr1,ncol,log10kf_var,log10kcat_var,Conc_prod_eq_sens_2Enz_background_eff,
#           abs="log10 (kf)",ord="log10 (kcat)",ncont=10,palette=pal,image=TRUE,scale="AUTO",labcex=1.25,sub="subtitle",subcex=2)
sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),
             expression(paste(eta,"=",10^-2,"/s")))

plot.new()
par(mfrow=c(1,1),oma=c(0,0,0,0),mai = c(0.75, 1, 0.75, 1))
addtxt<-list(l=0.035,h=0.975,txt=c("A","B"),srt = 0,font=2,col="black")

contour(w_kf_FD_sens_2Enz_background_eff[[1]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=1,lwd=1.5,lty=1,method ="flattest",xlim=c(0,1))
contour(w_kf_FD_sens_2Enz_background_eff[[3]],xaxt="n",yaxt="n",levels=c(0.9),labcex=1,col=2,lwd=1.5,lty=1,method ="flattest",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=1,lty=3,method ="edge",xlim=c(0,1),add=TRUE)
#contour(w_kf_FD_sens_2Enz_background_eff[[5+(l-1)*length(kcat1_set)]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=1,lty=5,method ="edge",xlim=c(0,1),add=TRUE)
title(main="",outer=FALSE,line=2.5,cex.main=1.5,font=2,xlab="log10 (kf)",ylab="log10 (kcat)",cex.lab=1.25)
legend(0.6,0.85,title="Degradation rate",legend=sublegend,bty="n",lty=c(1,1),col=c(1,2),ncol=1,cex=1)
axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
     labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1.25)
axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
     labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1.25)

for (l in 2:(length(eta_d_set))){
  for(s in 1:length(kcat1_set)){
    contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.9),lwd=1.5,labcex=1,col=s,lty=l+1,method = "flattest",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.999),labcex=1,col=s+1,lty=3,method = "edge",add=TRUE,xlim=c(0,1))
    #contour(w_kf_FD_sens_2Enz_background_eff[[(l-1)*length(kcat1_set)+s]],xaxt="n",yaxt="n",levels=c(0.99999),labcex=1,col=s+1,lty=5,method = "edge",add=TRUE,xlim=c(0,1))
    axis(1, at=seq(0,1,log10kf_var[3]/abs(log10kf_var[2]-log10kf_var[1])),
         labels=seq(log10kf_var[1],log10kf_var[2],log10kf_var[3]),cex.axis=1.25)
    axis(2, at=seq(0,1,log10kcat_var[3]/abs(log10kcat_var[2]-log10kcat_var[1])),
         labels=seq(log10kcat_var[1],log10kcat_var[2],log10kcat_var[3]),cex.axis=1.25)
    print((l-1)*length(kcat1_set)+s)
  }
}


rect(0.7, 0.65, 1.05, 1.05,
     col = "white", border = TRUE, lty = 1, lwd = 1)
#text(addtxt$l,addtxt$h,addtxt$txt[2],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
legend("topright",title="First enzyme efficiency",legend=c("Inefficient","  Perfect"),bty="n",lty=c(rep(1,3)),col=c(1,2),ncol=1,cex=0.7)


#title(main=bquote(atop("Variable degradation rates",
                    #   "and first enzyme efficiencies"))
    #  ,col.main="black",font.main=1,cex.main=1.1)

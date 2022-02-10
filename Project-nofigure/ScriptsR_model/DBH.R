setwd(dir="~")
setwd(dir="enzyme-evolution/ScriptsR_model")
source("AP11.Plotting_multiple_images.R")
library(Rmpfr)

library(lattice)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(plot.matrix)
library(fields)

jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

N_reso=100 #Resolution to build graphs
#Creation of two sets of parameters
KT=5*10^-5##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10Etot1_var<-c(-7,-1,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10Etot2_var<-c(-7,-1,2) #catalytic preferred to -1 for same reasons
kf<-10^6
kcat<-10^2
Se=c(10*KT)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(10^3) #An intermediate value seems reasonable to show general results
Etot1_var<-10^seq(log10Etot1_var[1],log10Etot1_var[2],length.out=N_reso)
Etot2_var<-10^seq(log10Etot2_var[1],log10Etot2_var[2],length.out=N_reso)
Ke<-10^-6

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list() #Equilibrium product flux
log10_P_Phi_eq<-list() #log10 of equilibrium product flux
kf_act<-c()
#for (l in 1:(length(kr_set))){
  #for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      #k_1=lr[l]*k2[j]
      kr=kr_set[1]
      for (j in 1:N_reso){
        E_tot_conc1<-Etot1_var[i]
        E_tot_conc2<-Etot2_var[j]
        #kf_act<-kf*(1-(E_tot_conc1+E_tot_conc2)/(E_tot_conc1+E_tot_conc2+Ke))
        kf_act<-kf*10^(-(E_tot_conc1+E_tot_conc2)/(5*10^-3))
        print(kf_act)
        a<- mpfr((Vm*kf_act+kf_act*E_tot_conc1*kcat*(1+Se/KT)),120)
        b<- mpfr((Vm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_tot_conc1*(KT+Se)),120)
        c<- mpfr(-Vm*Se*(kr+kcat),120)
        delta<-mpfr(b^2-4*a*c,120)
        S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
        ES_conc_eq=mpfr(E_tot_conc1*(kf_act*S_conc_eq)/(kr+kcat+kf_act*S_conc_eq),120)
        tab_P_Phi_eq[i,j]=as.numeric(kcat*ES_conc_eq)
      }
    }
    P_Phi_eq[[1]]<-tab_P_Phi_eq
    #log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10)
  #}
#}
    
w_BH<-list()
w_BH[[1]]<-P_Phi_eq[[1]]/max(P_Phi_eq[[1]])
    
multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10Etot1_var,log10Etot2_var,w_BH,
                 abs="log10 Etot1",ord="log10 Etot2",scale=c(0,1),lev=c(0.9,0.99,0.999,0.9999,0.99999),palette=pal,
             image=TRUE,pcex=1,subcex=1.5,labcex=1,axcex=1,colorkey=TRUE,globcex=0.5,legcex=1,contourlab=TRUE,meth="flattest",contcex=0.6)


###Regular scale

N_reso=100 #Resolution to build graphs
#Creation of two sets of parameters
KT=5*10^-5##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
#log10Etot1_var<-c(-7,-1,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
#log10Etot2_var<-c(-7,-1,2) #catalytic preferred to -1 for same reasons
kf<-10^6
kcat<-10^2
Se=c(10*KT)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(10^3) #An intermediate value seems reasonable to show general results
Etot1_var<-seq(10^-4,10^-2,length.out=N_reso)
Etot2_var<-seq(10^-4,10^-2,length.out=N_reso)
Ke<-10^-6

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list() #Equilibrium product flux
log10_P_Phi_eq<-list() #log10 of equilibrium product flux
kf_act<-c()
#for (l in 1:(length(kr_set))){
#for (s in 1:(length(Vm))){
for (i in 1:N_reso){
  #k_1=lr[l]*k2[j]
  kr=kr_set[1]
  for (j in 1:N_reso){
    E_tot_conc1<-Etot1_var[i]
    E_tot_conc2<-Etot2_var[j]
    #kf_act<-kf*(1-(E_tot_conc1+E_tot_conc2)/(E_tot_conc1+E_tot_conc2+Ke))
    kf_act<-kf*10^(-(E_tot_conc1+E_tot_conc2)/(5*10^-3))
    print(kf_act)
    a<- mpfr((Vm*kf_act+kf_act*E_tot_conc1*kcat*(1+Se/KT)),120)
    b<- mpfr((Vm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_tot_conc1*(KT+Se)),120)
    c<- mpfr(-Vm*Se*(kr+kcat),120)
    delta<-mpfr(b^2-4*a*c,120)
    S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
    ES_conc_eq=mpfr(E_tot_conc1*(kf_act*S_conc_eq)/(kr+kcat+kf_act*S_conc_eq),120)
    tab_P_Phi_eq[i,j]=as.numeric(kcat*ES_conc_eq)
  }
}
P_Phi_eq[[1]]<-tab_P_Phi_eq
#log10_P_Phi_eq[[s+length(Vm)*(l-1)]]<-log(P_Phi_eq[[s+length(Vm)*(l-1)]],10)
#}
#}

w_BH<-list()
w_BH[[1]]<-P_Phi_eq[[1]]/max(P_Phi_eq[[1]])

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,Etot1_var,Etot2_var,w_BH,
             abs="log10 Etot1",ord="log10 Etot2",scale=c(0,1),lev=c(0.9,0.99,0.999,0.9999,0.99999),palette=pal,
             image=TRUE,pcex=1,subcex=1.5,labcex=1,axcex=1,colorkey=TRUE,globcex=0.5,legcex=1,contourlab=TRUE,meth="flattest",contcex=0.6)


###Environonement stable à deux ressources (sans trade-off entre transporteurs)

N_reso=100 #Resolution to build graphs
#Creation of two sets of parameters
KT=5*10^-5##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10Etot1_var<-c(-7,-1,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10Etot2_var<-c(-7,-1,2) #catalytic preferred to -1 for same reasons
kf<-10^6
kcat<-10^2
Se=c(10*KT)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(10^3) #An intermediate value seems reasonable to show general results
Etot1_var<-10^seq(log10Etot1_var[1],log10Etot1_var[2],length.out=N_reso)
Etot2_var<-10^seq(log10Etot2_var[1],log10Etot2_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_2R_stable<-list() #Equilibrium product flux
kf_act<-c()
#for (l in 1:(length(kr_set))){
#for (s in 1:(length(Vm))){
for (i in 1:N_reso){
  #k_1=lr[l]*k2[j]
  kr=kr_set[1]
  for (j in 1:N_reso){
    E_tot_conc1<-Etot1_var[i]
    E_tot_conc2<-Etot2_var[j]
    #kf_act<-kf*(1-(E_tot_conc1+E_tot_conc2)/(E_tot_conc1+E_tot_conc2+Ke))
    kf_act<-kf*10^(-(E_tot_conc1+E_tot_conc2)/(5*10^-3))
    print(kf_act)
    a<- mpfr((Vm*kf_act+kf_act*(E_tot_conc1+E_tot_conc2)*kcat*(1+Se/KT)),120)
    b<- mpfr((Vm*(kr+kcat-Se*kf_act)+kf_act*kcat*(E_tot_conc1+E_tot_conc2)*(KT+Se)),120)
    c<- mpfr(-Vm*Se*(kr+kcat),120)
    delta<-mpfr(b^2-4*a*c,120)
    S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
    ES_conc_eq=mpfr((E_tot_conc1+E_tot_conc2)*(kf_act*S_conc_eq)/(kr+kcat+kf_act*S_conc_eq),120)
    tab_P_Phi_eq[i,j]=as.numeric(kcat*ES_conc_eq)
  }
}
P_Phi_eq_2R_stable[[1]]<-tab_P_Phi_eq


w_2R_stable<-list()
w_2R_stable[[1]]<-P_Phi_eq_2R_stable[[1]]/max(P_Phi_eq_2R_stable[[1]])

multiplePlot("","","","",KT,Vm[1],ncol=128,log10Etot1_var,log10Etot2_var,w_2R_stable,
             abs="log10 Etot1",ord="log10 Etot2",scale=c(0,1),lev=c(0.9,0.99,0.999,0.9999,0.99999),palette=pal,
             image=TRUE,pcex=1,subcex=1.5,labcex=1,axcex=1,colorkey=TRUE,globcex=0.5,legcex=1,contourlab=TRUE,meth="flattest",contcex=0.6)

###Environonement stable à deux ressources (en concentrations différentes)

N_reso=100 #Resolution to build graphs
#Creation of two sets of parameters
KT=5*10^-5##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10Etot1_var<-c(-7,-1,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10Etot2_var<-c(-7,-1,2) #catalytic preferred to -1 for same reasons
kf<-10^6
kcat<-10^2
Se1=c(10*KT)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
Se2=c(0.001*KT)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(10^3) #An intermediate value seems reasonable to show general results
Etot1_var<-10^seq(log10Etot1_var[1],log10Etot1_var[2],length.out=N_reso)
Etot2_var<-10^seq(log10Etot2_var[1],log10Etot2_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_2R_stable_diff<-list() #Equilibrium product flux
kf_act<-c()
#for (l in 1:(length(kr_set))){
#for (s in 1:(length(Vm))){
for (i in 1:N_reso){
  #k_1=lr[l]*k2[j]
  kr=kr_set[1]
  for (j in 1:N_reso){
    E_tot_conc1<-Etot1_var[i]
    E_tot_conc2<-Etot2_var[j]
    #kf_act<-kf*(1-(E_tot_conc1+E_tot_conc2)/(E_tot_conc1+E_tot_conc2+Ke))
    kf_act<-kf*10^(-(E_tot_conc1+E_tot_conc2)/(5*10^-3))
    print(kf_act)
    a1<- mpfr((Vm*kf_act+kf_act*(E_tot_conc1)*kcat*(1+Se1/KT)),120)
    b1<- mpfr((Vm*(kr+kcat-Se1*kf_act)+kf_act*kcat*(E_tot_conc1)*(KT+Se1)),120)
    c1<- mpfr(-Vm*Se1*(kr+kcat),120)
    delta1<-mpfr(b1^2-4*a1*c1,120)
    S_conc_eq1=mpfr((-b1+delta1^(1/2))/(2*a1),120)
    ES_conc_eq1=mpfr((E_tot_conc1)*(kf_act*S_conc_eq1)/(kr+kcat+kf_act*S_conc_eq1),120)
    
    a2<- mpfr((Vm*kf_act+kf_act*(E_tot_conc2)*kcat*(1+Se2/KT)),120)
    b2<- mpfr((Vm*(kr+kcat-Se2*kf_act)+kf_act*kcat*(E_tot_conc2)*(KT+Se2)),120)
    c2<- mpfr(-Vm*Se2*(kr+kcat),120)
    delta2<-mpfr(b2^2-4*a2*c2,120)
    S_conc_eq2=mpfr((-b2+delta2^(1/2))/(2*a2),120)
    ES_conc_eq2=mpfr((E_tot_conc2)*(kf_act*S_conc_eq2)/(kr+kcat+kf_act*S_conc_eq2),120)
    
    tab_P_Phi_eq[i,j]=as.numeric(kcat*ES_conc_eq1+kcat*ES_conc_eq2)
  }
}
P_Phi_eq_2R_stable_diff[[1]]<-tab_P_Phi_eq


w_2R_stable_diff<-list()
w_2R_stable_diff[[1]]<-P_Phi_eq_2R_stable_diff[[1]]/max(P_Phi_eq_2R_stable_diff[[1]])

multiplePlot("","","","",KT,Vm[1],ncol=128,log10Etot1_var,log10Etot2_var,w_2R_stable_diff,
             abs="log10 Etot1",ord="log10 Etot2",scale=c(0,1),lev=c(0.9,0.99,0.999,0.9999,0.99999),palette=pal,
             image=TRUE,pcex=1,subcex=1.5,labcex=1,axcex=1,colorkey=TRUE,globcex=0.5,legcex=1,contourlab=TRUE,meth="flattest",contcex=0.6)


setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R")
#Defining colors
ncol=128
pal<-colorRampPalette(c("hot pink","red","green"))(ncol)
#1.Variation of Se and facilitated diffusion transporters affinity

affinity=c("High affinity","Low affinity")
KT_set=c(0.002,0.1) #Molar, Michaelis-like constant for facilitated diffusion
Vm_set=c(4*10^-4,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
N_reso=101
#Creation of two sets of parameters
kr=10^3
log10kf_var<-c(0,10,2)
log10kcat_var<-c(-4,6,2)
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=c(0.01,0.5)#from low to high concentrations in the environment
#lr<-c(0,1)
E_tot_conc<-10^-3

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()
log10_P_Phi_eq<-list()

for (k in 1:(length(affinity))){
  KT<-KT_set[k]
  Vm<-Vm_set[k]
  for (s in 1:(length(Se))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        a<- mpfr((Vm*kf[i]+kf[i]*E_tot_conc*kcat[j]*(1+Se[s]/KT)),120)
        b<- mpfr((Vm*(kr+kcat[j]-Se[s]*kf[i])+kf[i]*kcat[j]*E_tot_conc*(KT+Se[s])),120)
        c<- mpfr(-Vm*Se[s]*(kr+kcat[j]),120)
        delta<-mpfr(b^2-4*a*c,120)
        S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),120)
        ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
        tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
      }
    }
    P_Phi_eq[[s+length(Se)*(k-1)]]<-tab_P_Phi_eq
    log10_P_Phi_eq[[s+length(Se)*(k-1)]]<-log(P_Phi_eq[[s+length(Se)*(k-1)]],10)
  }
}

palet=list(pal)
multiplePlot("","","[Se]=","M",affinity,Se,ncol,log10kf_var,log10kcat_var,P_Phi_eq,
             "Fitness landscape for the system facilitated diffusion + chemical reaction
             in relation to affinity transporters (kr=0)",
             "log10 k1","log10 k2",scale="AUTO",palette=palet,ncont=ncontour,image=TRUE)

#multiplePlot("","","[Se]=","M",affinity,Se,ncol,log10kf_var,log10kcat_var,log10_P_Phi_eq,
             #"Fitness landscape for the system facilitated diffusion + chemical reaction
             #in relation to affinity transporters (k-1=0)",
             #"log10 k1","log10 k2",palette=pal,ncont=5)

w_kf<-list()
for (a in 1:(length(affinity))){
  for (s in 1:(length(Se))){
    w_kf[[s+length(Se)*(a-1)]]<-round(P_Phi_eq[[s+length(Se)*(a-1)]]/max(P_Phi_eq[[s+length(Se)*(a-1)]]),digits=10)
  }
}

multiplePlot("","","[Se]=","M",affinity,Se,ncol,log10kf_var,log10kcat_var,w_kf,
             "Fitness landscape for the system facilitated diffusion + chemical reaction
             in relation to affinity transporters (kr=0)",
             "log10(kf)","log10(kcat)",scale="AUTO",palette=pal,lev=c(0.99,0.9999,0.999999))

#Comparing plots for different transporters

#Defining colors
ncol=128
pal<-colorRampPalette(c("blue","grey","yellow"))(ncol)
comp_flux_to_affinity<-list()
for (s in 1:(length(Se))){
  comp_flux_to_affinity[[s]]<-log10_P_Phi_eq[[s+length(Se)]]-log10_P_Phi_eq[[s]]
}

affinity_diff<-c("Yeast transporters extreme values")

multiplePlot("","","[Se]=","M",affinity_diff,Se,ncol,log10kf_var,log10kcat_var,comp_flux_to_affinity,
             "Fitness differences in relation to affinity ",
             "log10(kf)","log10(kcat)",scale=c(-1,1),palette=pal,ncont=ncontour)
#Blue: high affinity transporter moe efficient, yellow otherwise

#multiplePlot("","","[Se]=","M",affinity_diff,Se,ncol,log10k1_var,log10k2_var,comp_flux_to_affinity,
             #"Fitness differences in relation to affinity ",
             #"log10 k1","log10 k2",palette=pal,ncont=6,c(-0.8,0.8))

#2.Analysis with mechanistic values (not completely mechanistic however) and Vm variations
N_reso=101 #Resolution to build graphs
#Creation of two sets of parameters
KT=c(10^-1)##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=c(1,10^-4)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
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
        a<- mpfr((Vm[s]*kf[i]+kf[i]*E_tot_conc*kcat[j]*(1+Se[s]/KT)),120)
        b<- mpfr((Vm[s]*(kr+kcat[j]-Se[s]*kf[i])+kf[i]*kcat[j]*E_tot_conc*(KT+Se[s])),120)
        c<- mpfr(-Vm[s]*Se[s]*(kr+kcat[j]),120)
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

P_Phi_eq
palet=list(pal)
multiplePlot("KT=","M","Se=","M",KT,Se,ncol,log10kf_var,log10kcat_var,P_Phi_eq,
             "Fitness landscape for the system facilitated diffusion + chemical reaction
             in relation to affinity transporters (kr=0)",
             "log10(kf)","log10(kcat)",scale="AUTO",palette=palet,ncont=ncontour,image=TRUE)

w_kf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    w_kf[[s+length(Vm)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
  }
}

multiplePlot("KT=","M","[Se]=","M",KT,Se,ncol,log10kf_var,log10kcat_var,w_kf,
             "Fitness landscape for the system facilitated diffusion + chemical reaction
             in relation to affinity transporters (kr=0)",
             "log10(kf)","log10(kcat)",scale="AUTO",palette=palet,lev=c(0.99,0.9999,0.999999),image=TRUE)

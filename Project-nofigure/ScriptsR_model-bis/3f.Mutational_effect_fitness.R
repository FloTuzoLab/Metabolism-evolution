setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R") #To plot multiple gradient and color plots
source("AP2.RaphsonNewton.R") #To find equilibrium for 2 variables problems
require("Rmpfr") #to find equilibrium when very different values (dozens orders of magnitude) matter

###Plot configuration
options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncontour=10
#Color palette (128 colors)
ncol=128
jet.pal<-colorRampPalette(c("hot pink","red","green"))
palet <- jet.pal(ncol)

###Comparing flux when mutation occuring on one enzyme

##0.Calculating flux with all enzymes fixed at the same activity level

##With degradation in the pathway
N_reso=51
#Defining parameters
Se=10^-2
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
kr_set<-c(10^-1,10^5)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]
#Best possible values for the first reaction
kf1_set<-kf
kcat1_set<-kcat
kr1_set<-kr_set
#Ki<-10^-6
K_d<-10^-1
#eta_d<-c(10^-6,10^-3)
#eta_d_set<-c(10^-7,10^-4)

kf2_set<-kf
kcat2_set<-kcat
kr2_set<-kr_set##Testing with 2 values, but only one should be used later on in the paper
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc

error=10^-13
init<-c(10^-25,10^-25)

tab_P_Phi_eq <- matrix(nrow=N_reso, ncol=N_reso)
P_Phi_eq_id<-list() #Equilibrium product flux
MaxconcP<-0

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    for (i in 1:(N_reso)){
      Vm=Vm_set[s]
      eta_d<-Vm*0.1
      V_d=Vm
      for (j in 1:N_reso){
        print(paste(i,j))
        kf1<-kf1_set[i]
        kf2<-kf2_set[i]
        kcat1<-kcat1_set[j]
        kcat2<-kcat2_set[j]
        kr1=kr1_set[l]
        kr2=kr2_set[l]
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #eta_d<-eta_d_set[s+length(Vm_set)*(l-1)]
        #init<-c(10^-15,10^-15)
        #F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+V_d/(y+K_d))*y)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
      }
    }
    P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
  }
}

MaxconcP
#b.Analysing the results
##Ploting the equilibrium flux with regards to concentrations
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,P_Phi_eq_id,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",palette=palet,ncont=10) ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    w_kf[[s+length(Vm_set)*(l-1)]]<-round(P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]/max(P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]),digits=12)
  }
}

##Plotting the relative fitness
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,w_kf,
             "Relative fitness landscape for the system
             {facilitated diffusion+ 2 reactions}",
             "log10 kf","log10 kcat",palette=palet,scale=c(0,1),lev=c(0.5,
                                                                      0.99,
                                                                      0.9999,#below any Ne values
                                                                      0.999999#,#lower values for Ne
             ))#rounding to 8 digits, higher Ne values

##1.Mutation on kf

#a.Equilibrium calculus
N_reso=51
##Defining parameters
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=10^-2#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)

E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]

##Defining the whole variable set of parameters for both inside reactions
#delta_mut_kf=1
delta_mut_kf=0.1
#delta_mut_kcat=0.1
kf1_set<-kf
kcat1_set<-kcat
kr1_set<-kr_set
kf2_set<-10^seq(log10kf_var[1]-delta_mut_kf,log10kf_var[2]-delta_mut_kf,length.out=N_reso)
kcat2_set<-kcat
#kcat2_set<-10^seq(log10kcat_var[1]-delta_mut_kcat,log10kcat_var[2]-delta_mut_kcat,length.out=N_reso)
kr2_set<-kr_set
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc
#eta_d_set<-c(10^-7,10^-2,10^-4,10^-2)#degradation rate
##Error ceiling for the equilibrium finding
init<-c(10^-25,10^-25)
error=10^-13
MaxconcP<-0

##Finding the steady-state with numerical Raphson Newton method
tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_conc_P_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_mutkf<-list() #Equilibrium product flux
conc_P_eq_mutkf<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    for (i in 1:(N_reso)){
      Vm=Vm_set[s]
      eta_d=Vm*0.1
      V_d=Vm
      for (j in 1:N_reso){
        print(paste(i,j))
        kf1<-kf1_set[i]
        kf2<-kf2_set[i]
        kcat1<-kcat1_set[j]
        kcat2<-kcat2_set[j]
        kr1=kr1_set[l]
        kr2=kr2_set[l]
        Vm1<-kcat1*E1_conc
        KM1<-(kcat1+kr1)/kf1
        Vm2<-kcat2*E2_conc
        KM2<-(kcat2+kr2)/kf2
        T1<-KT+Se
        T2<-1+Se/KT
        #init<-c(10^-15,10^-15)
        #F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+V_d/(y+K_d))*y)
        F1e<-expression(Vm1*x/(KM1+x)-(Vm2/(KM2+y)+eta_d)*y)
        F2e<-expression(Vm*(Se-x)/(T1+x*T2)-Vm1*x/(KM1+x))
        Flist<-list(F1e,F2e)
        SSR<-RaphsonNewton(Flist,init,error)
        #print(paste(SSR$x,SSR$y))
        if(SSR$y>MaxconcP){
          MaxconcP<-SSR$y
        }
        tab_P_Phi_eq[i,j]=Vm2*SSR$y/(KM2+SSR$y)
        tab_conc_P_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_mutkf[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    conc_P_eq_mutkf[[s+length(Vm_set)*(l-1)]]<-tab_conc_P_eq
  }
}

MaxconcP
#
delta_P_Phi_eq_mutkf<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    delta_P_Phi_eq_mutkf[[s+length(Vm_set)*(l-1)]]<-round((P_Phi_eq_mutkf[[s+
    length(Vm_set)*(l-1)]]-P_Phi_eq_id[[s+length(Vm_set)*(l-1)]])/
    P_Phi_eq_id[[s+length(Vm_set)*(l-1)]],digits=15)
  }
}

multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,conc_P_eq_mutkf,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",palette=palet,lev=c(0.01,0.001,0.0001))

multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,P_Phi_eq_mutkf,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",palette=palet,ncont=10)

multiplePlot("kr=","","Vm=","M/s",kr2_set,Vm_set,ncol,log10kf_var,log10kcat_var,delta_P_Phi_eq_mutkf,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",scale=c(-1,0),palette=palet,lev=c(-10^-4,-10^-6,-10^-8))

##i.Mutation on kcat
#a.Equilibrium calculus
N_reso=101
##Defining parameters
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm_set=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=10^-2#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(10^-1,10^5) #An intermediate value seems reasonable to show general results
E_tot_conc=10^-3##c(10^-7,10^-4) Both values need to be studied as what matters most is kf[Etot]

##Defining the whole variable set of parameters for both inside reactions
delta_mut_kcat<-1
kf1_set<-kf
kcat1_set<-kcat
kr1_set<-kr_set
kf2_set<-kf#10^seq(log10kf_var[1]-delta_mut_kf,log10kf_var[2]-delta_mut_kf,length.out=N_reso)
kcat2_set<-10^seq(log10kcat_var[1]-delta_mut_kcat,log10kcat_var[2]-delta_mut_kcat,length.out=N_reso)
kr2_set<-kr_set
E1_conc<-E_tot_conc
E2_conc<-E_tot_conc
#eta_d<-10^-5#degradation rate
##Error ceiling for the equilibrium finding
error=10^-10

##Finding the steady-state with numerical Raphson Newton method
tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_conc_P_eq<-matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq_mutkcat<-list() #Equilibrium product flux
conc_P_eq_mutkcat<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    for (i in 1:(N_reso)){
      Vm=Vm_set[s]
      eta_d=Vm*10^-1
      for (j in 1:N_reso){
        print(paste(i,j))
        kf1<-kf1_set[i]
        kf2<-kf2_set[i]
        kcat1<-kcat1_set[j]
        kcat2<-kcat2_set[j]
        kr1=kr1_set[l]
        kr2=kr2_set[l]
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
        tab_conc_P_eq[i,j]=SSR$y
      }
    }
    P_Phi_eq_mutkcat[[s+length(Vm_set)*(l-1)]]<-tab_P_Phi_eq
    conc_P_eq_mutkcat[[s+length(Vm_set)*(l-1)]]<-tab_conc_P_eq
  }
}

#
delta_P_Phi_eq_mutkcat<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm_set))){
    delta_P_Phi_eq_mutkcat[[s+length(Vm_set)*(l-1)]]<-(P_Phi_eq_mutkcat[[s+length(Vm_set)*(l-1)]]-P_Phi_eq_id[[s+length(Vm_set)*(l-1)]])/P_Phi_eq_id[[s+length(Vm_set)*(l-1)]]
  }
}

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set,ncol,log10kf_var,log10kcat_var,delta_P_Phi_eq_mutkcat,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10 kf","log10 kcat",palette=palet,scale=c(-1,0.01),lev=c(-0.01,-0.0001,-0.00000001))


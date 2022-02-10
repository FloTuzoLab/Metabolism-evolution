setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R")
setwd(dir="Data")

###Plot configuration
options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
ncontour=5
palet<-colorRampPalette(c("hot pink","red","green"))(ncol)

###0.Drawing of fitness landscapes for passive diffusion

r_cell=1*10^-6#unit:m
SA<-4*pi*r_cell^2
V<-4/3*pi*r_cell^3

Se_high=10^-1#very high concentration
Se_low=10^-3#very high concentration
N_reso=101
Pd=10^-12#unit:m/s true glucose value:1e-12
E_tot_conc=10^-3#unit:M
log10kf_var<-c(0,10,2)
log10kcat_var<-c(-4,6,2)
kr<-10^3#intermediate level for reverse reaction constant
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
tab_P_Phi_eq_high <- matrix(nrow=N_reso , ncol=N_reso)
tab_P_Phi_eq_low <- matrix(nrow=N_reso , ncol=N_reso)
for (i in 1:N_reso){
  for (j in 1:N_reso){
    a<-mpfr(-Pd*SA/V*kf[i],120)
    b<-mpfr(Pd*SA/V*(Se_high*kf[i]-(kr+kcat[j]))-kcat[j]*kf[i]*E_tot_conc,120)
    c<-mpfr(Pd*SA/V*Se_high*(kr+kcat[j]),120)
    delta_pd<-mpfr(b^2-4*a*c,120)
    S_conc_eq=mpfr((-b-delta_pd^(1/2))/(2*a),120)
    ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
    tab_P_Phi_eq_high[i,j]=as.numeric(kcat[j]*ES_conc_eq)
  }
}

for (i in 1:N_reso){
  for (j in 1:N_reso){
    a<-mpfr(-Pd*SA/V*kf[i],120)
    b<-mpfr(Pd*SA/V*(Se_low*kf[i]-(kr+kcat[j]))-kcat[j]*kf[i]*E_tot_conc,120)
    c<-mpfr(Pd*SA/V*Se_low*(kr+kcat[j]),120)
    delta_pd<-mpfr(b^2-4*a*c,120)
    S_conc_eq=mpfr((-b-delta_pd^(1/2))/(2*a),120)
    ES_conc_eq=mpfr(E_tot_conc*(kf[i]*S_conc_eq)/(kr+kcat[j]+kf[i]*S_conc_eq),120)
    tab_P_Phi_eq_low[i,j]=as.numeric(kcat[j]*ES_conc_eq)
  }
}

P_Phi_eq<-list(tab_P_Phi_eq_low,tab_P_Phi_eq_high)
#Test results with image.plot
#image.plot(tab_P_Phi_eq)

#3D matrix plot of steady state flux
par(mfrow=c(1,1),pin=c(7,5))
x_log10kf<-seq(0,10,length=N_reso)
y_log10kcat<-seq(-4,6,length=N_reso)
fluxfacet <- (tab_P_Phi_eq_high[-1, -1] + tab_P_Phi_eq_high[-1, -N_reso] + tab_P_Phi_eq_high[-N_reso, -1] + tab_P_Phi_eq_high[-N_reso, -N_reso])/4
P_cut<-cut(fluxfacet,ncol)
persp(x_log10kf,y_log10kcat,xlab="log10 (kf)",ylab="log10 (kcat)",zlab="Equilibrium flux (M/s)",tab_P_Phi_eq_high,theta = -30, phi = 20, expand = 0.5,
      d=0.75,shade=0.5, col=palet[P_cut],border=NA,axes = T, ticktype="detailed",nticks = 5,
      box=T,cex.lab=0.75,cex.axis=0.4)
#2D comparison between high and low concentrations
par(mfrow=c(1,2),pin=c(4,5))
##Ploting the equilibrium flux with regards to concentrations
pal<-list(palet)
multiplePlot("kr=","","Se=","M",kr,c(Se_low,Se_high),ncol,log10kf_var,log10kcat_var,P_Phi_eq,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10,palette=pal,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_PD<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(c(Se_low,Se_high)))){
    w_kf_PD[[s+length(c(Se_low,Se_high))*(l-1)]]<-round(P_Phi_eq[[s+length(c(Se_low,Se_high))*(l-1)]]/max(P_Phi_eq[[s+length(c(Se_low,Se_high))*(l-1)]]),digits=10)
  }
} 
#Remark: the width of the plateau results from the low flux. Increasing it does reduce its width.

#Plotting the fitness
multiplePlot("kr=","","Se=","M",kr,c(Se_low,Se_high),ncol,log10kf_var,log10kcat_var,w_kf_PD,
             "Absolute fitness landscape for the system
{passive diffusion + chemical reaction}",
             "log10 kf","log10 kcat",palette=pal,image=TRUE,scale="AUTO",lev=c(0.9,0.99,0.9999,0.999999)) ##The unit and the value of Vm need to be checked and explained

##1.Optimal flux with facilitated diffusion

#I.Variation of Vm(flux thourgh membrane) and kr(reverse)
#1a.Analysis with mechanistic values (not completely mechanistic however) and Vm variations
N_reso=201 #Resolution to build graphs
#Creation of two sets of parameters
KT=10^-3##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-4.5,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=c(10^-2)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
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

#3D matrix plot of steady state flux
par(mfrow=c(1,1),pin=c(7,5))
x_log10kf<-seq(0,10,length=N_reso)
y_log10kcat<-seq(-4,6,length=N_reso)
fluxfacet <- (P_Phi_eq[[1]][-1, -1] + P_Phi_eq[[1]][-1, -N_reso] + P_Phi_eq[[1]][-N_reso, -1] + P_Phi_eq[[1]][-N_reso, -N_reso])/4
P_cut<-cut(fluxfacet,ncol)
persp(x_log10kf,y_log10kcat,xlab="log10 (kf)",ylab="log10 (kcat)",zlab="Equilibrium flux (M/s)",P_Phi_eq[[1]],theta = -30, phi = 15, expand = 0.6,
      d=0.75,shade=0.5, col=palet[P_cut],border=NA,axes = T, ticktype="detailed",nticks = 5,
      box=T,cex.lab=0.75,cex.axis=0.4)

##Ploting the equilibrium flux with regards to concentrations
Vm_set_plot=signif(c(10^-6,10^-4.5,10^-3),3)

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set_plot,ncol,log10kf_var,log10kcat_var,P_Phi_eq,
             "Absolute fitness landscape for the system
{facilitated diffusion + chemical reaction}",
             "log10 kf","log10 kcat",ncont=10,palette=pal,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_Vm<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    w_kf_Vm[[s+length(Vm)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
  }
}


##Plotting the relative fitness #Need to increase N_reso for elegant graphical representation ==201 or 251
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm_set_plot,ncol,log10kf_var,log10kcat_var,w_kf_Vm,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kf","log10 kcat",palette=pal,scale=c(0,1),lev=c(0.5,
                                                      0.99,
                                                      0.9999,#below any Ne values
                                                      0.999999#,#lower values for Ne
                                                      ),image=TRUE)#rounding to 8 digits, higher Ne values


##Plotting the contour for three different Vm
textlegend<-c("Vm=1e-06","Vm=1e-04.5","Vm=1e-03")

multiplePlot("kr=","","Vm=","M/s",kr,Vm_set_plot,ncol,log10kf_var,log10kcat_var,w_kf_Vm,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction} and many (Vm) values",
             "log10 kf","log10 kcat",
             palette=pal,lev=c(0.1,0.5,0.9),scale=c(0,1),image=FALSE,ncat=2,pcex=1,legpos="topright",legtext = textlegend,labcex=1)

#1b.Analysis with mechanistic values (not completely mechanistic however) and kr variations
N_reso=201 #Resolution to build graphs
#Creation of two sets of parameters
KT=10^-3##KT=c(0.001,0.01) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
log10kf_var<-c(0,10,2) #forward preferred to 1 as forward is already an emerging constant from diffusion and association
log10kcat_var<-c(-4,6,2) #catalytic preferred to -1 for same reasons
kf<-10^seq(log10kf_var[1],log10kf_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)
Se=c(10^-2)#c(10^-5,10^-3,0.1) #Environment concentration set to quasi saturation for transporter (tested for lower values but not very relevant as the cell would produce transporters with better affinities)
#lr<-c(1/100,100) #Used to set kr as a function of kcat
kr_set<-c(0,10^3,10^6) #An intermediate value seems reasonable to show general results
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

##Calculating the relative fitness with regards to the maximum achievable flux
w_kf_kr<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    w_kf_kr[[s+length(Vm)*(l-1)]]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
  }
}

##Plotting the contour for three different Vm
textlegend<-c("kr=0","kr=1e+03","kr=1e+06")

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10kf_var,log10kcat_var,w_kf_kr,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction} and many (Vm) values",
             "log10 kf","log10 kcat",
             palette=pal,lev=c(0.1,0.5,0.9),scale=c(0,1),image=FALSE,pcex=1,legpos="topright",legtext = textlegend,labcex=1)


#2.Analysis with phenomenological values kcat and KM
N_reso=251
KT=0.001 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se=10^-2
##Defining variables
log10KM_var<-c(-9,1,2)#with USI value that need to be changed for representation
KM<-10^seq(log10KM_var[1],log10KM_var[2],length.out=N_reso)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){ #Se can also be tested
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        kf=(kcat[j]+kr_set[l])/KM[i]
        if(kf/10^10>1){
          #k1=10^10
        #}
         tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          kr<-kr_set[l]
          #print(c("k1",k1))
          #print(c("k-1",k_1))
          a<- mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),60)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),60)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),60)
          delta<-mpfr(b^2-4*a*c,60)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),60)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),60)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_KM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_KM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_KM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_KM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_KM[[s+length(Vm)*(l-1)]]<-w_KM_tab
  }
}

##Plotting results including all enzymes data
#Importing the data
enzyme_data<-read.csv(file="BarEven11.txt",sep="\t",header=TRUE)
#Transforming the data
enzyme_data$log10kcat<-log10(enzyme_data$kcat)
enzyme_data$log10KM<-log10(enzyme_data$KM)
enzyme_data$log10kcat_KM<-log10(enzyme_data$kcat/enzyme_data$KM*10^6)#*10^6 to visualize in MicroMolar
KMkcatData<-as.data.frame(cbind(enzyme_data$log10KM,enzyme_data$log10kcat))

#KMkcatData<-as.data.frame(cbind(enzyme_data$log10KM+2,enzyme_data$log10kcat+2))#+2 to account for in vivo effects, but not relevant here as results are already significant
log10KM_var<-c(-3,7,2)#Rescaling to show KM in micromolar
##Ploting the equilibrium flux with regards to concentrations #N_reso=101
addtxt<-list(l=0.1,h=0.5,txt="Diffusion limit",srt = 90,font=3,col="red")
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10KM_var,log10kcat_var,P_Phi_eq,
             "Absolute fitness landscape for the system {facilitated diffusion
             + chemical reaction}",
             "log10(KM)","log10(kcat)",palette=pal,ncont=10,xyData=KMkcatData,TEXT_to_Add = addtxt,image=TRUE,scale="AUTO") ##The unit and the value of Vm need to be checked and explained


multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10KM_var,log10kcat_var,w_KM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 KM","log10 kcat",palette=pal,scale=c(0,1),lev=c(0.5,
                                                      0.99,
                                                      0.999,#below any Ne values
                                                      0.9999,#lower values for Ne
                                                      0.999999),TEXT_to_Add =addtxt,image=TRUE,scale="AUTO")#rounding to 8 digits, higher Ne values


##Including only enzymes for the few widely studied species (Bar-Even et al.,11)

#3.Plotting with kcat/KM as a variable

#Creation of kcat/KM set of parameters

log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-0.01
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

KMkcat_KmData<-as.data.frame(cbind(enzyme_data$log10kcat_KM,enzyme_data$log10kcat))
#KMkcat_KmData<-as.data.frame(cbind(enzyme_data$log10kcat_KM,enzyme_data$log10kcat+2))
addtxt<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=3,col="red")
Data_list<-list(KMkcat_KmData)
multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
             {facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",scale=c(0,1),
             xyData=Data_list,palette=pal,ncont=10,TEXT_to_Add = addtxt,image=TRUE,scale="AUTO")

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10kcatKM_var,log10kcat_var,P_Phi_eq,
             "Relative fitness landscape for the system
             {facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",palette=pal,scale=c(0,max(P_Phi_eq[[1]])),ncont=10,
             xyData=KMkcat_KmData,image=TRUE,scale="AUTO")


#4.Analyzing with kcat[E]tot

#Creation of the new set of parameter
log10E_tot_var<-c(-12,-2,2)
E_tot_conc<-10^seq(log10E_tot_var[1],log10E_tot_var[2],length.out=N_reso)
kcat=10^4 #Quite efficient catalytic property for the enzyme

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        kr<-kr_set[l]
        #kf=kcatKM[i]*(kr+kcat[j])/kcat[j]
        kf=kcatKM[i]*(1+kr/kcat)
        if(kf/10^10>1){
          #k1=10^10
          #}
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          #print(c("k1",k1))
          #print(c("k-1",k_1))
          a<- (Vm[s]*kf+kf*E_tot_conc[j]*kcat*(1+Se/KT))
          b<- (Vm[s]*(kr+kcat-Se*kf)+kf*kcat*E_tot_conc[j]*(KT+Se))
          c<- -Vm[s]*Se*(kr+kcat)
          delta<-b^2-4*a*c
          S_conc_eq=(-b+delta^(1/2))/(2*a)
          ES_conc_eq=E_tot_conc[j]*(kf*S_conc_eq)/(kr+kcat+kf*S_conc_eq)
          tab_P_Phi_eq[i,j]=kcat*ES_conc_eq
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_Etot_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_Etot_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_Etot_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_Etot_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_Etot_kcatKM[[s+length(Vm)*(l-1)]]<-w_Etot_kcatKM_tab
  }
}

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10kcatKM_var,log10E_tot_var,w_Etot_kcatKM,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 kcat/KM","log10 [Etot]",palette=pal,scale=c(0,1),lev=c(0.9,0.99,0.999,0.9999),image=TRUE,scale="AUTO")#rounding to 8 digits, higher Ne values

multiplePlot("kr=","","Vm=","M/s",kr_set,Vm,ncol,log10kcatKM_var,log10E_tot_var,P_Phi_eq,
             "Fitness landscape for the system facilitated diffusion
             + chemical reaction",
             "log10 kcat/KM","log10 [Etot]",palette=pal,ncont=10,image=TRUE,scale="AUTO")#rounding to 8 digits, higher Ne values

setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles")
source("AP1.Plotting_multiple_images.R")
source("5c.Enzyme_Evolutionary_Model_Results.R")
setwd(dir="Data")
require("Rmpfr")
  
  
##0.Importing the data
setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/Data")
enzyme_data<-read.csv(file="BarEven11.txt",sep="\t",header=TRUE)#,header=TRUE)
enzyme_data$log10kcat<-log10(enzyme_data$kcat)
enzyme_data$log10KM<-log10(enzyme_data$KM)
enzyme_data$log10kcat_KM<-log10(enzyme_data$kcat/enzyme_data$KM*10^6)

#Defining features of plots
ncol=128
ncontour=5
pal<-colorRampPalette(c("hot pink","red","green"))(ncol)

#1.Analysis with data and evolutionary results on phenomenological values kcat/KM and KM fitness landscapes
N_reso=100
KT=5*10^-5 #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se=10*KT
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
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),60)
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
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr))){
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

#load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/evo_model.rda")
###Need to open 5a. and to run it in order to produce results


#Results[["Ne3,bi01"]]

addtxt<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=3,col="red")

KMkcat_KmData<-as.data.frame(cbind(enzyme_data$log10kcat_KM,enzyme_data$log10kcat))
#KMkcat_KmData<-as.data.frame(cbind(enzyme_data$log10kcat_KM,enzyme_data$log10kcat+2))
#dev.on(width=5, height=4)
Data_list<-list(KMkcat_KmData)
palet<-list(pal)
multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",scale=c(0,1),lev=c(0.1,0.5,0.9),
             xyData=Data_list,palette=palet,TEXT_to_Add =addtxt,image=TRUE,col_data ="white",pcex=0.5)

multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,P_Phi_eq,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",ncont=10,
             xyData=Data_list,palette=palet,TEXT_to_Add =addtxt,scale="AUTO",image=TRUE)


Res<-list(Results[["Ne3,bi01"]])
multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",scale=c(0,1),lev=c(0.5,0.999,0.9999),
             xyData=Res,palette=palet,TEXT_to_Add =addtxt,image=TRUE)

Res<-list(Results[["Ne3,bi02"]])
multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",scale=c(0,1),lev=c(0.5,0.999,0.9999),
             xyData=Res,palette=palet,TEXT_to_Add =addtxt,image=TRUE)


Res<-list(Results[["Ne5,bi01"]])
multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",scale=c(0,1),lev=c(0.5,0.99999,0.999999),
             xyData=Res,palette=palet,TEXT_to_Add =addtxt,image=TRUE)

Res<-list(Results[["Ne5,bi02"]])
multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",scale=c(0,1),lev=c(0.5,0.99999,0.999999),
             xyData=Res,palette=palet,TEXT_to_Add =addtxt,image=TRUE)


Contour_results<-list(Results[["Ne3,bi01"]],Results[["Ne3,bi02"]],Results[["Ne5,bi01"]],Results[["Ne5,bi02"]])

ncol=64
pal<-colorRampPalette(c("white"))(ncol)
Ne_text<-c("Ne=1e+03","Ne=1e+05")
bias_text<-c("bias=-0.1","bias=-0.2")
legend_text<-c()
for (n in 1:length(Ne_text)){
  for (b in 1:length(bias_text)){
    legend_text<-c(legend_text,paste(Ne_text[n],",",bias_text[b]))
  }
}
palet<-list(pal)
multiplePlot("kr=","","Vm=","M/s",kr,Vm,ncol,log10kcatKM_var,log10kcat_var,w_kcatKM,
             "Relative fitness landscape for the system
{facilitated diffusion+ chemical reaction}",
             "log10 kcat/KM","log10 kcat",
             xyData=Contour_results,palette=palet,TEXT_to_Add =addtxt,lev=c(0.5,0.999,0.9999,0.99999,0.999999),scale=c(0,1),image=TRUE,ncat=2,pcex=1,labcex=1,legtext=legend_text,legpos="topleft")


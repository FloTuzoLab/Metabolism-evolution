setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/")
source("AP1.Plotting_multiple_images.R") #To plot multiple gradient and color plots
source("AP2.RaphsonNewton.R") #To find equilibrium for 2 variables problems
options(digits=20)

load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_equilibrium_time/eq100_4.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_equilibrium_time/eq1000_4.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_equilibrium_time/eq1000_5.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_equilibrium_time/eq10000_5.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_equilibrium_time/eq10000_6.rda")

#Testing the time to reach equilibrium
time_results
stime_results<-evo_results
mtime_results<-evo_results1000_5
ltime_results<-evo_results10000_6
P_Phi_ev_s<-list()
P_Phi_ev_m<-list()
P_Phi_ev_l<-list()

for(j in 1:10){
  P_Phi_ev_s[[j]]<-c(stime_results[[j]][[1]][3])
  P_Phi_ev_m[[j]]<-c(mtime_results[[j]][[1]][3])
  P_Phi_ev_l[[j]]<-c(ltime_results[[j]][[1]][3])
  for (i in 2:15){
    #print(P_Phi_ev)
    P_Phi_ev_s[[j]][i]<-stime_results[[j]][[i]][3]
    P_Phi_ev_m[[j]][i]<-mtime_results[[j]][[i]][3]
    P_Phi_ev_l[[j]][i]<-ltime_results[[j]][[i]][3]
  }
  #tab_P_Phi_ev[[j]]<-P_Phi_ev
}
moy_s<-c()
ec_s<-c()
moy_m<-c()
ec_m<-c()
moy_l<-c()
ec_l<-c()
for (i in 1:10){
  moy_s[i]=mean(P_Phi_ev_s[[i]])
  ec_s[i]<-sd(P_Phi_ev_s[[i]])
  moy_m[i]=mean(P_Phi_ev_m[[i]])
  ec_m[i]<-sd(P_Phi_ev_m[[i]])
  moy_l[i]=mean(P_Phi_ev_l[[i]])
  ec_l[i]<-sd(P_Phi_ev_l[[i]])
}
par(mfrow=c(3,2))
plot(moy_s,pch=19,col="blue",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=100,mu=1e-04")
plot(ec_s,pch=19,col="blue",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=100,mu=1e-04")
plot(moy_m,pch=19,col="blue",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=1000,mu=1e-05")
plot(ec_m,pch=19,col="blue",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(moy_l,pch=19,col="blue",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=10000,mu=1e-06")
plot(ec_l,pch=19,col="blue",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=10000,mu=1e-06")

mtime_results_quick<-evo_results1000_4
mtime_results_slow<-evo_results1000_5
ltime_results_quick<-evo_results10000_5
ltime_results_slow<-evo_results10000_6
P_Phi_ev_m_slow<-list()
P_Phi_ev_m_quick<-list()
P_Phi_ev_l_slow<-list()
P_Phi_ev_l_quick<-list()

for(j in 1:10){
  P_Phi_ev_m_slow[[j]]<-c(mtime_results_slow[[j]][[1]][3])
  P_Phi_ev_m_quick[[j]]<-c(mtime_results_quick[[j]][[1]][3])
  P_Phi_ev_l_slow[[j]]<-c(ltime_results_slow[[j]][[1]][3])
  P_Phi_ev_l_quick[[j]]<-c(ltime_results_quick[[j]][[1]][3])
  for (i in 2:15){
    #print(P_Phi_ev)
    P_Phi_ev_m_slow[[j]][i]<-mtime_results_slow[[j]][[i]][3]
    P_Phi_ev_m_quick[[j]][i]<-mtime_results_quick[[j]][[i]][3]
    P_Phi_ev_l_slow[[j]][i]<-ltime_results_slow[[j]][[i]][3]
    P_Phi_ev_l_quick[[j]][i]<-ltime_results_quick[[j]][[i]][3]
  }
  #tab_P_Phi_ev[[j]]<-P_Phi_ev
}

moy_m_slow<-c()
moy_m_quick<-c()
ec_m_slow<-c()
ec_m_quick<-c()

for (i in 6:10){
  moy_m_slow[i]=mean(P_Phi_ev_m_slow[[i]])
  ec_m_slow[i]<-sd(P_Phi_ev_m_slow[[i]])
  moy_m_quick[i]=mean(P_Phi_ev_m_quick[[i]])
  ec_m_quick[i]<-sd(P_Phi_ev_m_quick[[i]])
}
par(mfrow=c(2,2))
plot(moy_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=1000,mu=1e-05")
plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(moy_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=1000,mu=1e-04")
plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")


###
moy_l_slow<-c()
moy_l_quick<-c()
ec_l_slow<-c()
ec_l_quick<-c()

for (i in 6:10){
  moy_l_slow[i]=mean(P_Phi_ev_l_slow[[i]])
  ec_l_slow[i]<-sd(P_Phi_ev_l_slow[[i]])
  moy_l_quick[i]=mean(P_Phi_ev_l_quick[[i]])
  ec_l_quick[i]<-sd(P_Phi_ev_l_quick[[i]])
}
par(mfrow=c(2,2))
plot(moy_l_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=10000,mu=1e-06")
plot(ec_l_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=10000,mu=1e-06")
plot(moy_l_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=10000,mu=1e-05")
plot(ec_l_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=10000,mu=1e-05")

setwd(dir="~")
setwd(dir="Desktop/ScriptsR_modeles/Analysis/Evo_Results")

load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne10000_mu5_bi00.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne1000_mu4_bi01.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne10000_mu5_bi01.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne100000_mu6_bi01.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne1000_mu4_bi02.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne10000_mu5_bi02.rda")
load("/Users/florian/Desktop/ScriptsR_modeles/Analysis/Evo_Results/eq_Ne100000_mu6_bi02.rda")
#eq_Ne1000_mu4_bi01


P_Phi_ev_Ne10000_mu5_bi00<-list()
P_Phi_ev_Ne1000_mu4_bi01<-list()
P_Phi_ev_Ne1000_mu4_bi02<-list()
P_Phi_ev_Ne10000_mu5_bi01<-list()
P_Phi_ev_Ne10000_mu5_bi02<-list()
P_Phi_ev_Ne100000_mu6_bi01<-list()
P_Phi_ev_Ne100000_mu6_bi02<-list()

kcat_Ne10000_mu5_bi00<-list()
kcat_Ne1000_mu4_bi01<-list()
kcat_Ne1000_mu4_bi02<-list()
kcat_Ne10000_mu5_bi01<-list()
kcat_Ne10000_mu5_bi02<-list()
kcat_Ne100000_mu6_bi01<-list()
kcat_Ne100000_mu6_bi02<-list()

kf_Ne10000_mu5_bi00<-list()
kf_Ne1000_mu4_bi01<-list()
kf_Ne1000_mu4_bi02<-list()
kf_Ne10000_mu5_bi01<-list()
kf_Ne10000_mu5_bi02<-list()
kf_Ne100000_mu6_bi01<-list()
kf_Ne100000_mu6_bi02<-list()

for(j in 1:10){
  P_Phi_ev_Ne10000_mu5_bi00[[j]]<-c(eq_Ne10000_mu5_bi00[[j]][[1]][3])
  P_Phi_ev_Ne1000_mu4_bi01[[j]]<-c(eq_Ne1000_mu4_bi01[[j]][[1]][3])
  P_Phi_ev_Ne10000_mu5_bi01[[j]]<-c(eq_Ne10000_mu5_bi01[[j]][[1]][3])
  P_Phi_ev_Ne100000_mu6_bi01[[j]]<-c(eq_Ne100000_mu6_bi01[[j]][[1]][3])
  P_Phi_ev_Ne1000_mu4_bi02[[j]]<-c(eq_Ne1000_mu4_bi02[[j]][[1]][3])
  P_Phi_ev_Ne10000_mu5_bi02[[j]]<-c(eq_Ne10000_mu5_bi02[[j]][[1]][3])
  P_Phi_ev_Ne100000_mu6_bi02[[j]]<-c(eq_Ne100000_mu6_bi02[[j]][[1]][3])
  
  kcat_Ne10000_mu5_bi00[[j]]<-c(eq_Ne10000_mu5_bi00[[j]][[1]][1])
  kcat_Ne1000_mu4_bi01[[j]]<-c(eq_Ne1000_mu4_bi01[[j]][[1]][1])
  kcat_Ne10000_mu5_bi01[[j]]<-c(eq_Ne10000_mu5_bi01[[j]][[1]][1])
  kcat_Ne100000_mu6_bi01[[j]]<-c(eq_Ne100000_mu6_bi01[[j]][[1]][1])
  kcat_Ne1000_mu4_bi02[[j]]<-c(eq_Ne1000_mu4_bi02[[j]][[1]][1])
  kcat_Ne10000_mu5_bi02[[j]]<-c(eq_Ne10000_mu5_bi02[[j]][[1]][1])
  kcat_Ne100000_mu6_bi02[[j]]<-c(eq_Ne100000_mu6_bi02[[j]][[1]][1])
  
  kf_Ne10000_mu5_bi00[[j]]<-c(eq_Ne10000_mu5_bi00[[j]][[1]][2])
  kf_Ne1000_mu4_bi01[[j]]<-c(eq_Ne1000_mu4_bi01[[j]][[1]][2])
  kf_Ne10000_mu5_bi01[[j]]<-c(eq_Ne10000_mu5_bi01[[j]][[1]][2])
  kf_Ne100000_mu6_bi01[[j]]<-c(eq_Ne100000_mu6_bi01[[j]][[1]][2])
  kf_Ne1000_mu4_bi02[[j]]<-c(eq_Ne1000_mu4_bi02[[j]][[1]][2])
  kf_Ne10000_mu5_bi02[[j]]<-c(eq_Ne10000_mu5_bi02[[j]][[1]][2])
  kf_Ne100000_mu6_bi02[[j]]<-c(eq_Ne100000_mu6_bi02[[j]][[1]][2])
  
  for (i in 2:30){
    #print(P_Phi_ev)
    P_Phi_ev_Ne10000_mu5_bi00[[j]][i]<-eq_Ne10000_mu5_bi00[[j]][[i]][3]
    P_Phi_ev_Ne1000_mu4_bi01[[j]][i]<-eq_Ne1000_mu4_bi01[[j]][[i]][3]
    P_Phi_ev_Ne10000_mu5_bi01[[j]][i]<-eq_Ne10000_mu5_bi01[[j]][[i]][3]
    P_Phi_ev_Ne100000_mu6_bi01[[j]][i]<-eq_Ne100000_mu6_bi01[[j]][[i]][3]
    P_Phi_ev_Ne1000_mu4_bi02[[j]][i]<-eq_Ne1000_mu4_bi02[[j]][[i]][3]
    P_Phi_ev_Ne10000_mu5_bi02[[j]][i]<-eq_Ne10000_mu5_bi02[[j]][[i]][3]
    P_Phi_ev_Ne100000_mu6_bi02[[j]][i]<-eq_Ne100000_mu6_bi02[[j]][[i]][3]
    
    kcat_Ne10000_mu5_bi00[[j]][i]<-eq_Ne10000_mu5_bi00[[j]][[i]][1]
    kcat_Ne1000_mu4_bi01[[j]][i]<-eq_Ne1000_mu4_bi01[[j]][[i]][1]
    kcat_Ne10000_mu5_bi01[[j]][i]<-eq_Ne10000_mu5_bi01[[j]][[i]][1]
    kcat_Ne100000_mu6_bi01[[j]][i]<-eq_Ne100000_mu6_bi01[[j]][[i]][1]
    kcat_Ne1000_mu4_bi02[[j]][i]<-eq_Ne1000_mu4_bi02[[j]][[i]][1]
    kcat_Ne10000_mu5_bi02[[j]][i]<-eq_Ne10000_mu5_bi02[[j]][[i]][1]
    kcat_Ne100000_mu6_bi02[[j]][i]<-eq_Ne100000_mu6_bi02[[j]][[i]][1]
    
    kf_Ne10000_mu5_bi00[[j]][i]<-eq_Ne10000_mu5_bi00[[j]][[i]][2]
    kf_Ne1000_mu4_bi01[[j]][i]<-eq_Ne1000_mu4_bi01[[j]][[i]][2]
    kf_Ne10000_mu5_bi01[[j]][i]<-eq_Ne10000_mu5_bi01[[j]][[i]][2]
    kf_Ne100000_mu6_bi01[[j]][i]<-eq_Ne100000_mu6_bi01[[j]][[i]][2]
    kf_Ne1000_mu4_bi02[[j]][i]<-eq_Ne1000_mu4_bi02[[j]][[i]][2]
    kf_Ne10000_mu5_bi02[[j]][i]<-eq_Ne10000_mu5_bi02[[j]][[i]][2]
    kf_Ne100000_mu6_bi02[[j]][i]<-eq_Ne100000_mu6_bi02[[j]][[i]][2]
  }
  #tab_P_Phi_ev[[j]]<-P_Phi_ev
}

Fit_ev_Ne10000_mu5_bi00<-c()
Fit_ev_Ne1000_mu4_bi01<-c()
Fit_ev_Ne10000_mu5_bi01<-c()
Fit_ev_Ne100000_mu6_bi01<-c()
Fit_ev_Ne1000_mu4_bi02<-c()
Fit_ev_Ne10000_mu5_bi02<-c()
Fit_ev_Ne100000_mu6_bi02<-c()

kcatKM_Ne10000_mu5_bi00<-c()
kcatKM_Ne1000_mu4_bi01<-c()
kcatKM_Ne1000_mu4_bi02<-c()
kcatKM_Ne10000_mu5_bi01<-c()
kcatKM_Ne10000_mu5_bi02<-c()
kcatKM_Ne100000_mu6_bi01<-c()
kcatKM_Ne100000_mu6_bi02<-c()

for (i in 1:30){
  Fit_ev_Ne10000_mu5_bi00[i]=P_Phi_ev_Ne10000_mu5_bi00[[10]][i]/Max_val
  Fit_ev_Ne1000_mu4_bi01[i]=P_Phi_ev_Ne1000_mu4_bi01[[10]][i]/Max_val
  Fit_ev_Ne10000_mu5_bi01[i]=P_Phi_ev_Ne10000_mu5_bi01[[10]][i]/Max_val
  Fit_ev_Ne100000_mu6_bi01[i]=P_Phi_ev_Ne100000_mu6_bi01[[10]][i]/Max_val
  Fit_ev_Ne1000_mu4_bi02[i]=P_Phi_ev_Ne1000_mu4_bi02[[10]][i]/Max_val
  Fit_ev_Ne10000_mu5_bi02[i]=P_Phi_ev_Ne10000_mu5_bi02[[10]][i]/Max_val
  Fit_ev_Ne100000_mu6_bi02[i]=P_Phi_ev_Ne100000_mu6_bi02[[10]][i]/Max_val
  
  kcatKM_Ne10000_mu5_bi00[i]=kf_Ne10000_mu5_bi00[[10]][i]*kcat_Ne10000_mu5_bi00[[10]][i]/(kr_i+kcat_Ne10000_mu5_bi00[[10]][i])
  kcatKM_Ne1000_mu4_bi01[i]=kf_Ne1000_mu4_bi01[[10]][i]*kcat_Ne1000_mu4_bi01[[10]][i]/(kr_i+kcat_Ne1000_mu4_bi01[[10]][i])
  kcatKM_Ne1000_mu4_bi02[i]=kf_Ne1000_mu4_bi02[[10]][i]*kcat_Ne1000_mu4_bi02[[10]][i]/(kr_i+kcat_Ne1000_mu4_bi02[[10]][i])
  kcatKM_Ne10000_mu5_bi01[i]=kf_Ne10000_mu5_bi01[[10]][i]*kcat_Ne10000_mu5_bi01[[10]][i]/(kr_i+kcat_Ne10000_mu5_bi01[[10]][i])
  kcatKM_Ne10000_mu5_bi02[i]=kf_Ne10000_mu5_bi02[[10]][i]*kcat_Ne10000_mu5_bi02[[10]][i]/(kr_i+kcat_Ne10000_mu5_bi02[[10]][i])
  kcatKM_Ne100000_mu6_bi01[i]=kf_Ne100000_mu6_bi01[[10]][i]*kcat_Ne100000_mu6_bi01[[10]][i]/(kr_i+kcat_Ne100000_mu6_bi01[[10]][i])
  kcatKM_Ne100000_mu6_bi02[i]=kf_Ne100000_mu6_bi02[[10]][i]*kcat_Ne100000_mu6_bi02[[10]][i]/(kr_i+kcat_Ne100000_mu6_bi02[[10]][i])
}

Results<-list()
Results[["Ne3,bi01"]]<-data.frame(cbind(log10(kcatKM_Ne1000_mu4_bi01),log10(kcat_Ne1000_mu4_bi01[[10]])))
colnames(Results[["Ne3,bi01"]])<-c("kcat/KM","kcat")
Results[["Ne5,bi01"]]<-data.frame(cbind(log10(kcatKM_Ne100000_mu6_bi01),log10(kcat_Ne100000_mu6_bi01[[10]])))
colnames(Results[["Ne5,bi01"]])<-c("kcat/KM","kcat")

Results[["Ne3,bi02"]]<-data.frame(cbind(log10(kcatKM_Ne1000_mu4_bi02),log10(kcat_Ne1000_mu4_bi02[[10]])))
colnames(Results[["Ne3,bi02"]])<-c("kcat/KM","kcat")
Results[["Ne5,bi02"]]<-data.frame(cbind(log10(kcatKM_Ne100000_mu6_bi02),log10(kcat_Ne100000_mu6_bi02[[10]])))
colnames(Results[["Ne5,bi02"]])<-c("kcat/KM","kcat")

save(Results,file="evo_model.rda")

moy_P_Phi_ev_Ne10000_mu5_bi00<-c()
moy_P_Phi_ev_Ne1000_mu4_bi01<-c()
moy_P_Phi_ev_Ne10000_mu5_bi01<-c()
moy_P_Phi_ev_Ne100000_mu6_bi01<-c()
moy_P_Phi_ev_Ne1000_mu4_bi02<-c()
moy_P_Phi_ev_Ne10000_mu5_bi02<-c()
moy_P_Phi_ev_Ne100000_mu6_bi02<-c()

moy_kcat_Ne10000_mu5_bi00<-c()
moy_kcat_Ne1000_mu4_bi01<-c()
moy_kcat_Ne10000_mu5_bi01<-c()
moy_kcat_Ne100000_mu6_bi01<-c()
moy_kcat_Ne1000_mu4_bi02<-c()
moy_kcat_Ne10000_mu5_bi02<-c()
moy_kcat_Ne100000_mu6_bi02<-c()

moy_kf_Ne10000_mu5_bi00<-c()
moy_kf_Ne1000_mu4_bi01<-c()
moy_kf_Ne10000_mu5_bi01<-c()
moy_kf_Ne100000_mu6_bi01<-c()
moy_kf_Ne1000_mu4_bi02<-c()
moy_kf_Ne10000_mu5_bi02<-c()
moy_kf_Ne100000_mu6_bi02<-c()

#ec_m_slow<-c()
#ec_m_quick<-c()

for (i in 1:10){
  moy_P_Phi_ev_Ne10000_mu5_bi00[i]=mean(P_Phi_ev_Ne10000_mu5_bi00[[i]])
  moy_P_Phi_ev_Ne1000_mu4_bi01[i]=mean(P_Phi_ev_Ne1000_mu4_bi01[[i]])
  moy_P_Phi_ev_Ne10000_mu5_bi01[i]=mean(P_Phi_ev_Ne10000_mu5_bi01[[i]])
  moy_P_Phi_ev_Ne100000_mu6_bi01[i]=mean(P_Phi_ev_Ne100000_mu6_bi01[[i]])
  moy_P_Phi_ev_Ne1000_mu4_bi02[i]=mean(P_Phi_ev_Ne1000_mu4_bi02[[i]])
  moy_P_Phi_ev_Ne10000_mu5_bi02[i]=mean(P_Phi_ev_Ne10000_mu5_bi02[[i]])
  moy_P_Phi_ev_Ne100000_mu6_bi02[i]=mean(P_Phi_ev_Ne100000_mu6_bi02[[i]])
  #ec_m_slow[i]<-sd(P_Phi_ev_m_slow[[i]])
  moy_kcat_Ne10000_mu5_bi00[i]=mean(kcat_Ne10000_mu5_bi00[[i]])
  moy_kcat_Ne1000_mu4_bi01[i]=mean(kcat_Ne1000_mu4_bi01[[i]])
  moy_kcat_Ne10000_mu5_bi01[i]=mean(kcat_Ne10000_mu5_bi01[[i]])
  moy_kcat_Ne100000_mu6_bi01[i]=mean(kcat_Ne100000_mu6_bi01[[i]])
  moy_kcat_Ne1000_mu4_bi02[i]=mean(kcat_Ne1000_mu4_bi02[[i]])
  moy_kcat_Ne10000_mu5_bi02[i]=mean(kcat_Ne10000_mu5_bi02[[i]])
  moy_kcat_Ne100000_mu6_bi02[i]=mean(kcat_Ne100000_mu6_bi02[[i]])
  #ec_m_quick[i]<-sd(P_Phi_ev_m_quick[[i]])
  moy_kf_Ne10000_mu5_bi00[i]=mean(kf_Ne10000_mu5_bi00[[i]])
  moy_kf_Ne1000_mu4_bi01[i]=mean(kf_Ne1000_mu4_bi01[[i]])
  moy_kf_Ne10000_mu5_bi01[i]=mean(kf_Ne10000_mu5_bi01[[i]])
  moy_kf_Ne100000_mu6_bi01[i]=mean(kf_Ne100000_mu6_bi01[[i]])
  moy_kf_Ne1000_mu4_bi02[i]=mean(kf_Ne1000_mu4_bi02[[i]])
  moy_kf_Ne10000_mu5_bi02[i]=mean(kf_Ne10000_mu5_bi02[[i]])
  moy_kf_Ne100000_mu6_bi02[i]=mean(kf_Ne100000_mu6_bi02[[i]])
}



par(mfrow=c(2,2))

##Ne=10^4,bias=0 t=2*10^2
plot(moy_P_Phi_ev_Ne10000_mu5_bi00,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="Average flux(M/s)",main="Ne=10000,bias=0")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne10000_mu5_bi00),pch=19,col="green",xlab="Time (2e+02Ne)",ylab="Average kcat",main="Ne=10000,bias=0")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne10000_mu5_bi00),pch=19,col="green",xlab="Time (2e+02Ne)",ylab="Average kf",main="Ne=10000,bias=0")
plot(log10(moy_kcat_Ne10000_mu5_bi00)~log10(moy_kf_Ne10000_mu5_bi00),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=10000,bias=0")

##Ne=10^3,bias=-0.1
plot(moy_P_Phi_ev_Ne1000_mu4_bi01,pch=19,col="red",xlab="Time (5e+01Ne)",ylab="Average flux(M/s)",main="Ne=1000,bias=-0.1")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne1000_mu4_bi01),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kcat",main="Ne=1000,bias=-0.1")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne1000_mu4_bi01),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kf",main="Ne=1000,bias=-0.1")
plot(log10(moy_kcat_Ne1000_mu4_bi01)~log10(moy_kf_Ne1000_mu4_bi01),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=1000,bias=-0.1")

##Ne=10^4,bias=-0.1
plot(moy_P_Phi_ev_Ne10000_mu5_bi01,pch=19,col="red",xlab="Time (5e+01Ne)",ylab="Average flux(M/s)",main="Ne=10000,bias=-0.1")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne10000_mu5_bi01),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kcat",main="Ne=10000,bias=-0.1")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne10000_mu5_bi01),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kf",main="Ne=10000,bias=-0.1")
plot(log10(moy_kcat_Ne10000_mu5_bi01)~log10(moy_kf_Ne10000_mu5_bi01),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=10000,bias=-0.1")

##Ne=10^5,bias=-0.1
plot(moy_P_Phi_ev_Ne100000_mu6_bi01,pch=19,col="red",xlab="Time (5e+01Ne)",ylab="Average flux(M/s)",main="Ne=100000,bias=-0.1")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne100000_mu6_bi01),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kcat",main="Ne=100000,bias=-0.1")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne100000_mu6_bi01),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kf",main="Ne=100000,bias=-0.1")
plot(log10(moy_kcat_Ne100000_mu6_bi01)~log10(moy_kf_Ne100000_mu6_bi01),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=100000,bias=-0.1")

##Ne=10^3,bias=-0.2
plot(moy_P_Phi_ev_Ne1000_mu4_bi02,pch=19,col="red",xlab="Time (5e+01Ne)",ylab="Average flux(M/s)",main="Ne=1000,bias=-0.2")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne1000_mu4_bi02),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kcat",main="Ne=1000,bias=-0.2")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne1000_mu4_bi02),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kf",main="Ne=1000,bias=-0.2")
plot(log10(moy_kcat_Ne1000_mu4_bi02)~log10(moy_kf_Ne1000_mu4_bi02),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=1000,bias=-0.2")

##Ne=10^4,bias=-0.2
plot(moy_P_Phi_ev_Ne10000_mu5_bi02,pch=19,col="red",xlab="Time (5e+01Ne)",ylab="Average flux(M/s)",main="Ne=10000,bias=-0.2")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne10000_mu5_bi02),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kcat",main="Ne=10000,bias=-0.2")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne10000_mu5_bi02),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kf",main="Ne=10000,bias=-0.2")
plot(log10(moy_kcat_Ne10000_mu5_bi02)~log10(moy_kf_Ne10000_mu5_bi02),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=100000,bias=-0.2")

##Ne=10^5,bias=-0.2
plot(moy_P_Phi_ev_Ne100000_mu6_bi02,pch=19,col="red",xlab="Time (5e+01Ne)",ylab="Average flux(M/s)",main="Ne=100000,bias=-0.2")
#plot(ec_m_slow,pch=19,col="red",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-05")
plot(log10(moy_kcat_Ne100000_mu6_bi02),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kcat",main="Ne=100000,bias=-0.2")
#plot(ec_m_quick,pch=19,col="green",xlab="Time (2e+02Ne)",ylab="St. deviation flux(M/s)",main="Ne=1000,mu=1e-04")
plot(log10(moy_kf_Ne100000_mu6_bi02),pch=19,col="green",xlab="Time (5e+01Ne)",ylab="Average kf",main="Ne=100000,bias=-0.2")
plot(log10(moy_kcat_Ne100000_mu6_bi02)~log10(moy_kf_Ne100000_mu6_bi02),pch=19,col="blue",
     xlab="Average kf",ylab="Average kcat",main="Ne=100000,bias=-0.2")


##Ne=10^4,bias=0
par(mfrow=c(1,2),oma = c(0, 0, 0, 0))
plot(log10(kcat_Ne10000_mu5_bi00[[10]])~log10(kf_Ne10000_mu5_bi00[[10]]),
     col="blue",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne10000_mu5_bi00),nclass=6,density=50,cex.sub=0.5,
     col="blue",freq=FALSE,xlim=c(-8,-4),xlab="Fitness spread
from the maximal value (no bias)",main="",cex=0.5
)

##Ne=10^3,bias=-0.1
par(mfrow=c(2,2),oma = c(1, 0, 3, 0))
plot(log10(kcat_Ne1000_mu4_bi01[[10]])~log10(kf_Ne1000_mu4_bi01[[10]]),
     col="blue",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne1000_mu4_bi01),nclass=6,density=50,cex.sub=0.5,
     col="blue",freq=FALSE,xlim=c(-5,-2),xlab="Fitness spread
from the maximal value",main="",cex=0.5
)
mtext("Mutational bias = -0.1", outer = TRUE, cex = 1,font=1,side=1,line=-31)
##Ne=10^3,bias=-0.2
plot(log10(kcat_Ne1000_mu4_bi02[[10]])~log10(kf_Ne1000_mu4_bi02[[10]]),
     col="red",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne1000_mu4_bi02),nclass=6,density=50,cex.sub=0.5,
     col="red",freq=FALSE,xlim=c(-5,-2),xlab="Fitness spread
from the maximal value",main="",cex=0.5
)
mtext("Mutational bias = -0.2", outer = TRUE, cex = 1,font=1,side=1,line=-15)
mtext("Simulations for Ne=1000", outer = TRUE, cex = 1.25,font=2)

##Ne=10^4,bias=-0.1
par(mfrow=c(2,2),oma = c(1, 0, 3, 0))
plot(log10(kcat_Ne10000_mu5_bi01[[10]])~log10(kf_Ne10000_mu5_bi01[[10]]),
     col="blue",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne10000_mu5_bi01),nclass=6,density=50,cex.sub=0.5,
     col="blue",freq=FALSE,xlim=c(-7,-2),xlab="Fitness spread
from the maximal value",main="",cex=0.5
)
mtext("Mutational bias = -0.1", outer = TRUE, cex = 1,font=1,side=1,line=-31)
##Ne=10^4,bias=-0.2
plot(log10(kcat_Ne10000_mu5_bi02[[10]])~log10(kf_Ne10000_mu5_bi02[[10]]),
     col="red",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne10000_mu5_bi02),nclass=6,density=50,cex.sub=0.5,
     col="red",freq=FALSE,xlim=c(-7,-2),xlab="Fitness spread
from the maximal value",main="",cex=0.5
)
mtext("Mutational bias = -0.2", outer = TRUE, cex = 1,font=1,side=1,line=-15)
mtext("Simulations for Ne=10000", outer = TRUE, cex = 1.25,font=2)

##Ne=10^5,bias=-0.1
par(mfrow=c(2,2),oma = c(1, 0, 3, 0))
plot(log10(kcat_Ne100000_mu6_bi01[[10]])~log10(kf_Ne100000_mu6_bi01[[10]]),
     col="blue",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne100000_mu6_bi01),nclass=6,density=50,cex.sub=0.5,
     col="blue",freq=FALSE,xlim=c(-7,-2),xlab="Fitness spread
from the maximal value",main="",cex=0.5
)
mtext("Mutational bias = -0.1", outer = TRUE, cex = 1,font=1,side=1,line=-31)
##Ne=10^5,bias=-0.2
plot(log10(kcat_Ne100000_mu6_bi02[[10]])~log10(kf_Ne100000_mu6_bi02[[10]]),
     col="red",xlab="log10(kf)",ylab="log10(kcat)",pch=18)
hist(log10(1-Fit_ev_Ne100000_mu6_bi02),nclass=6,density=50,cex.sub=0.5,
     col="red",freq=FALSE,xlim=c(-7,-2),xlab="Fitness spread
from the maximal value",main="",cex=0.5
)
mtext("Mutational bias = -0.2", outer = TRUE, cex = 1,font=1,side=1,line=-15)
mtext("Simulations for Ne=100000", outer = TRUE, cex = 1.25,font=2)

#Fitness comparison between conditions
par(mfrow=c(2,1))
boxplot(log10(1-Fit_ev_Ne1000_mu4_bi01),log10(1-Fit_ev_Ne10000_mu5_bi01),log10(1-Fit_ev_Ne100000_mu6_bi01),
        xlab="",col="blue")
axis(1, at=1:3, labels=c("1e+03","1e+04","1e+05"))
mtext("Mutational bias = -0.1", outer = TRUE, cex = 1,font=1,side=1,line=-27)
boxplot(log10(1-Fit_ev_Ne1000_mu4_bi02),log10(1-Fit_ev_Ne10000_mu5_bi02),log10(1-Fit_ev_Ne100000_mu6_bi02),
        xlab="",col="red")
axis(1, at=1:3, labels=c("1e+03","1e+04","1e+05"))
mtext("Mutational bias = -0.2", outer = TRUE, cex = 1,font=1,side=1,line=-13)
mtext("Fitness spread vs. Ne
for different mutational bias", outer = TRUE, cex = 1.25,font=2)

par(mfrow=c(1,1))
plot(log10(kcat_Ne10000_mu5_bi00[[10]])~log10(kcatKM_Ne10000_mu5_bi00),
     col="green",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Enzyme activity when no mutational bias
for Ne=1e+04")

par(mfrow=c(2,3))

plot(log10(kcat_Ne1000_mu4_bi01[[10]])~log10(kcatKM_Ne1000_mu4_bi01),
     col="orange",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Ne=1e+03")

plot(log10(kcat_Ne10000_mu5_bi01[[10]])~log10(kcatKM_Ne10000_mu5_bi01),
     col="orange",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Ne=1e+04")

plot(log10(kcat_Ne100000_mu6_bi01[[10]])~log10(kcatKM_Ne100000_mu6_bi01),
     col="orange",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Ne=1e+05")

mtext("Mutational bias = -0.1", outer = TRUE, cex = 1,font=1,side=1,line=-36)

plot(log10(kcat_Ne1000_mu4_bi02[[10]])~log10(kcatKM_Ne1000_mu4_bi02),
     col="red",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Ne=1e+03")

plot(log10(kcat_Ne10000_mu5_bi02[[10]])~log10(kcatKM_Ne10000_mu5_bi02),
     col="red",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Ne=1e+04")

plot(log10(kcat_Ne100000_mu6_bi02[[10]])~log10(kcatKM_Ne100000_mu6_bi02),
     col="red",xlab="log10(kcat/KM)",ylab="log10(kcat)",pch=18,main="Ne=1e+05")

mtext("Mutational bias = -0.2", outer = TRUE, cex = 1,font=1,side=1,line=-18)
mtext("Enzyme activity", outer = TRUE, cex = 1.25,font=2,line=1)

summary(lm(log10(kcat_Ne1000_mu4_bi01[[10]])~log10(kcatKM_Ne1000_mu4_bi01)))
summary(lm(log10(kcat_Ne10000_mu5_bi01[[10]])~log10(kcatKM_Ne10000_mu5_bi01)))


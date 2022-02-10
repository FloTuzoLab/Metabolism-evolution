library(MASS)
mut_corr=TRUE
simul=7
flux="H"

setwd(dir="~")
#setwd(dir="Desktop/ScriptsR_modeles/")
#source("AP1.Plotting_multiple_images.R") #To plot multiple gradient and color plots
#source("AP2.RaphsonNewton.R") #To find equilibrium for 2 variables problems
options(digits=15)
setwd(dir="Desktop/enzyme-evolution/Evo_Results")

Ne_set=c(10^2,10^3,10^4,10^5)
N_rep=30
mut_bias_set<-c(-0.1,-0.2,"increasing")#<-c(-0.1)
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

kcat_i<-10^-2
kf_i<-10^2
kr_i<-10^3
E_tot_conc=10^-3

KT_set=c(10^-5,10^-1)
Vm_set=c(10^-6,10^-3)

if(flux=="W"){
  KT=KT_set[1]
  Vm=Vm_set[1]
}
if(flux=="H"){
  KT=KT_set[2]
  Vm=Vm_set[2]
}

Se=10*KT

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

Max_val<-fit(10^10,10^0,10^6,10^-3)
evo_gen_results<-list()
evo_results<-list()

#G<-5*10^2*Ne_set[1]#used for mutational bias
#for(l in 1:length(mut_bias_set)){
 # evo_results[[l]]<-list()
  #mut_bias<-mut_bias_set[l]
  #for (n in 1:length(Ne_set)){
#for (g in 1:length(G_set)){

mut_bias=as.character(param_set$bias[simul])
if(mut_bias=="increasing"){
  mut_bias_coeff=-0.05
}else{
  mut_bias<-as.numeric(mut_bias)
}

  Ne<-param_set$Ne_set[simul]
  G<-5*10^2*Ne##2.5*10^2#In principle 10^3
  mu<-10^-1/Ne
  Sigma_mut<-0.3
  Rho_mut=0
  if(mut_corr){
    Rho_mut=0.5
  }
  Sigma_mut_corr1 <- matrix(c(Sigma_mut^2,Sigma_mut^2,Sigma_mut^2,Sigma_mut^2),2,2)#VarCov Matrix for mutations affecitng enzyme efficiencies
  Sigma_mut_corr2 <- matrix(c(Sigma_mut^2,Rho_mut*Sigma_mut^2,Rho_mut*Sigma_mut^2,Sigma_mut^2),2,2)#VarCov Matrix for mutations affecitng enzyme efficiencies
  nam<-paste("eq_",flux,"_Ne",Ne,"_mu",-log10(mu),"_bi0",-mut_bias*10,sep="")
  #G<-G_set[g]
  #evo_results[[g]]<-list()
  for (g in 1:10){
    evo_results[[g]]<-list()
  }
  
    for(r in 1:N_rep){
      print(r)
      g=1
      #print(l)
    #print(r)
      
    #List of (Ne) cells
      Cell_tab<-list()
    #Cell_tab_curr<-list()
      Cell_tab[["kcat"]]<-c()
      Cell_tab[["kf"]]<-c()
      Cell_tab[["kr"]]<-c()
      Cell_tab[["fit"]]<-c()
      Cell_tab[["w"]]<-c()
      Cell_tab[["N"]]<-c()
    #Resident phenotype
      for (i in 1:1){
        Cell_tab[["kcat"]]<-c(Cell_tab[["kcat"]],kcat_i)
        Cell_tab[["kf"]]<-c(Cell_tab[["kf"]],kf_i)
        Cell_tab[["kr"]]<-c(Cell_tab[["kr"]],kr_i)
        Cell_tab[["fit"]]<-c(Cell_tab[["fit"]],1)
        Cell_tab[["w"]]<-c(Cell_tab[["w"]],1)
        Cell_tab[["N"]]<-c(Cell_tab[["N"]],Ne)
      }

      
      for (t in 1:G){
        ##Absolute fitness calculation
        for (i in 1:length(Cell_tab[["N"]])){
          Cell_tab[["fit"]][i]<-fit(Cell_tab[["kf"]][i],Cell_tab[["kr"]][i],Cell_tab[["kcat"]][i],E_tot_conc)
        }
        if(Cell_tab[["fit"]][i]>10^-3){
          print(paste(Cell_tab[["fit"]][i],Cell_tab[["kr"]][i],Cell_tab[["kf"]][i],Cell_tab[["kcat"]][i]))
        }
        #Mean absolute fitness
        w_m<-weighted.mean(Cell_tab[["fit"]],Cell_tab[["N"]])
      ##Relative fitness calculation
        for (i in 1:length(Cell_tab[["N"]])){
          Cell_tab[["w"]][i]<-Cell_tab[["fit"]][i]/w_m
        }
        #print(length(Cell_tab[["w"]]))
        ##Wright-Fisher sampling
        Cell_tab[["N"]]<-as.vector(rmultinom(1, Ne, prob=Cell_tab[["w"]]*Cell_tab[["N"]]))
        i=1
        while (i < length(Cell_tab[["N"]])+1){
          if(Cell_tab[["N"]][i]==0){
            Cell_tab[["kcat"]]<-Cell_tab[["kcat"]][-i]
            Cell_tab[["kf"]]<-Cell_tab[["kf"]][-i]
            Cell_tab[["kr"]]<-Cell_tab[["kr"]][-i]
            Cell_tab[["fit"]]<-Cell_tab[["fit"]][-i]
            Cell_tab[["w"]]<-Cell_tab[["w"]][-i]
            Cell_tab[["N"]]<-Cell_tab[["N"]][-i]
          }
          else{
            i=i+1
          }
        }
        mut_table<-c()
        for (i in 1:length(Cell_tab[["N"]])){
          mut_table<-c(mut_table,rbinom(1,Cell_tab[["N"]][i],mu))
        }
        pop_Size_prev<-length(Cell_tab[["N"]])
        for (i in 1:pop_Size_prev){
          if(mut_table[i]>0){
            Cell_tab[["N"]][i]=Cell_tab[["N"]][i]-mut_table[i]
            for (m in 1:mut_table[i]){
              if(mut_bias=="increasing"){
                mut_delta_keff<-mvrnorm(n = 1, c(mut_bias_coeff*log10(Cell_tab[["kcat"]][i]/10^-4),
                                                   mut_bias_coeff*log10(Cell_tab[["kr"]][i])), Sigma_mut_corr)
                #mut_delta_kcat<-rnorm(1,mean=-0.2+mut_bias_coeff*log10(Cell_tab[["kcat"]][i]/10^-4),sd=0.30)
                #mut_delta_kf<-rnorm(1,mean=-0.2+mut_bias_coeff*log10(Cell_tab[["kf"]][i]),sd=0.30)
                Cell_tab[["kcat"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kcat"]][i])+mut_delta_keff[1])
                Cell_tab[["kr"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kr"]][i])+mut_delta_keff[2])
              }
              else{
                mut_delta_1keff<-mvrnorm(n = 1, rep(mut_bias,2), Sigma_mut_corr1)
                mut_delta_2keff<-mvrnorm(n = 1, rep(mut_bias,2), Sigma_mut_corr2)
                Cell_tab[["kcat"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kcat"]][i])+mut_delta_2keff[2])
                Cell_tab[["kr"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kr"]][i])+mut_delta_1keff[1]-mut_delta_2keff[1])
                Cell_tab[["kf"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kf"]][i])+mut_delta_1keff[2])
                #print(mut_delta_1keff[1])
              }
              
              #if(Cell_tab[["kcat"]][length(Cell_tab[["N"]])]>10^6){
               # Cell_tab[["kcat"]][length(Cell_tab[["N"]])]<-10^6
              #}
              if(Cell_tab[["kf"]][length(Cell_tab[["N"]])]>10^10){
                Cell_tab[["kf"]][length(Cell_tab[["N"]])]<-10^10
              }
              if(Cell_tab[["kcat"]][length(Cell_tab[["N"]])]>10^6){
                Cell_tab[["kcat"]][length(Cell_tab[["N"]])]<-10^6
              }
              Cell_tab[["fit"]][length(Cell_tab[["N"]])+1]<-1
              Cell_tab[["w"]][length(Cell_tab[["N"]])+1]<-1
              Cell_tab[["N"]][length(Cell_tab[["N"]])+1]<-1
            }
          }
        }
        #if(t%%10^5==0){
         # print(paste(Cell_tab[["kcat"]],Cell_tab[["kf"]]))
          #print(t)
        #}
        if(t==G/10*g){
          print(paste(Cell_tab[["kcat"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],Cell_tab[["kf"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],Cell_tab[["kr"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],Cell_tab[["fit"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))]))
          evo_results[[g]][[r]]=c(Cell_tab[["kcat"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],
                                  Cell_tab[["kf"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],
                                  Cell_tab[["kr"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],
                                  Cell_tab[["fit"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))])
          print(g)
          g=g+1
        }
      }
      #evo_results[[l]][[r]]
    #print(m)
    }
#}
  #}
#}
  if(mut_bias=="increasing"){
    save(simul,evo_results,file=paste("eq_",flux,"_Ne",Ne,"_mu",-log10(mu),"_bi0",3,"3PMutCorr.rda",sep=""))
  }
  if(mut_bias!="increasing"){
    save(simul,evo_results,file=paste("eq_",flux,"_Ne",Ne,"_mu",-log10(mu),"_bi0",-mut_bias*10,"3PMutCorr.rda",sep=""))
  }



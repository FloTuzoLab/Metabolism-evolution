install.packages("Rmpfr")
require("Rmpfr")
#mpfr allows for more precision, needed for estimations of flux with low permeability

###TESTING MODELS

##0.Estimating maximal size of cell that can be sustained by passive diffusion
P<-10^-6
Se<-10^-1
R<-(P*Se*(4*pi)*6.02*10^8/(5*10^4)/(4/3*pi)^1.88)^(1/3.64)
V<-4/3*pi*R^3
#R around 0.1 to 1 micrometer

###1.Simulation of the dynamics to test the analytical resolution accuracy for facilitated diffusion

##a.Simulation procedure and results
#Particular conditions to test models
KT=0.005 #Molar, Michaelis-like constant for facilitated diffusion
Vm=10^-3 #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3
Se=0.001
kf=10^10
kr=10^-3
kcat=10^6
E_tot_conc=10^-3
dt=0.0000001#timestep in seconds
N_steps=500000##will further need to be replace by a stoping condition of the algorithm, when concentrations nearby equilibrium
## Initializing current concentrations - representing the initial concentrations
cES=0
cEf=E_tot_conc
cSf=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
##Initalization of the process
P=0
for (i in 1:N_steps){
  cSf_c<-cSf
  cEf_c<-cEf
  cES_c<-cES
  FD_dt=Vm*(Se-cSf_c)/(KT+(Se+cSf_c)+Se*cSf_c/KT) ##Facilitated diffusion
  CR1_dt=kf*cSf_c*cEf_c-(kr)*cES_c
  CR2_dt=kcat*cES_c
  cSf=cSf_c+(FD_dt-CR1_dt)*dt 
  cEf=cEf_c+(CR2_dt-CR1_dt)*dt
  cES=cES_c-(CR2_dt-CR1_dt)*dt
}
cP=kcat*cES
cEf_cES<-cEf/cES
#Results
#Optimized parameters and low concentrations
#cP=0.000166667
#With high values of parameters : KT=0.01,Vm=10^-3,Se=0.2,kf=10^10,kr=10^3,kcat=10^6,E_tot_conc=10^-3
#cP=0.0009523809 cEf/cES=1049999 with dt=0.0000001 and N_steps=500000
#With low values of parameters : KT=0.1,Vm=10^-6,Se=0.2,kf=10^2,kr=10^-2,kcat=10^-4,E_tot_conc=10^-6
#cP=9.994855e-11 cEf/cES=0.0005147541 with dt=0.1 and N_steps=30000000 (very slow, explains the sight difference)

##b.Analytical determination quasi-equilibrium and results

a<- (Vm*kf+kf*E_tot_conc*kcat*(1+Se/KT))
b<- (Vm*(kr+kcat-Se*kf)+kf*kcat*E_tot_conc*(KT+Se))
c<- -Vm*Se*(kr+kcat)
delta<-b^2-4*a*c
S_conc_eq=(-b+delta^(1/2))/(2*a)

ES_conc_eq=E_tot_conc*(kf*S_conc_eq)/(kr+kcat+kf*S_conc_eq)
P_Phi_eq=kcat*ES_conc_eq

#Results
#With high values of parameters
#P_Phi_eq=0.0009523809 (E_tot_conc-ES_conc_eq)/ES_conc_eq=1049999
#With low values of parameters
#P_Phi_eq=9.99495e-11 (E_tot_conc-ES_conc_eq)/ES_conc_eq=0.0005052272

###2.Simulation of the dynamics to test the analytical resolution accuracy for passive diffusion

##a.Simulation procedure and results
#Fixing the size of a cell
r_cell=2.5*10^-6#unit:m
SA<-4*pi*r_cell^2
V<-4/3*pi*r_cell^3

Se=0.1
kf=10^10
kr=10^6
kcat=10^3
Pd=10^-7#unit:m/s not posible to test with values too low because of numerical problems above 10^20
#Pd 5 orders of magnitude above glucose/fructose values
E_tot_conc=10^-3
dt=0.00000001 #timestep in seconds
N_steps=5000000##will further need to be replace by a stoping condition of the algorithm, when concentrations nearby equilibrium
## Initializing current concentrations - representing the initial concentrations
cES=0
cEf=E_tot_conc
cSf=0 ##No substrate at the beginning, as it "comes from" the outside of the cell.
##Initalization of the process
P=0

for (i in 1:N_steps){
  cSf_c<-cSf
  cEf_c<-cEf
  cES_c<-cES
  PD_dt=Pd*SA/V*(Se-cSf_c) ##Facilitated diffusion
  CR1_dt=kf*cSf_c*cEf_c-(kr)*cES_c
  CR2_dt=kcat*cES_c
  cSf=cSf_c+(PD_dt-CR1_dt)*dt 
  cEf=cEf_c+(CR2_dt-CR1_dt)*dt
  cES=cES_c-(CR2_dt-CR1_dt)*dt
}
cP=kcat*cES


#Results
#Comparing with facilitated diffusion using true values for cell permeability
#For high concentrations (0.2M)
#cP=2.4e-07 More than 2 orders of magnitude below facilitated diffusion
#For low concentrations (10^-3M)
#cP=1.2e-09 More than 5 orders of magnitude below facilitated diffusion
#With high values of parameters : r_cell=2.5*10^-6,Pd=10^-7,Se=0.1,kf=10^10,kr=10^3,kcat=10^6,E_tot_conc=10^-3
#cP=(0.01199999986) 0.01199985411(more accuracy in simulations) with dt=0.0000001,N_steps=500000
#With low values of parameters : r_cell=2.5*10^-6,Pd=10^-10,Se=0.1,kf=10^2,kr=10^-2,kcat=10^-4,E_tot_conc=10^-6
#cP=9.989910107e-11 with dt=0.1,N_steps=10000000

##b.Analytical determination quasi-equilibrium
a<--mpfr(Pd*SA/V*kf,120)
b<-mpfr(Pd*SA/V*(Se*kf-(kr+kcat))-kcat*kf*E_tot_conc,120)
c<-mpfr(Pd*SA/V*Se*(kr+kcat),120)
delta_pd<-mpfr(b^2-4*a*c,120)
S_conc_eq=mpfr((-b-delta_pd^(1/2))/(2*a),120)
ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat+kf*S_conc_eq),120)
P_Phi_eq=mpfr(kcat*ES_conc_eq,120)
P_Phi_eq

#Results
#For low concentrations (10^-3M)
#cP=1.1999999998582581069570827355765e-9 More than 5 orders of magnitude below facilitated diffusion
#With high values of parameters
#P_Phi_eq=0.011999854107058456031037862183631 (slight difference, but better precision for analytical)
#With low values of parameters
#P_Phi_eq=9.989910107e-11 (no difference here)





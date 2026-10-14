# Run from Structural Roles of Asymptomatic Infection and a Downstream Bacterial State in a Cholera Model: Rscript simulations_paper1.R
library(deSolve)
dir.create('figures',showWarnings=FALSE); dir.create('tables',showWarnings=FALSE)
dir.create('validation',showWarnings=FALSE)
pars_base <- c(
  Lambda = 20,
  beta = 0.45,
  K = 1e5,
  theta = 0.25,
  mu_S = 0.00004,
  mu_A = 0.00003,
  mu_R = 0.00004,
  delta = 0.005,
  gamma_A = 1 / 5,
  gamma_S = 1 / 7,
  xi_A = 0.1,
  xi_S = 0.3,
  eta = 0.2,
  rho = 0.05,
  eta_nv = 0.15
)
p <- pars_base[!names(pars_base)%in%c('mu_S','mu_A','mu_R')]
p <- c(p,mu=unname(pars_base['mu_S']))
R0 <- function(p) with(as.list(p),beta*Lambda/(mu*K*(eta+rho))*((1-theta)*xi_A/(mu+gamma_A)+theta*xi_S/(mu+delta+gamma_S)))
rhs <- function(t,y,p) with(as.list(c(y,p)),{
 inc <- beta*S*Bv/(K+Bv)
 list(c(Lambda-inc-mu*S,(1-theta)*inc-(mu+gamma_A)*IA,
 theta*inc-(mu+delta+gamma_S)*IS,gamma_A*IA+gamma_S*IS-mu*R,
 xi_A*IA+xi_S*IS-(eta+rho)*Bv,rho*Bv-eta_nv*Bnv))
})
endemic <- function(p) with(as.list(p),{
 aA<-mu+gamma_A;aS<-mu+delta+gamma_S;cB<-eta+rho
 Q<-(1-theta)*xi_A/aA+theta*xi_S/aS
 B<-mu*K*(R0(p)-1)/(beta+mu)
 S<-(cB*K+Q*Lambda)/(Q*(beta+mu))
 IA<-(1-theta)*cB*B/(aA*Q);IS<-theta*cB*B/(aS*Q)
 c(S=S,IA=IA,IS=IS,R=(gamma_A*IA+gamma_S*IS)/mu,Bv=B,Bnv=rho*B/eta_nv)
})
y0<-c(S=unname(p['Lambda']/p['mu']-10),IA=5,IS=5,R=0,Bv=100,Bnv=0)
times<-0:1095
run<-function(par,rt=1e-9,at=1e-10) as.data.frame(ode(y0,times,rhs,par,method='lsoda',rtol=rt,atol=at))
baseline<-run(p)
low<-p;low['beta']<-p['beta']*0.8/R0(p)
extinction<-run(low)
oldcontrol<-p;oldcontrol['beta']<-0.12;oldcontrol['eta']<-0.45
write.csv(data.frame(parameter=names(p),value=unname(p)), 'tables/parameters.csv',row.names=FALSE)
write.csv(baseline,'tables/baseline_trajectories.csv',row.names=FALSE)
write.csv(extinction,'tables/subthreshold_trajectories.csv',row.names=FALSE)
E<-endemic(p)
stopifnot(max(abs(rhs(0,E,p)[[1]]))<1e-8)
checks<-c(equilibrium_residual=max(abs(rhs(0,E,p)[[1]])),minimum_state=min(as.matrix(baseline[,-1])),R0_baseline=R0(p),R0_subthreshold=R0(low),R0_previous_control_under_name2=R0(oldcontrol))
# Check the NGM explicitly against the analytic expression.
aA<-p['mu']+p['gamma_A'];aS<-p['mu']+p['delta']+p['gamma_S'];cB<-p['eta']+p['rho'];S0<-p['Lambda']/p['mu'];q<-p['beta']*S0/p['K']
F<-matrix(c(0,0,(1-p['theta'])*q,0,0,p['theta']*q,0,0,0),3,byrow=TRUE)
V<-matrix(c(aA,0,0,0,aS,0,-p['xi_A'],-p['xi_S'],cB),3,byrow=TRUE)
checks<-c(checks,NGM_error=abs(max(Mod(eigen(F%*%solve(V))$values))-R0(p)))
# Natural population balance at an arbitrary state.
checks<-c(checks,population_balance_error=abs(sum(rhs(0,y0,p)[[1]][1:4])-(p['Lambda']-p['mu']*sum(y0[1:4])-p['delta']*y0['IS'])))
strict<-run(p,1e-11,1e-12)
checks<-c(checks,tolerance_scaled_difference=max(abs(as.matrix(strict[,-1]-baseline[,-1]))/(1+abs(as.matrix(strict[,-1])))))
# Exact decoupling check: eta_nv changes no transmission equations.
alt<-p;alt['eta_nv']<-2*p['eta_nv'];alternate<-run(alt)
checks<-c(checks,downstream_scaled_difference=max(abs(as.matrix(alternate[,2:6]-baseline[,2:6]))/(1+abs(as.matrix(baseline[,2:6])))))
# Local numerical stability at the illustrated endemic equilibrium (not a global proof).
jac<-function(y,par){sapply(seq_along(y),function(k){h<-1e-5*max(1,abs(y[k]));yp<-ym<-y;yp[k]<-yp[k]+h;ym[k]<-ym[k]-h;(rhs(0,yp,par)[[1]]-rhs(0,ym,par)[[1]])/(2*h)})}
checks<-c(checks,endemic_max_real_eigenvalue=max(Re(eigen(jac(E,p))$values)))
# Compute sensitivity elasticties numerically, check key exact identities.
el<-sapply(names(p),function(n){pp<-pm<-p;pp[n]<-p[n]*(1+1e-5);pm[n]<-p[n]*(1-1e-5);(R0(pp)-R0(pm))/(2e-5*R0(p))})
stopifnot(abs(el['beta']-1)<1e-7,abs(el['K']+1)<1e-7,abs(el['eta_nv'])<1e-7,checks['NGM_error']<1e-8,min(as.matrix(baseline[,-1])) > -1e-7)
write.csv(data.frame(parameter=names(el),elasticity=unname(el)),'tables/elasticities.csv',row.names=FALSE)
write.csv(data.frame(check=names(checks),value=unname(checks)),'validation/numerical_checks.csv',row.names=FALSE)
write.csv(data.frame(compartment=names(E),equilibrium=unname(E)),'tables/endemic_equilibrium.csv',row.names=FALSE)
summary_one<-function(s,par,label){inc<-with(s,par['theta']*par['beta']*S*Bv/(par['K']+Bv));data.frame(scenario=label,R0=R0(par),peak_prevalence=max(s$IS),peak_day=s$time[which.max(s$IS)],symptomatic_person_days=sum(diff(s$time)*(head(s$IS,-1)+tail(s$IS,-1))/2),cumulative_symptomatic_incidence=sum(diff(s$time)*(head(inc,-1)+tail(inc,-1))/2),final_IS=tail(s$IS,1))}
summary<-rbind(summary_one(baseline,p,'baseline'),summary_one(extinction,low,'subthreshold'))
write.csv(summary,'tables/scenario_summary.csv',row.names=FALSE)
# Publication figures: population counts and environmental concentrations on separate axes.
pdf('figures/fig1_local_trajectories.pdf',width=8,height=6.5)
par(mfrow=c(3,2),mar=c(3.3,4,1.6,1),mgp=c(2.3,.7,0))
for(n in names(y0)){plot(baseline$time,baseline[[n]],type='l',col='#2166AC',lwd=1.8,xlab='Time (days)',ylab=if(n%in%c('Bv','Bnv')) 'Concentration (model units)' else 'Persons',main=n);abline(h=E[n],lty=3,col='#666666')}
dev.off()
pdf('figures/fig2_threshold.pdf',width=8,height=3.8)
par(mfrow=c(1,2),mar=c(4,4,2,1))
for(n in c('IS','Bv')){plot(baseline$time,pmax(baseline[[n]],1e-10),type='l',log='y',ylim=c(1e-8,max(baseline[[n]])*1.2),col='#2166AC',lwd=2,xlab='Time (days)',ylab=if(n=='IS') 'Symptomatic prevalence (persons)' else 'Concentration (model units)',main=n);lines(extinction$time,pmax(extinction[[n]],1e-10),col='#B2182B',lwd=2);legend('topright',c(sprintf('Baseline: R0 = %.3f',R0(p)),'Below threshold: R0 = 0.800'),col=c('#2166AC','#B2182B'),lty=1,bty='n',cex=.72)}
dev.off()
pdf('figures/fig3_elasticities.pdf',width=7,height=4.5)
par(mar=c(4,5,1,1));o<-order(el);barplot(el[o],names.arg=names(el)[o],horiz=TRUE,las=1,col=ifelse(el[o]>=0,'#2166AC','#B2182B'),xlab='Normalized sensitivity of R0',cex.names=.8);abline(v=0)
dev.off()
# Exact equilibrium branch: no time-series inference of asymptotic convergence.
grid<-seq(.5,2,length.out=151)
branch<-sapply(grid,function(r){pr<-p;pr['beta']<-p['beta']*r/R0(p);if(r<=1) 0 else endemic(pr)['IS']})
write.csv(data.frame(R0=grid,IS_equilibrium=branch),'tables/equilibrium_branch.csv',row.names=FALSE)
pdf('figures/fig4_equilibrium_branch.pdf',width=6.5,height=4)
plot(grid,branch,type='l',lwd=2,col='#2166AC',xlab='Basic reproduction number R0',ylab='Endemic symptomatic prevalence (persons)');abline(v=1,lty=3);segments(1,0,2,0,lty=2,col='#B2182B');legend('topleft',c('Nonnegative equilibrium branch','Unstable disease-free branch (R0 > 1)'),col=c('#2166AC','#B2182B'),lty=c(1,2),bty='n',cex=.8)
dev.off()
sink('validation/sessionInfo.txt');print(sessionInfo());sink()
# Numerical results are exported to files; keep the console quiet.

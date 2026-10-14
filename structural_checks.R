# Run from Structural Roles of Asymptomatic Infection and a Downstream Bacterial State in a Cholera Model: Rscript structural_checks_paper1.R
# Exact structural comparisons for the local model; no surveillance fitting.
source("simulations.R")
upstream <- c("S","IA","IS","R","Bv")
metric <- function(a,b) max(abs(as.matrix(a)-as.matrix(b))/(1+abs(as.matrix(b))))
integrate <- function(initial,fun,parameters) {
  as.data.frame(ode(initial,times,fun,parameters,method="lsoda",
                    rtol=1e-10,atol=1e-11))
}
reference <- integrate(y0,rhs,p)
rhs_reduced <- function(t,y,parameters) {
  full <- c(y,Bnv=0)
  list(rhs(t,full,parameters)[[1]][1:5])
}
reduced <- integrate(y0[upstream],rhs_reduced,p)
p_down <- p
p_down["eta_nv"] <- 2*p["eta_nv"]
changed_down <- integrate(y0,rhs,p_down)
initial_down <- y0
initial_down["Bnv"] <- 10000
changed_initial <- integrate(initial_down,rhs,p)
p_pair <- p
p_pair["eta"] <- p["eta"]+p["rho"]/2
p_pair["rho"] <- p["rho"]/2
same_sum <- integrate(y0,rhs,p_pair)
kappa <- 10
p_scale <- p
p_scale[c("K","xi_A","xi_S")] <- kappa*p[c("K","xi_A","xi_S")]
initial_scale <- y0
initial_scale[c("Bv","Bnv")] <- kappa*y0[c("Bv","Bnv")]
scaled <- integrate(initial_scale,rhs,p_scale)

# All-symptomatic limiting model, retaining the same infectious removal.
rhs_one <- function(t,y,parameters) with(as.list(c(y,parameters)),{
  inc <- beta*S*Bv/(K+Bv)
  list(c(Lambda-inc-mu*S,inc-(mu+delta+gamma_S)*IS,
         gamma_S*IS-mu*R,xi_S*IS-(eta+rho)*Bv))
})
initial_one <- c(S=unname(y0["S"]),IS=unname(y0["IA"]+y0["IS"]),
                 R=unname(y0["R"]),Bv=unname(y0["Bv"]))
one <- integrate(initial_one,rhs_one,p)
p_one <- p
p_one["theta"] <- 1
theta_grid <- seq(0,1,length.out=101)
threshold_grid <- data.frame(theta=theta_grid,
 R0=sapply(theta_grid,function(th){z<-p;z["theta"]<-th;R0(z)}))
write.csv(threshold_grid,"tables/hierarchy_thresholds.csv",row.names=FALSE)
write.csv(data.frame(model=c("H1_all_symptomatic","H2_split_reduced","H3_split_full"),
 R0=c(R0(p_one),R0(p),R0(p)),
 initial_total_infected=rep(sum(y0[c("IA","IS")]),3),
 infectious_loss=rep(sum(p[c("eta","rho")]),3)),
 "tables/hierarchy_summary.csv",row.names=FALSE)
outputs <- list(H1=one,H2=reduced,H3=reference,
 H3_changed_removal=changed_down,H3_changed_initial=changed_initial,
 same_total_removal=same_sum,concentration_scaled=scaled)
for(n in names(outputs)) write.csv(outputs[[n]],
 paste0("tables/structural_",n,".csv"),row.names=FALSE)
write.csv(data.frame(experiment=c("reference","same_total_removal","scaled_concentration"),
 eta=c(p["eta"],p_pair["eta"],p_scale["eta"]),
 rho=c(p["rho"],p_pair["rho"],p_scale["rho"]),
 eta_nv=c(p["eta_nv"],p_pair["eta_nv"],p_scale["eta_nv"]),
 K=c(p["K"],p_pair["K"],p_scale["K"]),concentration_scale=c(1,1,kappa)),
 "tables/identifiability_transformations.csv",row.names=FALSE)

# Explicit split next-generation decomposition.
Fsplit <- matrix(c(0,0,(1-p["theta"])*q,0,0,p["theta"]*q,
                   p["xi_A"],p["xi_S"],0),3,byrow=TRUE)
Vsplit <- diag(as.numeric(c(aA,aS,cB)))
split_radius <- max(Mod(eigen(Fsplit%*%solve(Vsplit))$values))
structural_checks <- c(
 reduced_upstream_difference=metric(reduced[upstream],reference[upstream]),
 downstream_removal_difference=metric(changed_down[upstream],reference[upstream]),
 downstream_initial_difference=metric(changed_initial[upstream],reference[upstream]),
 equal_sum_upstream_difference=metric(same_sum[upstream],reference[upstream]),
 concentration_scale_human_difference=metric(scaled[c("S","IA","IS","R")],
                                                reference[c("S","IA","IS","R")]),
 concentration_scale_environment_difference=metric(scaled[c("Bv","Bnv")]/kappa,
                                                   reference[c("Bv","Bnv")]),
 concentration_scale_R0_error=unname(abs(R0(p_scale)-R0(p))),
 split_convention_error=unname(abs(split_radius^2-R0(p))),
 linearization_decomposition_error=max(abs((Fsplit-Vsplit)-(F-V)))
)
stopifnot(all(is.finite(structural_checks)),max(structural_checks)<1e-6)

# Analytic Jacobian: no trajectory convergence inference.
analytic_jacobian <- function(y,parameters) with(as.list(c(y,parameters)),{
  gS <- beta*Bv/(K+Bv)
  gB <- beta*S*K/(K+Bv)^2
  matrix(c(-(mu+gS),0,0,0,-gB,0,
   (1-theta)*gS,-(mu+gamma_A),0,0,(1-theta)*gB,0,
   theta*gS,0,-(mu+delta+gamma_S),0,theta*gB,0,
   0,gamma_A,gamma_S,-mu,0,0,
   0,xi_A,xi_S,0,-(eta+rho),0,
   0,0,0,0,rho,-eta_nv),6,byrow=TRUE)
})
core <- c(1,2,3,5)
R0_grid <- sort(unique(c(1+10^seq(-6,-2,length.out=61),
                        seq(1.01,2,length.out=151))))
stability <- do.call(rbind,lapply(R0_grid,function(r){
 z<-p; z["beta"]<-p["beta"]*r/R0(p)
 J<-analytic_jacobian(endemic(z),z)
 ev_full<-eigen(J,only.values=TRUE)$values
 ev_core<-eigen(J[core,core],only.values=TRUE)$values
 data.frame(R0=r,beta=unname(z["beta"]),
            core_max_real=max(Re(ev_core)),full_max_real=max(Re(ev_full)))
}))
write.csv(stability,"tables/endemic_eigenvalues.csv",row.names=FALSE)
J <- analytic_jacobian(E,p)
structural_checks <- c(structural_checks,
 analytic_jacobian_difference=max(abs(J-jac(E,p)))/(1+max(abs(J))))
stopifnot(tail(structural_checks,1)<1e-7)
write.csv(data.frame(check=names(structural_checks),value=unname(structural_checks)),
 "validation/structural_checks.csv",row.names=FALSE)
# Record empirical sign checks rather than assuming a theorem over the range.
write.csv(data.frame(sampled_points=nrow(stability),
 min_R0=min(stability$R0),max_R0=max(stability$R0),
 sampled_core_stable=all(stability$core_max_real<0),
 sampled_full_stable=all(stability$full_max_real<0)),
 "validation/eigenvalue_sampling_summary.csv",row.names=FALSE)

pdf("figures/fig5_structural_hierarchy.pdf",width=8.5,height=7.2)
par(mfrow=c(2,2),mar=c(4.2,4.5,2.3,1),mgp=c(2.9,.75,0))
ix <- reference$time<=120
plot(reference$time[ix],reference$Bv[ix],type="l",lwd=2.4,col="#2166AC",
 xlab="Time (days)",ylab="Infectious concentration (model units)",
 main="(a) Reduced and full transmission")
lines(reduced$time[ix],reduced$Bv[ix],col="#B2182B",lty=2,lwd=1.8)
lines(changed_down$time[ix],changed_down$Bv[ix],col="#1B7837",lty=3,lwd=1.8)
lines(changed_initial$time[ix],changed_initial$Bv[ix],col="#762A83",lty=4,lwd=1.8)
legend("topright",c("H3 full","H2 reduced","H3: doubled eta_nv","H3: changed Bnv(0)"),
 col=c("#2166AC","#B2182B","#1B7837","#762A83"),lty=1:4,bty="n",cex=.68)
limits <- range(reference$Bnv[ix],changed_down$Bnv[ix],changed_initial$Bnv[ix])
plot(reference$time[ix],reference$Bnv[ix],type="l",ylim=limits,lwd=2,col="#2166AC",
 xlab="Time (days)",ylab="Downstream concentration (model units)",
 main="(b) Distinct downstream outputs")
lines(changed_down$time[ix],changed_down$Bnv[ix],col="#1B7837",lty=3,lwd=2)
lines(changed_initial$time[ix],changed_initial$Bnv[ix],col="#762A83",lty=4,lwd=2)
legend("topright",c("Reference H3","Doubled eta_nv","Changed Bnv(0)"),
 col=c("#2166AC","#1B7837","#762A83"),lty=c(1,3,4),bty="n",cex=.7)
plot(threshold_grid$theta,threshold_grid$R0,type="l",lwd=2,col="#2166AC",
 xlab="Symptomatic fraction theta",ylab="Human-generation R0",
 main="(c) Symptom allocation and threshold")
abline(h=R0(p_one),lty=2,col="#B2182B")
points(unname(p["theta"]),R0(p),pch=16,col="#2166AC")
legend("topright",c("H2 / H3","H1: all symptomatic","Illustrative baseline"),
 col=c("#2166AC","#B2182B","#2166AC"),lty=c(1,2,NA),
 pch=c(NA,NA,16),bty="n",cex=.7)
plot(reference$time[ix],reference$IS[ix],type="l",
 ylim=range(reference$IS[ix],one$IS[ix]),lwd=2,col="#2166AC",
 xlab="Time (days)",ylab="Symptomatic prevalence (persons)",
 main="(d) Matched initial total infections")
lines(one$time[ix],one$IS[ix],lwd=2,col="#B2182B",lty=2)
legend("topright",c("H2 / H3: split infection","H1: all symptomatic"),
 col=c("#2166AC","#B2182B"),lty=c(1,2),bty="n",cex=.7)
invisible(dev.off())

pdf("figures/fig6_endemic_eigenvalues.pdf",width=8.5,height=4)
par(mfrow=c(1,2),mar=c(4.5,5,2,1),mgp=c(3.1,.75,0))
near <- stability$R0<=1.01
draw_spectrum <- function(rows,near_threshold=FALSE) {
 x<-if(near_threshold) stability$R0[rows]-1 else stability$R0[rows]
 limits<-range(0,stability$core_max_real[rows],stability$full_max_real[rows])
 plot(x,stability$core_max_real[rows],type="l",lwd=2,col="#2166AC",
 log=if(near_threshold) "x" else "",ylim=limits,
 xlab=if(near_threshold) "Distance above threshold (R0 - 1)" else "Human-generation R0",
 ylab="Largest real eigenvalue (per day)",
 main=if(near_threshold) "(a) Near the threshold" else "(b) Sampled endemic branch")
 lines(x,stability$full_max_real[rows],col="#B2182B",lty=2,lwd=2)
 abline(h=0,col="#555555",lty=3)
 legend("bottomleft",c("Transmission core","Full system"),
 col=c("#2166AC","#B2182B"),lty=c(1,2),bty="n",cex=.75)
}
draw_spectrum(which(near),TRUE)
draw_spectrum(seq_len(nrow(stability)))
invisible(dev.off())

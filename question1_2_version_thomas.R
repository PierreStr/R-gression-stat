library(stats)
library(ggplot2)
library(pracma)

# Approximation de la constante d'Euler-Mascheroni par une somme de série
n.tilt <- 1000000
gamma_approx <- sum(1 / seq(1, n.tilt)) - log(n.tilt)

set.seed("311")
# Fonction pour initialiser les listes avec des valeurs aléatoires de distribution Gumbel
initialisation <- function(n) {
  unif.distrib <- runif(n, min = 0, max = 1)
  gumb.distrib.initial <- vector("numeric", length = n)
  for (i in 1:n){
    gumb.distrib.initial[i] <- -3 * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + 2
  }
  return(gumb.distrib.initial)
}


# Fonction de log-vraisemblance
log_vraissemblance <- function(params, L) {
  beta <- params[1]
  mu <- params[2]
  n <- length(L)
  result <- sum( (mu - L)/beta - exp((mu - L)/beta))
  result<-result-n*log(beta)
  return(-result)
}

# Initialisation des listes de différentes tailles
l_250 <- numeric(250)
l_500 <- numeric(500)
l_750 <- numeric(750)
l_1000 <- numeric(1000)
l_1250 <- numeric(1250)

# Nombre de simulations
M <- 100

# MSE pour beta
MSE_250 <- numeric(M)
MSE_500 <- numeric(M)
MSE_750 <- numeric(M)
MSE_1000 <- numeric(M)
MSE_1250 <- numeric(M)

#MSE de mu
MSEmu_250 <- numeric(M)
MSEmu_500 <- numeric(M)
MSEmu_750 <- numeric(M)
MSEmu_1000 <- numeric(M)
MSEmu_1250 <- numeric(M)


# Biais pour beta
B_250 <- numeric(M)
B_500 <- numeric(M)
B_750 <- numeric(M)
B_1000 <- numeric(M)
B_1250 <- numeric(M)

#biais de mu
Bmu_250 <- numeric(M)
Bmu_500 <- numeric(M)
Bmu_750 <- numeric(M)
Bmu_1000 <- numeric(M)
Bmu_1250 <- numeric(M)

# Variance de beta
V_250 <- numeric(M)
V_500 <- numeric(M)
V_750 <- numeric(M)
V_1000 <- numeric(M)
V_1250 <- numeric(M)


# Variance de mu
Vmu_250 <- numeric(M)
Vmu_500 <- numeric(M)
Vmu_750 <- numeric(M)
Vmu_1000 <- numeric(M)
Vmu_1250 <- numeric(M)

#valeur des estimateurs de beta vraissemblance

beta_hat_250<- numeric(M)
beta_hat_500<- numeric(M)
beta_hat_750<- numeric(M)
beta_hat_1000<- numeric(M)
beta_hat_1250<- numeric(M)

#valeur des estimateurs de mu vraiasemblance

mu_hat_250<- numeric(M)
mu_hat_500<- numeric(M)
mu_hat_750<- numeric(M)
mu_hat_1000<- numeric(M)
mu_hat_1250<- numeric(M)


# valeur des estimateurs mu des moments

mu_mom_250<- numeric(M)
mu_mom_500<- numeric(M)
mu_mom_750<- numeric(M)
mu_mom_1000<- numeric(M)
mu_mom_1250<- numeric(M)

# valeur des estimateurs beta des moments

beta_mom_250<- numeric(M)
beta_mom_500<- numeric(M)
beta_mom_750<- numeric(M)
beta_mom_1000<- numeric(M)
beta_mom_1250<- numeric(M)


# Boucle de simulations
for (i in 1:M) {
  l_250 <- initialisation(250)
  l_500 <- initialisation(500)
  l_750 <- initialisation(750)
  l_1000 <- initialisation(1000)
  l_1250 <- initialisation(1250)
  
  res_250 <- optim(c(4, 3), log_vraissemblance, L=l_250, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_500 <- optim(c(4, 3), log_vraissemblance, L=l_500, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_750 <- optim(c(4, 3), log_vraissemblance, L=l_750, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_1000 <- optim(c(4, 3), log_vraissemblance, L=l_1000, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_1250 <- optim(c(4, 3), log_vraissemblance, L=l_1250, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  
  MSE_250[i] <- (res_250$par[1] - 3)^2
  MSE_500[i] <- (res_500$par[1] - 3)^2
  MSE_750[i] <- (res_750$par[1] - 3)^2
  MSE_1000[i] <- (res_1000$par[1] - 3)^2
  MSE_1250[i] <- (res_1250$par[1] - 3)^2
  
  MSEmu_250[i] <- (res_250$par[2] - 2)^2
  MSEmu_500[i] <- (res_500$par[2] - 2)^2
  MSEmu_750[i] <- (res_750$par[2] - 2)^2
  MSEmu_1000[i] <- (res_1000$par[2] - 2)^2
  MSEmu_1250[i] <- (res_1250$par[2] - 2)^2
  
  B_250[i] <- (res_250$par[1] - 3)
  B_500[i] <- (res_500$par[1] - 3)
  B_750[i] <- (res_750$par[1] - 3)
  B_1000[i] <- (res_1000$par[1] - 3)
  B_1250[i] <- (res_1250$par[1] - 3)
  
  Bmu_250[i] <- (res_250$par[2] - 2)
  Bmu_500[i] <- (res_500$par[2] - 2)
  Bmu_750[i] <- (res_750$par[2] - 2)
  Bmu_1000[i] <- (res_1000$par[2] - 2)
  Bmu_1250[i] <- (res_1250$par[2] - 2)
  
  
  beta_hat_250[i]<- res_250$par[1]
  beta_hat_500[i]<- res_500$par[1]
  beta_hat_750[i]<- res_750$par[1]
  beta_hat_1000[i]<- res_1000$par[1]
  beta_hat_1250[i]<-res_1250$par[1]
  
  mu_hat_250[i]<- res_250$par[2]
  mu_hat_500[i]<- res_500$par[2]
  mu_hat_750[i]<- res_750$par[2]
  mu_hat_1000[i]<- res_1000$par[2]
  mu_hat_1250[i]<-res_1250$par[2]
  
  
  
  # calcul des estimateurs des moments
  
  beta_mom_250[i]<-sqrt((6/pi^2)*sum(mean(l_250^2)-mean(l_250)^2))
  beta_mom_500[i]<-sqrt((6/pi^2)*sum(mean(l_500^2)-mean(l_500)^2))
  beta_mom_750[i]<-sqrt((6/pi^2)*sum(mean(l_750^2)-mean(l_750)^2))
  beta_mom_1000[i]<-sqrt((6/pi^2)*sum(mean(l_1000^2)-mean(l_1000)^2))
  beta_mom_1250[i]<-sqrt((6/pi^2)*sum(mean(l_1250^2)-mean(l_1250)^2))
  
  mu_mom_250[i]<-mean(l_250) - gamma_approx*beta_mom_250[i]
  mu_mom_500[i]<-mean(l_500) - gamma_approx*beta_mom_500[i]
  mu_mom_750[i]<-mean(l_750) - gamma_approx*beta_mom_750[i]
  mu_mom_1000[i]<-mean(l_1000) - gamma_approx*beta_mom_1000[i]
  mu_mom_1250[i]<-mean(l_1250) - gamma_approx*beta_mom_1250[i]
  
  
}






abcis <- c(250, 500, 750, 1000, 1250)

#MSE,biais et variance moyenne pour beta likly
           
MSE_tot <- c(mean(MSE_250), mean(MSE_500), mean(MSE_750), mean(MSE_1000), mean(MSE_1250))
B_tot <- c(mean(B_250), mean(B_500), mean(B_750), mean(B_1000), mean(B_1250))
B_tot <-abs(B_tot)
V_tot<-c(var(beta_hat_250),var(beta_hat_500),var(beta_hat_750),var(beta_hat_1000),var(beta_hat_1250))

#MSE,biais et variance moyenne pour mu likly

MSE_tot_mu <- c(mean(MSEmu_250), mean(MSEmu_500), mean(MSEmu_750), mean(MSEmu_1000), mean(MSEmu_1250))
B_tot_mu <- c(mean(Bmu_250), mean(Bmu_500), mean(Bmu_750), mean(Bmu_1000), mean(Bmu_1250))
B_tot_mu <- abs(B_tot_mu)
V_tot_mu<-c(var(mu_hat_250),var(mu_hat_500),var(mu_hat_750),var(mu_hat_1000),var(mu_hat_1250))



# on calcul le biais des estimateurs de moment

B_tot_mom_mu <-c(mean(mu_mom_250 -2 ),mean(mu_mom_500 -2 ),mean(mu_mom_750 -2 ),mean(mu_mom_1000 -2 ),mean(mu_mom_1250 -2 ))
B_tot_mom_beta <-c(mean(beta_mom_250 -3 ),mean(beta_mom_500 -3 ),mean(beta_mom_750 -3 ),mean(beta_mom_1000 -3 ),mean(beta_mom_1250 -3 ))
B_tot_mom_beta<-abs(B_tot_mom_beta)
B_tot_mom_mu<-abs(B_tot_mom_mu)


# on calcul la MSE des estimateurs de moment
MSE_tot_mom_mu <-c(mean((mu_mom_250 -2)^2 ),mean((mu_mom_500 -2)^2 ),mean((mu_mom_750 -2)^2 ),mean((mu_mom_1000 -2)^2 ),mean((mu_mom_1250 -2 )^2))
MSE_tot_mom_beta <-c(mean((beta_mom_250 -3)^2 ),mean((beta_mom_500 -3)^2 ),mean((beta_mom_750 -3)^2 ),mean((beta_mom_1000 -3)^2 ),mean((beta_mom_1250 -3)^2 ))


# on calcul la variance des estimateurs de moment

Var_tot_mom_mu <-c(var(mu_mom_250),var(mu_mom_500),var(mu_mom_750),var(mu_mom_1000),var(mu_mom_1250))
Var_tot_mom_beta <-c(var(beta_mom_250),var(beta_mom_500),var(beta_mom_750),var(beta_mom_1000),var(beta_mom_1250))


# Tracé des MSE pour beta
plot(abcis, MSE_tot, type="b", col="blue", xlab="Taille de l'échantillon", ylab="MSE", main="MSE de beta en fonction de la taille de l'échantillon")
lines(abcis, MSE_tot_mom_beta, col="red")
legend("topright", legend = c("Likley", "moment"), col = c("blue", "red"), lty = 1)


# Tracé des MSE pour mu
plot(abcis, MSE_tot_mu, type="b", col="blue", xlab="Taille de l'échantillon", ylab="MSE", main="MSE de mu en fonction de la taille de l'échantillon")
lines(abcis, MSE_tot_mom_beta, col="red")
legend("topright", legend = c("Likley", "moment"), col = c("blue", "red"), lty = 1)




# Tracé des biais pour beta
plot(abcis, B_tot, type="b", col="blue", xlab="Taille de l'échantillon", ylab="Biais", main="Biais de beta en fonction de la taille de l'échantillon")
lines(abcis, B_tot_mom_mu, col="red")
legend("topright", legend = c("likely", "moment"), col = c("blue", "red"), lty = 1)


# Tracé des biais pour mu
plot(abcis, B_tot_mu, type="b", col="blue", xlab="Taille de l'échantillon", ylab="Biais", main="Biais de mu en fonction de la taille de l'échantillon")
lines(abcis, B_tot_mom_beta, col="red")
legend("topright", legend = c("likely", "moment"), col = c("blue", "red"), lty = 1)


# Tracé de la variance pour beta
plot(abcis, V_tot, type="b", col="blue", xlab="Taille de l'échantillon", ylab="Variance", main="Variance de beta en fonction de la taille de l'échantillon")
lines(abcis, Var_tot_mom_beta, col="red")
legend("topright", legend = c("likely", "moment"), col = c("blue", "red"), lty = 1)


# Tracé de la variance pour mu
plot(abcis, V_tot_mu, type="b", col="blue", xlab="Taille de l'échantillon", ylab="Variance", main="Variance de mu en fonction de la taille de l'échantillon")
lines(abcis, Var_tot_mom_mu, col="red")
legend("topright", legend = c("likely", "moment"), col = c("blue", "red"), lty = 1)



B_V_tot <-c(B_tot[1]^2 + V_tot[1],B_tot[2]^2 + V_tot[2],B_tot[3]^2 + V_tot[3],B_tot[4]^2 + V_tot[4],B_tot[5]^2 + V_tot[5])
plot(abcis, B_V_tot, type="b", col="blue", xlab="Taille de l'échantillon", ylab="MSE", main="biais de beta carré + variance de beta")



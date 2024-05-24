library(stats)
library(ggplot2)
library(pracma)

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
  result <- sum(log(1/beta) + (mu - L)/beta - exp((mu - L)/beta))
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
  
  B_250[i] <- (res_250$par[1] - 3)^2
  B_500[i] <- (res_500$par[1] - 3)^2
  B_750[i] <- (res_750$par[1] - 3)^2
  B_1000[i] <- (res_1000$par[1] - 3)^2
  B_1250[i] <- (res_1250$par[1] - 3)^2
  
  Bmu_250[i] <- (res_250$par[2] - 2)^2
  Bmu_500[i] <- (res_500$par[2] - 2)^2
  Bmu_750[i] <- (res_750$par[2] - 2)^2
  Bmu_1000[i] <- (res_1000$par[2] - 2)^2
  Bmu_1250[i] <- (res_1250$par[2] - 2)^2
  
  
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
  
  
}



abcis <- c(250, 500, 750, 1000, 1250)
MSE_tot <- c(mean(MSE_250), mean(MSE_500), mean(MSE_750), mean(MSE_1000), mean(MSE_1250))
B_tot <- c(mean(B_250), mean(B_500), mean(B_750), mean(B_1000), mean(B_1250))
V_tot<-c(var(beta_hat_250),var(beta_hat_500),var(beta_hat_750),var(beta_hat_1000),var(beta_hat_1250))

# Tracé des MSE
plot(abcis, MSE_tot, type="b", col="blue", xlab="Taille de l'échantillon", ylab="MSE", main="MSE en fonction de la taille de l'échantillon")

# Tracé des biais
plot(abcis, B_tot, type="b", col="red", xlab="Taille de l'échantillon", ylab="Biais", main="Biais en fonction de la taille de l'échantillon")

# Tracé de la variance 
plot(abcis, V_tot, type="b", col="red", xlab="Taille de l'échantillon", ylab="Variance", main="Variance en fonction de la taille de l'échantillon")

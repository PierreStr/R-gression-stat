#Partie 1.2 du projet du cours "Probabiités et analyse statistiques"


#Tout d'abord, commencons par implémenter certains paramètres importants


# Approximation de la constante d'Euler-Mascheroni par une somme de série
n.tilt <- 1000000
gamma_approx <- sum(1 / seq(1, n.tilt)) - log(n.tilt)






#1) Création artificielle d'une distribution de Gumbel à l'aide de la méthode proposé dans l'énoncé
n = 200
unif.distrib <- runif(n, min = 0, max = 1) #échantillon d'une distribution uniforme entre 0 et 1
gumb.distrib.betha.mu <- seq(1, n)
betha <- 1
mu <- 0

for (i in 1:n){
  gumb.distrib.betha.mu[i] <- -betha * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
}

#Généralisons en construisant une fonction
gumb.distrib <- function(n, betha, mu) {
  unif.distrib <- runif(n, min = 0, max = 1)
  gumb.distrib.initial <- seq(1, n)
  for (i in 1:n){
    gumb.distrib.initial[i] <- -betha * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
  }
  return(gumb.distrib.initial)
}



#Nous avons 200 valeurs aléatoires provenant d'une distribution de Gumbel, nous allons calculer l'estimateur des 
#moments et l'estimateur du maximum de vraisemblance sur base de cet échantillon

#1.1 On calcule betha.hat et mu.hat avec les valeurs de mu et betha, sans utiliser les valeurs issues de notre 
#échantillon aléatoires
density.function <- function(x) {
  betha^-1 * exp(-exp(-(x - mu)/betha)) * exp(-(x - mu)/betha)
}

first.moment.function <- function(x) {
  x * density.function(x)
}

second.moment.function <- function(x) {
  x^2 * density.function(x)
}

first.moment <- integrate(first.moment.function, lower = -500, upper = 500)
second.moment <- integrate(second.moment.function, lower = -500, upper = 500)

mu.theo <- first.moment$value - (sqrt(6/pi^2) * sqrt(second.moment$value - (first.moment$value * first.moment$value)) * gamma_approx)
betha.theo <- sqrt(6/pi^2) * sqrt(second.moment$value - (first.moment$value * first.moment$value))







#Calculons les deux premiers moments sur base de notre échantillon de 200 valeurs aléatoires 

first.moment.hat <- sum(gumb.distrib.betha.mu)/ (n-1)
second.moment.hat <- sum(gumb.distrib.betha.mu^2)/(n-1)

#pour comparer
first.moment$value
second.moment$value

mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
betha.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))

#Généralisons cette méthode avec une fonction
moment.method.function <- function(betha, mu, data) {
  first.moment.hat <- sum(gumb.distrib.betha.mu)/ (n-1)
  second.moment.hat <- sum(gumb.distrib.betha.mu^2)/(n-1)
  mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
  betha.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
  return(c(betha.hat, mu.hat))
}

#Appliquons maintenant la méthode du maximum de vraisemblance 

log.likelihood.function <- function(betha, mu, data) {
  n <- length(data)
  stack <- numeric(n)
  for (k in 1:n) {
    stack[k] <- -log(betha) - (data[k] - mu) / betha - exp(-(data[k] - mu) / betha)
  }
  return(sum(stack))
}

neg.log.likelihood <- function(params, data) {
  betha <- params[1]
  mu <- params[2]
  return(-log.likelihood.function(betha, mu, data))
}

#Utilisons la fonction optim pour obtenir les paramètres idéaux
betha.initial <- 1  #Valeur initiale pour beta
mu.initial <- 0   #Valeur initiale pour mu
initial.params <- c(betha.initial, mu.initial)
data <- gumb.distrib.betha.mu     

opt.result <- optim(par = initial.params, fn = neg.log.likelihood, data = data)
opt.params <- opt.result$par
opt.params



#2) Effectuons cette opération m = 100 fois

#RAPPEL DES FONCTIONS
#Généralisons en construisant une fonction
gumb.distrib <- function(n, betha, mu) {
  unif.distrib <- runif(n, min = 0, max = 1)
  gumb.distrib.initial <- seq(1, n)
  for (i in 1:n){
    gumb.distrib.initial[i] <- -betha * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
  }
  return(gumb.distrib.initial)
}
moment.method.function <- function(betha, mu, data) {
    first.moment.hat <- sum(data)/ (n-1)
    second.moment.hat <- sum(data^2)/(n-1)
    mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
    betha.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
    return(c(betha.hat, mu.hat))
}
log.likelihood.function <- function(betha, mu, data) {
  n <- length(data)
  stack <- numeric(n)
  for (k in 1:n) {
    stack[k] <- -log(betha) - (data[k] - mu) / betha - exp(-(data[k] - mu) / betha)
  }
  return(sum(stack))
}
neg.log.likelihood <- function(params, data) {
  betha <- params[1]
  mu <- params[2]
  return(-log.likelihood.function(betha, mu, data))
}



m <- 100
betha.M <- seq(1:m)
mu.M <- seq(1:m)
betha.MV <- seq(1:m)
mu.MV <- seq(1:m)

for (l in 1:m) {
  data.use <- gumb.distrib(200, 1, 0)
  betha.M[l] <- moment.method.function(1, 0, data.use)[1]
  mu.M[l] <- moment.method.function(3, 2, data.use)[2]
  opt.result <- optim(par = initial.params, fn = neg.log.likelihood, data = data.use)
  opt.params <- opt.result$par
  betha.MV[l] <- opt.params[1]
  mu.MV[l] <- opt.params[2]
}

theta.M <- c(betha.M, mu.M)
theta.MV <- c(betha.MV, mu.MV)
theta.M <- matrix(nrow = 2, ncol = 100)
theta.M[1, ] <- betha.M
theta.M[2, ] <- mu.M

moment.method.function(1, 0, data.use)




#Créons les histogrammes et boxplot
boxplot(x = betha.M)

png("boxplot.png")
boxplot(theta.M, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs",
        names = c("betha_hat", "mu_hat"))

boxplot(betha.M, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
boxplot(mu.M, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
boxplot(betha.MV, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
boxplot(mu.MV, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
?hist
?boxplot


# Sauvegarder l'histogramme de betha_hat comme une image PNG
png("hist_betha_hat.png")
hist(theta.M[1,], main = "Histogramme de betha_hat", xlab = "betha_hat", ylab = "Fréquence")
dev.off()

# Sauvegarder l'histogramme de mu_hat comme une image PNG
png("hist_mu_hat.png")
hist(theta.M[2,], main = "Histogramme de mu_hat", xlab = "mu_hat", ylab = "Fréquence")
dev.off()


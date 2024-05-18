#Code de la partie 1.2 du projet du cours "Probabilités et analyse statistique"

#Commençons par implémenter des éléments et fonctions utile pour la suite 


# Approximation de la constante d'Euler-Mascheroni par une somme de série
n.tilt <- 1000000
gamma_approx <- sum(1 / seq(1, n.tilt)) - log(n.tilt)


#a)

#Utilisons la méthode proposée dans l'énoncé afin de générer un échantillon de valeurs issues d'une loi de Gumbel
n = 200
betha <- 3
mu <- 2
theta <- c(betha, mu)

unif.distrib <- runif(n, min = 0, max = 1)
gumb.distrib.betha.mu <- vector("numeric", length = n)
for (i in 1:n){
  gumb.distrib.betha.mu[i] <- -betha * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
}


#Maintenant que nous avons trouvé comment obtenir des valeurs aléatoires de notre loi de Gumbel, généralisons cela
#à l'aide d'une fonction 

rv.generator <- function(n, betha, mu) {
  unif.distrib <- runif(n, min = 0, max = 1)
  gumb.distrib.initial <- vector("numeric", length = n)
  for (i in 1:n){
    gumb.distrib.initial[i] <- -betha * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
  }
  return(gumb.distrib.initial)
}


#En utilisant cet échantillon, nous allons calculer l'estimateur des moements et l'estimateur du maximum de vraisemblance de 
#notre distribution

first.moment.hat <- sum(gumb.distrib.betha.mu)/ (n)
second.moment.hat <- sum(gumb.distrib.betha.mu^2)/(n)
#REMARQUE : Faut-il diviser par n ou plutôt par n-1 ? 

mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
betha.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
theta.hat <- c(betha.hat, mu.hat)

#Comme pour la première partie, généralisons ce que nous venons de faire avec des fonctions 

moment.est <- function(betha, mu, data) {
  n <- length(data)
  first.moment.hat <- sum(data)/ (n-1)
  second.moment.hat <- sum(data^2)/(n-1)
  mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
  betha.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
  return(c(betha.hat, mu.hat))
}

#Appliquons maintenant la méthode du maximum de vraisemblance 

log.likelihood.est <- function(theta, data) {
  n <- length(data)
  var.inter <- numeric(n)
  for (k in 1:n) {
    var.inter[k] <- -log(betha) - (data[k] - mu) / betha - exp(-(data[k] - mu) / betha)
  }
  return(sum(var.inter))
}

log.likelihood.est <- function(betha, mu, data) {
  n <- length(data)
  var.inter <- numeric(n)
  for (k in 1:n) {
    var.inter[k] <- -log(betha) - (data[k] - mu) / betha - exp(-(data[k] - mu) / betha)
  }
  return(sum(var.inter))
}

log.likelihood.est(3, 2, data = rv.generator(200, 1, 0))

neg.log.likelihood.est <- function(params, data) {
  betha <- params[1]
  mu <- params[2]
  return(-log.likelihood.est(betha, mu, data))
}


#Utilisons la fonction optim pour obtenir les paramètres idéaux
betha.initial <- 1  
mu.initial <- 0
initial.params <- c(betha.initial, mu.initial)
data <- gumb.distrib.betha.mu     

opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data)
opt.params <- opt.result$par
opt.params


#b) Effectuons cette opération m = 100 fois pour obtenir 100 valeurs des estimateurs théta selon chacune des deux méthodes

m <- 100
betha.M <- vector("numeric", length = m)
mu.M <- vector("numeric", length = m)
betha.MV <- vector("numeric", length = m)
mu.MV <- vector("numeric", length = m)

for (l in 1:m) {
  data.use <- rv.generator(200, betha, mu)
  betha.M[l] <- moment.est(betha, mu, data.use)[1]
  mu.M[l] <- moment.est(betha, mu, data.use)[2]
  opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data.use)
  opt.params <- opt.result$par
  betha.MV[l] <- opt.params[1]
  mu.MV[l] <- opt.params[2]
}

#Deux manière de stocker les données
theta.M <- c(betha.M, mu.M)
theta.MV <- c(betha.MV, mu.MV)

theta.M <- matrix(nrow = 2, ncol = 100)
theta.MV <- matrix(nrow = 2, ncol = 100)
theta.M[1, ] <- betha.M
theta.M[2, ] <- mu.M
theta.MV[1, ] <- betha.MV
theta.MV[2, ] <- mu.MV

#Créons maintenant les histogrammes et boxplot

boxplot(betha.M, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
boxplot(mu.M, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
boxplot(betha.MV, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")
boxplot(mu.MV, main = "Boxplot des estimateurs", xlab = "Paramètres", ylab = "Estimateurs")

png("hist_betha_hat.png")
hist(theta.M[1,], main = "Histogramme de betha_hat", xlab = "betha_hat", ylab = "Fréquence")
dev.off()

png("hist_mu_hat.png")
hist(theta.M[2,], main = "Histogramme de mu_hat", xlab = "mu_hat", ylab = "Fréquence")
dev.off()

png("hist_betha_hat.png")
hist(theta.MV[1,], main = "Histogramme de betha_hat", xlab = "betha_hat", ylab = "Fréquence")
dev.off()

png("hist_mu_hat.png")
hist(theta.MV[2,], main = "Histogramme de mu_hat", xlab = "mu_hat", ylab = "Fréquence")
dev.off()


#c) Calculons le biais, la variance et le MSE 

#Calculons quelques valeurs nécessaire pour la suite 
esp.bheta.hat.M <- sum(betha.M)/m
esp.mu.hat.M <- sum(mu.M)/m
esp.bheta.hat.MV <- sum(betha.MV)/m
esp.mu.hat.MV <- sum(mu.MV)/m

esp.betha.hat.squared.M <- sum(betha.M^2)/m
esp.mu.hat.squared.M <- sum(mu.M^2)/m
esp.betha.hat.squared.MV <- sum(betha.MV^2)/m
esp.mu.hat.squared.MV <- sum(mu.MV^2)/m

bias.M <- c(esp.bheta.hat.M - betha, esp.mu.hat.M- mu)
bias.MV <- c(esp.bheta.hat.MV - betha, esp.mu.hat.MV - mu)

var.M <- c(esp.betha.hat.squared.M - (esp.bheta.hat.M)^2, esp.mu.hat.squared.M - (esp.mu.hat.M)^2)
var.MV <- c(esp.betha.hat.squared.MV - (esp.bheta.hat.MV)^2, esp.mu.hat.squared.MV - (esp.mu.hat.MV)^2)

MSE.M <- (bias.M)^2 + var.M
MSE.MV <- (bias.MV)^2 + var.MV

#IL FAUT VERIFIER QUE LES CALCULS SONT CORRECTS ET COHERENTS ET PUIS FAIRE DES CONCLUSIONS DESSUS
#Les valeurs concernant mu ne devraient pas être plus basse dans la stratégie MV que la stratégie M ?
#Est-ce qu'on fait des graphes ?

#d)

#Créons une formule générale pour les calculs effectués en c) afin de faciliter les choses


theta.M.est <- function(n, m, betha, mu) {
  betha.M <- vector("numeric", length = m)
  mu.M <- vector("numeric", length = m)
  for (l in 1:m) {
    data.use <- rv.generator(n, betha, mu)
    betha.M[l] <- moment.est(betha, mu, data.use)[1]
    mu.M[l] <- moment.est(betha, mu, data.use)[2]
  }
  theta.M[1, ] <- betha.M
  theta.M[2, ] <- mu.M
  return(theta.M)
}

theta.MV.est <- function(n, m, betha, mu) {
  betha.MV <- vector("numeric", length = m)
  mu.MV <- vector("numeric", length = m)
  for (l in 1:m) {
    data.use <- rv.generator(n, betha, mu)
    opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data.use)
    opt.params <- opt.result$par
    betha.MV[l] <- opt.params[1]
    mu.MV[l] <- opt.params[2]
  }
  theta.MV[1, ] <- betha.MV
  theta.MV[2, ] <- mu.MV
  return(theta.MV)
}

#ATTENTION GROSSES INCOHERENTES POUR DES VALEURS HAUTES TELS QUE BHETA = 7 
aa <-   theta.M.est(2000, 100, betha = 7, mu = 0)
ab <- theta.M.est(200, 100, betha = 7, mu = 0)
ba <- theta.MV.est(200, 100, betha = 7, mu = 2)
bb <- theta.MV.est(2000, 100, betha = 7, mu = 2)
bc <- theta.MV.est(20000, 100, betha = 7, mu = 2)

rv.generator(2000, 1, 0)


bias.var.MSE <- function(n, m, betha, mu) {
  esp.bheta.hat.M <- sum(theta.M.est(n, m, betha, mu)[1])/m
  esp.mu.hat.M <- sum(theta.M.est(n, m, betha, mu)[2])/m
  esp.bheta.hat.MV <- sum(theta.MV.est(n, m, betha, mu)[1])/m
  esp.mu.hat.MV <- sum(theta.MV.est(n, m, betha, mu)[2])/m
  esp.betha.hat.squared.M <- sum(theta.M.est(n, m, betha, mu)[1]^2)/m
  esp.mu.hat.squared.M <- sum(theta.M.est(n, m, betha, mu)[2]^2)/m
  esp.betha.hat.squared.MV <- sum(theta.MV.est(n, m, betha, mu)[1]^2)/m
  esp.mu.hat.squared.MV <- sum(theta.MV.est(n, m, betha, mu)[2]^2)/m
  
  bias.M <- c(esp.bheta.hat.M - betha, esp.mu.hat.M- mu)
  bias.MV <- c(esp.bheta.hat.MV - betha, esp.mu.hat.MV - mu)
  
  var.M <- c(esp.betha.hat.squared.M - (esp.bheta.hat.M)^2, esp.mu.hat.squared.M - (esp.mu.hat.M)^2)
  var.MV <- c(esp.betha.hat.squared.MV - (esp.bheta.hat.MV)^2, esp.mu.hat.squared.MV - (esp.mu.hat.MV)^2)
  
  MSE.M <- (bias.M)^2 + var.M
  MSE.MV <- (bias.MV)^2 + var.MV
  
  result <- matrix(nrow = 2, ncol = 6)
  result.M <- c(bias.M, var.M, MSE.M)
  result.MV <- c(bias.MV, var.MV, MSE.MV)
  result[1,] <- result.M
  result[2,] <- result.MV
  return(result)
}

test <- bias.var.MSE(200, 100, 3, 2)

bias.M
var.M
MSE.M

test



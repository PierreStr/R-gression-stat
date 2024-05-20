#Code de la partie 1.2 du projet du cours "Probabilités et analyse statistique"

#Commençons par implémenter des éléments et fonctions utile pour la suite 


# Approximation de la constante d'Euler-Mascheroni par une somme de série
n.tilt <- 1000000
gamma_approx <- sum(1 / seq(1, n.tilt)) - log(n.tilt)

#Importation d'un package
library(ggplot2)

#Implémentons une seed
set.seed('313')

#a)

#Utilisons la méthode proposée dans l'énoncé afin de générer un échantillon de valeurs issues d'une loi de Gumbel
n = 200
beta.choosen <- 3
mu.choosen <- 2
theta.choosen <- c(beta.choosen, mu.choosen)

unif.distrib <- runif(n, min = 0, max = 1)
gumb.distrib.beta.mu <- vector("numeric", length = n)
for (i in 1:n){
  gumb.distrib.beta.mu[i] <- -beta.choosen * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu.choosen
}

#Maintenant que nous avons trouvé comment obtenir des valeurs aléatoires de notre loi de Gumbel, généralisons cela
#à l'aide d'une fonction 

rv.generator <- function(n, beta, mu) {
  unif.distrib <- runif(n, min = 0, max = 1)
  gumb.distrib.initial <- vector("numeric", length = n)
  for (i in 1:n){
    gumb.distrib.initial[i] <- -beta * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
  }
  return(gumb.distrib.initial)
}


#En utilisant cet échantillon, nous allons calculer l'estimateur des moements et l'estimateur du maximum de vraisemblance de 
#notre distribution

first.moment.hat <- sum(gumb.distrib.beta.mu)/ (n)
second.moment.hat <- sum(gumb.distrib.beta.mu^2)/(n)

mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
beta.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
theta.hat <- c(beta.hat, mu.hat)

#Comme pour la première partie, généralisons ce que nous venons de faire avec des fonctions 

moment.est <- function(beta, mu, data) {
  n <- length(data)
  first.moment.hat <- sum(data)/ (n)
  second.moment.hat <- sum(data^2)/(n)
  mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
  beta.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
  return(c(beta.hat, mu.hat))
}

moment.est(3, 2, data = rv.generator(200, 3, 2))
#Appliquons maintenant la méthode du maximum de vraisemblance 

log.likelihood.est <- function(theta, data) {
  n <- length(data)
  var.inter <- numeric(n)
  for (k in 1:n) {
    var.inter[k] <- -log(beta) - (data[k] - mu) / beta - exp(-(data[k] - mu) / beta)
  }
  return(sum(var.inter))
}

log.likelihood.est <- function(beta, mu, data) {
  n <- length(data)
  var.inter <- numeric(n)
  for (k in 1:n) {
    var.inter[k] <- -log(beta) - (data[k] - mu) / beta - exp(-(data[k] - mu) / beta)
  }
  return(sum(var.inter))
}

neg.log.likelihood.est <- function(params, data) {
  beta <- params[1]
  mu <- params[2]
  return(-log.likelihood.est(beta, mu, data))
}


#Utilisons la fonction optim pour obtenir les paramètres idéaux
beta.initial <- 1  
mu.initial <- 0
initial.params <- c(beta.initial, mu.initial)
data <- gumb.distrib.beta.mu     

opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data)
opt.params <- opt.result$par
opt.params


#b) Effectuons cette opération m = 100 fois pour obtenir 100 valeurs des estimateurs théta selon chacune des deux méthodes

m <- 100
beta.M <- vector("numeric", length = m)
mu.M <- vector("numeric", length = m)
beta.MV <- vector("numeric", length = m)
mu.MV <- vector("numeric", length = m)

for (l in 1:m) {
  data.use <- rv.generator(200, beta.choosen, mu.choosen)
  beta.M[l] <- moment.est(beta.choosen, mu.choosen, data.use)[1]
  mu.M[l] <- moment.est(beta.choosen, mu.choosen, data.use)[2]
  opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data.use)
  opt.params <- opt.result$par
  beta.MV[l] <- opt.params[1]
  mu.MV[l] <- opt.params[2]
}

#Deux manière de stocker les données

theta.M <- matrix(nrow = 2, ncol = 100)
theta.MV <- matrix(nrow = 2, ncol = 100)
theta.M[1, ] <- beta.M
theta.M[2, ] <- mu.M
theta.MV[1, ] <- beta.MV
theta.MV[2, ] <- mu.MV

#Créons maintenant les histogrammes et boxplot
png("boxplot_beta_hat_M.png")
boxplot(beta.M, main = "Boxplot des estimateurs par la méthode des moments pour beta", xlab = "Paramètres", ylab = "Estimateurs")
dev.off()

png("boxplot_mu_hat_M.png")
boxplot(mu.M, main = "Boxplot des estimateurs par la méthode des moments pour mu", xlab = "Paramètres", ylab = "Estimateurs")
dev.off()

png("boxplot_beta_hat_MV.png")
boxplot(beta.MV, main = "Boxplot des estimateurs par la méthode du MV pour beta", xlab = "Paramètres", ylab = "Estimateurs")
dev.off()

png("boxplot_mu_hat_MV.png")
boxplot(mu.MV, main = "Boxplot des estimateurs par la méthode du MV pour mu", xlab = "Paramètres", ylab = "Estimateurs")
dev.off()

png("hist_beta_hat_M.png")
hist(theta.M[1,], main = "Histogramme de β_hat", xlab = 'valeurs de β_hat', ylab = "Fréquence")
dev.off()

png("hist_mu_hat_M.png")
hist(theta.M[2,], main = "Histogramme de μ_hat", xlab = "valeurs de μ_hat", ylab = "Fréquence")
dev.off()

png("hist_beta_hat_MV.png")
hist(theta.MV[1,], main = "Histogramme de β_hat", xlab = "valeurs de β_hat", ylab = "Fréquence")
dev.off()

png("hist_mu_hat_MV.png")
hist(theta.MV[2,], main = "Histogramme de μ_hat", xlab = "valeurs de μ_hat", ylab = "Fréquence")
dev.off()


#c) Calculons le biais, la variance et le MSE 

#Calculons quelques valeurs nécessaire pour la suite 
esp.beta.hat.M <- sum(beta.M)/m
esp.mu.hat.M <- sum(mu.M)/m
esp.beta.hat.MV <- sum(beta.MV)/m
esp.mu.hat.MV <- sum(mu.MV)/m

esp.beta.hat.squared.M <- sum(beta.M^2)/m
esp.mu.hat.squared.M <- sum(mu.M^2)/m
esp.beta.hat.squared.MV <- sum(beta.MV^2)/m
esp.mu.hat.squared.MV <- sum(mu.MV^2)/m

bias.M <- c(esp.beta.hat.M - beta.choosen, esp.mu.hat.M- mu.choosen)
bias.MV <- c(esp.beta.hat.MV - beta.choosen, esp.mu.hat.MV - mu.choosen)

var.M <- c(esp.beta.hat.squared.M - (esp.beta.hat.M)^2, esp.mu.hat.squared.M - (esp.mu.hat.M)^2)
var.MV <- c(esp.beta.hat.squared.MV - (esp.beta.hat.MV)^2, esp.mu.hat.squared.MV - (esp.mu.hat.MV)^2)

MSE.M <- (bias.M)^2 + var.M
MSE.MV <- (bias.MV)^2 + var.MV

#d)

#Créons une formule générale pour les calculs effectués en c) afin de faciliter les choses


theta.M.est <- function(n, m, beta, mu) {
  beta.M <- vector("numeric", length = m)
  mu.M <- vector("numeric", length = m)
  for (l in 1:m) {
    data.use <- rv.generator(n, beta, mu)
    beta.M[l] <- moment.est(beta, mu, data.use)[1]
    mu.M[l] <- moment.est(beta, mu, data.use)[2]
  }
  theta.M[1, ] <- beta.M
  theta.M[2, ] <- mu.M
  return(theta.M)
}

theta.MV.est <- function(n, m, beta, mu) {
  beta.MV <- vector("numeric", length = m)
  mu.MV <- vector("numeric", length = m)
  initial.params <- c(beta, mu)
  for (l in 1:m) {
    data.use <- rv.generator(n, beta, mu)
    opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data.use, method = "BFGS")
    opt.params <- opt.result$par
    beta.MV[l] <- opt.params[1]
    mu.MV[l] <- opt.params[2]
  }
  theta.MV[1, ] <- beta.MV
  theta.MV[2, ] <- mu.MV
  return(theta.MV)
}


bias.var.MSE <- function(n, m, beta, mu) {
  esp.beta.hat.M.temp <- sum(theta.M.est(n, m, beta, mu)[1, ])/m
  esp.mu.hat.M.temp <- sum(theta.M.est(n, m, beta, mu)[2, ])/m
  esp.beta.hat.MV.temp <- sum(theta.MV.est(n, m, beta, mu)[1, ])/m
  esp.mu.hat.MV.temp <- sum(theta.MV.est(n, m, beta, mu)[2, ])/m
  esp.beta.hat.squared.M.temp <- sum(theta.M.est(n, m, beta, mu)[1, ]^2)/m
  esp.mu.hat.squared.M.temp <- sum(theta.M.est(n, m, beta, mu)[2, ]^2)/m
  esp.beta.hat.squared.MV.temp <- sum(theta.MV.est(n, m, beta, mu)[1, ]^2)/m
  esp.mu.hat.squared.MV.temp <- sum(theta.MV.est(n, m, beta, mu)[2, ]^2)/m
  
  bias.M.temp <- c(esp.beta.hat.M.temp - beta, esp.mu.hat.M.temp- mu)
  bias.MV.temp <- c(esp.beta.hat.MV.temp - beta, esp.mu.hat.MV.temp - mu)
  
  var.M.temp <- c(esp.beta.hat.squared.M.temp - (esp.beta.hat.M.temp)^2, esp.mu.hat.squared.M.temp - (esp.mu.hat.M.temp)^2)
  var.MV.temp <- c(esp.beta.hat.squared.MV.temp - (esp.beta.hat.MV.temp)^2, esp.mu.hat.squared.MV.temp - (esp.mu.hat.MV.temp)^2)
  
  MSE.M.temp <- (bias.M.temp)^2 + var.M.temp
  MSE.MV.temp <- (bias.MV.temp)^2 + var.MV.temp
  
  result <- matrix(nrow = 2, ncol = 6)
  result.M <- c(bias.M.temp, var.M.temp, MSE.M.temp)
  result.MV <- c(bias.MV.temp, var.MV.temp, MSE.MV.temp)
  result[1,] <- result.M
  result[2,] <- result.MV
  return(result)
}

n.values <- c(250, 500, 750, 1000, 1250)

bias.n.M.beta <- vector("numeric", length = 5)
bias.n.MV.beta <- vector("numeric", length = 5)
var.n.M.beta <- vector("numeric", length = 5)
var.n.MV.beta <- vector("numeric", length = 5)
MSE.n.M.beta <- vector("numeric", length = 5)
MSE.n.MV.beta <- vector("numeric", length = 5)
bias.n.M.mu <- vector("numeric", length = 5)
bias.n.MV.mu <- vector("numeric", length = 5)
var.n.M.mu <- vector("numeric", length = 5)
var.n.MV.mu <- vector("numeric", length = 5)
MSE.n.M.mu <- vector("numeric", length = 5)
MSE.n.MV.mu <- vector("numeric", length = 5)

for (i in 1:5) {
  var.int <- bias.var.MSE(n.values[i], 100, 3, 2) #Matrice 2 lignes 6 colonnes 
  bias.n.M.beta[i] <- var.int[1, 1]
  var.n.M.beta[i] <- var.int[1, 3]
  MSE.n.M.beta[i] <- var.int[1, 5]
  bias.n.MV.beta[i] <- var.int[2, 1]
  var.n.MV.beta[i] <- var.int[2, 3]
  MSE.n.MV.beta[i] <- var.int[2, 5]
  bias.n.M.mu[i] <- var.int[1, 2]
  var.n.M.mu[i] <- var.int[1, 4]
  MSE.n.M.mu[i] <- var.int[1, 6]
  bias.n.MV.mu[i] <- var.int[2, 2]
  var.n.MV.mu[i] <- var.int[2, 4]
  MSE.n.MV.mu[i] <- var.int[2, 6]
}

# Calculer les estimations pour chaque taille d'échantillon et stocker les résultats dans des dataframes ou des vecteurs

# Créer un dataframe pour les résultats
results_df.beta <- data.frame(
  n = rep(n.values, each = 2),  # Répéter chaque valeur de n deux fois (une fois pour chaque estimateur)
  estimateur = rep(c("M", "MV"), length(n.values)),  # Indiquer l'estimateur
  biais = abs(c(bias.n.M.beta, bias.n.MV.beta)),  # Les valeurs des biais pour chaque taille d'échantillon et chaque estimateur
  variance = abs(c(var.n.M.beta, var.n.MV.beta)),  # Les valeurs des variances pour chaque taille d'échantillon et chaque estimateur
  MSE = abs(c(MSE.n.M.beta, MSE.n.MV.beta))  # Les valeurs des MSE pour chaque taille d'échantillon et chaque estimateur
)

# Tracer les graphiques
# Graphique pour les biais
ggplot(results_df.beta, aes(x = n, y = biais, color = estimateur, linetype = estimateur)) +
  geom_line() +
  labs(title = "Évolution des biais de beta en fonction de la taille de l'échantillon",
       x = "Taille de l'échantillon (n)",
       y = "Biais") +
  theme_minimal()

# Graphique pour les variances
ggplot(results_df.beta, aes(x = n, y = variance, color = estimateur, linetype = estimateur)) +
  geom_line() +
  labs(title = "Évolution des variances de beta en fonction de la taille de l'échantillon",
       x = "Taille de l'échantillon (n)",
       y = "Variance") +
  theme_minimal()

# Graphique pour les MSE
ggplot(results_df.beta, aes(x = n, y = MSE, color = estimateur, linetype = estimateur)) +
  geom_line() +
  labs(title = "Évolution des MSE de beta en fonction de la taille de l'échantillon",
       x = "Taille de l'échantillon (n)",
       y = "MSE") +
  theme_minimal()



# Créer un dataframe pour les résultats
results_df.mu <- data.frame(
  n = rep(n.values, each = 2),  # Répéter chaque valeur de n deux fois (une fois pour chaque estimateur)
  estimateur = rep(c("M", "MV"), length(n.values)),  # Indiquer l'estimateur
  biais = c(bias.n.M.mu, bias.n.MV.mu),  # Les valeurs des biais pour chaque taille d'échantillon et chaque estimateur
  variance = c(var.n.M.mu, var.n.MV.mu),  # Les valeurs des variances pour chaque taille d'échantillon et chaque estimateur
  MSE = c(MSE.n.M.mu, MSE.n.MV.mu)  # Les valeurs des MSE pour chaque taille d'échantillon et chaque estimateur
)

# Tracer les graphiques
# Graphique pour les biais
ggplot(results_df.mu, aes(x = n, y = biais, color = estimateur, linetype = estimateur)) +
  geom_line() +
  labs(title = "Évolution des biais de mu en fonction de la taille de l'échantillon",
       x = "Taille de l'échantillon (n)",
       y = "Biais") +
  theme_minimal()

# Graphique pour les variances
ggplot(results_df.mu, aes(x = n, y = variance, color = estimateur, linetype = estimateur)) +
  geom_line() +
  labs(title = "Évolution des variances de mu en fonction de la taille de l'échantillon",
       x = "Taille de l'échantillon (n)",
       y = "Variance") +
  theme_minimal()

# Graphique pour les MSE
ggplot(results_df.mu, aes(x = n, y = MSE, color = estimateur, linetype = estimateur)) +
  geom_line() +
  labs(title = "Évolution des MSE de mu en fonction de la taille de l'échantillon",
       x = "Taille de l'échantillon (n)",
       y = "MSE") +
  theme_minimal()


#e)
bias.var.MSE(200, 100, 3, 2)
bias.var.MSE(1000, 100, 3, 2)

# PARTIE 1.1 : Estimation #######################################################

# Chargement des packages
library(stats)
library(ggplot2)
library(pracma)

# Fixer la seed
set.seed("311")

# Chargement des données
dataproject <- read.csv("data_project2024.txt", sep=",")

# Approximation de la constante d'Euler-Mascheroni par une somme de série
n.tilt <- 1000000
gamma_approx <- sum(1 / seq(1, n.tilt)) - log(n.tilt)


# f) ---------------------------------------------------------------


Yi <- dataproject$Y
Y.barre <- mean(Yi)

# Implémentation de l'opposé de la fonction du log de vraissemblance 
neg_log_likelihood<-function(param){
  
  beta <- param[1]
  
  if(beta <= 0){return (NaN)}  # vérifier que beta n'est pas négatif
  
  mu <- param[2]
  n<-length(Yi)
  
  resultat <- -n*log(beta)
  
  for (i in 1:n) {
    tmp<- (mu - Yi[i])/beta
    resultat <- resultat + tmp-exp(tmp)
  }
  
  return (-resultat)}

# Optimisation de la fonction définie plus haut
initial_param<-c(3,5)  

result <- optim(par = initial_param, fn = neg_log_likelihood, method = "BFGS")
opt <- result$par
opt

# Vérification des estimateurs avec les équations données (la valeur retournée doit se rapprocher 0)
# Pour mu
verif.mu<-function(beta.hat, mu.hat){
  
  result <- -beta.hat*log((1/length(Yi))*sum(exp(-(Yi/beta.hat)))) - mu.hat
  return(result)}

#Pour beta
verif.beta<-function(beta.hat, mu.hat){
  
  result <- Y.barre - (sum(Yi*exp(-(Yi/beta.hat))))/(sum(exp(-(Yi/beta.hat)))) - beta.hat
  return(result)}


verifmu <- verif.mu(0.3610547, 2.8413821)
verifbeta <- verif.beta(0.3610547, 2.8413821)
verifmu
verifbeta


# PARTIE 1.2 : Simulations #######################################################


# a) ---------------------------------------------------------------


# Choix des variables à fixer
n = 200
beta.choosen <- 3
mu.choosen <- 2


# Implémentation de la fonction donnant dans une liste n valeurs provenant d'une distribution de Gumbel
rv.generator <- function(n, beta, mu) {
  unif.distrib <- runif(n, min = 0, max = 1)
  gumb.distrib.initial <- vector("numeric", length = n)
  for (i in 1:n){
    gumb.distrib.initial[i] <- -beta * log(-log(unif.distrib[i], base = exp(1)), base = exp(1)) + mu
  }
  return(gumb.distrib.initial)
}

# Création de notre premier échantillons contenant 200 valeurs
gumb.distrib.beta.mu <- rv.generator(n, beta.choosen, mu.choosen)


# Sur base de cet échantillon, on calcule l'estimateur de beta et mu avec la méthode des moments

moment.est <- function(beta, mu, data) {
  n <- length(data)
  first.moment.hat <- sum(data)/ (n)
  second.moment.hat <- sum(data^2)/(n)
  mu.hat <- first.moment.hat - (sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat)) * gamma_approx)
  beta.hat <- sqrt(6/pi^2) * sqrt(second.moment.hat - (first.moment.hat * first.moment.hat))
  return(c(beta.hat, mu.hat))
}

moment.est(3, 2, gumb.distrib.beta.mu) # = (beta_hat, mu_hat)


# De même pour la méthode de vraisemblance

# Implémentation de la fonction du log de vraisemblance dépendante des beta et mu donnés
log.likelihood.est <- function(beta, mu, data) {
  n <- length(data)
  var.inter <- numeric(n)
  for (k in 1:n) {
    var.inter[k] <- -log(beta) - (data[k] - mu) / beta - exp(-(data[k] - mu) / beta)
  }
  return(sum(var.inter))
}

# Return l'opposé du log de vraissemblance (pour minimiser l'opposé et avoir le max)
neg.log.likelihood.est <- function(params, data) {
  beta <- params[1]
  mu <- params[2]
  return(-log.likelihood.est(beta, mu, data))
}

# Utilisation de la fonction optim pour obtenir les paramètres idéaux
beta.initial <- 1  
mu.initial <- 0
initial.params <- c(beta.initial, mu.initial)
data <- gumb.distrib.beta.mu     

opt.result <- optim(par = initial.params, fn = neg.log.likelihood.est, data = data)
opt.params <- opt.result$par
opt.params


# b) ---------------------------------------------------------------


# On effectue m = 100 fois les calculs précédents

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

# Stockage des données 
theta.M <- matrix(nrow = 2, ncol = 100)
theta.MV <- matrix(nrow = 2, ncol = 100)
theta.M[1, ] <- beta.M
theta.M[2, ] <- mu.M
theta.MV[1, ] <- beta.MV
theta.MV[2, ] <- mu.MV

# Création des histogrammes et boxplots
boxplot(beta.M, main = "Boxplot des estimateurs \n par la méthode des moments pour beta", xlab = "Paramètres", ylab = "Estimateurs")

boxplot(mu.M, main = "Boxplot des estimateurs \n par la méthode des moments pour mu", xlab = "Paramètres", ylab = "Estimateurs")

boxplot(beta.MV, main = "Boxplot des estimateurs \n par la méthode du MV pour beta", xlab = "Paramètres", ylab = "Estimateurs")

boxplot(mu.MV, main = "Boxplot des estimateurs \n par la méthode du MV pour mu", xlab = "Paramètres", ylab = "Estimateurs")

hist(theta.M[1,], main = "Histogramme des estimateurs \n par la méthode des moments pour beta", xlab = "valeurs de l'estimateur de beta", ylab = "Fréquence")

hist(theta.M[2,], main = "Histogramme des estimateurs \n par la méthode des moments pour mu", xlab = "valeurs de l'estimateur de mu", ylab = "Fréquence")

hist(theta.MV[1,], main = "Histogramme des estimateurs \n par la méthode du MV pour beta", xlab = "valeurs de l'estimateur de beta", ylab = "Fréquence")

hist(theta.MV[2,], main = "Histogramme des estimateurs \n par la méthode du MV pour mu", xlab = "valeurs de l'estimateur de mu", ylab = "Fréquence")


# c) ---------------------------------------------------------------

# Calculs des biais, des variances et des MSE pour les 2 méthodes

# Calculs de valeurs utiles
esp.beta.hat.M <- sum(beta.M)/m
esp.mu.hat.M <- sum(mu.M)/m
esp.beta.hat.MV <- sum(beta.MV)/m
esp.mu.hat.MV <- sum(mu.MV)/m

esp.beta.hat.squared.M <- sum(beta.M^2)/m
esp.mu.hat.squared.M <- sum(mu.M^2)/m
esp.beta.hat.squared.MV <- sum(beta.MV^2)/m
esp.mu.hat.squared.MV <- sum(mu.MV^2)/m

# Calcul des biais
bias.M <- c(esp.beta.hat.M - beta.choosen, esp.mu.hat.M- mu.choosen)
bias.MV <- c(esp.beta.hat.MV - beta.choosen, esp.mu.hat.MV - mu.choosen)

# Calcul des variances
var.M <- c(esp.beta.hat.squared.M - (esp.beta.hat.M)^2, esp.mu.hat.squared.M - (esp.mu.hat.M)^2)
var.MV <- c(esp.beta.hat.squared.MV - (esp.beta.hat.MV)^2, esp.mu.hat.squared.MV - (esp.mu.hat.MV)^2)

# Calcul des MSE
MSE.M <- (bias.M)^2 + var.M
MSE.MV <- (bias.MV)^2 + var.MV

# Affichage des résultats
bias.M     # = (biais de beta_hat, biais de mu_hat)
var.M      # = ...
MSE.M
      
bias.MV
var.MV
MSE.MV


# d) ---------------------------------------------------------------


# Listes des échantillons suivant une loi de Gumbel
l_250 <- numeric(250)
l_500 <- numeric(500)
l_750 <- numeric(750)
l_1000 <- numeric(1000)
l_1250 <- numeric(1250)

# Valeurs des estimateurs mu des moments
mu_mom_250<- numeric(M)
mu_mom_500<- numeric(M)
mu_mom_750<- numeric(M)
mu_mom_1000<- numeric(M)
mu_mom_1250<- numeric(M)

# Valeur des estimateurs beta des moments
beta_mom_250<- numeric(M)
beta_mom_500<- numeric(M)
beta_mom_750<- numeric(M)
beta_mom_1000<- numeric(M)
beta_mom_1250<- numeric(M)

# Biais pour beta
B_250 <- numeric(M)
B_500 <- numeric(M)
B_750 <- numeric(M)
B_1000 <- numeric(M)
B_1250 <- numeric(M)

# Biais de mu
Bmu_250 <- numeric(M)
Bmu_500 <- numeric(M)
Bmu_750 <- numeric(M)
Bmu_1000 <- numeric(M)
Bmu_1250 <- numeric(M)

# Variances de beta
V_250 <- numeric(M)
V_500 <- numeric(M)
V_750 <- numeric(M)
V_1000 <- numeric(M)
V_1250 <- numeric(M)

# Variances de mu
Vmu_250 <- numeric(M)
Vmu_500 <- numeric(M)
Vmu_750 <- numeric(M)
Vmu_1000 <- numeric(M)
Vmu_1250 <- numeric(M)

# Valeurs des estimateurs de beta vraisemblance
beta_hat_250<- numeric(M)
beta_hat_500<- numeric(M)
beta_hat_750<- numeric(M)
beta_hat_1000<- numeric(M)
beta_hat_1250<- numeric(M)

# Valeurs des estimateurs de mu vraisemblance
mu_hat_250<- numeric(M)
mu_hat_500<- numeric(M)
mu_hat_750<- numeric(M)
mu_hat_1000<- numeric(M)
mu_hat_1250<- numeric(M)

# MSE pour beta
MSE_250 <- numeric(M)
MSE_500 <- numeric(M)
MSE_750 <- numeric(M)
MSE_1000 <- numeric(M)
MSE_1250 <- numeric(M)

# MSE de mu
MSEmu_250 <- numeric(M)
MSEmu_500 <- numeric(M)
MSEmu_750 <- numeric(M)
MSEmu_1000 <- numeric(M)
MSEmu_1250 <- numeric(M)


# Boucle des simulations
for (i in 1:m) {
  l_250 <- rv.generator(250, beta.choosen, mu.choosen)
  l_500 <- rv.generator(500, beta.choosen, mu.choosen)
  l_750 <- rv.generator(750, beta.choosen, mu.choosen)
  l_1000 <- rv.generator(1000, beta.choosen, mu.choosen)
  l_1250 <- rv.generator(1250, beta.choosen, mu.choosen)
  
  initial.params <- c(4, 3)
  
  # Les optimisations sont bornées par [0.1, 5] pour l'estimateur de beta et [0, 5] pour celui de mu
  res_250 <- optim(initial.params, neg.log.likelihood.est, data=l_250, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_500 <- optim(initial.params, neg.log.likelihood.est, data=l_500, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_750 <- optim(initial.params, neg.log.likelihood.est, data=l_750, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_1000 <- optim(initial.params, neg.log.likelihood.est, data=l_1000, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
  res_1250 <- optim(initial.params, neg.log.likelihood.est, data=l_1250, method="L-BFGS-B", lower=c(0.1, 0), upper=c(5, 5))
 
  # Calculs des estimateurs de vraisemblance 
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
  
  # MSE et biais de beta et mu pour la méthode de vraisemblance
  MSE_250[i] <- (res_250$par[1] - beta.choosen)^2
  MSE_500[i] <- (res_500$par[1] - beta.choosen)^2
  MSE_750[i] <- (res_750$par[1] - beta.choosen)^2
  MSE_1000[i] <- (res_1000$par[1] - beta.choosen)^2
  MSE_1250[i] <- (res_1250$par[1] - beta.choosen)^2
  
  MSEmu_250[i] <- (res_250$par[2] - mu.choosen)^2
  MSEmu_500[i] <- (res_500$par[2] - mu.choosen)^2
  MSEmu_750[i] <- (res_750$par[2] - mu.choosen)^2
  MSEmu_1000[i] <- (res_1000$par[2] - mu.choosen)^2
  MSEmu_1250[i] <- (res_1250$par[2] - mu.choosen)^2
  
  B_250[i] <- (res_250$par[1] - beta.choosen)
  B_500[i] <- (res_500$par[1] - beta.choosen)
  B_750[i] <- (res_750$par[1] - beta.choosen)
  B_1000[i] <- (res_1000$par[1] - beta.choosen)
  B_1250[i] <- (res_1250$par[1] - beta.choosen)
  
  Bmu_250[i] <- (res_250$par[2] - mu.choosen)
  Bmu_500[i] <- (res_500$par[2] - mu.choosen)
  Bmu_750[i] <- (res_750$par[2] - mu.choosen)
  Bmu_1000[i] <- (res_1000$par[2] - mu.choosen)
  Bmu_1250[i] <- (res_1250$par[2] - mu.choosen)
  
  
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


# MSE, biais et variance moyenne pour l'estimateur de beta par la méthode de vraisemblance 
MSE_tot <- c(mean(MSE_250), mean(MSE_500), mean(MSE_750), mean(MSE_1000), mean(MSE_1250))
B_tot <- c(mean(B_250), mean(B_500), mean(B_750), mean(B_1000), mean(B_1250))
V_tot<-c(var(beta_hat_250),var(beta_hat_500),var(beta_hat_750),var(beta_hat_1000),var(beta_hat_1250))

# MSE, biais et variance moyenne pour l'estimateur de mu par la méthode de vraisemblance 
MSE_tot_mu <- c(mean(MSEmu_250), mean(MSEmu_500), mean(MSEmu_750), mean(MSEmu_1000), mean(MSEmu_1250))
B_tot_mu <- c(mean(Bmu_250), mean(Bmu_500), mean(Bmu_750), mean(Bmu_1000), mean(Bmu_1250))
V_tot_mu<-c(var(mu_hat_250),var(mu_hat_500),var(mu_hat_750),var(mu_hat_1000),var(mu_hat_1250))


# on calcul le biais des estimateurs de moment (beta.choosen = 3 et mu.choosen = 2)
B_tot_mom_mu <-c(mean(mu_mom_250 -2 ),mean(mu_mom_500 -2 ),mean(mu_mom_750 -2 ),mean(mu_mom_1000 -2 ),mean(mu_mom_1250 -2 ))
B_tot_mom_beta <-c(mean(beta_mom_250 -3 ),mean(beta_mom_500 -3 ),mean(beta_mom_750 -3 ),mean(beta_mom_1000 -3 ),mean(beta_mom_1250 -3 ))

# on calcul la MSE des estimateurs de moment
MSE_tot_mom_mu <-c(mean((mu_mom_250 -2)^2 ),mean((mu_mom_500 -2)^2 ),mean((mu_mom_750 -2)^2 ),mean((mu_mom_1000 -2)^2 ),mean((mu_mom_1250 -2 )^2))
MSE_tot_mom_beta <-c(mean((beta_mom_250 -3)^2 ),mean((beta_mom_500 -3)^2 ),mean((beta_mom_750 -3)^2 ),mean((beta_mom_1000 -3)^2 ),mean((beta_mom_1250 -3)^2 ))

# on calcul la variance des estimateurs de moment
Var_tot_mom_mu <-c(var(mu_mom_250),var(mu_mom_500),var(mu_mom_750),var(mu_mom_1000),var(mu_mom_1250))
Var_tot_mom_beta <-c(var(beta_mom_250),var(beta_mom_500),var(beta_mom_750),var(beta_mom_1000),var(beta_mom_1250))


# On prend la valeur absolue des biais pour une meilleure analyse
B_tot_mu <- abs(B_tot_mu)
B_tot <-abs(B_tot)
B_tot_mom_mu<-abs(B_tot_mom_mu)
B_tot_mom_beta<-abs(B_tot_mom_beta)



# Construction des graphiques
abcis <- c(250, 500, 750, 1000, 1250)


# Fonction pour tracer les courbes
trace_courbes <- function(x, y1, y2, ylab, main) {
  
  # Définir les limites de l'axe y pour inclure toutes les valeurs
  ylim_range <- range(c(y1, y2), na.rm = TRUE)
  
  # Tracer la première courbe
  plot(x, y1, type = "b", col = "blue", ylim = ylim_range,
       xlab = "Taille de l'échantillon", ylab = ylab, main = main)
  
  # Ajouter la deuxième courbe
  lines(x, y2, type = "b", col = "red")
  
  # Ajouter une légende
  legend("topright", legend = c("maximum de vraisemblance", "méthode des moments"),
         col = c("blue", "red"), lty = 1, pch = 1)
}

# Tracé des biais pour beta
trace_courbes(abcis, B_tot, B_tot_mom_beta, "Biais", "Biais de l'estimateur de beta\nen fonction de la taille de l'échantillon")

# Tracé des biais pour mu
trace_courbes(abcis, B_tot_mu, B_tot_mom_mu, "Biais", "Biais de l'estimateur de mu\nen fonction de la taille de l'échantillon")

# Tracé de la variance pour beta
trace_courbes(abcis, V_tot, Var_tot_mom_beta, "Variance", "Variance de l'estimateur de beta\nen fonction de la taille de l'échantillon")

# Tracé de la variance pour mu
trace_courbes(abcis, V_tot_mu, Var_tot_mom_mu, "Variance", "Variance de l'estimateur de mu\nen fonction de la taille de l'échantillon")

# Tracé des MSE pour beta
trace_courbes(abcis, MSE_tot, MSE_tot_mom_beta, "MSE", "MSE de l'estimateur de beta\nen fonction de la taille de l'échantillon")

# Tracé des MSE pour mu
trace_courbes(abcis, MSE_tot_mu, MSE_tot_mom_mu, "MSE", "MSE de l'estimateur de mu\nen fonction de la taille de l'échantillon")


# PARTIE 2 : Régression #######################################################



# a) ---------------------------------------------------------------


#Importation des données "data_project2024"
dataproject <- read.csv("data_project2024.txt", sep=",")

#si True alors on transforme Y en log(y-2), si False on ne change rien
toImprove <- F

#Données séparées et leur moyenne
Xi <- dataproject$X
Yi <- dataproject$Y

if(toImprove){Yi <- log(Yi-2)}

X_barre <- mean(Xi)
Y_barre <- mean(Yi)

#Calculs des Sxx Syy et Sxy
Sxx <- sum((Xi - X_barre)^2)
Sxy <- sum((Yi - Y_barre)*(Xi - X_barre))
Syy <- sum((Yi - Y_barre)^2)

#Calculs des coefficients de la régression linéaire
beta_1_chapeau <- Sxy / Sxx
beta_0_chapeau <- Y_barre - beta_1_chapeau * X_barre

#Affichage des résultats
beta_1_chapeau
beta_0_chapeau

#Vérification des résultats
mydata <- data.frame(Xi, Yi)
coefficients_XY <- lm(Yi ~ Xi, data = mydata)
coefficients_XY

#Affichage du graphique de la régression linéaire
plot(x = Xi, y = Yi,
     ylab = "Temps maximal entre 2 pannes (en mois)",
     xlab = "Temps moyen d'utilisation par jour (en heures)",
     pch = 20)

#Ajout de la ligne de régression
abline(a = beta_0_chapeau, b = beta_1_chapeau, lty = 1, col = "red", lwd = 2)


# b) ---------------------------------------------------------------


#On réitère les opérations avec W
Wi <- dataproject$W
W_barre <- mean(Wi)

Sxw <- sum((Wi - W_barre)*(Xi - X_barre))
Sww <- sum((Wi - W_barre)^2)

#Calculs des coefficients de la régression linéaire (W en fonction de X)
meilleur.beta.1.chapeau <- Sxw / Sxx
meilleur.beta.0.chapeau <- W_barre - meilleur.beta.1.chapeau * X_barre

#Affichage des résultats
meilleur.beta.0.chapeau
meilleur.beta.1.chapeau

#Vérification des résultats
mydata2 <- data.frame(Xi, Wi)
coefficients_XW <- lm(Wi ~ Xi, data = mydata2)
coefficients_XW

#Affichage du graphique de la régression linéaire
plot(x = Xi, y = Wi,
     ylab = "Fonction du temps maximal entre 2 pannes",
     xlab = "Temps moyen d'utilisation par jour (en heures)",
     pch = 20) 

#Ajout de la ligne de régression
abline(a = meilleur.beta.0.chapeau, b = meilleur.beta.1.chapeau, lty = 1, col = "red", lwd = 2)


# c) ---------------------------------------------------------------


#calcul de la p valeur pour le testing

#On calcule sigma chapeau
sigma.chapeau.carre <- (Sww - meilleur.beta.1.chapeau * Sxw)/(250 - 2)
sigma.chapeau <- sqrt(sigma.chapeau.carre)

#On calcule maintenant T.obs
T.obs <- (meilleur.beta.1.chapeau/sigma.chapeau) * sqrt(Sxx)
abs.T.obs <- abs(T.obs)

#On peut maintenant calculer la probabilité que T.0 soit supérieure à T.obs
p.value <- 2 * pt(abs.T.obs, df = 248, lower.tail = FALSE ) 
p.value


# d) ---------------------------------------------------------------


#Définition de la fonction extrapolant le temps max entre 2 pannes en fonction du nbr moyen d'heures d'utilisation
extrapolation <- function(nbr_heures) {meilleur.beta.0.chapeau + meilleur.beta.1.chapeau*nbr_heures}
extrapolation(5)

#Calculs des bornes inférieure et supérieure de l'intervalle
borne.gauche <- function(alpha, nbr_heures) {
  extrapolation(nbr_heures) - (qt(p=alpha/2, df=248, lower.tail = FALSE)*sigma.chapeau*sqrt(1 + 1/250 + (((nbr_heures - X_barre)^2))/Sxx))
}
borne.droite <- function(alpha, nbr_heures) {
  extrapolation(nbr_heures) + (qt(p=alpha/2, df=248, lower.tail = FALSE)*sigma.chapeau*sqrt(1 + 1/250 + (((nbr_heures - X_barre)^2))/Sxx))
}

borne.gauche(0.05, 5)
borne.droite(0.05, 5)


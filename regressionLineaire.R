#Partie 2: Régression ----

#(a)----------------------

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



#(b)----------------------

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


#(c)----------------------


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



#(d)------------------------


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





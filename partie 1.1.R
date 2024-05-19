# Partie 1.1 exercice f

# On charge les données
dataproject <- read.csv("data_project2024.txt", sep=",")

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
param<-c(3,5)  

result <- optim(par = param, fn = neg_log_likelihood, method = "BFGS")
opt <- result$par
opt

# Vérification des estimateurs avec les équations données
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

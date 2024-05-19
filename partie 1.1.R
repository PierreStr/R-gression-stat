# Partie 1.1 exercice f

# On charge les données
dataproject <- read.csv("data_project2024.txt", sep=",")

Yi <- dataproject$Y

# Fonction de vraissemblance 

log_likelihood<-function(param){

  beta <- param[1]
  
  if (beta <= 0) {
    print("m'enfin")
    return(0)
  }
  
  mu <- param[2]
  n<-length(Yi)
  
  resultat <- -n*log(beta)
  
  for (i in 1:n) {
    tmp<- (mu - Yi[i])/beta
    resultat <- resultat + exp(tmp-exp(tmp))
  }
  
  return (-resultat)}


param<-c(0.00000001,6000000000)

x<-log_likelihood(param)  
x

lower_bounds <- c(0.00000002, -Inf)

result <- optim(par = param, fn = log_likelihood, method = "L-BFGS-B", lower = lower_bounds)
opt <- result$par
opt
result






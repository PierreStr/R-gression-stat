# Partie 1.1 exercice f

# On charge les données
dataproject <- read.csv("data_project2024.txt", sep=",")

Yi <- dataproject$Y

# Fonction de vraissemblance 

likelyhood<-function(param, Yi){
  n<-length(Yi)
  result<-1
  beta<-param[1]
  mu<-param[2]
  for (i in 1:n) {
    result<-(1/beta)*result*exp(((mu -Yi[i])/beta)-exp((mu -Yi[i])/beta))
  }
  return (-result)}


param<- c(1,2)  
x<-likelyhood(param,Yi)  
print(x)

opt.result <- optim(par = param, fn = likelyhood, Yi = Yi,method = "BFGS")
opt.params <- opt.result$par
opt.params
print(opt.result$par)


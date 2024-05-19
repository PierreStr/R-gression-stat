# Partie 1.1 exercice f

# On charge les données
dataproject <- read.csv("data_project2024.txt", sep=",")

Yi <- dataproject$Y

# Fonction de vraissemblance 

likelyhood<-function(param){
  
  result<-1
  beta <- param[1]
  mu <- param[2]
  dataproject <- read.csv("data_project2024.txt", sep=",")
  
  Yi <- dataproject$Y
  n<-length(Yi)

  
  for (i in 1:n) {
    tmp<-(mu -Yi[i])/beta
    result<-(1/beta)*result*exp(tmp-exp(tmp))
    print(result)
  }
  
  return (-result)}


param<-c(1,0) 
x<-likelyhood(param)  
cprint(x)

result <- optim(par = param, fn = likelyhood,method = "BFGS")
opt <- result$par
print(opt)
opt.params
print(opt.result$par)


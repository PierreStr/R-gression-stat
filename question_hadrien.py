import numpy as np
import pylab as pl
from scipy import stats as sts
from matplotlib import pyplot as plt
from scipy import optimize


l_250 = np.zeros(250)
l_500 = np.zeros(500)
l_750 = np.zeros(750)
l_1000 = np.zeros(1000)
l_1250 = np.zeros(1250)
def initialisation(l_250,l_500,l_750,l_1000,l_1250):
    for i in range(250):
        l_250[i] = np.random.gumbel(2,3)
    for i in range(500):
        l_500[i] = np.random.gumbel(2,3)
    for i in range(750):
        l_750[i] = np.random.gumbel(2,3)
    for i in range(1000):
        l_1000[i] = np.random.gumbel(2,3)
    for i in range(1250):
        l_1250[i] = np.random.gumbel(2,3)
    return l_250,l_500,l_750,l_1000,l_1250
def log_vraissemb(x,L):
    n = len(L)
    beta = x[0]
    mu = x[1]
    result = 0
    for i in range(n):
        tmp = (mu-L[i])/beta
        result+= np.log(1/beta)+tmp-np.exp(tmp)
    return -result
M = 100

# MSE pour beta
MSE_250 = np.zeros(M)
MSE_500 = np.zeros(M)
MSE_750 = np.zeros(M)
MSE_1000 = np.zeros(M)
MSE_1250 = np.zeros(M)


#biais pour beta
B_250 = np.zeros(M)
B_500 = np.zeros(M)
B_750 = np.zeros(M)
B_1000 = np.zeros(M)
B_1250 = np.zeros(M)


#variance de beta
V_250 = np.zeros(M)
V_500 = np.zeros(M)
V_750 = np.zeros(M)
V_1000 = np.zeros(M)
V_1250 = np.zeros(M)

l_250,l_500,l_750,l_1000,l_1250=initialisation(l_250,l_500,l_750,l_1000,l_1250)
for i in range(M):

    l_250,l_500,l_750,l_1000,l_1250=initialisation(l_250,l_500,l_750,l_1000,l_1250)
    res_250 = optimize.minimize(log_vraissemb,[4,3],args=l_250,method="L-BFGS-B",bounds=[(0.1,5),(0,5)])
    res_500 = optimize.minimize(log_vraissemb, [4, 3], args=l_500, method="L-BFGS-B",bounds=[(0.1,5),(0,5)])
    res_750 = optimize.minimize(log_vraissemb, [4, 3], args=l_750, method="L-BFGS-B",bounds=[(0.1,5),(0,5)])
    res_1000 = optimize.minimize(log_vraissemb, [4, 3], args=l_1000, method="L-BFGS-B",bounds=[(0.1,5),(0,5)])
    res_1250 = optimize.minimize(log_vraissemb, [4, 3], args=l_1250, method="L-BFGS-B",bounds=[(0.1,5),(0,5)])

    MSE_250[i] = (res_250.x[0]-3)**2
    MSE_500[i] = (res_500.x[0] - 3) ** 2
    MSE_750[i] = (res_750.x[0] - 3) ** 2
    MSE_1000[i] = (res_1000.x[0] - 3) ** 2
    MSE_1250[i] = (res_1250.x[0] - 3) ** 2

    B_250[i] = (res_250.x[0] - 3) ** 2
    B_500[i] = (res_500.x[0] - 3) ** 2
    B_750[i] = (res_750.x[0] - 3) ** 2
    B_1000[i] = (res_1000.x[0] - 3) ** 2
    B_1250[i] = (res_1250.x[0] - 3) ** 2

    print("M ========== ",i)
abcis = [250,500,750,1000,1250]
print(MSE_250)
MSE_tot =[np.mean(MSE_250),np.mean(MSE_500),np.mean(MSE_750),np.mean(MSE_1000),np.mean(MSE_1250)]
B_tot = [np.mean(B_250),np.mean(B_500),np.mean(B_750),np.mean(B_1000),np.mean(B_1250)]
pl.plot(abcis,MSE_tot)

plt.show()

plt.plot(abcis,B_tot)
plt.show()

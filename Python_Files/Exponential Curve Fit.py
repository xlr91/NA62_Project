import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

def exponential(x, lamb, A, B, C):
    return A*(1-np.exp(-lamb*(x-C)))+B

def f(x, A, B, C):
    return -A/(x-B) + C 
def g(x, a0, a1, a2, a3):
    return a0 + a1/x + a2/(x**2) + a3/(x**3)

def rich(x, F, n, m, e):
    return F*np.sqrt(2-(2/n)*(np.sqrt(1+(m**2/x**2)))) + e

#Rich = np.vectorize(rich)

def exponential_curve_fit(x, y):  
    popt = curve_fit(rich, x, y)[0]
    xnew = np.linspace(np.min(x), np.max(x), 10000)
    ynew = rich(xnew, *popt)
    fig = plt.figure(figsize=[10, 5])
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xnew, ynew, 'b', lw=1)
    ax.scatter(x, y, s=10, c='b')
    ax.set_ylim([0, 250])
    ax.set_xlim([0, 70000])
    ax.set_title('Radius of Ring vs Particle Momentum (MC)')
    plt.show()
    return popt

#‎⁨Macintosh HD⁩ ▸ ⁨Users⁩ ▸ ⁨emirmuhammad⁩ ▸ ⁨Desktop⁩ ▸ ⁨NA62Files⁩ ▸ ⁨Python_Files⁩
namep = "MCARingp.txt"
namer = "MCARingr.txt"
fp = open(namep, "r")
fr = open(namer, "r")

pdata = []
rdata = []


for x in fp:
    pdata.append(float(x[:-1]))
for y in fr:
    rdata.append(float(y[:-1]))



fp.close()
fr.close()


while min(rdata) == 0.0:
    I = rdata.index(min(rdata))
    del pdata[I]
    del rdata[I]
while max(rdata) > 250:
    I = rdata.index(max(rdata))
    del pdata[I]
    del rdata[I]  


p = np.asarray(pdata)
r = np.asarray(rdata)
#exponential_curve_fit(p, r)

#px = np.linspace()



xt = np.linspace(15000, 35000, 1000)
#yt = rich(xt, 1.77027651e+04, 1.00005792e+00, 1.34279662e+02)
yt2 = rich(xt, 17000, 1.000061, 139.570,8)
yt1 = rich(xt, 17000, 1.000061, 139.570,-2)



#plt.ylim(0, 250)
plt.ylabel('Radius (mm)')
plt.xlabel('Momentum (MeV)')
plt.title("Radius of Ring vs Particle Momentum (MC)")
plt.plot(xt, yt2, 'r', label='r(p, 8)')
plt.scatter(pdata, rdata, s=1, label='Events',)
plt.plot(xt, yt1, color='k', label='r(p, -2)')
plt.legend()
plt.show()







if __name__ != '__main__':
    xt = np.linspace(0, 10, 100)
    yt = exponential(xt, 2, 3, 4)+1
    exponential_curve_fit(xt, yt)


if __name__ != '__main__':
    
    f = open("MCARingp.txt", "r")
    #for x in f:
        #print(f.readline())
    
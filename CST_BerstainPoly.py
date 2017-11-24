# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 10:02:30 2016

@author: YinCheang
"""
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import leastsq
import pandasNACA0012 as naca
import math
import os

class CST_Berstain():
    N1  =   0
    N2  =   0
    Ve  =   0
    def __init__(self, N1, N2):
        self.N1     =   N1
        self.N2     =   N2
        
    def CFunc(self,x):
        C   =   (x**self.N1)*(1.0-x)**self.N2
        return C
    
    def SFunc(self,x,n,r):
#        ns  =   np.arange(1,n+1)
#        ns  =   np.prod(ns)
#        rs  =   np.arange(1,r+1)
#        rs  =   np.prod(rs)
#        ks  =   np.arange(1,(n-r)+1)
#        ks  =   np.prod(ks)
#        K   =   ns/(rs*ks)
        n   =   n-1
        K   =   math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        S   =   K*(x**r)*(1.0-x)**(n-r)
        return S
    
    def Dime(self,n,x):
        S1  =   []
        for r in range(n):
            S1.append(self.CFunc(x)*self.SFunc(x,n,r))
        return S1
    
    def solve(self,n,x,z):
        n   =   n+1
        S2  =   []
        b   =   []
        zg  =   0.001260
        for r in range(n):
            S2.append(self.Dime(n,x[r]))
            b.append(z[r]-x[r]*zg)
        S2  =   np.array(S2)
        b   =   np.array(b)
        self.Ve     =   np.linalg.solve(S2,b)
        return self.Ve
    
    def vecfit(self,n,x,vector):
        n   =   n+1
        z   =   0
        zg  =   0.001260
        for r in range(n):
            z   =  z + vector[r]*self.SFunc(x,n,r)*self.CFunc(x)
        z   =   z+zg
        return z

def UpandLow(x,y,half_length):
    length  =   80
    Upx     =   []
    Upy     =   []
    Lowx    =   []
    Lowy    =   []
    for i in range(46):
        Upx.append(x[i])
        Upy.append(y[i])
    for i in range(81-46):
        Lowx.append(x[46+i])
        Lowy.append(y[46+i])
    return Upx,Upy,Lowx,Lowy
    
def fit_fun(p,x):
    f   =   np.poly1d(p)
    return f(x)
def residuals_func(p, y, x):  
    ret     =   fit_fun(p, x) - y  
    return ret  

def write_file(generation,fitness,x,y):
   # flag    =   os.path.exists('C:\Users\YinCheang\code\StandardStripAnalysis\%f' %generation)
   # if flag == False:
   #     os.makedirs('C:\Users\YinCheang\code\StandardStripAnalysis\%f' %generation)
   # 
    fi      =   file('C:\Users\YinCheang\code\StandardStripAnalysis\%f\%f' %(generation,fitness),'w+')
    li      =   ["%ffoil\n"%fitness]
    fi.writelines(li)
    for i in range(len(x)):
        li = ["     %f    %f\n"%(x[i],y[i])]
        fi.writelines(li)
    fi.close()
    return 1    
    
def write_vec(vector):
    fi  =   file('INPUT.DAT','w+')
    for i in range(len(vector)):
        li  =   ["Vec%f   \t%f\n"%(i,vector[i])]
        fi.writelines(li)
    fi.close()
    return 1

def main():
    n               =   5
    m               =   20
    length          =   81
    temx,temz       =   naca.naca0012(length)        
    Ux,Uz,Lx,Lz     =   UpandLow(temx,temz,((length-1)/2))
    Uz = np.array(Uz)
    Lz = np.array(Lz)
    Uz              =   Uz
    Lz              =   Lz

    
    p_init          =   np.random.randn(m)
    plsq            =   leastsq(residuals_func, p_init, args=(Uz, Ux)) 
    plsq2           =   leastsq(residuals_func, p_init, args=(Lz, Lx))
    
    test            =   CST_Berstain(0.5,1)
    x               =   [(test.N1+i)/(test.N1+test.N2+n) for i in range(n+1)]
    
    z               =   fit_fun(plsq[0], x)
    z2              =   fit_fun(plsq2[0],x)
    
    vector          =   test.solve(n,x,z)
    vector2         =   test.solve(n,x,z2)
    zans            =   []
    zans2           =   []
    
#    plt.plot(x , z,'bo', label='fitted curve') 

#    plt.legend()
    plt.plot(Ux,Uz,'ro') 
    plt.plot(Lx,Lz,'ro')
        
    x               =   Ux    
    for i in range(len(x)):
        zans.append(test.vecfit(n,x[i],vector))
        zans2.append(test.vecfit(n,x[i],vector2))
    
    plt.plot(x,zans,'g', label='fitted curve',linewidth=3.0) 
    plt.plot(x,zans2,'g',linewidth=3.0)
    plt.ylabel(r'$Y$',fontsize=24)
    plt.xlabel(r'$X$',fontsize=24) 
    plt.legend()
    plt.show()  
        
    write_vec(np.concatenate((vector,vector2)))
    print vector
    print vector2

    return 1
    
def run_main(vector):
    n       =       5
    pointn  =       81
    x       =       []
    z       =       []
    tempx   =       np.linspace(0,1,pointn) 
    fit     =       CST_Berstain(0.5,1)    
    
    for i in range(pointn-1,-1,-1):
        x.append(tempx[i])
        z.append(fit.vecfit(n,tempx[i],vector[0:6]))
    for i in range(1,pointn):
        x.append(tempx[i])
        z.append(fit.vecfit(n,tempx[i],vector[6:12]))
        
    plt.plot(x,z,'-*')
#    write_file(generation,num,x,z)
    return 1

if __name__ == '__main__':
#    vector = [0.28068683, 0.44616852,  0.28168954,  0.35790707,  0.29068174,  0.6281203,-0.09617941, 0.03389814,0.09662825,0.20311073,0.19948541,  0.51931562]
#    run_main(vector)
    main()
    #[68.47304513   11.40601003   31.39704246   63.19292648  19.37812644 206.52604391]
    #[-26.33561019  -4.38944935 -12.07369813 -24.31427484   7.45956332 -79.46003125]
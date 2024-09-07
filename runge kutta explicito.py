# -*- coding: utf-8 -*-
from pylab import *

def runge_kutta_explicito(f,t0,tf,N,y0,A,b,c):
    h = (tf - t0)/float(N) # Paso de malla
    y = zeros(N+1) # Aproximaciond de la solucion
    t = zeros(N+1) # Paso temporal
    y[0] = y0 # Condicion inicial
    
    s = len(b)
    
    for k in range(0,N):
        t[k+1] = t[k] + h        
        ti = t[k] + c*h
        
        yi = zeros(s)
        yi[0] = y[k]
        
        for j in range(1,s):
            yi[j] = y[k] + h*A[j,:j]@f(ti[:j],yi[:j])
            
        y[k+1] = y[k] + h*b@f(ti,yi)
        
    return (t,y)
 
def f(t,y):
    return y


A = array([[0,0,0,0],
          [0.5,0,0,0],
          [0,0.5,0,0],
          [0,0,1,0]])
b = array([1/6,1/3,1/3,1/6])
c = array([0,1/2,1/2,1])

t0 = 0
tf = 2
N = 50
y0 = 1

# (t,y) = runge_kutta_explicito(f,t0,tf,N,y0,A,b,c)

def runge_kutta_explicito_sistema(f,t0,tf,N,y0,A,b,c):
    h = (tf - t0)/float(N) # Paso de malla
    m = len(y0)
    y = zeros([m,N+1]) # Aproximaciond de la solucion
    t = zeros(N+1) # Paso temporal
    y[:,0] = y0 # Condicion inicial
    
    s = len(b)
    
    for k in range(0,N):
        t[k+1] = t[k] + h        
        ti = t[k] + c*h
        
        yi = zeros([m,s])
        yi[:,0] = y[:,k]
        
        for j in range(1,s):
            print(A[j,:j])
            print(yi[:,:j])
            yi[:,j] = y[:,k] + h*A[j,:j]@f(ti[:j],yi[:,:j])
            
        y[:,k+1] = y[:,k] + h*b@f(ti,yi)
        
    return (t,y)
 
def f(t,y):
    return array([y,2*y])


A = array([[0,0,0,0],
          [0.5,0,0,0],
          [0,0.5,0,0],
          [0,0,1,0]])
b = array([1/6,1/3,1/3,1/6])
c = array([0,1/2,1/2,1])

t0 = 0
tf = 2
N = 50
y0 = array([1])

(t,y) = runge_kutta_explicito(f,t0,tf,N,y0,A,b,c)





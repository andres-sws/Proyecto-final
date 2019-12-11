#!/usr/bin/python

#Libraries to be used
import numpy as np#To use np.exp() and np.linspace()
from scipy.integrate import odeint#to use odeint()
import matplotlib.pyplot as plt#to use plt.plot() and plt.show()

#Numerical constants
d_default = 8.64
mu_1_default = 4.95*1e2
mu_2_default = 4.95*1e-2
v_s_default = 0.12
v_d_default = 1.23
w_default = 0.001
k_1_default = 2.19*1e-4
k_2_default = 6.12*1e-5
k_3_default = 0.997148
k_4_default = 6.79*1e-2
c1 = 2.0
c2 = 10.5
c3 = 2.7
t1 = 1988
t2 = 2100
t3 = 2265
s1 = 21
s2 = 96
s3 = 57

#Initial values for the variables
p_initial = 1.00
sigma_s_initial = 2.01
sigma_d_initial = 2.23
alpha_s_initial = 2.20
alpha_d_initial = 2.26
t_initial = 1000
t_final = 5000
#Functions and equations of the model

def f(t):
    #An function that approximates the CO2 emission between the years 2000 and 5000
    return c1*np.exp(-((t-t1)**2)/(s1**2))+c2*np.exp(-((t-t2)**2)/(s2**2))+c3*np.exp(-((t-t3)**2)/(s3**2))

def p_t(p_s,p,t,funky = f,d = d_default,mu_1 = mu_1_default):
    #The time variation of the partial pressure of CO_2 in the atmosphere
    return (p_s-p)/d + funky(t)/mu_1 

def sigma_s_t(sigma_d,sigma_s,p_s,p,v_s = v_s_default, w = w_default,k_1 = k_1_default,mu_2 = mu_2_default,d = d_default):
    #The time variation of the C concentration in shallow oceans.
    return (1/v_s)*(w*(sigma_d-sigma_s)-k_1-mu_2*(p_s-p)/d)

def sigma_d_t(sigma_d,sigma_s,v_d = v_d_default,k_1 = k_1_default,w = w_default):
    #The time variation of the C concentration in deep oceans.
    return (1/v_d)*(k_1-w*(sigma_d-sigma_s))

def alpha_s_t(alpha_d,alpha_s,v_s = v_s_default,k_2 = k_2_default,w = w_default):
    #The time variation of the alkalinity of shallow oceans
    return (1/v_s)*(w*(alpha_d-alpha_s)-k_2)

def alpha_d_t(alpha_d,alpha_s,v_d = v_d_default,k_2 = k_2_default,w = w_default):
    #The time variation of the alkalinity of deep oceans
    return (1/v_d)*(k_2-w*(alpha_d-alpha_s))

#Complementary ocean equilibrium equations
def h_s(alpha_s,sigma_s,k_3 = k_3_default):
    #Complementary shallow ocean equilibrium equation h_s
    return (sigma_s-((sigma_s**2)-k_3*alpha_s*(2*sigma_s-alpha_s))**(0.5))/k_3

def c_s(alpha_s,h_s):
    #Complementary shallow ocean equilibrium equation c_s
    return (alpha_s-h_s)/2

def p_s(h_s,c_s,k_4 = k_4_default):
    #Complementary shallow ocean equilibrium equation p_s
    return k_4*((h_s**2)/c_s)

#Model compacted as a function
def model(state,t):
    #This is a function that allows using odeint for the system
    #"state" is a list with the initial conditions of the system ordered as below, t is a number.
    p_status = state[0]
    sigma_s_status = state[1]
    sigma_d_status = state[2]
    alpha_s_status = state[3]
    alpha_d_status = state[4]

    h_s_status = h_s(alpha_s_status,sigma_s_status)
    c_s_status = c_s(alpha_s_status,h_s_status)
    p_s_status = p_s(h_s_status,c_s_status)

    p_t_status = p_t(p_s_status,p_status,t)
    sigma_s_t_status = sigma_s_t(sigma_d_status,sigma_s_status,p_s_status,p_status)
    sigma_d_t_status = sigma_d_t(sigma_d_status,sigma_s_status)
    alpha_s_t_status = alpha_s_t(alpha_d_status,alpha_s_status)
    alpha_d_t_status = alpha_d_t(alpha_d_status,alpha_s_status)

    #returns a vector of the time derivatives of each variable
    return [p_t_status,sigma_s_t_status,sigma_d_t_status,alpha_s_status,alpha_d_t_status]

#Initial state
state_initial = [p_initial,sigma_s_initial,sigma_d_initial,alpha_s_initial,alpha_d_initial]

#Time interval
t_range = np.linspace(t_initial,t_final,t_final-t_initial)

#First graph
plt.plot(t_range,f(t_range))
plt.show()

#solving the system
sol = odeint(model,state_initial,t_range)

p_sol = sol[:,0]
sigma_s_sol = sol[:,1]
sigma_d_sol = sol[:,2]
alpha_s_sol = sol[:,3]
alpha_d_sol = sol[:,4]

#NOT WORKING, GOES TO INFINITY
plt.plot(t_range,p_sol[:t_final])
plt.show()

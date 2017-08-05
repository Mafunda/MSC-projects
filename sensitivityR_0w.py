from __future__ import division
from sympy import*
#import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
#from matplotlib import pyplot as plt
from fractions import Fraction
import csv
#symbols

N_v,alpha_v,alpha_c,alpha_h,Lambda_s, beta_h, mu_s,mu_v,mu_c,mu_h,Lambda_H,beta_H,mu_H,alpha_H,delta_H,V_0 = symbols('N_v alpha_v alpha_c alpha_h Lambda_s beta_h mu_s mu_v mu_c mu_h Lambda_H beta_H mu_H alpha_H delta_H V_0')

#R_0w = (N_v*alpha_v*alpha_c*Lambda_s*beta_h)/(mu_s*(alpha_v + mu_v)*(alpha_c + mu_c)*(alpha_h + mu_h));

#V_h = (mu_s/beta_h)*(R_0w -1)

#R_0b = (Lambda_H*beta_H*alpha_h*mu_s*(R_0w-1))/(beta_h*mu_H*alpha_H*V_0*(mu_H + delta_H))

#V_H = (mu_H*V_0*(R_0b-1))/(beta_H + mu_H) 


A = [N_v,alpha_v,alpha_c,alpha_h,Lambda_s, beta_h, mu_s,mu_v,mu_c,mu_h,Lambda_H,beta_H,mu_H,alpha_H,delta_H,V_0]

J = [650,0.09,0.00045,0.09,50000,0.00135,0.02,0.15,0.02,1.5,0.034,0.03,0.000002,0.24,0.000006,2e8]

A_0 = zip(A,J)

R_0w = (N_v*alpha_v*alpha_c*Lambda_s*beta_h)/(mu_s*(alpha_v + mu_v)*(alpha_c + mu_c)*(alpha_h + mu_h));

V_h = (mu_s/beta_h)*(((N_v*alpha_v*alpha_c*Lambda_s*beta_h)/(mu_s*(alpha_v + mu_v)*(alpha_c + mu_c)*(alpha_h + mu_h)))-1);

R_0b = (Lambda_H*beta_H*alpha_h*mu_s*((N_v*alpha_v*alpha_c*Lambda_s*beta_h)/(mu_s*(alpha_v + mu_v)*(alpha_c + mu_c)*(alpha_h + mu_h))-1))/(beta_h*mu_H*alpha_H*V_0*(mu_H + delta_H));

V_H = (mu_H*V_0*(R_0b-1))/(beta_H + mu_H) 

R_0w0 = []
V_h0 = []
R_0b0 = []
V_H0 = []

for i in A:
	x35 = zip(A,J)
	y_1 = diff(R_0w, i)*(i/R_0w)
	A1 = y_1.subs(x35)
	R_0w0.append(A1)

print R_0w0;

for i in A:
	x35 = zip(A,J)
	y_1 = diff(V_h, i)*(i/V_h)
	A1 = y_1.subs(x35)
	V_h0.append(A1)

print V_h0;

for i in A:
	x35 = zip(A,J)
	y_1 = diff(R_0b, i)*(i/R_0b)
	A1 = y_1.subs(x35)
	R_0b0.append(A1)

print R_0b0;

for i in A:
	x35 = zip(A,J)
	y_1 = diff(V_H, i)*(i/V_H)
	A1 = y_1.subs(x35)
	V_H0.append(A1)

print V_H0;

B_sy = ['$\\N_v$','$\\alpha_v$','$\\alpha_c$','$\\alpha_h$','$\\Lambda_s$', '$\\beta_h$', '$\\mu_s$','$\\mu_v$','$\\mu_c$','$\\mu_h$','$\\Lambda_H$','$\\beta_H$','$\\mu_H$','$\\alpha_H$','$\\delta_H$','$V_0$']
	
with open("SensitivityAnalysis.csv", "w") as out_file:
	writer = csv.DictWriter(out_file, fieldnames = ['Number','&','Parameter','&','$R_{0w}$','&', '$V_{h}$','&','$R_{0b}$','&', '$V_{H}$','\\'])
    	writer.writeheader()
	for i in range(len(B_sy)):
		out_string = ""
		out_string += str(i+1)
		out_string += "," +str("&")
		out_string += "," +str(B_sy[i])
		out_string += "," +str("&")
		out_string += "," +str(R_0w0[i])
		out_string += "," +str("&")
		out_string += "," +str(V_h0[i])
		out_string += "," +str("&")
		out_string += "," +str(R_0b0[i])
		out_string += "," +str("&")
		out_string += "," +str(V_H0[i])
		out_string += "," +str("\\")
		out_string += "\n"
		out_file.write(out_string)

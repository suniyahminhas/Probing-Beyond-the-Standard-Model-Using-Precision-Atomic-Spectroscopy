from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot
from scipy.misc import imread
from PIL import Image, ImageDraw



#attempt to make a constraint plot
h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge
M=0.5109989461 #mass of an electron Mev

r=Symbol("r", real=True, positive=True)
m_b_values= numpy.logspace(-4,6, num=30)

def pert_for_function((N,L,Z)):
    #this is what would be in the integral (i,e functions to integrate)
    data_matrix = numpy.zeros((len(m_b_values),1))
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    return data_matrix
    
def constraint(a,b,error_bar):
    k_new= pert_for_function(a)-pert_for_function(b)
    error_E_bar= error_bar*h/(e*k)
    return error_E_bar/k_new


pyplot.figure()
img = imread("fade.png")
pyplot.imshow(img,zorder=0, extent= (10**-6, 10**3, 10**-16, 1), aspect='auto')
pyplot.show()
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),(constraint((5,4,79),(15,5,79),10000))/197, color='darkgreen', linewidth=4.0, linestyle='--') 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),(constraint((4,1,79),(12,3,79),400000))/197, color='blue', linewidth=4.0, linestyle='--') 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),(constraint((7,0,79),(16,1,79),2*10**6))/197, color='purple', linewidth=4.0, linestyle='--') 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,0,1),56), color='red', linewidth=4.0) 

#pyplot.legend(loc=2, fontsize = 20)
pyplot.xlabel('$m_\phi (MeV)$', fontsize=40)
#pyplot.yticks(numpy.array[(10**-17, 10**-14, 10**-8, 10**-5)])
pyplot.ylabel('$g_{e}g_{n}$', fontsize=40)
pyplot.rc('xtick', labelsize=30) 
pyplot.rc('ytick', labelsize=30) 
pyplot.yscale('log')
pyplot.ylim(10**-16, 1)
pyplot.xlim(10**-6, 10**3)
pyplot.xscale('log')
pyplot.xticks([10**-6,10**-3,1,10**3], fontsize=30)
pyplot.yticks([10**0,10**-4,10**-8, 10**-12, 10**-16], fontsize=30)

pyplot.subplots_adjust(hspace=0.5, wspace=0.5)

pyplot.show()

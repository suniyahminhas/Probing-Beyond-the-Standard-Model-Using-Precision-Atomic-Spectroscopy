from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot
import matplotlib.axis as ax
from scipy.misc import imread
from PIL import Image, ImageDraw
import RadialFunctions as f
R=numpy.arange(0,30,0.3)
r=Symbol("r", real=True, positive=True)
Z=Symbol("Z", real=True, positive=True)

print R_nl(4,3,r,Z)
def R_4_3(r,Z):
    return numpy.sqrt(35)*Z**(9/2)*r**3*numpy.exp(-Z*r/4)/26880
    
print(f.R_43(0))

pyplot.subplot(1,3,1)
pyplot.plot(R, f.R_40(R), linewidth=3.0, label=0, color='green')
pyplot.plot(R, f.R_41(R), linewidth=3.0,linestyle='--', label=1,color='green')
pyplot.plot(R, f.R_42(R), linewidth=3.0,linestyle=':', label=2,color='green')
pyplot.plot(R, f.R_43(R), linewidth=3.0,linestyle='-.', label=3,color='green')
pyplot.xlabel('$r (a_0)$', fontsize=35, labelpad=-5)
pyplot.xticks([0,10,20,30], fontsize=25)
pyplot.ylabel('$R_{4l}$', fontsize=35, labelpad=-5)
pyplot.title('l', fontsize=25)
pyplot.legend(fontsize=20)
pyplot.yticks(fontsize=25)

pyplot.subplot(1,3,2)
pyplot.plot(R, f.R_10(R), linewidth=3.0, label=1,color='purple')
pyplot.plot(R, f.R_20(R), linewidth=3.0, label=2,color='purple',linestyle='--')
pyplot.plot(R, f.R_30(R), linewidth=3.0, label=3, color='purple',linestyle='-.')
pyplot.title('n', fontsize=25)
pyplot.xlabel('$r (a_0)$', fontsize=35, labelpad=-5)
pyplot.ylabel('$R_{n0}$', fontsize=35, labelpad=-5)
pyplot.xticks([0,10,20,30], fontsize=25)
pyplot.legend(fontsize=20)
pyplot.yticks(fontsize=25)


pyplot.subplot(1,3,3)
pyplot.plot(R,  R_4_3(R,1), linewidth=3.0, label=1,color='blue')
pyplot.plot(R, R_4_3(R,2), linewidth=3.0, label=2, color='blue',linestyle='--')
pyplot.plot(R, R_4_3(R,3), linewidth=3.0, label=3, color='blue',linestyle='-.')
pyplot.xlabel('$r (a_0)$', fontsize=35, labelpad=-5)
pyplot.ylabel('$R_{10}$', fontsize=35, labelpad=-5)
pyplot.xticks([0,10,20,30], fontsize=25)
pyplot.legend(fontsize=20)
pyplot.title('Z', fontsize=25)
pyplot.yticks(fontsize=25)

pyplot.subplots_adjust(hspace=0.3, wspace=0.3)

pyplot.show()
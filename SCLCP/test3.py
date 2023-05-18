import math

import mpmath.functions.expintegrals
import numpy as np
# yxb = "I want to fly, whcih means I have wings!"
# list = yxb.split()
# print(list)
# for i in list:
#     if ',' in i:
#         list.remove(i)
#         list.insert(3,i.split(',')[0])
# print(list)

# array1 = np.ones((2,3))
# array2 = np.zeros((3,4))
# array3 = np.array([10])
# array4 = np.array([20])
# array = np.array([array1, array2, array3, array4],dtype=object)
# # print(array[0])
# # print(array[2])
# print(array2+10)
#import math
import scipy.sparse.csgraph
import sympy
from scipy.constants import Boltzmann,elementary_charge,epsilon_0,R, N_A
delta_voltage = 0.009
#print(Boltzmann)
#print(elementary_charge/(Boltzmann*300))
from sympy.abc import x
from sympy import exp
# a = math.exp(-elementary_charge/(Boltzmann*300))  # 1.5875937562011973e-17
# print(a)
#print(b)
# k = 31.6
# v = 13.06
#print(-elementary_charge/(Boltzmann*300)) #-38.681727071833606
# y = sympy.solve(k*exp(a*(x-0.5)/100)+v*exp(a*(x-0.1)/100)-k-v,x)
# print(y)
# y = sympy.solve(exp(-1 * x * 29) - 0.1, x)
# if y > 0:
#     print(y)



# i = 1
# a = float(input("数字："))
# while i <= 3:
#     i += 1
#     if a >= 0:
#         break
#
#     b = i+1
# print(b)
# import input_parameters as ip
# temp = ip.Condition().temperature
# permittivity = 60 * epsilon_0
# phi_0 = 0.484
# phi_bulk = 0
# c_0 = 26.05*math.exp(-elementary_charge/(Boltzmann*temp*10)*(phi_0-phi_bulk))
# print(c_0)
# rou_0 = elementary_charge*c_0-elementary_charge*26.05
# #rou_1 = elementary_charge*c_
# efs = math.sqrt(2*Boltzmann*temp/permittivity*abs(c_0-26.05))
# phi_1 = phi_0 - 0.01*efs - 0.01**2/(2*permittivity)*(rou_0)
# print(phi_1)

# test = 4.9E-4*exp(-0.226*elementary_charge/(Boltzmann*300))
# print(test)
# y = 4%4
# # # print(y)string
# string_b = []
# string_a = ["a","cat","is","eating","."]
# for i in string_a:
#     b = " "+i
#     string_b.append(b)
# for j in string_b:
#     print("{}\n".format(j))

cond = 8.5E-4 #电导率
diff = 2.83E-8 #扩散系数
concen = Boltzmann*300*cond/(elementary_charge**2*diff)*10**(-21) #浓度
# print(concen)

voltage = 0.487
permittivity = 5.743
concentration = 28.685
xy = math.sqrt(2*permittivity*epsilon_0*voltage/(elementary_charge*concentration)*1E-9)
# print(xy)

s1 = Boltzmann*300/elementary_charge*1E-3
s2 = 0.321
s3 = (s1/s2)**(1/500)
#print(s3)

# a = math.exp(-elementary_charge/(Boltzmann*300))  # 1.5875937562011973e-17
#print(-elementary_charge/(10*Boltzmann*300)) #-38.681727071833606 #乘以beta=10 —— -3.8681727071833603

k = 0.1/(math.exp(0.32*elementary_charge/(Boltzmann*300))*0.1+1)
g= Boltzmann*300/elementary_charge
#print(g)

yin = 1/(exp(elementary_charge/(Boltzmann*300))*0.1+1)
bin = 1/(exp(5.2*elementary_charge/(Boltzmann*300)))
guding = elementary_charge/(Boltzmann*300)
print(bin)
# print(guding)
print(k)

yin1 = 1/(math.exp(elementary_charge/(Boltzmann*300)*(0.8))+10)
print(yin1)
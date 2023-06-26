# 用户需要输入的参数
import numpy as np
from scipy.constants import epsilon_0
from math import exp


class Condition():
    temperature = 300
    #delta_c = 0.001
    cm2nm2 = 1E14
    total_step = 1000
    delta_x = 10/total_step
    alpha = 0.1
    beta = 10
'''
体相锂离子浓度单位：/nm3
扩散系数单位：cm2/s
'''


#定义氧化物正极材料的各种性质
def licoo2():
    array = np.ones(5)
    z_i = 1
    c_bulk = 24.53
    surface = 1
    permittivity = 12.9 #2018-Space-charge layers in all-solid-state batteries; important or negligible
    diffusivity = 2E-10 #2018-Space-charge layers in all-solid-state batteries; important or negligible

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    return array


#定义硫化物正极的各种性质
def litis2():
    array = np.ones(5)
    z_i = 1
    c_bulk =1
    surface = 1
    permittivity =1
    diffusivity =1

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    return array


#定义氧化物涂层的性质
def linbo3():
    array = np.ones(5)
    z_i = 1
    c_bulk =1
    surface = 1
    permittivity =1
    diffusivity =1

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    return array


#定义氧化物固态电解质材料的各种性质
def li7la3zr2o12():
    array = np.ones(7)
    z_i = 1
    c_bulk = 26.047 #晶胞
    vacancy_e = 4.388
    #23.445 #26.05
    surface = 1
    permittivity = 5.743
    diffusivity = 4E-9  # 2016-Data mining of molecular dynamics data reveals Li diffusion characteristics in garnet Li7La3Zr2O12
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    #array[6] = c_bulk_minimum
    return array


def lilatio():
    array = np.ones(7)
    z_i = 1
    c_bulk = 5.632 #晶胞
    #22.41 #5.63
    surface = 1
    permittivity = 10.048
    diffusivity = 7.2E-9 #7Li NMR diffusion studies in micrometre-space for perovskite-type Li0. 33La0. 55TiO3 (LLTO) influenced by grain boundaries
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def liocl():
    #离子电导率：
    array = np.ones(7)
    z_i = 1
    c_bulk = 53.603 #晶胞
    #50#11.97 #53.6
    surface = 1
    permittivity = 8.39 #Concentration of Charge Carriers, Migration, and Stability in Li3OCl Solid Electrolytes
    diffusivity = 1.87E-8 #Multi-ion Conduction in Li3OCl Glass Electrolytes
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def litipo():
    #离子电导率：
    array = np.ones(7)
    z_i = 1
    c_bulk = 3.898 #晶胞
    #0.4246
    surface = 1
    permittivity = 7.651
    diffusivity = 6.93E-9 #2003-LiTi2 (PO4) 3 with NASICON-type structure as lithium-storage materials
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def li2po2n():
    #离子电导率：
    array = np.ones(7)
    z_i = 1
    c_bulk = 35.190 #晶胞
    #0.061 #35.19
    surface = 1
    permittivity = 5.601
    diffusivity = 1.7E-8 #Dielectric properties, conductivity and Li+ ion motion in LiPON thin films
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


#定义硫化物固态电解质的各种性质
def li10gep2s12():
    array = np.ones(7)
    z_i = 1
    c_bulk = 21.414 # 晶胞
    #公式推导 24.729
    surface = 1
    permittivity = 8.89
    diffusivity = 7.83E-8 #Molecular dynamics simulations of lithium superionic conductor Li10GeP2S12 using a machine learning potential
    voltage =0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def li3ps4():
    array = np.ones(7)
    z_i = 1
    c_bulk = 19.316 #晶胞
    #13.588 #公式推导 #19.32
    surface = 1
    permittivity = 9.104
    diffusivity = 1.9E-9 #The Journal of Physical Chemistry C, 2019, 123(16): 10280-10290.
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def li6ps5cl():
    array = np.ones(7)
    z_i = 1
    c_bulk = 23.219 #晶胞
    #29.689  #公式推导 #7.218
    surface = 1
    permittivity = 38
    diffusivity = 2.5E-8 #Structure and Diffusion Pathways in Li6PS5Cl Argyrodite from Neutron Diffraction, Pair-Distribution Function Analysis, and NMR
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


#卤化物电解质
def li3incl6():
    array = np.ones(7)
    z_i = 1
    c_bulk = 13.434 #晶胞
    #29.689  #公式推导 #7.218
    surface = 1
    permittivity = 5.876
    diffusivity = 2.177E-9 #Investigation of the effect of F-doping on the solid-electrolyte property of Li3InCl6
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def li3ycl6():
    array = np.ones(7)
    z_i = 1
    c_bulk = 14.104  # 晶胞
    # 29.689  #公式推导 #7.218
    surface = 1
    permittivity = 5.939
    diffusivity = 1E-9
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array
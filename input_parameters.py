# input parameter
import numpy as np
from constants import epsilon_0
from math import exp


class Condition():
    temperature = 300 #Kelvin temperature
    cm2nm2 = 1E14  #Conversion coefficient between cm2 and nm2
    cm3nm3 = 1E21  #Conversion coefficient between cm3 and nm3
    total_step = 10000   #Total calculation step
    cathode = "oco"   #oca, oco, sca
    se = "li3incl6"    #SSE: li7la3zr2o12, lilatio, li3ocl, litipo, li2po2n; li10gep2s12, li3ps4, li6ps5cl; li3incl6, li3ycl6
    se_type = "halide" #SSE type: oxide, sulfide, halide or other
    picture = 0 #If you want to plot concentration versus range, you can set this parameter to 1. The function still needs to be further improved and is not open for the time being.

'''
Li-ion bulk concentration：/nm3
Diffusivity：cm2/s
'''
"""
z_i: the number of charge carriers, such as +1 for lithium-ions and −1 for electrons
c_bulk: the mobile lithium-ion concentration, which can be obtained through the Nernst-Einstein equation
surface: the the cross-sectional area of the interface and is set to 1 nm2
permittivity: the relative permittivity of the SSE, which is calculated by DFT
diffusivity: the diffusivity of lithium-ion in SSE, which can be obtained by DFT calculation or directly from experiments
voltage: bulk electric-potential in SSE, defaults to 0 V
"""

#oxide cathode property
def oca():
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


#sulfide cathode property
def sca():
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


#oxide coating property
def oco():
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


#oxide sse property
def li7la3zr2o12():
    array = np.ones(7)
    z_i = 1
    c_bulk = 10.085
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
    c_bulk = 16.136
    surface = 1
    permittivity = 10.048
    diffusivity = 1.7E-8 #7Li NMR diffusion studies in micrometre-space for perovskite-type Li0. 33La0. 55TiO3 (LLTO) influenced by grain boundaries
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def li3ocl():
    array = np.ones(7)
    z_i = 1
    c_bulk = 8.238
    surface = 1
    permittivity = 8.39 #Concentration of Charge Carriers, Migration, and Stability in Li3OCl Solid Electrolytes
    diffusivity = 1.9E-8 #Multi-ion Conduction in Li3OCl Glass Electrolytes
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def litipo():
    array = np.ones(7)
    z_i = 1
    c_bulk = 23.385
    surface = 1
    permittivity = 7.651
    diffusivity = 6.9E-9 #2003-LiTi2 (PO4) 3 with NASICON-type structure as lithium-storage materials
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


def li2po2n():
    array = np.ones(7)
    z_i = 1
    c_bulk = 6.169
    surface = 1
    permittivity = 5.601
    diffusivity = 1.7E-9 #Dielectric properties, conductivity and Li+ ion motion in LiPON thin films
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


#sulfide sse property
def li10gep2s12():
    array = np.ones(7)
    z_i = 1
    c_bulk = 6.454
    surface = 1
    permittivity = 8.89
    diffusivity = 3E-7 #Molecular dynamics simulations of lithium superionic conductor Li10GeP2S12 using a machine learning potential
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
    c_bulk = 18.441
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
    c_bulk = 22.475
    surface = 1
    permittivity = 7.218
    diffusivity = 2.5E-8 #Structure and Diffusion Pathways in Li6PS5Cl Argyrodite from Neutron Diffraction, Pair-Distribution Function Analysis, and NMR
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage

    return array


#halide sse property
def li3incl6():
    array = np.ones(7)
    z_i = 1
    c_bulk = 57.941
    surface = 1
    permittivity = 5.876
    diffusivity = 2.2E-9 #Investigation of the effect of F-doping on the solid-electrolyte property of Li3InCl6
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
    c_bulk = 82.291
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

#Space-charge layers in all-solid-state batteries; important or negligible?
def llzo():
    array = np.ones(7)
    z_i = 1
    c_bulk = 25.65
    surface = 1
    permittivity = 60
    diffusivity = 4E-9
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage
    return array

def latp():
    array = np.ones(7)
    z_i = 1
    c_bulk = 5.37
    surface = 1
    permittivity = 15
    diffusivity = 3E-9
    voltage = 0

    array[0] = z_i
    array[1] = c_bulk
    array[2] = surface
    array[3] = permittivity * epsilon_0
    array[4] = diffusivity
    array[5] = voltage
    return array
import numpy as np
import pandas as pd
import input_parameters as ip
import matplotlib.pyplot as plt
import sys
import csv
import math
import input_mcfss
from constants import Boltzmann, R, elementary_charge, epsilon_0, N_A, faraday


def get_concentration_li(concentration_bulk, z_i, potential, potential_bulk):
    '''
    The carrier concentration at different locations is solved based on Boltzmann distribution
    :param concentration_bulk: Carrier bulk concentration （nm-3）
    :param z_i: The charge of the carrier
    :param potential: Potential real-time sites, the initial value by | MCFSS |
    :param potential_bulk: The potential for bulk, let's call it 0. How to find the potential of SE.
    :return: Interfacial carrier concentration
    '''
    boltzmann_distribution = 1/(math.exp(z_i*elementary_charge/(Boltzmann*temperature)*(potential - potential_bulk)))
    concentration_x = concentration_bulk*boltzmann_distribution
    return concentration_x


def get_electrical_field_strength(permittivity, concentration_x_last, concentration_bulk, z_i, c_im, potential_bulk, potential):
    '''
    If the initial potential is too small, it may not be possible to open the root.
    :param permittivity: permittivity of SSE
    :param concentration_x_last: The last calculated carrier concentration
    :param concentration_bulk: The carrier concentration in bulk
    :param z_i: The charge number of the carrier
    :param c_im: The concentration of electron
    :param potential_bulk: Electric-potential in bulk
    :param potential: Interfacial electric-potential
    :return: Interfacial electric field intensity
    '''
    electrical_field_strength = sgn*math.sqrt((2*Boltzmann*temperature/permittivity)
                                              *abs((concentration_x_last-concentration_bulk)
                                                +z_i*elementary_charge*c_im/(Boltzmann*temperature) * (potential_bulk- potential)))
    return electrical_field_strength


def get_factor_alpha():
    '''
    Regulatory factor in calculating the electric-potential near the bulk
    :return: factor alpha
    '''
    if cathode == "oca":
        if se_type == "oxide":
            a = 0.05
        elif se_type == "sulfide":
            a = 0.5
        else:
            a = 0.05
    elif cathode == "sca":
        if se_type == "oxide":
            a = 0.05
        elif se_type == "sulfide":
            a = 0.05
        else:
            a = 0.05
    elif cathode == "oco":
        if se_type == "oxide":
            a = 0.05
        elif se_type == "sulfide":
            a = 0.05
        else:
            a = 0.05
    else:
        a = 1
    # if se == "li2po2n":
    #     a = 0.01
    return a


def get_delta_x(charge_density_scl, charge_density):
    '''
    :param charge_density_scl: Total charge density in SCL
    :param charge_density: The charge density of the carrier
    :return: Distance of different calculation steps
    '''
    delta_charge_density_scl = charge_density_scl/total_step
    delta_x = abs(delta_charge_density_scl/charge_density)
    return delta_x


def get_potential(potential_last,delta_x,electrical_field_strength,permittivity,charge_density,delta_charge_density):
    '''
    :param potential_last: The last calculated electric-potential
    :param delta_x: distances for different computation steps
    :param electrical_field_strength: Interfacial electric field intensity
    :param permittivity: Permittivity of SSE
    :param charge_density: Charge_density of carrier
    :param delta_charge_density: distances for different charge_density of carrier
    :return: Electric-potential
    '''
    potential_next = potential_last - delta_x*electrical_field_strength-delta_x**2/(2*permittivity)*\
                     (charge_density+delta_charge_density/3)
    return potential_next


def get_resistance_scl(z_i, surface, diff, delta_x_array, c_array):
    '''
    :param z_i: The charge number of the carrier
    :param surface: Surface area (1 nm2)
    :param diff: Diffusivity of SSE
    :param delta_x_array: An array of distances for different computation steps
    :param delta_c_array: An array of distances for different computation concentration
    :return: Space charge layer resistance
    '''
    r_scl = (Boltzmann*temperature)/(surface*z_i**2*elementary_charge**2*diff)*sum(delta_x_array/c_array)
    return r_scl


if __name__ == '__main__':
    cathode = ip.Condition().cathode #licoo2: OCA; litis2: SCA; linbo3: OCO
    se = ip.Condition().se      #SSE: li7la3zr2o12, lilatio, li3ocl, litipo, li2po2n; li10gep2s12, li3ps4, li6ps5cl; li3incl6, li3ycl6
    se_type = ip.Condition().se_type #SSE type: oxide, sulfide and halide
    print("Input Completion".center(30,'*'))
    print("Cathode：{}".format(cathode))
    print("SEs：{}".format(se))
    print("Start Calculating".center(30, '*'))

    mod = sys.modules["input_parameters"]
    array1 = mod.__dict__[cathode]()  # An array of properties of cathode or coating materials
    array2 = mod.__dict__[se]()  # An array of solid-state-electrolyte materials

    total_step = ip.Condition().total_step
    temperature = ip.Condition().temperature
    cm2nm2 = ip.Condition().cm2nm2
    z_i_cathode = array1[0]
    z_i_se = array2[0]
    z_i_im = -1
    c_bulk_cathode = array1[1] #nm-3
    c_bulk_se = array2[1] #nm-3
    c_e_se = c_bulk_se #nm-3
    surface_cathode = array1[2] #nm2
    surface_se = array2[2] #nm2
    permittivity_cathode = array1[3]*1E-9 #F/nm
    permittivity_se = array2[3]*1E-9 #F/nm
    diffusivity_cathode = array1[4]*1E14 #nm2/s
    diffusivity_se = array2[4]*1E14 #nm2/s
    se_vol = abs(array2[5]) #V

    if cathode not in ["oca", "oco", "sca"]:
        print("cathode error!")
        exit()
    if se_type not in ["oxide", "sulfide", "halide"]:
        print("se type error!")
        exit()

    potential_list = input_mcfss.potential_list

    '''test'''
    #potential_list = [-0.4,-0.5,-0.6,-0.7,-0.8]
    #potential_list = [-0.2, -0.3]
    #potential_list = [2.5, 3]
    del array1, array2

    rou_scl = 0
    resistance_list = []
    count = 0
    output_dic = {}
    con_dic = {}
    phi_end = abs(Boltzmann*temperature/elementary_charge*math.log(1+1E-7))
    for cathode_vol in potential_list:
        phi_0 = abs(cathode_vol)
        count += 1
        step = 0
        if phi_0 > se_vol:
            sgn = 1
        else:
            sgn = -1
        c_list = []
        delta_c_list = []
        delta_x_list = []
        rou_list = []
        efs_list = []
        phi_list = [phi_0]
        len_scl = 0
        len_scl_list = []
        dx = 0
        sym = 0
        tp = 0
        alpha = get_factor_alpha()
        phi_ratio = alpha * (phi_end / phi_list[0]) ** (1 / total_step)
        # deby_length = math.sqrt(permittivity_se*Boltzmann*temperature/(2*elementary_charge**2*c_bulk_se))
        # print(deby_length)
        while True:
            step += 1
            phi_x = phi_list[-1]
            if cathode == "oco":
                phi_x = phi_x*0.01

            c_x_se = get_concentration_li(concentration_bulk=c_bulk_se, z_i = z_i_se, potential=phi_x*alpha, potential_bulk= se_vol)
            c_list.append(c_x_se)
            if step > 1:
                delta_c = abs(c_list[-1] - c_list[-2])
            else:
                delta_c = 0
            delta_c_list.append(delta_c)

            rou_x = (z_i_se * c_x_se + z_i_im * c_e_se) * elementary_charge  #(<0)
            rou_list.append(rou_x)
            if step != 1:
                delta_rou = rou_list[-1] - rou_list[-2]  #(>0)
            else:
                delta_rou = 0

            efs = get_electrical_field_strength(permittivity = permittivity_se, concentration_x_last = c_x_se,
                                                concentration_bulk = c_bulk_se, z_i = z_i_im, c_im = c_e_se,
                                                potential_bulk= se_vol, potential=phi_list[-1])
            efs_list.append(efs)
            if step == 1:
                rou_scl = -efs_list[0]*permittivity_se
            else:
                rou_scl = rou_scl

            if step <= 0.5*total_step:
                dx = get_delta_x(charge_density_scl=rou_scl, charge_density=rou_list[-1])
                phi_x = get_potential(potential_last=phi_list[-1],delta_x=dx,electrical_field_strength=efs_list[-1],
                                      permittivity=permittivity_se,charge_density=rou_list[-1],delta_charge_density=delta_rou/3)
            else:
                factor_b = (phi_end / phi_list[-1]) ** (1 / (total_step - step))
                phi_x = factor_b * phi_list[-1] * alpha
                dx = get_delta_x(charge_density_scl=rou_scl, charge_density=rou_list[-1])

            '''test'''
            # if step == 10:
                # print(phi_x / phi_list[-1])
                # print(phi_ratio)
                # print(phi_list[1]/phi_list[0])

                # print((- 2 / permittivity_se * (phi_x - phi_list[-1] * (rou_list[-1] + delta_rou / 3))))
                # print(efs ** 2)
                # print(permittivity_se / (rou_list[-1] + delta_rou / 3))
                # print(dx)
                # print(-dx ** 2 / (2 * permittivity_se) * (rou_list[-1] + delta_rou / 3))
                # print(-dx * efs)
                # print(phi_x)
                # print((phi_end / phi_list[-1]) ** (1 / (total_step - step)))
            '''test over'''

            delta_x_list.append(dx)
            phi_list.append(phi_x)
            len_scl += abs(dx)
            len_scl_list.append(len_scl)

            '''test'''
            #print(step)
            # #print(phi_end)
            #print(c_list[-1])
            # #print(delta_x_list)
            # #print(rou_x)
            # #print(electrical_field)
            #print(phi_x)
            # # print()
            '''test over'''

            if phi_list[-1] < 1E-7 or step >= total_step-1:
                break
        del delta_c_list[0]

        #del delta_c_list[-1]
        c_se_array = np.array(c_list)
        dx_array = np.array(delta_x_list)
        resistance_se = get_resistance_scl(z_i=z_i_se, surface=surface_se, diff=diffusivity_se, delta_x_array=dx_array,
                                           c_array=c_se_array)

        print("No.{}".format(count).center(30,"-"))
        print("step:{}".format(step))
        #print("turning point:{}".format(tp))
        print("phi_0:{}".format(phi_0))
        #print("delta_x_list: {}".format(delta_x_list[-1]))
        print("resistance:{}".format(resistance_se))
        print("scl length:{}".format(len_scl))
        #print("concentration:{}".format(c_list))
        output_dic.setdefault('index', []).append(count)
        output_dic.setdefault('phi_0', []).append(phi_0)
        output_dic.setdefault('scl_resistance', []).append(resistance_se)
        con_dic.setdefault('concentration',[]).append(c_list)
        con_dic.setdefault('length',[]).append(len_scl_list)

    filename = 'scl_resistance_{0}_{1}.csv'.format(cathode,se)
    df = pd.DataFrame(output_dic)
    df.to_csv(filename, index=False)

    print("Calculation Completion".center(30, '*'))
    if ip.Condition().picture == 1:
        filename2 = 'concentration_{0}_{1}.csv'.format(cathode, se)
        df2 = pd.DataFrame(con_dic)
        df2.to_csv(filename2, index=False)
    else:
        pass
    print("Start Drawing".center(30, '*'))
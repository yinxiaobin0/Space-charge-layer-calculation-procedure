import numpy as np
import pandas as pd
import input_parameters as ip
import matplotlib.pyplot as plt
import sys
import math
from scipy.constants import Boltzmann, R, elementary_charge, epsilon_0, N_A


def get_concentration_li(concentration_bulk, z_i, potential, potential_bulk, beta=10, alpha=0.1):
    '''
    基于玻尔兹曼分布求解不同位置的载流子浓度
    :param concentration_bulk: 载流子体相的浓度 （nm-3）
    :param z_i: 载流子的电荷数
    :param potential: 实时位点处的电势，初值由|MCFSS|提供
    :param potential_bulk: 体相处的电势，此处取0。后续可拓展怎么求固态电解质体相的电势
    :param alpha: 0.1 实际上并没有什么影响
    :param beta: 10 由于我们主要对比的是相对值，因此用该系数同比地缩放指数项地变化，默认为10。为了针对不同扩散系数物质，对于扩散系数快的物质，beta下降，指数项的变化会更明显；反之，则不明显
    :return:
    '''
    fermi_distribution = 1/(math.exp(z_i*elementary_charge/(beta*Boltzmann*temperature)*(potential - potential_bulk))+alpha)+0.001
    concentration_x = concentration_bulk*fermi_distribution
    return concentration_x
    #0.001：加入这个尾项是为了防止因电势过大而导致急剧变化，一般超过三个数量级后边可忽略


def get_electrical_field_strength(permittivity, concentration_x_last, concentration_bulk, z_i, c_im, potential_bulk, potential):
    '''

    :param permittivity:
    :param concentration_x_last:
    :param concentration_bulk:
    :param z_i: 同上
    :param c_im:
    :param potential_bulk:
    :param potential:
    :return:
    '''
    electrical_field_strength = sgn*math.sqrt((2*Boltzmann*temperature/permittivity)*abs((concentration_x_last-concentration_bulk)+z_i*elementary_charge*c_im/(Boltzmann*temperature) * (potential_bulk- potential)))
    return electrical_field_strength


def get_factor_b(step_x):
    if step_x < 0.5*total_step:
        b = (phi_end/phi_0)**(1/(0.5*total_step))
    else:
        b = (phi_end/phi_x)**(1/(total_step-step_x))
    return b


def get_delta_x(charge_density_scl, charge_density):
    '''

    :param charge_density_scl:
    :param charge_density:
    :return:
    '''
    delta_charge_density_scl = charge_density_scl/total_step
    delta_x = abs(delta_charge_density_scl/charge_density)
    return delta_x


def get_potential(potential_last,delta_x,electrical_field_strength,permittivity,charge_density,delta_charge_density):
    '''

    :param potential_last:
    :param delta_x:
    :param electrical_field_strength:
    :param permittivity:
    :param charge_density:
    :param delta_charge_density:
    :return:
    '''
    potential_next = potential_last - delta_x*electrical_field_strength-delta_x**2/(2*permittivity)*\
                     (charge_density+delta_charge_density/3)
    return potential_next


def get_resistance_scl(z_i, surface, diff, delta_x_array, delta_c_array):
    '''

    :param z_i:
    :param surface:
    :param diff:
    :param delta_x_array:
    :param delta_c_array:
    :return:
    '''
    r_scl = (Boltzmann*temperature)/(surface*z_i**2*elementary_charge**2*diff)*sum(delta_x_array/delta_c_array)
    return r_scl


if __name__ == '__main__':
    cathode = "linbo3" #oxide or sulfide or coating
    #cathode_vol = abs(-0.45) #此处修改正极的MCFSSS
    se = "li3ycl6"
    se_type = "sulfide"
    print("Input Completion".center(30,'*'))
    print("Cathode：{}".format(cathode))
    print("SEs：{}".format(se))
    print("Start Calculating".center(30, '*'))

    mod = sys.modules["input_parameters"]
    array1 = mod.__dict__[cathode]()  # 集合了正极材料性质的数组
    array2 = mod.__dict__[se]()  # 集合了固态电解质材料的数组

    total_step = ip.Condition().total_step
    temperature = ip.Condition().temperature
    alpha = ip.Condition().alpha
    beta = ip.Condition().beta
    #delta_x = ip.Condition().delta_x
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

    '''硫化物正极'''
    #potential_list = [-0.241,-0.225,-0.205,-0.211,-0.299,-0.266,-0.219,-0.123,-0.26,-0.239,-0.193,-0.23,-0.242,-0.153,-0.321]
    '''氧化物正极'''
    #potential_list = [-0.551,-0.57,-0.487,-0.498,-0.678,-0.718,-0.775,-0.814,-0.4,-0.403,-0.438,-0.691]
    '''氧化物涂层'''
    #potential_list = [2.669, 2.375, 2.523, 3.273]
    '''能带-氧化物正极'''
    #potential_list = [-0.7639, -1.3979, -1.7768, -1.2131, -0.0747, -0.3073, -0.3853, -0.384, -1.7929, -1.8312, -1.1405, -2.0677]
    '''能带-硫化物正极'''
    # potential_list = [-]
    '''能带-氧化物涂层'''
    potential_list = [0.75141, 1.62629, 1.0014, 1.79498]
    '''测试用'''
    #potential_list = [-0.4,-0.403,-0.438]
    del array1, array2

    if c_bulk_se <= 10:
        sp = 0.5
    else:
        sp =1
    '''调控影响因子'''
    if cathode == "licoo2":
        if se_type == "sulfide":
            if 1E-8 < diffusivity_se:
                beta = 2
            elif 1E-9 < diffusivity_se < 1E-8:
                beta = 7
            else:
                beta = 10
        elif se_type == "halide":
            if 1E-8 < diffusivity_se:
                beta = 10
            elif 1E-9 < diffusivity_se < 1E-8:
                beta = 12
            else:
                beta = 15
            sp = 0.1
        else:
            if 1E-8 < diffusivity_se:
                beta = 10
            elif 1E-9 < diffusivity_se < 1E-8:
                beta = 12
            else:
                beta = 15
    elif cathode == "litis2":
        if 1E-8 < diffusivity_se:
            beta = 10
        elif 1E-9 < diffusivity_se < 1E-8:
            beta = 12
        else:
            beta = 15
        if se == "li3ps4":
            sp = 0.5
        if se_type == "halide":
            sp = 0.1

    elif cathode == "linbo3":
        beta = 10
        if se_type  == "halide":
            sp =0.1
        #sp = 0.5
    else:
        print("Input error:cathode!")

    # if "o" in cathode and "s" in se:
    #     if se == "li6ps5cl":
    #         alpha = 10
    #     elif se == "li10gep2s12":
    #         alpha = 100
    #     else:
    #         alpha = 1
    # elif "o" in cathode and "o" in se:
    #     alpha = 1
    # elif "s" in cathode and "s" in se:
    #     alpha = 1
    #     if se == "li3ps4":
    #         alpha = 0.1
    # elif "s" in cathode and "o" in se:
    #     alpha = 1
    # else:
    #     alpha = 1
    rou_scl = 0
    resistance_list = []
    count = 0
    phi_end = Boltzmann*temperature/elementary_charge*math.log(1+1E-3)
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
        c_ref = c_bulk_se / (1 + alpha)

        while True:
            step += 1
            phi_x = phi_list[-1]
            if cathode == "linbo3":
                phi_x = phi_x*0.1

            c_x_se = get_concentration_li(concentration_bulk=c_bulk_se, z_i = z_i_se, potential=phi_x, potential_bulk= se_vol, alpha= alpha, beta=beta)
            c_list.append(c_x_se)
            if step > 1:
                delta_c = abs(c_list[-1] - c_list[-2])
            else:
                delta_c = 0
            delta_c_list.append(delta_c)

            rou_x = (z_i_se * c_x_se + z_i_im * c_e_se) * elementary_charge #固态电解质中的电荷密度(<0)
            rou_list.append(rou_x)
            if step != 1:
                delta_rou = rou_list[-1] - rou_list[-2] #(>0)
            else:
                delta_rou = 0

            efs = get_electrical_field_strength(permittivity = permittivity_se, concentration_x_last = c_x_se, concentration_bulk = c_bulk_se, z_i = -1, c_im = c_e_se, potential_bulk= se_vol, potential=phi_list[-1])
            efs_list.append(efs)
            if step == 1:
                rou_scl = -efs_list[0]*permittivity_se
            else:
                rou_scl = rou_scl

            dx = get_delta_x(charge_density_scl=rou_scl, charge_density=rou_list[-1])
            #dx =0.01
            phi_x = phi_list[-1] - dx * efs - dx **2 / (2*permittivity_se) * (rou_list[-1] + delta_rou/3)
            #factor_b = get_factor_b(step_x = step)

            # if step == 1:
            #     delta_x = abs((charge_density/total_step)/rou_x)
            #     phi_x = phi_list[-1] - delta_x * efs - delta_x **2 / (2*permittivity_se) * (rou_x+delta_rou/3)
            # elif step < 0.5*total_step:
            #     delta_x = abs((charge_density /total_step) / rou_x)
            #     phi_x = phi_list[-1] - delta_x * efs - delta_x ** 2 / (2 * permittivity_se) * (
            #                 rou_x + delta_rou / 3)
            # else:
            #     phi_x = factor_b*phi_list[-1] #单位：V
            #     delta_x = permittivity_se/(rou_x+delta_rou/3)*(-efs_list[-1]+math.sqrt(abs(efs_list[-1]**2-2/permittivity_se*(phi_x-phi_list[-1])*(rou_x+delta_rou/3))))
            #     '''delta_x 单位：nm'''
            delta_x_list.append(dx)
            phi_list.append(phi_x)
            len_scl += abs(dx)

            '''测试用'''
            #print(step)
            # #print(phi_end)
            #print(c_list[-1])
            # #print(delta_x_list)
            # #print(rou_x)
            # #print(electrical_field)
            #print(phi_x)
            # # print()
            '''测试用完毕'''

            # if 1-c_list[-1]/c_ref > 0.995: #or step > 10000:
            #     break
            if phi_list[-1] < 0:
                break
        del delta_c_list[0]

        #del delta_c_list[-1]
        c_se_array = np.array(delta_c_list)
        dx_array = np.array(delta_x_list)
        #print(c_se_array)
        #print(phi_list)
        resistance_se = get_resistance_scl(z_i=z_i_se, surface=surface_se, diff=diffusivity_se, delta_x_array=0.01, delta_c_array=c_se_array)*sp

        print("No.{}".format(count).center(30,"-"))
        print("step:{}".format(step))
        print("phi_0:{}".format(phi_0))
        print("delta_x_list: {}".format(delta_x_list[-1]))
        print("resistance:{}".format(resistance_se))
        print("scl length:{}".format(len_scl))
        print("concentration:{}".format(c_list[-1]))
        #resistance_list.append(resistance_se)

    # test = pd.DataFrame({'concentration':c_list,'rou':rou_list,'phi':phi_list})
    # test.to_csv('results',encoding='utf8')


    #resistance_cathode = get_resistance_scl(array1[0], array1[2], array1[4], delta_c_list_cathode)
    #resistance_total = resistance_cathode+resistance_se


    #print("concentration:{}".format(c_list))
    #print("phi:{}".format(phi_list))
    #print("delta_c:{}".format(delta_c_list))
    #print("efs:{}".format(efs_list))

    print("Calculation Completion".center(30,'*'))
    #print("Start Drawing".center(30, '*'))
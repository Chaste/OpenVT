from multiprocessing.pool import Pool
import os

end_time = 10000
# base_runner = f"/home/chaste/build/projects/OpenVT/test/Test02MonolayerGrowth_VT -end_time {end_time}"
base_runner = f"/home/chaste/build/projects/OpenVT/test/Test02MonolayerGrowth -end_time {end_time}"
sample_rate = 10000
cut_off_length = 2.0
# beta_vect =  [0.0, 0.7, 0.80, 0.90, 0.9200]
# gamma_vect = [0.0, 0.9, 0.95, 0.99, 0.9975]


# beta_vect =  [0.0, 0.968654, 0.984427, 0.999778, 1.000000]
# gamma_vect = [0.0, 0.677758, 0.812121, 0.967786, 1.000000]


beta_vect =  [0.0, 0.971458, 0.984134, 0.9917,   0.993896]
gamma_vect = [0.0, 0.677817, 0.800954, 0.890508, 0.902354]

# ------------------------------------------------------------------------------
# This code block is to run simulations sequentially
# ------------------------------------------------------------------------------
# for beta_ind, beta_parameter in enumerate(beta_vect):
#     for gamma_ind, gamma_parameter in enumerate(gamma_vect):
#         if beta_ind + gamma_ind > 0:
# 
#             output_name = f"LongTests_7/MonolayerGrowth_Beta{beta_ind}_Gamma{gamma_ind}_EndTime{end_time}"
# 
#             model_runner = base_runner + f" -cut_off_length {cut_off_length} -p_beta {beta_parameter} -p_gamma {gamma_parameter} -output_name {output_name} -sample_rate {sample_rate}"
# 
#             print(model_runner)
#             os.system(model_runner)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# This code block is to run simulations in parallel
# ------------------------------------------------------------------------------
# def run_simulation(params):
#     beta_ind, beta_parameter, gamma_ind, gamma_parameter = params
    
#     output_name = f"LongTests_Calibrted_V5/MonolayerGrowth_Beta{beta_ind}_Gamma{gamma_ind}_EndTime{end_time}"
#     force_law = 'quadratic'

#     model_runner = base_runner + f" -p_force_law {force_law} -cut_off_length {cut_off_length} -p_beta {beta_parameter} -p_gamma {gamma_parameter} -output_name {output_name} -sample_rate {sample_rate} -random_seed {0}"

#     print(model_runner)
#     os.system(model_runner)

# tasks = []
# for beta_ind, beta_parameter in enumerate(beta_vect):
#     for gamma_ind, gamma_parameter in enumerate(gamma_vect):
#         if beta_ind>0 + gamma_ind > 0:
#             tasks.append((beta_ind, beta_parameter, gamma_ind, gamma_parameter))

# if __name__ == '__main__':
#     with Pool(processes=2) as pool:
#         pool.map(run_simulation, tasks)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# This code block is to run a single simulation
# ------------------------------------------------------------------------------
# beta_ind = 0
# beta_parameter = -1
# gamma_ind = 0
# gamma_parameter = -1
# sample_rate = 500
# index = 0
# # output_name = f"NoInhibition_weak/MonolayerGrowth_Beta{beta_ind}_Gamma{gamma_ind}_EndTime{end_time}"
# output_name = f"NoInhibition_weak/MonolayerGrowth_Beta_Index{index}"


# model_runner = base_runner + f" -cut_off_length {cut_off_length} -p_beta {beta_parameter} -p_gamma {gamma_parameter} -output_name {output_name} -sample_rate {sample_rate}"

# print(model_runner)
# os.system(model_runner)
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# This code block is to run simulations in parallel to produce the histograms with no inhibition
# ------------------------------------------------------------------------------
def run_simulation(params):
    index, beta_parameter, gamma_parameter = params
    
    force_law = "quadratic" # "log"; # "quadratic"; # "linear";

    output_name = f"New/NoInhibition_strong_100/{force_law}/MonolayerGrowth_Beta_Index{index}"

    model_runner = base_runner + f" -p_force_law {force_law} -cut_off_length {cut_off_length} -p_beta {beta_parameter} -p_gamma {gamma_parameter} -output_name {output_name} -sample_rate {sample_rate} -random_seed {index}"

    print(model_runner)
    os.system(model_runner)

tasks = []
beta_parameter = 0.0
gamma_parameter = 0.0
sample_rate = 0.25/0.002

for index in range(0, 100, 1):
    tasks.append((index, beta_parameter, gamma_parameter))

if __name__ == '__main__':
    with Pool(processes=4) as pool:
        pool.map(run_simulation, tasks)

# ------------------------------------------------------------------------------









# beta_ind = 1
# beta_parameter = 0.0 #0.971458
# gamma_ind = 1
# gamma_parameter = 0.0 #0.677817
# sample_rate = 0.5/0.002
# index = 0

# force_law = "quadratic" # "log"; # "quadratic"; # "linear";

# output_name = f"TestNewForceLaws/{force_law}/"

# model_runner = base_runner + f" -p_force_law {force_law} -cut_off_length {cut_off_length} -p_beta {beta_parameter} -p_gamma {gamma_parameter} -output_name {output_name} -sample_rate {sample_rate} -random_seed {0}"

# print(model_runner)
# os.system(model_runner)
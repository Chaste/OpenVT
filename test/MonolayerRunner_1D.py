from multiprocessing.pool import Pool
import os

end_time = 10000

p_quadratic_spring_stiffness = 155.0

base_runner = f"/home/chaste/build/projects/OpenVtMonolayer/test/TestMonolayerGrowth1dExamples -p_quadratic_spring_stiffness {p_quadratic_spring_stiffness}"

model_runner = base_runner

print(model_runner)
os.system(model_runner)
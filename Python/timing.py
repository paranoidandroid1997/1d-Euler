import subprocess
import os
import time

# List of commands to execute

command = "./slugEuler1d"
inputs = [f"./inputs/input{num}.txt" for num in range(26, 30)]
inputs += [f"./inputs/input{num}.txt" for num in range(38, 42)]
inputs += [f"./inputs/input{num}.txt" for num in range(42, 46)]

os.chdir("../Fortran") 
for i in inputs:
    os.system(f"cp {i} ./slug.init")


    print(f"Timing input: {i}")
    # Execute the command and capture its output

    num_meas = 5
    total_elapsed_time = 0.0
    for _ in range(num_meas):
        start_time = time.time()
        result = subprocess.run([command], stdout=subprocess.DEVNULL)
        end_time = time.time()

        elapsed_time = end_time - start_time
        total_elapsed_time += elapsed_time
    
    
    # Print the timing information
    print(f"Elapsed time: {total_elapsed_time / num_meas}")


os.chdir("../Python") 
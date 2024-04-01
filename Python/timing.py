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

    order = 0
    ic = ""

    # Open the text file
    with open('./slug.init', 'r') as file:
        # Read the lines
        lines = file.readlines()

        # Iterate through the lines
        for line in lines:
            # Check if the line contains 'sim_order'
            if 'sim_order' in line:
                # Split the line by whitespace
                parts = line.split()
                # Get the second part which should be the number
                order = int(parts[1])
            elif "IC_type" in line:
                # Split the line by whitespace
                parts = line.split()
                # Get the second part which should be the ic
                ic = parts[1]



    print(f"Timing input: {i}")
    print(f"Order: {order}")
    print(f"IC: {ic}")
    # Execute the command and capture its output

    num_meas = 10
    total_elapsed_time = 0.0
    for _ in range(num_meas):
        start_time = time.time()
        result = subprocess.run([command], stdout=subprocess.DEVNULL)
        end_time = time.time()

        elapsed_time = end_time - start_time
        total_elapsed_time += elapsed_time
    
    
    # Print the timing information
    print(f"Elapsed time: {total_elapsed_time / num_meas}")
    print()


os.chdir("../Python") 
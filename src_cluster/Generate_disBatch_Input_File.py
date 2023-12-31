##################################################################
## Generate input file for disBatch
## Run holoQUADS circuits in parallel
##################################################################

import numpy as np
import os

def generate_input_file(input_index, task_file):
    '''Generate corresponding folders and input files based on chemical potential'''
    folder_name = "Sample" + "{}".format(input_index) + "/"
    task_file.write("cd " + folder_name \
        + " &&  julia --threads=1 Disentangler.jl" + " &> disentangler_" \
        + "{}".format(input_index) + ".log" + "\n")
    
def main():
    sample_seed=1
    lower_bound=(sample_seed-1)*500+1
    upper_bound=sample_seed*500+1
    sample_list = np.arange(lower_bound, upper_bound, 1)
    # location = os.path.dirname(os.path.realpath(__file__))

    submit_file = open("N100", "a")
    for tmp in sample_list:
        generate_input_file(tmp, submit_file)
    submit_file.close()    

main()
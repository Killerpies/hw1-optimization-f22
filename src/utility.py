import csv
import math
import numpy as np
import os
import random


def get_start(file_name):
    """
    Get starting position
    @return np.array of 2 values
    """
    if os.path.exists('data'):
        full_file_name = os.path.join('data', file_name)
        if os.path.exists(full_file_name):
            # Get existing data and return
            # Just delete file if you want to create new random data
            with open(full_file_name, 'rt') as fin:
                reader = csv.reader(fin, delimiter=',')
                terms = []
                for line in reader:
                    assert len(line) == 2, "Require two values for starting point"
                    return np.array([float(val) for val in line])

        # Create random data and write to file for next time
        xs = 10.0*(random.random()-0.5) # Random start in -5, 5 range
        ys = 10.0*(random.random()-0.5) # Random start in -5, 5 range
        start = np.array([xs, ys])

        print(f"Save starting position to {file_name}!")
        with open(full_file_name, 'wt') as fout:
            writer = csv.writer(fout, delimiter=',')
            writer.writerow(start)
        return start

    else:
        # No data folder, so just generate random terms
        print("Generate random starting point ...")
        xs = 10.0*(random.random()-0.5) # Random start in -5, 5 range
        ys = 10.0*(random.random()-0.5) # Random start in -5, 5 range
        start = np.array([xs, ys])
        return start

def get_terms(file_name, n):

    if os.path.exists('data'):
        full_file_name = os.path.join('data', file_name)
        if os.path.exists(full_file_name):
            # Get existing data and return
            # Just delete file if you want to create new random data
            with open(full_file_name, 'rt') as fin:
                reader = csv.reader(fin, delimiter=',')
                terms = []
                for line in reader:
                    terms.append([float(val) for val in line])

                if len(terms) >= n:
                    print(f"Return partial of {n} terms from existing file!")
                    return terms[:n]
                else:
                    print(f"Return {len(terms)} terms as defined in existing file!")
                    return terms

        # Create random data and write to file for next time
        print(f"Generate {n} terms with mostly random data and save to {file_name}!")
        terms = get_random_terms(n)
        with open(full_file_name, 'wt') as fout:
            writer = csv.writer(fout, delimiter=',')
            for line in terms:
                writer.writerow(line)
            return terms

    else:
        # No data folder, so just generate random terms
        print(f"Generate {n} random terms ")
        return get_random_terms(n)


def get_random_terms(n=4):
    """
    Generate the term values, some logic to keep somewhat reasonable for problem
    But the values really don't impact the math or code you need to write

    """
    if n is None:
        n = random.randint(2,7)

    terms = []
    for i in range(n):

        if i == 0:
            amp = random.choice([-1, 1])*10.0*(random.random() - 0.5)  # constant offset in -5., 5 range
            wx = 0.  # First is constant offset (w=0)
            wy = 0.  # First is constant offset (w=0)
        elif i == 1:
            # Keep main basin in normal bounds
            amp = random.choice([-1, 1])*10.0 + 3*(random.random() - 0.5)  # constant offset in -5., 5 range
            wx = math.pi/15 # First is single cup across x range
            wy = math.pi/15

        else:
            amp = random.choice([-1, 1])*8.0*(1. + 0.25*(random.random() - 0.5))/i
            wx =  i*0.05 + i*1.0*random.random()
            wy =  i*0.05 + i*1.0*random.random()

        xc = 4.0*(random.random() - 0.5)  # x centroid in -2, 2 range
        yc = 4.0*(random.random() - 0.5)  # y centroid in -2, 2 range

        terms.append([amp, wx, wy, xc, yc])
    return terms

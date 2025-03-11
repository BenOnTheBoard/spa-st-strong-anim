import os
import random
import shutil
from itertools import product

import numpy as np
from tqdm import tqdm

from instanceGenerator import SPASTIG

######
# Panel:
DENSITY_STEPS = 8
FILES_PER_DENSITY_PAIR = 100
DIRECTORY = "examples/gen/"
######

print("Directory cleanup...")
if os.path.isdir(DIRECTORY):
    shutil.rmtree(DIRECTORY)
os.mkdir(DIRECTORY)
print("Cleanup done.")

density_levels = np.arange(0, 1, 1 / DENSITY_STEPS)
density_pairs = list(product(density_levels, density_levels))

for densities in tqdm(density_pairs):
    students = random.randint(1, 4)
    projects = random.randint(1, 16)
    lecturers = random.randint(1, projects)

    for k in range(FILES_PER_DENSITY_PAIR):
        S = SPASTIG(
            students=students,
            projects=projects,
            lecturers=lecturers,
            pref_list_length_lb=0,
            pref_list_length_ub=projects,
            student_tie_density=densities[0],
            lecturer_tie_density=densities[1],
        )
        filename = f"{DIRECTORY}{int(densities[0] * DENSITY_STEPS)}_{int(densities[1] * DENSITY_STEPS)}_{k}.txt"
        S.write_instance_with_ties(filename)

from itertools import product

import numpy as np
import random
import os
import shutil

from instanceGenerator import SPASTIG

######
# Panel:
DENSITY_STEPS = 5
FILES_PER_DENSITY_PAIR = 5
DIRECTORY = "examples/gen/"
######

if os.path.isdir(DIRECTORY):
    shutil.rmtree(DIRECTORY)
os.mkdir(DIRECTORY)

density_levels = np.arange(0, 1, 1 / DENSITY_STEPS)
density_pairs = product(density_levels, density_levels)

for densities in density_pairs:
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
        print(filename)
        S.write_instance_with_ties(filename)

import os
import random
import shutil
from itertools import product

import numpy as np
from tqdm import tqdm

from instanceGenerator import SPASTIG

######
# Panel:
DENSITY_STEPS = 6
FILES_PER_DENSITY_PAIR = 5000
DIRECTORY = "examples/gen/"

MAX_STUDENTS = 5
MAX_PROJECTS = 6
######

print("Directory cleanup...")
if os.path.isdir(DIRECTORY):
    shutil.rmtree(DIRECTORY)
os.mkdir(DIRECTORY)
print("Cleanup done.")

step_size = 1 / DENSITY_STEPS
density_levels = np.arange(0, 1 + step_size / 2, step_size)
density_pairs = list(product(density_levels, density_levels))

for densities in tqdm(density_pairs):
    for k in range(FILES_PER_DENSITY_PAIR):
        students = random.randint(2, MAX_STUDENTS)
        projects = random.randint(2, MAX_PROJECTS)
        lecturers = random.randint(2, projects)

        S = SPASTIG(
            students=students,
            projects=projects,
            lecturers=lecturers,
            pref_list_length_lb=projects,
            pref_list_length_ub=projects,
            student_tie_density=densities[0],
            lecturer_tie_density=densities[1],
        )

        filename = f"{DIRECTORY}{int(densities[0] * DENSITY_STEPS)}_{int(densities[1] * DENSITY_STEPS)}_{k}.txt"
        S.write_instance_with_ties(filename)

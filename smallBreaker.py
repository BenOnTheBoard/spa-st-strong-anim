import os
import random
import shutil

from tqdm import tqdm

from bruteforce import STSMBruteForce
from instanceGenerator import SPASTIG
from spaststrong import SPAST_STRONG

######
# Panel:
DIRECTORY = "examples/small_breakers/"
MAX_TRIALS = 1000
densities = (1, 1)
######

print("Directory cleanup...")
if os.path.isdir(DIRECTORY):
    shutil.rmtree(DIRECTORY)
os.mkdir(DIRECTORY)
print("Cleanup done.")


students = 4
projects = random.randint(3, 4)
lecturers = random.randint(1, projects)

for k in tqdm(range(MAX_TRIALS)):
    filename = f"{DIRECTORY}{k}_breaker.txt"

    instance_ssm_list = []
    while not instance_ssm_list:
        S = SPASTIG(
            students=students,
            projects=projects,
            lecturers=lecturers,
            pref_list_length_lb=0,
            pref_list_length_ub=projects,
            student_tie_density=densities[0],
            lecturer_tie_density=densities[1],
        )
        S.write_instance_with_ties(filename)

        bruteforcer = STSMBruteForce(filename)
        bruteforcer.choose()
        instance_ssm_list = bruteforcer.get_ssm_list()

    strong_solver = SPAST_STRONG(filename)
    matching = strong_solver.run()

    if matching not in instance_ssm_list:
        break
    else:
        os.remove(filename)

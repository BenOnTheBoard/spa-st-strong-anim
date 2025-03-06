import os

from bruteforce import STSMBruteForce
from spaststrong import SPAST_STRONG


def pprint_G(G):
    for k, v in G.items():
        print(f"{k} :\t {v}")


def break_line():
    print(f"\n{'---' * 10}")


def test(filename, verbose=False):
    bruteforcer = STSMBruteForce(filename)
    bruteforcer.choose(1)
    bf_ssm_list = bruteforcer.get_ssm_list()

    strong_solver = SPAST_STRONG(filename)
    solver_matching = strong_solver.run()

    if verbose:
        break_line()
        print("Bruteforced results:")
        bf_ssm_list = bruteforcer.get_ssm_list()
        if bf_ssm_list:
            for matching in bf_ssm_list:
                print(matching)
        else:
            print("There are no stable matchings for this instance.")
        break_line()
        print("Solver result:")
        print(solver_matching)
        break_line()

    if not bf_ssm_list:
        return True
    elif solver_matching in bf_ssm_list:
        return True
    else:
        return False


def all_tests():
    for subdir, _, files in os.walk("examples"):
        for filename in files:
            if filename.endswith(".txt"):
                filepath = subdir + os.sep + filename
                if test(filepath):
                    print(f"{filename}:\tpass\t")
                else:
                    print(f"{filename}:\t\tfail")


def single_test(filepath):
    if test(filepath):
        print(f"{filepath}:\tpass\t")
    else:
        print(f"{filepath}:\t\tfail")


single_test("examples/K instances/K33.txt")

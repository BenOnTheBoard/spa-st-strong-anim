import os

from bruteforce import STSMBruteForce
from spaststrong import SPAST_STRONG
from SPAST_Strong_Anim import SPAST_STRONG_ANIM


def pprint_G(G):
    for k, v in G.items():
        print(f"{k} :\t {v}")


def break_line():
    print(f"\n{'---' * 10}")


def test(filename, verbose=False):
    bruteforcer = STSMBruteForce(filename)
    bruteforcer.choose()
    bf_ssm_list = bruteforcer.get_ssm_list()

    strong_solver = SPAST_STRONG(filename)
    solver_matching = strong_solver.run()

    if verbose:
        break_line()
        print("Bruteforced results:")
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


def all_tests(verbose=False):
    wins = 0
    seen = 0
    for subdir, _, files in os.walk("examples"):
        print(f"{subdir}:")
        for filename in files:
            if filename.endswith(".txt"):
                filepath = subdir + os.sep + filename
                seen += 1
                if test(filepath, verbose=verbose):
                    wins += 1
                    # print(f"\t{filename}:\tpass\t\t{wins}/{seen}")
                else:
                    print(f"\t{filename}:\t\tfail\t{wins}/{seen}")
                if verbose:
                    print("\n\n")


def single_test(filepath, verbose=True):
    if test(filepath, verbose=verbose):
        print(f"{filepath}:\tpass\t")
        return True
    else:
        print(f"{filepath}:\t\tfail")
        return False


# all_tests()
# single_test("examples/misc/5_5_30.txt")

# cache abuse trick
f = "examples/misc/5530simple.txt"
if not single_test(f, verbose=False):
    SPAST_STRONG_ANIM(f)

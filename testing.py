import os

from bruteforce import STSMBruteForce
from spaststrong import SPAST_STRONG


def pprint_G(G):
    for k, v in G.items():
        print(f"{k} :\t {v}")


def break_line():
    print(f"\n{'---' * 10}")


def test(filename, verbose=False):
    bruteforcer = STSMBruteForce(filename=filename)
    bruteforcer.choose()
    bf_ssm_list = bruteforcer.get_ssm_list()

    strong_solver = SPAST_STRONG(filename=filename)
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
                if test(filepath):
                    wins += 1
                    if verbose:
                        print(f"\t{filename}:\tpass\t\t{wins}/{seen}")
                else:
                    if verbose:
                        print(f"\t{filename}:\t\tfail\t{wins}/{seen}")
        if verbose:
            print("\n\n")

    print(f"Final score:\t{wins}/{seen},\t{100 * (wins / seen):.2f}")


def single_test(filepath, verbose=True):
    if test(filepath, verbose=verbose):
        print(f"{filepath}:\tpass\t")
        return True
    else:
        print(f"{filepath}:\t\tfail")
        return False


# all_tests(verbose=True)
single_test("examples/small_breakers/no_del_2bii.txt")
# single_test("examples/misc/5530simple.txt", verbose=False)

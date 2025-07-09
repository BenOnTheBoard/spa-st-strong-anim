import random
import time
import argparse
from multiprocessing import Process, Value, Event, Lock

from bruteforce import STSMBruteForce
from instanceGenerator import SPASTIG
from spaststrong import SPAST_STRONG


def check_condition_student(solver):
    for si in solver.sp:
        if len(solver.G[si]["bound"]) >= 2:
            return False
    return True


def check_condition_lecturer(solver):
    for lk in solver.lp:
        lk_instance_info = solver.lp[lk]

        M_lk = set()
        for s, p in solver.M.items():
            if p in lk_instance_info["projects"]:
                M_lk.add(s)

        if solver.G[lk]["replete"]:
            if lk_instance_info["cap"] > len(M_lk):
                return False
        else:
            if solver.lquota(lk) > len(M_lk):
                return False
    return True


def check_condition_project(solver):
    for pj in solver.plc:
        lk = solver.plc[pj]["lec"]
        lk_instance_info = solver.lp[lk]
        # lk_alg_info = solver.G[lk]
        pj_instance_info = solver.plc[pj]
        pj_alg_info = solver.G[pj]

        M_pj = set()
        M_lk = set()
        for s, p in solver.M.items():
            if p == pj:
                M_pj.add(s)
            if p in lk_instance_info["projects"]:
                M_lk.add(s)

        if pj_alg_info["replete"]:
            if lk_instance_info["cap"] > len(M_lk):
                if pj_instance_info["cap"] > len(M_pj):
                    return False

    return True


def list_pref(lec_list, s_prime, s):
    """Return true when s_prime > s on Lk"""
    for tie in lec_list:
        if s in tie:
            return False
        if s_prime in tie:
            return True
    raise ValueError("Neither student is on the list.")


def list_indif(lec_list, s_prime, s):
    """Return true when s_prime = s on Lk"""
    for tie in lec_list:
        if s in tie and s_prime in tie:
            return True
    return False


def capacity_validity(solver):
    project_apps = {pj: 0 for pj in solver.plc.keys()}
    lecturer_apps = {lk: 0 for lk in solver.lp.keys()}
    for si_info in solver.sp.values():
        for pj in si_info["list_rank"].keys():
            project_apps[pj] += 1
            lk = solver.plc[pj]["lec"]
            lecturer_apps[lk] += 1

    for pj, apps in project_apps.items():
        if solver.plc[pj]["cap"] > apps:
            return False

    for pj, apps in project_apps.items():
        if solver.plc[pj]["cap"] > apps:
            return False

    return True


def worker(shared_counter, max_trials, found_event, lock, process_id):
    print(f"[Process {process_id}] Started.")
    while True:
        if found_event.is_set():
            print(f"[Process {process_id}] Exiting due to: found_event set.")
            return

        with lock:
            if shared_counter.value >= max_trials:
                print(f"[Process {process_id}] Exiting due to: max trials reached")
                return
            shared_counter.value += 1
            trial_num = shared_counter.value

        filename = f"test_{process_id}.txt"

        valid = False
        while not valid:
            densities = (random.uniform(0, 1), random.uniform(0, 1))
            students = random.randint(1, 5)
            projects = random.randint(1, 4)
            lecturers = random.randint(1, projects)

            S = SPASTIG(
                students=students,
                projects=projects,
                lecturers=lecturers,
                pref_list_length_lb=1,
                pref_list_length_ub=projects,
                student_tie_density=densities[0],
                lecturer_tie_density=densities[1],
            )
            S.instance_generator_with_ties()

            solver = SPAST_STRONG(generator=S)
            valid = capacity_validity(solver)

        bruteforcer = STSMBruteForce(generator=S)
        bruteforcer.choose()
        instance_ssm_list = bruteforcer.get_ssm_list()
        exists_ssm = bool(instance_ssm_list)

        solver.run()

        if exists_ssm and solver.M not in instance_ssm_list:
            print(f"[Process {process_id}] reports an algorithm error at {trial_num}")
            found_event.set()
            return

        lemma_s = check_condition_student(solver)
        lemma_l = check_condition_lecturer(solver)
        lemma_p = check_condition_project(solver)
        conditions = (lemma_s, lemma_l, lemma_p)

        is_ssm_by_lemmas = all(conditions)

        if exists_ssm and not is_ssm_by_lemmas:
            print(f"[Process {process_id}] got one at trial {trial_num}")
            print(f"Exists SSM: {exists_ssm}")
            print(f"Condition tuple: {conditions}")
            found_event.set()
            return


def test_single_file(filename):
    bruteforcer = STSMBruteForce(filename=filename)
    bruteforcer.choose()
    instance_ssm_list = bruteforcer.get_ssm_list()
    exists_ssm = bool(instance_ssm_list)

    solver = SPAST_STRONG(filename=filename)
    solver.run()

    lemma_s = check_condition_student(solver)
    lemma_l = check_condition_lecturer(solver)
    lemma_p = check_condition_project(solver)
    conditions = (lemma_s, lemma_l, lemma_p)

    is_ssm_by_lemmas = all(conditions)

    if exists_ssm == is_ssm_by_lemmas:
        print("Pass")
    else:
        print("Fail with:")
        print(f"\tSSM:\t{exists_ssm}")
        print(f"\tLemmas:\t{conditions}")

    for s in solver.sp:
        print(f"{s}   {solver.G[s]}")
    for p in solver.plc:
        print(f"{p} {solver.plc[p]['cap']} {solver.G[p]}")
    for L in solver.lp:
        print(f"{L} {solver.lp[L]['cap']} {solver.G[L]}")


def main():
    parser = argparse.ArgumentParser(
        description="Parallel stable matching instance checker."
    )
    parser.add_argument("trials", type=int, help="Total number of trials to run.")
    parser.add_argument("num_processes", type=int, help="Number of parallel processes.")
    args = parser.parse_args()

    max_trials = args.trials
    num_processes = args.num_processes

    shared_counter = Value("i", 0)
    found_event = Event()
    lock = Lock()
    processes = []

    for i in range(num_processes):
        p = Process(
            target=worker, args=(shared_counter, max_trials, found_event, lock, i)
        )
        p.start()
        processes.append(p)

    while not found_event.is_set() and shared_counter.value < max_trials:
        time.sleep(1)
        print(f"[Main] Trials completed: {shared_counter.value}")

    for p in processes:
        p.join()

    if not found_event.is_set():
        print(f"[Main] No instance found in {max_trials}")


if __name__ == "__main__":
    # test_single_file("test_0.txt")
    main()

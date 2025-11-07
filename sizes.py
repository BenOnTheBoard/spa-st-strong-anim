import random
import time
import argparse
from multiprocessing import Process, Value, Event, Lock

from bruteforce import STSMBruteForce
from instanceGenerator import SPASTIG

def capacity_validity(generator):
    project_apps = {pj: 0 for pj in generator.plc.keys()}
    lecturer_apps = {lk: 0 for lk in generator.lp.keys()}
    for si_info in generator.sp.values():
        applied_lecturers = set()
        for pj in si_info[0]:
            project_apps[pj] += 1
            lk = generator.plc[pj][1]
            applied_lecturers.add(lk)
        
        for a_lk in applied_lecturers:
            lecturer_apps[a_lk] += 1

    for pj, apps in project_apps.items():
        if generator.plc[pj][0] > apps:
            return False

    for lk, apps in lecturer_apps.items():
        if generator.lp[lk][0] > apps:
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

        instance_ssm_list = []
        while len(instance_ssm_list) < 2:
            valid = False
            while not valid:
                densities = (random.uniform(0, 1), random.uniform(0, 1))
                students = random.randint(2, 4)
                projects = random.randint(2, 3)
                lecturers = 2

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
                valid = capacity_validity(S)

            bruteforcer = STSMBruteForce(generator=S)
            bruteforcer.choose()
            instance_ssm_list = bruteforcer.get_ssm_list()

        size_list = []
        for ssm in instance_ssm_list:
            size = sum(1 for pj in ssm.values() if pj != "")
            size_list.append(size)

        if not all(elt == size_list[0] for elt in size_list):
            print(f"[Process {process_id}] got one at trial {trial_num}")
            print(f"[Process {process_id}] claims {instance_ssm_list}.")
            S.write_instance_with_ties("sizes_example.txt")
            found_event.set()
            return


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
    main()

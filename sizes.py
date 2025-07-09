import random
import time
import argparse
from multiprocessing import Process, Value, Event, Lock

from bruteforce import STSMBruteForce
from instanceGenerator import SPASTIG


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
        instance_ssm_list = []
        while not instance_ssm_list or len(instance_ssm_list) < 2:
            densities = (random.uniform(0, 1), random.uniform(0, 1))
            students = random.randint(2, 5)
            projects = random.randint(2, 7)
            lecturers = random.randint(2, projects)

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
            S.write_instance_with_ties(filename)

            bruteforcer = STSMBruteForce(filename=filename)
            bruteforcer.choose()
            instance_ssm_list = bruteforcer.get_ssm_list()

        size_list = []
        for ssm in instance_ssm_list:
            size = sum(1 for pj in ssm.values() if pj != "")
            size_list.append(size)

        if not all(elt == size_list[0] for elt in size_list):
            print(f"[Process {process_id}] got one at trial {trial_num}")
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

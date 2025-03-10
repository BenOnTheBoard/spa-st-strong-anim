import math
import numpy as np
import random


class SPASTInstanceGenerator:
    def __init__(
        self,
        num_students,
        lower_bound,
        upper_bound,
        num_projects,
        num_lecturers,
        num_dimensions=5,
        force_project_capacity=0,
        force_lecturer_capacity=0,
    ) -> None:
        assert lower_bound <= upper_bound, (
            "Lower bound must be less than or equal to upper bound."
        )
        assert upper_bound <= num_projects, (
            "Upper bound must be less than or equal to the number of projects."
        )

        self._num_students = num_students
        self._num_projects = num_projects
        self._num_lecturers = num_lecturers

        self._force_project_capacity = force_project_capacity
        self._force_lecturer_capacity = force_lecturer_capacity
        self._total_project_capacity = int(math.ceil(1.1 * self._num_students))

        self._li = lower_bound  # lower bound of student preference list
        self._lj = upper_bound  # upper bound of student preference list

        self._reset_instance()

        self._num_dimesions = num_dimensions
        self.to_project_string = lambda x: f"p{x + 1}"

    def _reset_instance(self):
        self._sp = {
            f"s{i}": [] for i in range(1, self._num_students + 1)
        }  # student -> [project preferences]
        self._plc = {
            f"p{i}": [1, ""] for i in range(1, self._num_projects + 1)
        }  # project -> [capacity, lecturer]
        self._lp = {
            f"l{i}": [0, [], 0, 0] for i in range(1, self._num_lecturers + 1)
        }  # lecturer -> [capacity, project preferences, max of all c_j, sum of all c_j]

    def _assign_project_lecturer(self, project, lecturer):
        self._plc[project][1] = lecturer
        self._lp[lecturer][1].append(project)
        self._lp[lecturer][3] += self._plc[project][0]  # track sum of all c_j
        if self._plc[project][0] > self._lp[lecturer][2]:  # track max of all c_j
            self._lp[lecturer][2] = self._plc[project][0]

    def _sample_points(self, num_points: int):
        return np.random.uniform(0, 1, (num_points, self._num_dimesions))

    def _distance_function(self, points, point):
        return np.linalg.norm(points - point, axis=1)

    def _get_ordered_list(self, points_list, idx, length=None, reverse=False):
        return list(
            map(
                self.to_project_string,
                np.argsort(
                    self._distance_function(self._project_points, points_list[idx])
                )[:: -1 if reverse else 1][:length],
            )
        )

    def _generate_projects(self):
        """
        Generates projects for the SPA-P problem.
        """
        project_list = list(self._plc.keys())
        if self._force_project_capacity:
            for project in self._plc:
                self._plc[project][0] = self._force_project_capacity
        else:
            # randomly assign remaining project capacities
            for _ in range(self._total_project_capacity - self._num_projects):
                self._plc[random.choice(project_list)][0] += 1

    def _generate_students(self):
        for i in range(self._num_students):
            self._sp[f"s{i + 1}"] = self._get_ordered_list(
                self._student_points, i, random.randint(self._li, self._lj)
            )

    def _generate_lecturers(self):
        lecturer_list = list(self._lp.keys())

        upper_bound_lecturers = self._num_projects // self._num_lecturers
        project_list = list(self._plc.keys())

        for lecturer in self._lp:
            num_projects = random.randint(1, upper_bound_lecturers)
            for _ in range(num_projects):
                p = random.choice(project_list)
                project_list.remove(p)
                self._assign_project_lecturer(p, lecturer)

        # while some projects are unassigned
        while project_list:
            p = random.choice(project_list)
            project_list.remove(p)
            lecturer = random.choice(lecturer_list)
            self._assign_project_lecturer(p, lecturer)

        # decide ordered preference and capacity
        for i, lecturer in enumerate(self._lp):
            ordered_project_list = self._get_ordered_list(self._lecturer_points, i)
            self._lp[lecturer][1] = [
                p for p in ordered_project_list if p in self._lp[lecturer][1]
            ]

            if self._force_lecturer_capacity:
                self._lp[lecturer][0] = self._force_lecturer_capacity
            else:
                self._lp[lecturer][0] = random.randint(
                    self._lp[lecturer][2], self._lp[lecturer][3]
                )

    def sample_all_points(self):
        self._student_points = self._sample_points(self._num_students)
        self._project_points = self._sample_points(self._num_projects)
        self._lecturer_points = self._sample_points(self._num_lecturers)

    def generate_instance(self) -> None:
        """
        Generates a random instance for the SPA-P problem.
        Stores details in self._sp, self._plc, self._lp.
        """
        self.sample_all_points()
        self._reset_instance()
        self._generate_projects()
        self._generate_students()
        self._generate_lecturers()

    def write_instance_to_file(self, filename: str) -> None:
        """
        Writes instances to filename specified.
        """
        if filename.endswith(".txt"):
            delim = " "
        elif filename.endswith(".csv"):
            delim = ","

        with open(filename, "w") as f:
            f.write(
                delim.join(
                    map(
                        str,
                        [self._num_students, self._num_projects, self._num_lecturers],
                    )
                )
                + "\n"
            )

            # student index, preferences
            for student in self._sp:
                f.write(
                    delim.join(
                        map(
                            str,
                            [
                                student[1:],
                                delim.join([p[1:] for p in self._sp[student]]),
                            ],
                        )
                    )
                    + "\n"
                )

            # project index, capacity, lecturer
            for project in self._plc:
                f.write(
                    delim.join(
                        map(
                            str,
                            [
                                project[1:],
                                self._plc[project][0],
                                self._plc[project][1][1:],
                            ],
                        )
                    )
                    + "\n"
                )

            # lecturer index, capacity, projects
            for lecturer in self._lp:
                f.write(
                    delim.join(
                        map(
                            str,
                            [
                                lecturer[1:],
                                self._lp[lecturer][0],
                                delim.join([p[1:] for p in self._lp[lecturer][1]]),
                            ],
                        )
                    )
                    + "\n"
                )

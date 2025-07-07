from copy import deepcopy

import networkx as nx

from readFile import SPASTFileReader


class SPAST_STRONG:
    def __init__(self, filename):
        self.filename = filename

        r = SPASTFileReader(self.filename)
        r.read_file()

        self.num_students = r.students
        self.num_projects = r.projects
        self.num_lecturers = r.lecturers

        self.sp = r.sp
        self.sp_no_tie_deletions = r.sp_no_tie_deletions
        self.plc = r.plc
        self.lp = r.lp

        self.og_sp = deepcopy(self.sp)
        self.og_plc = deepcopy(self.plc)
        self.og_lp = deepcopy(self.lp)

        self.unassigned_and_non_empty_list = list(self.sp.keys())
        self.G = {}
        for student in self.sp.keys():
            self.G[student] = {"projects": set(), "bound": set(), "unbound": set()}
        for project in self.plc.keys():
            self.G[project] = {
                "students": set(),
                "replete": False,
                "rejected": [],
                "num_bound_edges": set(),
                "revised_quota": 0,
            }
        for lecturer in self.lp.keys():
            self.G[lecturer] = {
                "students": set(),
                "projects": set(),
                "replete": False,
                "num_bound_edges": set(),
                "revised_quota": 0,
            }

        self.M = {}
        self.build_Gr = False
        self.usable_Mr = False
        self.Mr_prefigure_flow = {}
        self.max_flow = {}
        self.Zp = set()
        self.Zs = set()

    def pquota(self, project):
        return min(self.plc[project]["cap"], len(self.G[project]["students"]))

    def lquota(self, lecturer):
        alpha_k = sum(
            [self.pquota(project) for project in self.G[lecturer]["projects"]]
        )
        return min(self.lp[lecturer]["cap"], alpha_k)

    def add_edge_to_G(self, student, project, lecturer):
        self.G[student]["projects"].add(project)
        self.G[project]["students"].add(student)
        self.G[lecturer]["students"].add(student)
        self.G[lecturer]["projects"].add(project)

    def remove_edge_from_G(self, student, project, lecturer):
        if project in self.G[student]["projects"]:
            self.G[student]["projects"].remove(project)
            self.G[project]["students"].remove(student)
            self.G[project]["rejected"].append(student)

            if len(self.G[project]["students"]) == 0:
                self.G[lecturer]["projects"].remove(project)

            if (
                len(
                    self.G[student]["projects"].intersection(
                        self.G[lecturer]["projects"]
                    )
                )
                == 0
            ):
                self.G[lecturer]["students"].remove(student)

    def delete(self, student, project, lecturer):
        if project in self.sp_no_tie_deletions[student]:
            p_rank = self.sp[student]["list_rank"][project]
            self.sp[student]["list"][p_rank[0]][p_rank[1]] = "dp"
            self.sp_no_tie_deletions[student].remove(project)

        self.remove_edge_from_G(student, project, lecturer)

        if (
            len(self.G[student]["projects"]) == 0
            and student not in self.unassigned_and_non_empty_list
            and len(self.sp_no_tie_deletions[student]) > 0
        ):
            self.unassigned_and_non_empty_list.append(student)

    def p_dominated_students(self, project):
        qpj = self.pquota(project)
        Lkj_students = self.plc[project]["list"]
        Lkj_tail_index = self.plc[project]["tail_idx"]
        Gpj = self.G[project]["students"]
        dominated_index = None
        dominated_students = []
        count = 0
        for i in range(Lkj_tail_index + 1):
            assigneees = Gpj.intersection(Lkj_students[i])
            count += len(assigneees)
            if count >= qpj:
                self.plc[project]["tail_idx"] = i
                dominated_index = i + 1
                dominated_students = self.plc[project]["list"][dominated_index:]

                return dominated_index, dominated_students

    def l_dominated_students(self, lecturer):
        Lk_students = self.lp[lecturer]["list"]
        qlk = self.lquota(lecturer)
        projects_inG = self.G[lecturer]["projects"]
        p_reduced_quotas = {p: self.pquota(p) for p in projects_inG}
        Lk_tail_index = self.lp[lecturer]["tail_idx"]
        dominated_index = None
        dominated_students = []

        count = 0
        for i in range(Lk_tail_index + 1):
            for s in Lk_students[i]:
                for p in self.G[s]["projects"].intersection(projects_inG):
                    if p_reduced_quotas[p] > 0:
                        count += 1
                        p_reduced_quotas[p] -= 1
                    if count >= qlk:
                        self.lp[lecturer]["tail_idx"] = i
                        dominated_index = i + 1
                        dominated_students = Lk_students[dominated_index:]
                        return dominated_index, dominated_students

    def next_tie_student_head(self, student):
        s_preference = self.sp[student]["list"]
        head_index = self.sp[student]["head_idx"]
        pref_list_length = self.sp[student]["list_len"]
        if head_index < pref_list_length:
            self.sp[student]["head_idx"] += 1
            return s_preference[head_index]
        return None

    def while_loop(self):
        while self.unassigned_and_non_empty_list:
            student = self.unassigned_and_non_empty_list.pop(0)
            tie_at_head = self.next_tie_student_head(student)
            if tie_at_head is not None:
                for project in tie_at_head:
                    if project == "dp":
                        continue
                    lecturer = self.plc[project]["lec"]
                    self.add_edge_to_G(student, project, lecturer)

                    if self.pquota(project) == self.plc[project]["cap"]:
                        self.G[project]["replete"] = True
                        dominated_index, dominated_students = self.p_dominated_students(
                            project
                        )

                        if len(dominated_students) > 0:
                            self.plc[project]["list"] = self.plc[project]["list"][
                                :dominated_index
                            ]

                        for tie in dominated_students:
                            for st in tie:
                                self.delete(st, project, lecturer)

                    if self.lquota(lecturer) == self.lp[lecturer]["cap"]:
                        self.G[lecturer]["replete"] = True
                        dominated_index, dominated_students = self.l_dominated_students(
                            lecturer
                        )

                        self.lp[lecturer]["list"] = self.lp[lecturer]["list"][
                            :dominated_index
                        ]
                        p_k = self.lp[lecturer]["projects"]
                        for tie in dominated_students:
                            for st in tie:
                                a_t = self.sp_no_tie_deletions[st]
                                common_projects = p_k.intersection(a_t)
                                for pu in common_projects:
                                    self.delete(st, pu, lecturer)

            if (
                len(self.G[student]["projects"]) != 0
                and student in self.unassigned_and_non_empty_list
            ):
                self.unassigned_and_non_empty_list.remove(student)

            if (
                len(self.G[student]["projects"]) == 0
                and student not in self.unassigned_and_non_empty_list
                and len(self.sp_no_tie_deletions[student]) > 0
            ):
                self.unassigned_and_non_empty_list.append(student)

    def is_bound(self, student, project, lecturer):
        p_tail_idx = self.plc[project]["tail_idx"]
        l_tail_idx = self.lp[lecturer]["tail_idx"]

        qpj = self.pquota(project)
        qlk = self.lquota(lecturer)

        p_boolean = (len(self.G[project]["students"]) == qpj) or (
            student not in self.plc[project]["list"][p_tail_idx]
        )
        l_boolean = (
            sum([self.pquota(project) for project in self.G[lecturer]["projects"]])
            == qlk
        ) or student not in self.lp[lecturer]["list"][l_tail_idx]

        return p_boolean and l_boolean

    def reset_bound(self):
        for project in self.plc:
            self.G[project]["num_bound_edges"], self.G[project]["revised_quota"] = 0, 0
        for lecturer in self.lp:
            self.G[lecturer]["num_bound_edges"], self.G[lecturer]["revised_quota"] = (
                0,
                0,
            )

    def update_bound_unbound(self):
        self.build_Gr = False
        self.reset_bound()
        for s in self.sp:
            Gs = self.G[s]["projects"]
            bound_projects = set()
            for p in Gs:
                L = self.plc[p]["lec"]
                if self.is_bound(s, p, L):
                    bound_projects.add(p)
                    self.G[p]["num_bound_edges"] += 1
                    self.G[L]["num_bound_edges"] += 1
            unbound_projects = Gs.difference(bound_projects)
            self.G[s]["bound"] = bound_projects
            self.G[s]["unbound"] = unbound_projects
            if len(bound_projects) == 0 and len(unbound_projects) > 0:
                self.build_Gr = True

    def revised_pquota(self, project):
        qpj = self.pquota(project)
        pbound_edges = self.G[project]["num_bound_edges"]
        return qpj - pbound_edges

    def revised_lquota(self, lecturer):
        qlk = self.lquota(lecturer)
        lbound_edges = self.G[lecturer]["num_bound_edges"]
        return qlk - lbound_edges

    def update_revised_quota(self):
        for project in self.plc:
            self.G[project]["revised_quota"] = self.revised_pquota(project)
        for lecturer in self.lp:
            self.G[lecturer]["revised_quota"] = self.revised_lquota(lecturer)

    def buildGr(self):
        Gr = nx.DiGraph()
        for si in self.sp:
            if len(self.G[si]["bound"]) == 0 and len(self.G[si]["unbound"]) > 0:
                Gr.add_edge("s", si, capacity=1)
                for pj in self.G[si]["unbound"]:
                    Gr.add_edge(si, pj, capacity=1)

        for lk in self.lp:
            if self.G[lk]["revised_quota"] > 0:
                for pj in self.G[lk]["projects"]:
                    if self.G[pj]["revised_quota"] > 0:
                        Gr.add_edge(pj, lk, capacity=self.G[pj]["revised_quota"])
                Gr.add_edge(lk, "t", capacity=self.G[lk]["revised_quota"])

        max_flow = nx.max_flow_min_cost(Gr, "s", "t")
        return max_flow

    def unhappy_projects(self):
        Gr_students = set(self.max_flow["s"].keys())
        typeII_Us = set([])
        Up = set([])
        for student in Gr_students:
            if self.max_flow["s"][student] == 1:
                adjacent_projects = set(self.max_flow[student].keys())
                assigned_project = None
                for project in adjacent_projects:
                    if self.max_flow[student][project] == 1:
                        assigned_project = project
                        break

                lecturer = self.plc[assigned_project]["lec"]
                lecturer_projects = set(self.G[lecturer]["projects"])
                intersect = adjacent_projects.intersection(lecturer_projects)
                intersect.remove(assigned_project)

                for project in intersect:
                    if (
                        self.max_flow[project][lecturer]
                        < self.G[project]["revised_quota"]
                    ):
                        typeII_Us.add(student)
                        Up.add(project)

        return Up, typeII_Us

    def criticalset_projects(self, Up):
        unexplored_projects = set([p for p in Up])
        explored_projects, visited_students = set(), set()
        while unexplored_projects:
            project = unexplored_projects.pop()
            lk = self.plc[project]["lec"]
            explored_projects.add(project)
            unmatched_svertices = set([])

            for s in self.max_flow["s"]:
                if (
                    s not in visited_students
                    and project in self.max_flow[s]
                    and self.max_flow[s][project] == 0
                ):
                    unmatched_svertices.add(s)

            for s in unmatched_svertices:
                visited_students.add(s)
                for p in self.max_flow[s]:
                    if (
                        self.max_flow[s][p] == 1
                        and self.plc[p]["lec"] == lk
                        and p not in explored_projects
                    ):
                        unexplored_projects.add(p)

        return explored_projects

    def Zp_deletions(self, typeII_Us):
        for si in typeII_Us:
            nsi_crit = []
            for project in self.max_flow[si].keys():
                if project in self.Zp:
                    nsi_crit.append(project)

            for px in nsi_crit:
                px_lec = self.plc[px]["lec"]
                for py in nsi_crit:
                    py_lec = self.plc[py]["lec"]
                    if px == py:
                        continue
                    if px_lec == py_lec:
                        self.delete(si, px, px_lec)
                        self.delete(si, py, py_lec)

    def unhappy_students(self):
        Gr_students = set(self.max_flow["s"].keys())
        Us = set(
            [student for student in Gr_students if self.max_flow["s"][student] == 0]
        )
        return Us

    def criticalset_students(self, Us):
        unexplored_students = set([s for s in Us])
        explored_students, visited_projects, explored_lecturers = set(), set(), set()
        while unexplored_students:
            student = unexplored_students.pop()
            explored_students.add(student)
            projects_to_explore = set(
                [p for p in self.max_flow[student] if self.max_flow[student][p] == 0]
            )
            new_projects = set()

            for p in projects_to_explore:
                lecturer = self.plc[p]["lec"]
                if (
                    lecturer not in explored_lecturers
                    and self.max_flow[p][lecturer] < self.G[p]["revised_quota"]
                ):
                    explored_lecturers.add(lecturer)
                    PknPr = self.G[lecturer]["projects"].intersection(
                        set(self.max_flow.keys())
                    )
                    saturated_projects = set(
                        [
                            pj
                            for pj in PknPr
                            if self.max_flow[pj][lecturer]
                            == self.G[pj]["revised_quota"]
                        ]
                    )
                    new_projects.update(saturated_projects)

            projects_to_explore.update(new_projects)

            for project in projects_to_explore:
                if project not in visited_projects:
                    visited_projects.add(project)
                    all_matched_students = set()
                    Gr_students = self.G[project]["students"].intersection(
                        set(self.max_flow.keys())
                    )
                    for sr in Gr_students:
                        if self.max_flow[sr][project] == 1:
                            all_matched_students.add(sr)
                    for st in all_matched_students:
                        if st not in explored_students:
                            unexplored_students.add(st)

        return explored_students

    def Zs_deletions(self):
        N_Zs = set()
        for student in self.Zs:
            N_Zs.update(set(self.max_flow[student].keys()))

        for project in N_Zs:
            lecturer = self.plc[project]["lec"]
            Lkj_students = self.plc[project]["list"]
            Lkj_tail_index = self.plc[project]["tail_idx"]
            Lkj_tail = Lkj_students[Lkj_tail_index]
            for st in Lkj_tail:
                self.delete(st, project, lecturer)

            self.plc[project]["list"] = self.plc[project]["list"][:Lkj_tail_index]
            self.plc[project]["tail_idx"] -= 1

    def inner_repeat(self):
        self.usable_Mr = False

        self.Zs, self.Zp = set([None]), set([None])
        while self.Zs.union(self.Zp):
            self.Zs, self.Zp = set(), set()
            self.while_loop()
            self.update_bound_unbound()

            if self.build_Gr:
                self.update_revised_quota()
                self.max_flow = self.buildGr()
                self.Mr_prefigure_flow = deepcopy(self.max_flow)

                ### student ###
                Us = self.unhappy_students()
                self.Zs = self.criticalset_students(Us)
                self.Zs_deletions()

    def most_preferred_reject(self, project):
        rejects = self.G[project]["rejected"]

        for tie in self.og_plc[project]["list"]:
            for student in tie:
                if student in rejects:
                    return student
        raise ValueError("most_preferred_reject found nobody")

    def tail_no_better(self, lecturer, reject):
        lecturer_tail = self.get_lecturer_tail(lecturer)
        found_reject = False
        found_tail = False
        for tie in self.og_lp[lecturer]["list"]:
            for student in tie:
                if student in lecturer_tail:
                    found_tail = True
            if reject in tie:
                found_reject = True

            if found_reject:
                return True
            if found_tail:
                return False

    def tail_worse(self, lecturer, reject):
        lecturer_tail = self.get_lecturer_tail(lecturer)
        found_reject = False
        found_tail = False
        for tie in self.og_lp[lecturer]["list"]:
            for student in tie:
                if student in lecturer_tail:
                    found_tail = True
            if reject in tie:
                found_reject = True

            if found_tail:
                return False
            if found_reject:
                return True

    def get_lecturer_tail(self, lecturer):
        info = self.lp[lecturer]
        return info["list"][info["tail_idx"]]

    def project_neighbourhood_comparison(self, student, project, neighbourhood):
        """
        1 - student prefers project to something in neighbourhood
        0 - student is indifferent between project and something in neighbourhood
        -1 - neither
        """
        student_ranking = self.sp[student]["list_rank"]
        project_rank = student_ranking[project][0]
        indifference = False
        strict_preference = False

        for pn in neighbourhood:
            pn_rank = student_ranking[pn][0]
            if pn_rank > project_rank:
                strict_preference = True
            elif pn_rank == project_rank:
                indifference = True

        if strict_preference and indifference:
            raise ValueError("Oh the horror")

        if strict_preference:
            return 1
        elif indifference:
            return 0
        else:
            return -1

    def get_acceptable_offered_projects(self, lecturer, student):
        results = []
        Pk = self.og_lp[lecturer]["projects"]
        for tie in self.og_sp[student]["list"]:
            for project in tie:
                if project in Pk:
                    results.append(project)
        return results

    def repletion_tail_wipe(self, lk):
        for st in self.get_lecturer_tail(lk):
            for pu in self.get_acceptable_offered_projects(lk, st):
                self.delete(st, pu, lk)

    def repletion_deletions(self):
        for pj, pj_info in self.plc.items():
            if self.G[pj]["replete"] and self.pquota(pj) < self.plc[pj]["cap"]:
                lk = pj_info["lec"]
                sr = self.most_preferred_reject(pj)
                nsr = self.G[sr]["projects"]

                if len(nsr) == 0:
                    if self.tail_no_better(lk, sr):
                        self.repletion_tail_wipe(lk)
                else:
                    comp = self.project_neighbourhood_comparison(sr, pj, nsr)
                    if comp == 1:
                        if self.tail_no_better(lk, sr):
                            self.repletion_tail_wipe(lk)

                    if comp == 0:
                        if self.tail_worse(lk, sr):
                            self.repletion_tail_wipe(lk)

    def get_feasible_matching(self):
        Gf = nx.DiGraph()
        for si in self.sp:
            Gf.add_edge("s", si, capacity=1)

            if self.G[si]["bound"]:
                for pj in self.G[si]["bound"]:
                    Gf.add_edge(si, pj, capacity=1)

            elif self.G[si]["unbound"]:
                for pj in self.G[si]["unbound"]:
                    Gf.add_edge(si, pj, capacity=1)

        for lk in self.lp:
            for pj in self.G[lk]["projects"]:
                Gf.add_edge(pj, lk, capacity=self.plc[pj]["cap"])
            Gf.add_edge(lk, "t", capacity=self.lp[lk]["cap"])

        self.max_flow = nx.max_flow_min_cost(Gf, "s", "t")
        self.max_flow_to_feasible_matching()

    def max_flow_to_feasible_matching(self):
        self.M = {s: "" for s in self.sp}
        for s in self.sp:
            for p, flow in self.max_flow[s].items():
                if flow == 1:
                    self.M[s] = p

    def get_prefigure_project(self, student):
        if student in self.Mr_prefigure_flow:
            for k, v in self.Mr_prefigure_flow[student].items():
                if v == 1:
                    return k
        return None

    def run(self):
        while self.unassigned_and_non_empty_list:
            self.inner_repeat()  # lines 2 - 26
            self.repletion_deletions()  # lines 27 - 34

        self.update_bound_unbound()
        self.get_feasible_matching()

        return self.M

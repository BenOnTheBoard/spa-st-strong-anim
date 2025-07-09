#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 23:20:17 2019

@author: sofiat
"""

from copy import deepcopy

from readinputSPAST import READSPAST

from genReader import generator_to_dicts, generator_to_other_bruteforce_args


class STSMBruteForce:
    def __init__(self, filename=None, generator=None):
        if filename is not None:
            self.filename = filename

            r = READSPAST(self.filename)
            r.read_file()

            self.sp = r.sp
            self.sp_no_tie = r.sp_no_tie
            self.sp_no_tie_deletions = r.sp_no_tie_deletions
            self.plc = r.plc
            self.lp = r.lp

            self.lp_rank = r.lp_rank
            self.proj_rank = r.proj_rank

        elif generator is not None:
            self.sp, self.sp_no_tie_deletions, self.plc, self.lp = generator_to_dicts(
                generator
            )
            self.sp_no_tie, self.lp_rank, self.proj_rank = (
                generator_to_other_bruteforce_args(generator)
            )

        else:
            raise ValueError("Neither initialisation method supplied.")

        self.num_students = len(self.sp)
        self.num_projects = len(self.plc)
        self.num_lecturers = len(self.lp)

        self.sp_copy = deepcopy(self.sp)
        self.lp_copy = deepcopy(self.lp)

        self.M = {s: "" for s in self.sp}
        self.project_wstcounter = {
            "p" + str(i): [0, []] for i in range(1, len(self.plc) + 1)
        }
        self.lecturer_wstcounter = {
            "l" + str(i): [0, []] for i in range(1, len(self.lp) + 1)
        }

        self.blocking_pair = False
        self.ssm_list = []

    def blockingpair_1bi(self, student, project, lecturer):
        if self.plc[project][1] > 0 and self.lp[lecturer][0] > 0:
            return True
        return False

    def blockingpair_1bii(self, student, project, lecturer):
        if self.plc[project][1] > 0 and self.lp[lecturer][0] == 0:
            proj_in_M = self.M[student]
            if proj_in_M != "" and self.plc[proj_in_M][0] == lecturer:
                return True
            lec_worst_pointer = self.lecturer_wstcounter[lecturer][0]
            student_rank_Lk = self.lp_rank[lecturer][student]
            if student_rank_Lk <= lec_worst_pointer:
                return True
        return False

    def blockingpair_1biii(self, student, project, lecturer):
        if self.plc[project][1] == 0:
            proj_worst_pointer = self.project_wstcounter[project][0]
            student_rank_Lkj = self.proj_rank[project][student]
            if student_rank_Lkj <= proj_worst_pointer:
                return True
        return False

    def blockingpair_2bi(self, student, project, lecturer):
        if self.plc[project][1] > 0 and self.lp[lecturer][0] > 0:
            return True
        return False

    def blockingpair_2bii(self, student, project, lecturer):
        if self.plc[project][1] > 0 and self.lp[lecturer][0] == 0:
            lec_worst_pointer = self.lecturer_wstcounter[lecturer][0]
            student_rank_Lk = self.lp_rank[lecturer][student]
            if student_rank_Lk < lec_worst_pointer:
                return True
        return False

    def blockingpair_2biii(self, student, project, lecturer):
        if self.plc[project][1] == 0:
            proj_worst_pointer = self.project_wstcounter[project][0]
            student_rank_Lkj = self.proj_rank[project][student]
            if student_rank_Lkj < proj_worst_pointer:
                return True
        return False

    def check_stability(self):
        self.blocking_pair = False
        for student in self.sp:
            preferred_projects, indifference = [], []
            if self.M[student] == "":
                preferred_projects = self.sp_no_tie[student]
            else:
                matched_project = self.M[student]
                rank_matched_project = self.sp_copy[student][2][matched_project][0]
                A_si = self.sp_copy[student][1]
                preferred_projects = [
                    s for tie in A_si[:rank_matched_project] for s in tie
                ]
                indifference = A_si[rank_matched_project][:]
                indifference.remove(matched_project)

            for project in preferred_projects:
                lecturer = self.plc[project][0]
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_1bi(
                        student, project, lecturer
                    )
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_1bii(
                        student, project, lecturer
                    )
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_1biii(
                        student, project, lecturer
                    )
                if self.blocking_pair:
                    break

            for project in indifference:
                lecturer = self.plc[project][0]
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_2bi(
                        student, project, lecturer
                    )
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_2bii(
                        student, project, lecturer
                    )
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_2biii(
                        student, project, lecturer
                    )
                if self.blocking_pair:
                    break

            if self.blocking_pair:
                break

    def choose(self, i=1):
        if i > self.students:
            for project in self.plc:
                if self.project_wstcounter[project][1] != []:
                    self.project_wstcounter[project][0] = max(
                        self.project_wstcounter[project][1]
                    )
            for lecturer in self.lp:
                if self.lecturer_wstcounter[lecturer][1] != []:
                    self.lecturer_wstcounter[lecturer][0] = max(
                        self.lecturer_wstcounter[lecturer][1]
                    )

            self.check_stability()
            if self.blocking_pair is False:
                self.ssm_list.append(deepcopy(self.M))

        else:
            student = "s" + str(i)
            for project in self.sp_no_tie[student]:
                lecturer = self.plc[project][0]
                if self.plc[project][1] > 0 and self.lp[lecturer][0] > 0:
                    self.M[student] = project

                    self.plc[project][1] -= 1
                    self.lp[lecturer][0] -= 1

                    student_rank_Lk = self.lp_rank[lecturer][student]
                    student_rank_Lkj = self.proj_rank[project][student]
                    self.project_wstcounter[project][1].append(student_rank_Lkj)
                    self.lecturer_wstcounter[lecturer][1].append(student_rank_Lk)

                    self.choose(i + 1)

                    self.M[student] = ""
                    student_rank_Lk = self.lp_rank[lecturer][student]
                    student_rank_Lkj = self.proj_rank[project][student]

                    self.project_wstcounter[project][1].remove(student_rank_Lkj)
                    self.lecturer_wstcounter[lecturer][1].remove(student_rank_Lk)

                    self.plc[project][1] += 1
                    self.lp[lecturer][0] += 1
            self.choose(i + 1)

    def get_ssm_list(self):
        return self.ssm_list

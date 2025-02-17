#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:47:18 2023

@author: sofiat
"""


class SPASTFileReader:
    def __init__(self, filename):
        self.filename = filename
        self.students = 0
        self.projects = 0
        self.lecturers = 0  # assume number of lecturers <= number of projects
        self.sp = {}
        self.sp_no_tie = {}
        self.sp_no_tie_deletions = {}
        self.plc = {}
        self.lp = {}
        self.lp_rank = {}
        self.proj_rank = {}

    def read_file(self):  # reads the SPAT instance
        with open(self.filename) as t:
            t = t.readlines()
        entry1 = t[0].strip().split(" ")
        self.students, self.projects, self.lecturers = (
            int(entry1[0]),
            int(entry1[1]),
            int(entry1[2]),
        )

        # -------------------------------------------------------------------------------------------------------------------
        #  we build the student's dictionary
        for i in range(1, self.students + 1):
            entry = t[i].rstrip(" \n").split(" ")
            student = "s" + str(entry[0])
            preferencelist = []
            self.sp_no_tie[student] = []
            self.sp_no_tie_deletions[student] = (
                set()
            )  # because removing from a set is constant time
            for k in entry[1:]:
                if k[0] == "(":
                    split_tie = k.rstrip(")").lstrip("(").split(":")
                    tie = ["p" + str(j) for j in split_tie]
                    preferencelist.append(tie)
                    self.sp_no_tie[student].extend(tie)
                else:
                    preferencelist.append(["p" + str(k)])
                    self.sp_no_tie[student].append("p" + str(k))
            self.sp_no_tie_deletions[student] = set(
                self.sp_no_tie[student]
            )  # because removing from a set is constant time
            length = len(preferencelist)
            # ---------------------------------------------------------------------------------------------------------------------
            rank = {}  # store the index of each project on each student's preference list, for ease of deletion later on in the allocation
            for i in range(length):
                orderedpreference = preferencelist[i]
                count_tie = 0
                for p in orderedpreference:
                    rank[p] = (i, count_tie)
                    count_tie += 1
            self.sp[student] = {
                "list": preferencelist,
                "list_rank": rank,
                "list_len": length,
                "head_idx": 0,
            }
            # the last entry is to keep track of when a student has an empty preference list - to be incremented and compared with length

        # # -------------------------------------------------------------------------------------------------------------------
        # #  we build the projects's dictionary

        for i in range(self.students + 1, self.students + 1 + self.projects):
            entry = t[i].strip().split(" ")
            # project = [lecturer, project_capacity_yet_to_be_filled, full(project) = False, keep track of students that was rejected from project]
            # length of the preferred students for p_j according to l_k to be appended when we have more information..
            self.plc["p" + str(entry[0])] = {
                "cap": int(entry[1]),
                "lec": f"l{entry[2]}",
            }
        # # -------------------------------------------------------------------------------------------------------------------

        # # -------------------------------------------------------------------------------------------------------------------
        #  we build the lecturer's dictionary
        for i in range(
            self.students + 1 + self.projects,
            self.students + 1 + self.projects + self.lecturers,
        ):
            entry = t[i].strip().split(" ")
            lecturer = "l" + str(entry[0])
            self.lp_rank[lecturer] = {}  # stores rank of each student in L_k
            lecturerList = []
            rank = 1
            for k in entry[2:]:
                if k[0] == "(":
                    split_tie = k.rstrip(")").lstrip("(").split(":")
                    tie = ["s" + str(j) for j in split_tie]
                    lecturerList.append(tie)
                    for s in tie:
                        self.lp_rank[lecturer][s] = rank
                        # self.lp[lecturer]["list_rank"][s] = rank
                else:
                    s = "s" + str(k)
                    lecturerList.append([s])
                    self.lp_rank[lecturer][s] = rank
                    # self.lp[lecturer]["list_rank"][s] = rank

                rank += 1
            length = len(lecturerList)
            self.lp[lecturer] = {
                "cap": int(entry[1]),
                "list": lecturerList,
                "list_len": length,
                "tail_idx": length - 1,
            }
            # -------------------------------------------------------------------------------------------------------------------
            #  another useful dictionary is created here and attached to the lecturer's dictionary - L_k_j
            #  the lecturer's ranked preference list according to each project they offer
            # the students preferred by this lecturer   ------ [['s7'], ['s4', 's1'], ['s2'], ['s3', 's5', 's6']] for lecturer 1
            lec_projects = set()
            for project in self.plc:
                if self.plc[project]["lec"] == lecturer:
                    lec_projects.add(project)
                    Lkj = []
                    self.proj_rank[project] = {}  # stores rank of each student in L_k_j
                    rank = 1
                    for tie in lecturerList:  # Is it the case that a student will not be on the ranked preference of a lecturer - same project in common
                        student_tie = []
                        for student in tie:
                            if project in self.sp_no_tie_deletions[student]:
                                student_tie.append(student)
                                self.proj_rank[project][student] = rank
                        if len(student_tie) == 0:
                            continue
                        else:
                            Lkj.append(student_tie)
                        rank += 1
                    self.plc[project]["list"] = Lkj
                    self.plc[project]["list_len"] = len(Lkj)
                    self.plc[project]["tail_idx"] = len(Lkj) - 1
            self.lp[lecturer]["projects"] = lec_projects


# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------
# s = SPASTFileReader("caldam.txt")
# # s = SPASTFileReader("ex7.txt")
# s.read_file()

# for i in s.sp:
#     print(f"{i} :::> {s.sp[i]}")
#     #print(f"{i} :::> {s.sp[i]} :::> {s.sp_no_tie_deletions[i]}")
# print()
# for p in s.plc:
#     print(f"{p} :::> {s.plc[p]}")
# print()
# for l in s.lp:
#     print(f"{l} :::> {s.lp[l]}")


"""  
=== Output for SPA-ST instance in caldam.txt

s1 :::> {'list': [['p1'], ['p6']], 'list_rank': {'p1': (0, 0), 'p6': (1, 0)}, 'list_len': 2, 'head_idx': 0}
s2 :::> {'list': [['p1'], ['p2']], 'list_rank': {'p1': (0, 0), 'p2': (1, 0)}, 'list_len': 2, 'head_idx': 0}
s3 :::> {'list': [['p1'], ['p4']], 'list_rank': {'p1': (0, 0), 'p4': (1, 0)}, 'list_len': 2, 'head_idx': 0}
s4 :::> {'list': [['p2'], ['p5', 'p6']], 'list_rank': {'p2': (0, 0), 'p5': (1, 0), 'p6': (1, 1)}, 'list_len': 2, 'head_idx': 0}
s5 :::> {'list': [['p2', 'p3']], 'list_rank': {'p2': (0, 0), 'p3': (0, 1)}, 'list_len': 1, 'head_idx': 0}
s6 :::> {'list': [['p2', 'p4']], 'list_rank': {'p2': (0, 0), 'p4': (0, 1)}, 'list_len': 1, 'head_idx': 0}
s7 :::> {'list': [['p3'], ['p1']], 'list_rank': {'p3': (0, 0), 'p1': (1, 0)}, 'list_len': 2, 'head_idx': 0}
s8 :::> {'list': [['p5'], ['p1']], 'list_rank': {'p5': (0, 0), 'p1': (1, 0)}, 'list_len': 2, 'head_idx': 0}

p1 :::> {'cap': 2, 'lec': 'l1', 'list': [['s8'], ['s7'], ['s1', 's2', 's3']], 'list_len': 3, 'tail_idx': 2}
p2 :::> {'cap': 2, 'lec': 'l1', 'list': [['s2'], ['s4', 's5'], ['s6']], 'list_len': 3, 'tail_idx': 2}
p3 :::> {'cap': 1, 'lec': 'l2', 'list': [['s5'], ['s7']], 'list_len': 2, 'tail_idx': 1}
p4 :::> {'cap': 1, 'lec': 'l2', 'list': [['s6'], ['s3']], 'list_len': 2, 'tail_idx': 1}
p5 :::> {'cap': 1, 'lec': 'l3', 'list': [['s4'], ['s8']], 'list_len': 2, 'tail_idx': 1}
p6 :::> {'cap': 2, 'lec': 'l3', 'list': [['s1', 's4']], 'list_len': 1, 'tail_idx': 0}

l1 :::> {'cap': 3, 'list': [['s8'], ['s7'], ['s1', 's2', 's3'], ['s4', 's5'], ['s6']], 'list_len': 5, 'tail_idx': 4, 'projects': {'p1', 'p2'}}
l2 :::> {'cap': 2, 'list': [['s6'], ['s5'], ['s7', 's3']], 'list_len': 3, 'tail_idx': 2, 'projects': {'p4', 'p3'}}
l3 :::> {'cap': 3, 'list': [['s1', 's4'], ['s8']], 'list_len': 2, 'tail_idx': 1, 'projects': {'p6', 'p5'}}
"""

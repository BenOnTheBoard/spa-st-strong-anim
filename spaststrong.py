#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last edited on Sat Feb  1 20:20:49 2020

@author: sofiat
"""


from readinputSPAST import READSPAST
from hk import HopcroftKarp
from copy import deepcopy


class SPASTSTRONG:
    def __init__(self, input):

        self.filename = input
        r = READSPAST()
        r.read_file(self.filename)

        self.students = r.students
        self.projects = r.projects
        self.lecturers = r.lecturers

        self.sp = r.sp
        self.sp_copy = r.sp_copy
        self.sp_no_tie = r.sp_no_tie
        self.sp_no_tie_deletions = r.sp_no_tie_deletions
        self.plc = r.plc
        self.lp = r.lp
        self.lp_copy = r.lp_copy
        

        # -------------------------------------------------------------------------------------------------------------------
        self.unassigned_empty_list = [] # keeps track of unassigned students who has an empty list
        # =======================================================================
        # initialising data structure for each student, project, and lecturer in G
        # =======================================================================
        self.G = {} #provisional assignment graph
        for student in r.sp:
            self.unassigned_empty_list.append(student)
            self.G[student] = [set(), set(), set()] # [set of p_j's adjacent to s_i, bound p_j's, unbound p_j's]
        for project in r.plc:
            self.G[project] = [set(), r.plc[project][1], False] # [students assigned to p_j, remnant of c_j --- I don't think I care about remnant of cj or dk]
        for lecturer in r.lp_copy:
            self.G[lecturer] = [set(), set(), r.lp_copy[lecturer][0], False] # [students assigned to l_k, non-empty p_j's in G offered by l_k, remnant of d_k, replete]
        
        self.M = {}
        self.M_w_proj_lec = {}
        self.Zs = set()
        self.Zp = set()
        # -------------------------------------------------------------------------------------------------------------------
        self.not_assigned = 0
        self.assigned = 0
        self.student_check = False
        self.project_check = False
        self.lecturer_check = False
        self.blocking_pair = False
        self.full_projects = set()
        self.found_stsm = 'N'

#        print('===================================')
#        for k,v in self.lp.items():
#            print(k, '::>', v)
#        print('===================================')
#        for k,v in self.plc.items():
#            print(k, '::>', v)
#        print('===================================')
        
    def pquota(self, project): # q_j = min(c_j, d_G(p_j))
        return min(self.plc[project][1], len(self.G[project][0]))
        
    
    def lquota(self, lecturer): # q_k = min(d_k, alpha_k)
        alpha_k = sum([self.pquota(project) for project in self.G[lecturer][1] if len(self.G[project][0]) != 0]) # can the if statement here be removed if self.G[lecturer][1] is updated?
        return min(self.lp[lecturer][0], alpha_k)
        
    # =======================================================================
    # add edge (s_i, p_j) to G
    # =======================================================================    
    def add_edge_to_G(self, student, project, lecturer):
        self.G[student][0].add(project) 
        self.G[project][0].add(student)
        #self.G[project][1] -= 1  # reduce c_j
        # dk is reduced for each student even if the student is already assigned to a project offered by lk        
        self.G[lecturer][1].add(project)
        self.G[lecturer][0].add(student)
        #self.G[lecturer][2] -= 1  # reduce d_k        
        
        # if student not in self.G[lecturer][0]:
        #     self.G[lecturer][0].add(student)
        #     self.G[lecturer][2] -= 1  # reduce d_k


    # =======================================================================
    # remove edge (s_i, p_j) from G
    # =======================================================================       
    def remove_edge_from_G(self, student, project, lecturer):
        if project in self.G[student][0]:
            self.G[student][0].remove(project) 
            # discard is fine for bound and unbound edges
            self.G[student][1].discard(project)
            self.G[student][2].discard(project)            
            
            self.G[project][0].remove(student)
            #self.G[project][1] += 1  # increment c_j             
            # if the project becomes an isolated vertex
            if len(self.G[project][0]) == 0:
                self.G[lecturer][1].remove(project)
            
            # if in G, student no longer has any project in common w/ lecturer
            if len(self.G[student][0].intersection(self.G[lecturer][1])) == 0:
                self.G[lecturer][0].remove(student)        
                #self.G[lecturer][2] += 1  # increment d_k'
    
  
    # =======================================================================
    # delete (s_i, p_j) from A(s_i)   ---! and not L_k^j
    # =======================================================================
    def delete(self, student, project, lecturer):       
        # replace project with dp in the student's preference list
        # also remove it from sp_no_tie_deletions
        if project in self.sp_no_tie_deletions[student]:
            # get the rank of pj on the student's list
            # say (2, 0): 2 is the position of p_j's tie T in A(s_i)
            # and 0 is the position of p_j in T
            p_rank = self.sp[student][2][project]
            self.sp[student][1][p_rank[0]][p_rank[1]] = 'dp'  # we replace the project in this position by a dummy project
            self.sp_no_tie_deletions[student].remove(project)
            self.plc[project][3].append(student)  # keep track of students who were rejected from pj
        
        # remove the edge from G
        self.remove_edge_from_G(student, project, lecturer)
            
        # if student is not adjacent to an edge in G and has a non-empty list
        # add her to the list of unassigned students
        if len(self.G[student][0]) == 0 and student not in self.unassigned_empty_list and len(self.sp_no_tie_deletions[student]) > 0:
            self.unassigned_empty_list.append(student) 
       
        
    # =======================================================================
    # find dominated students in L_k^j
    # =======================================================================    
    def p_dominated_students(self, project):
        """
        param: p_j
        return: starting index of dominated students in Lkj as well as the students
        """
        lecturer = self.plc[project][0]
        cj = self.plc[project][1]
        Lkj_students = self.lp[lecturer][2][project]  # students who chose p_j according to Lk
        Lkj_tail_index = self.plc[project][5] # self.plc[project][4] is length of list while self.plc[project][5] is initial tail index
        Gpj = self.G[project][0]
        dominated_index = None
        dominated_students = []
        count = 0 # this will increment the number of edges adjacent to students who are better than dominated students in Lkj
        for i in range(Lkj_tail_index+1):
            assigneees = Gpj.intersection(Lkj_students[i])
            count += len(assigneees)
            if count >= cj:
                self.plc[project][5] = i # the updated tail of Lkj is now the tie in position i, again 0-based
                dominated_index = i+1
                dominated_students = self.lp[lecturer][2][project][dominated_index:]
                # print('L_k_j dominated index is: ', dominated_index)
                # print('L_k_j dominated student is: ', dominated_students)
                return dominated_index, dominated_students
  
    
    # =======================================================================
    # find dominated students in L_k
    # =======================================================================
    def l_dominated_students(self, lecturer):
        """
        param: l_k
        return: starting index of dominated students in Lk as well as those students
        """
        Lk_students = self.lp[lecturer][1]  # students in L_k  
        dk = self.lp[lecturer][0]
        #Glk = self.G[lecturer][0]
        projects_inG = self.G[lecturer][1]
        p_reduced_quotas = {p:self.pquota(p) for p in projects_inG}
#        print(p_reduced_quotas)
        Lk_tail_index = self.lp[lecturer][5]
        dominated_index = None
        dominated_students = []
        count = 0
        for i in range(Lk_tail_index+1):
            for s in Lk_students[i]:
                for p in self.G[s][0].intersection(projects_inG):
                    if p_reduced_quotas[p] > 0:
                        count += 1
                        p_reduced_quotas[p] -= 1
                    # print(lecturer, s, p, count, i)
                    if count >= dk:
                        self.lp[lecturer][5] = i
                        dominated_index = i+1
                        dominated_students = Lk_students[dominated_index:]
                        # print('L_k dominated index is: ', dominated_index)
                        # print('L_k dominated student is: ', dominated_students)
                        return dominated_index, dominated_students
   
                
        
        
    
    # =======================================================================
    # get the next set of projects tied together at the head of a student's preference list
    # ======================================================================= 
    def next_tie_student_head(self, student):
        s_preference = self.sp[student][1]  # the projected preference list for this student.. this changes during the allocation process.
            
        # self.sp[student][3] points to the tie at the head of s_i's list 
        # if tie pointer is not up to length of pref list --- length is not 0-based!
        head_index = self.sp[student][3] # 0-based
        pref_list_length = self.sp[student][0] # !0-based
        if  head_index < pref_list_length:
            self.sp[student][3] += 1  # we increment head_index pointer --> moves inward by 1
            return s_preference[head_index] # projects at the head of the list ::: could be length 1 or more
        return None
             
    # =======================================================================
    # while loop that constructs G from students preference lists
    # =======================================================================    
    def while_loop(self):
        while self.unassigned_empty_list:            
            student = self.unassigned_empty_list.pop(0)  
            #print('unassigned students: ' , self.unassigned_empty_list)
            tie_at_head = self.next_tie_student_head(student)
            if tie_at_head is not None: # a None value here implies that the student has an empty preference list
                for project in tie_at_head:
                    if project == 'dp': continue
                    else:                        
                        lecturer = self.plc[project][0]
                        # add the edge (student, project) to G 
                        self.add_edge_to_G(student, project, lecturer)
                        # ----------- if the quota of pj in G  is equal to cj ----------
                        if self.pquota(project) == self.plc[project][1]:  
                            self.G[project][2] = True   # set replete(p_j) to True
                            dominated_index, dominated_students = self.p_dominated_students(project) # finds dominated students and the starting index on L_k^j
                            # we delete dominated students from L_k^j by replacing L_k_j with non-dominated students
                            if len(dominated_students) > 0:
                                self.lp[lecturer][2][project] = self.lp[lecturer][2][project] [:dominated_index] 
                                #print('remaining L_k^j: ', self.lp[lecturer][2][project])
                            #print()
                            # for each dominated student, delete (student, project)                          
                            for tie in dominated_students:
                                for st in tie:
                                    self.delete(st, project, lecturer)
                                    # print('line10 delete', st, project, lecturer)

                        # ----------- if the quota of lecturer is equal to dk -----------
                        if self.lquota(lecturer) == self.lp[lecturer][0]:
                            # print(lecturer, " replete")
                            self.G[lecturer][3] = True # set replete(l_k) to True 
                            dominated_index, dominated_students = self.l_dominated_students(lecturer) # finds dominated students and the starting index on l_k
                            
                            # print(lecturer, dominated_index, dominated_students, self.G[lecturer])
                            # print(self.lp[lecturer])
                            # print()
#                                for k,v in self.G.items():
#                                    print(k, '::>', v)
                            
                            # we delete dominated students from l_k by replacing l_k with non-dominated students
                            self.lp[lecturer][1] = self.lp[lecturer][1][:dominated_index]
                            p_k = set([i for i in self.lp[lecturer][2].keys()])  # all the projects that lecturer is offering
                            for tie in dominated_students:
                                for st in tie:
                                    a_t = set(self.sp_no_tie_deletions[st])  # the student's altered preference list without ties..
                                    common_projects = p_k.intersection(a_t)
                                    for pu in common_projects:
                                        self.delete(st, pu, lecturer)
                                        # print('lne 14 delete', st,pu,lecturer)
    #                        print('***********************************')
                        # self.update_project_tail() # !important --- do not delete.
    #                        print('***********************************')
    #                        print(lecturer, self.lp[lecturer][1])
    #                        for k,v in s.G.items():
    #                            print('\t', k, '::>', v)
    #                        print('***********************************')

            #########################################################################################################################################
            if len(self.G[student][0]) != 0 and student in self.unassigned_empty_list:
                self.unassigned_empty_list.remove(student)
                
            # !* if the current student is unassigned in the matching, with a non-empty preference list, we re-add the student to the unassigned list
            if len(self.G[student][0]) == 0 and student not in self.unassigned_empty_list and len(self.sp_no_tie_deletions[student]) > 0:  #--- caught in an infinite while loop for tie-9 (check later!**)
            #if self.M[student] == set() and self.sp[student][3] < self.sp[student][0]:
                self.unassigned_empty_list.append(student)
            #########################################################################################################################################
    
    
    # =======================================================================
    # update lecturer tail pointer 
    # =======================================================================
    def update_lecturer_tail(self):
        for lecturer in self.lp:
            alpha_k = sum([self.pquota(project) for project in self.G[lecturer][1] if len(self.G[project][0]) != 0]) # can the if statement here be removed if self.G[lecturer][1] is updated?
            d_k = self.lp[lecturer][0]
            # we only update the tail if alpha_k >= d_k is full or oversubscribed
#            if self.G[lecturer][2] <= 0: # caught in an infinite loop with because tail could be not be updated and it forced the critical set to be non-empty
            if alpha_k >= d_k:
                # find the pointer of the worst assigned tie on L_k
                Lk_students = self.lp[lecturer][1]  # students in L_k        
                found = False
                tail_index = self.lp[lecturer][4]
#                print(lecturer, self.lp[lecturer])
#                for k,v in self.G.items():
#                    print('\t', k, '::>', v)
#                print('***********************************')
                while not found and tail_index >= 0:   # stop as soon as a worst tie containing a student assigned to l_k is found
                    worst_student_tie = Lk_students[tail_index]  # tied students at the tail of L_k^j
                    for find_worst_student in worst_student_tie:
                        if find_worst_student in self.G[lecturer][0]:
                            found = True                    
                            break   # exit the for loop as soon as a worst student assigned to p_j is found
                    # if no worst assigned student is found, we decrement the project worst counter and repeat the process again..
                    if found is False:
                        tail_index -= 1
                        
                if found:     
                    self.lp[lecturer][4] = tail_index
                    self.lp[lecturer][1] = self.lp[lecturer][1][:tail_index+1]
                
    # =======================================================================
    # update project tail pointer - why do we still need this? find confidence
    # =======================================================================
    def update_project_tail(self):
        for project in self.plc:
            # tail update can only be done if project is full or oversubscribed
            if self.G[project][1] <= 0:
                lecturer = self.plc[project][0]
                Lkj_students = self.lp[lecturer][2][project][:]  # students who chose p_j according to Lk -
#                print(project, Lkj_students, self.plc[project])
                found = False
                tail_index = self.plc[project][4]
                while not found and tail_index >= 0:   # exit the loop as soon as a worst student assigned to p_j is found
                    worst_student_tie = Lkj_students[tail_index]  # tied students at the tail of L_k^j
                    for find_worst_student in worst_student_tie:
                        if find_worst_student in self.G[project][0]:
                            found = True                    
                            break   # exit the for loop as soon as a worst student assigned to p_j is found
                    # if no worst assigned student is found, we decrement the project worst counter and repeat the process again..
                    if not found:
                        tail_index -= 1  
                
                if found:
                    self.plc[project][4] = tail_index
                    self.lp[lecturer][2][project] = self.lp[lecturer][2][project][:tail_index+1]
                    

    # =======================================================================
    # is the edge a lower rank edge? returns boolean
    # =======================================================================
    def is_lower_rank_edge(self, student, project, lecturer):
        
        l_preflist = self.lp[lecturer][1][:]
        alpha_k = sum([self.pquota(project) for project in self.G[lecturer][1] if len(self.G[project][0]) != 0]) # can the if statement here be removed if self.G[lecturer][1] is updated?
        d_k = self.lp[lecturer][0]
#        if project == 'p2' and student == 's6':
#            print(l_preflist, alpha_k)
        if student in l_preflist[-1] and  alpha_k > d_k:
            return True
        
        return False
    
    
    # =======================================================================
    # is the edge bound? returns boolean
    # =======================================================================
    def is_bound(self, student, project, lecturer):   
        p_preflist = self.lp[lecturer][2][project][:]
        l_preflist = self.lp[lecturer][1][:]     
        # (number of edges adjacent to p_j in G is <= cj) or s_i is not in the tail of L_k_j (or both)
        p_boolean = len(self.G[project][0]) <= self.plc[project][1] or student not in p_preflist[-1]
#        if project == 'p3' and student == 's1':
#            print(project, p_preflist, p_boolean, self.G[project])
        # sum_{pj in Pk} (quota of pj) <= dk or (s_i is not in the tail of L_k) (or both)
        l_boolean = sum([self.pquota(project) for project in self.G[lecturer][1]]) <= self.lp[lecturer][0] or student not in l_preflist[-1]
        # if project == 'p2' and student == 's4':
        #     print(p_boolean, l_boolean, l_preflist[-1])
#            print(self.is_lower_rank_edge(student, project, 'l2'))
        return p_boolean and l_boolean
            
    # =======================================================================
    # update bound and unbound projects (edges) for each assigned student in self.G
    # ======================================================================= 
    def update_bound_unbound(self):
        for s in self.sp:
            Gs = self.G[s][0]
            bound_projects = set()
            for p in Gs:
                l = self.plc[p][0]
                if self.is_bound(s, p, l):
                    bound_projects.add(p)
            unbound_projects = Gs.difference(bound_projects)
            self.G[s][1] = bound_projects
            self.G[s][2] = unbound_projects
                    
 
filename = "caldam.txt"
I = SPASTSTRONG(filename)
I.while_loop()
I.update_bound_unbound()
for s in I.sp:
    print(f"{s} ;;;> {I.G[s]}")
print()
for p in I.plc:
    print(f"{p} ;;;> {I.G[p]}")
print()
for l in I.lp:
    print(f"{l} ;;;> {I.G[l]}")

"""
    # =======================================================================
    # lower rank edges: returns set of lower rank edges for the lecturer
    # =======================================================================
    def lower_rank_edges(self, graph, lecturer):
        l_preflist = self.lp[lecturer][1][:]
        #students at the tail of lecturer that are also in Gr
        
        tail = graph[lecturer][0].intersection(l_preflist[-1])
        alpha_k = sum([self.pquota(project) for project in self.G[lecturer][1] if len(self.G[project][0]) != 0]) # can the if statement here be removed if self.G[lecturer][1] is updated?
        d_k = self.lp[lecturer][0]
#        dG_k, d_k = len(self.G[lecturer][0]), self.lp[lecturer][0]
#        print(lecturer, graph[lecturer], alpha_k, d_k)
#        print()
        if alpha_k > d_k:            
            # projects offered by l_k in Gr that are assigned to a student in her tail
            lower_rank_projects = set()
            for s in tail:
                for p in graph[s][0].intersection(graph[lecturer][1]):                   
                    lower_rank_projects.add(p)
            #lower_rank_projects = set([p for p in graph[s][0].intersection(graph[lecturer][1]) for s in tail])
            #print(lower_rank_projects)        
            return lower_rank_projects
        
    
    # =======================================================================
    # form reduced assignment graph Gr
    # =======================================================================
    def reduced_assignment_graph(self):
       
        # to form Gr, first delete isolated student, project and lecturer vertices in G
#        print(self.G)
        Gr = deepcopy(self.G)
        Gr = {k:v for k,v in Gr.items() if v[0] != set()} 
        students_in_Gr = set([s for s in self.sp if s in Gr])
        projects_in_Gr = set([p for p in self.plc if p in Gr])
        lecturers_in_Gr = set([l for l in self.lp if l in Gr])
        # make sure last elt for each project and lecturer in Gr is their true quota
        for p in projects_in_Gr:
            Gr[p].append(self.pquota(p))
            # Gr[p] = [Gr(p_j), remnant of c_j, q_j]
#            cj = self.plc[p][1]            
#            Gr[p][1] = min(len(Gr[p][0]), cj)
        for l in lecturers_in_Gr:
            Gr[l].append(self.lquota(l))
            # Gr[l] = [Gr(l_k), P_k n P_r, remnant of d_k, q_k]
#            dk = self.lp[l][0]                
#            Gr[l][2] = min(len(Gr[l][0]), dk)

#        print('(s1,p3)', self.is_bound('s1','p3','l2'))
#        print('*************************************************************')
        for s in students_in_Gr:
            
            bound_projects = set([p for p in self.G[s][0] if self.is_bound(s, p, self.plc[p][0]) is True])
            unbound_projects = self.G[s][0].difference(bound_projects)
            #print('---------------------------------')
#            print(s, 'bp:', bound_projects, ' up:', unbound_projects)
            # keep track of bound and unbound edges for s
            # update keeps projects that changed from bound to unbound *** do no use .update
            # or better still empty the bound and unbound projects at each iteration
            self.G[s][1] = bound_projects
            self.G[s][2] = unbound_projects
#            print('before: ', s, Gr)
            # for p_j in bound_projects, remove (s_i, p_j) from Gr and reduce quota of p_j in Gr by one
            for bp in bound_projects:
                
                bl = self.plc[bp][0] # lecturer that offers bp
#                print('::::::::>', s, bp, bl, Gr)
                
                Gr[s][0].remove(bp)
                Gr[bp][0].remove(s)
                Gr[bp][2] -= 1
                Gr[bl][3] -= 1
#                print(s, Gr[bp], Gr[bl])
#                print()
                if s in Gr[bl][0]:
                    Gr[bl][0].remove(s)
                    
                #if quota of p_j in Gr is reduced to 0 or p_j becomes an
                # isolated vertex, remove p_j from Gr                    
                if Gr[bp][2] <= 0 and Gr[bp][0] == set():
                    del Gr[bp]
                    Gr[bl][1].remove(bp)
                    projects_in_Gr.remove(bp)     
                    
                #if quota of l_k in Gr is reduced to 0 or l_k becomes an
                # isolated vertex, remove l_k from Gr 
                if Gr[bl][3] <= 0 and Gr[bl][0] == set() and Gr[bl][1] == set():
                    del Gr[bl]
                    lecturers_in_Gr.remove(bl)
                        
            # if bound_projects is non-empty, then we remove all other unbound edges
            # incident to s_i in Gr and we delete the isolated vertex s_i
            # we achieve this by deleting the key s from the dictionary Gr
            if len(bound_projects) > 0:
#                print(s, Gr[s])
                del Gr[s]
                # students_in_Gr.remove(s) --- cannot delete elt from a set during iteration
                # and we can remove the unbound edges from Gr without reducing quotas
                for up in unbound_projects:
                    ul = self.plc[up][0] # lecturer that offers up
                    Gr[up][0].remove(s)
                    
                    #if quota of p_j in Gr is reduced to 0 or p_j becomes an
                    # isolated vertex, remove p_j from Gr                    
                    if Gr[up][2] <= 0 or Gr[up][0] == set():
                        del Gr[up]
                        projects_in_Gr.remove(up)                    
                        Gr[ul][1].remove(up)
                        
                    Gr[ul][0].discard(s)
                    #if quota of l_k in Gr is reduced to 0 or l_k becomes an
                    # isolated vertex, remove l_k from Gr
                    if Gr[ul][3] <= 0 and Gr[ul][0] == set() and Gr[ul][1] == set():
                        del Gr[ul]
                        lecturers_in_Gr.remove(ul)
#        print('***********************************')                      
#        for k,v in Gr.items():
#            print(k, '::>', v)
#            print('***********************************')         
        
        dummy_subgraph = {}
        # ------- extend cloned_Gr with dummy students if needed
        for l in lecturers_in_Gr:
            l_Grquota = Gr[l][3] # Gr[l][3] = q_k^*
            ps_Grquota = sum([Gr[p][2] for p in Gr[l][1]])
            s_adj_l = len(Gr[l][0]) # no of students adjacent to l in Gr
            n1 = s_adj_l - l_Grquota
            n2 = ps_Grquota - l_Grquota
#            n = sum([Gr[p][2] for p in Gr[l][1]]) - l_Grquota # Gr[l][3] = q_k^*
            # if students are more than l's quota in Gr and n > 0, proceed with dummies

            if n1 > 0 and n2 > 0:
                existing_dummies = len(dummy_subgraph)                
                dummy_students = set(['sd'+str(i) for i in range(existing_dummies+1, existing_dummies+n2+1)])
                students_in_Gr.update(dummy_students)
                lower_rank_projects = self.lower_rank_edges(Gr, l)
                Gr[l][0].update(dummy_students)
                for dummy in dummy_students:
                    Gr[dummy] = [set([p for p in lower_rank_projects]), set(), set()]
                    dummy_subgraph[dummy] = [] # just a placeholder for each dummy students in the dummy subgraph
                for project in lower_rank_projects:
                    Gr[project][0].update(dummy_students)
                
        # we clone each p_j q_j^* times and populate the edges to form cloned_Gr
        # we find maximum_matching on cloned_Gr.... 
        # but not before we make sure each dummy students are assigned
        
        cloned_Gr = {}
        for s in students_in_Gr.intersection(Gr):
            s_cloned = set([str(i)+':'+p for p in Gr[s][0] for i in range(1, min(len(Gr[p][0]), Gr[p][2])+1)])
            cloned_Gr[s] = s_cloned
            if s in dummy_subgraph:
                dummy_subgraph[s] = cloned_Gr[s]
                    
#        print(cloned_Gr)
        return dummy_subgraph, cloned_Gr, Gr
    
    
    # =======================================================================
    # find critical set of students
    # =======================================================================s    
    def critical_set(self, graph, matching, exposed_students):
        unvisited_students = deepcopy(exposed_students) # starts with all unassigned students (left set)
        explored_students = set()
        visited_projects = set()        

        #we find all the students in graph reachable from a student in exposed_students via an alternating path
        while unvisited_students:
            s = unvisited_students.pop()
            if s not in explored_students:
                explored_students.add(s)
                #follow all unvisited projects adjacent to s via an unmatched edge                
                unmatched_p_vertices = set([p for p in graph[s] if p not in visited_projects])
                visited_projects.update(unmatched_p_vertices)
                # next we follow all unvisited/unexplored students 
                # adjacent to a project in unmatched_p_vertices via a matched edge
                for p in unmatched_p_vertices:                    
                    if p in matching: 
                        matched_s_vertex = matching[p]
                        if matched_s_vertex not in explored_students.union(unvisited_students):
                            unvisited_students.add(matched_s_vertex)
        return explored_students
    
    # =======================================================================
    # inner repeat until loop -- terminates when self.Z is empty
    # =======================================================================           
    def inner_repeat(self):
        self.Zs, self.Zp = set([None]), set([None])
        while self.Zs.union(self.Zp) != set():
            self.Zs, self.Zp = set(), set()
            self.while_loop() # execute while loop            
            dummy_Gr, cloned_Gr, Gr = self.reduced_assignment_graph() # form reduced assignment graph  
            self.update_project_tail()
            self.update_lecturer_tail() 
            
#            for k,v in self.G.items():
#                print(k, '::>', v)
            dummy_Gr_copy, initial_matching, exposed_dummies, exposed_dprojects = HopcroftKarp(dummy_Gr).maximum_matching() # pre assign dummy students
            cloned_Gr_copy, max_matching, exposed_students, exposed_cloned_projects = HopcroftKarp(cloned_Gr, initial_matching).maximum_matching() # subject to above, augment cloned_Gr to a maximum matching            
            self.Zp = self.critical_set(cloned_Gr_copy, max_matching, exposed_cloned_projects) # find critical set of projects
            
#            print('Gr', Gr)
#            print('cloned Gr copy:', cloned_Gr_copy)
#            print('critical projects', self.Zp)
#            
#            print('max_matching', max_matching)
            #print('********************************')
            #neighbour_Z = set([p for p in Gr[s][0] for s in self.Z]) #** this raises KeyError
            for cp in self.Zp: # where cp is a cloned project
                pu = cp.split(':')[1]
                l_pu = self.plc[pu][0]
                Lkj = self.lp[l_pu][2][pu]
                Lkj_tail_index = self.plc[pu][4]
#                print('worst pointer:', Lkj_tail_index)
                found = False
                for i in range(Lkj_tail_index+1):
                    for s in Lkj[i]:
                        if s in Gr[pu][0]:
                            found = True
                            self.plc[pu][4] = i-1 # what if this becomes negative
                            # well, that means every tie has been deleted, and no student can apply to the project
                            break
                    if found:
                        break
                Lkj_tail_index = self.plc[pu][4]      
                pu_tail = self.lp[l_pu][2][pu][Lkj_tail_index+1:]
#                print(pu_tail)
                # delete every student worse than this index in Lkj
                self.lp[l_pu][2][pu] = self.lp[l_pu][2][pu][:Lkj_tail_index+1]
                for tie in pu_tail:
                    for st in tie:
#                        print(st,pu)
                        # if st is adjacent to pu in Gr and in cloned_Gr_copy, delete the edge
                        if st in Gr:
                            Gr[pu][0].discard(st)
                            Gr[st][0].discard(pu)
                            
                            cloned_Gr_copy[cp].discard(st)
                            cloned_Gr_copy[st].discard(cp)
                        if st in max_matching and max_matching[st] == cp:
                            del max_matching[st]
                            del max_matching[cp]
                        self.delete(st, pu, l_pu)
#                        print('line 20 delete', st,pu,l_pu)
            for student in cloned_Gr_copy:
                if student in self.sp:
                    if cloned_Gr_copy[student] != set() and student not in max_matching:
                        exposed_students.add(student)
            self.Zs = self.critical_set(cloned_Gr_copy, max_matching, exposed_students) # find critical set of students         
#            print('critical students', exposed_students, self.Zs)
            neighbour_Z = set()
            for s in self.Zs:
                for p in Gr[s][0]:
                    neighbour_Z.add(p)
            #print('N(Z)', neighbour_Z)
#            self.update_project_tail() --- not needed because reduced assignment graph will update and no deletions have been afterwards
            for pu in neighbour_Z:                
                l_pu = self.plc[pu][0]
                Lkj = self.lp[l_pu][2][pu]
                Lkj_tail_index = self.plc[pu][4]
#                print('worst pointer:', Lkj_tail_index)
                found = False
                for i in range(Lkj_tail_index+1):
                    for s in Lkj[i]:
                        if s in Gr[pu][0]:
                            found = True
                            self.plc[pu][4] = i-1 # what if this becomes negative
                            # well, that means every tie has been deleted, and no student can apply to the project
                            break
                    if found:
                        break
                Lkj_tail_index = self.plc[pu][4]      
                pu_tail = self.lp[l_pu][2][pu][Lkj_tail_index+1:]
#                print(pu_tail)
                # delete every student worse than this index in Lkj
                self.lp[l_pu][2][pu] = self.lp[l_pu][2][pu][:Lkj_tail_index+1]
                for tie in pu_tail:
                    for st in tie:
#                        print(st,pu)
                        self.delete(st, pu, l_pu)
#                        print('line 20 delete', st,pu,l_pu)
   
#                pu_tail = self.lp[l_pu][2][pu][pu_tail_pointer][:] # Q:what can go wrong for pu_tail_pointer to not point to pu's tail?
#                Q:what can go wrong for pu_tail_pointer to not point to pu's tail?
#                A: if the current tail_index of some project does not contain some student in Gr, because cj > 0
#                as a result, when the stored index is used to delete students at the tail, we don't delete the correct set of students
#                I guess you figured the answer out the hard way! :D :D :D


                    
    # =======================================================================
    # outer repeat until loop -- terminates when every unassigned student has an empty list
    # ======================================================================= 
    def outer_repeat(self):
        while self.unassigned_empty_list:
            self.inner_repeat()       
            self.update_project_tail()
            self.update_lecturer_tail()
            # ** extra deletions (lines 22 - 29 of Algorithm SPA-ST-strong)
            for pj in self.plc:
                deleted_students = self.plc[pj][3]
                cj_remnant = self.G[pj][1]
                if deleted_students and cj_remnant > 0:                    
                    lk = self.plc[pj][0] # lecturer who offers pj
                    sr = self.plc[pj][3][-1] # most preferred student rejected from pj
                    lk_tail_pointer = self.lp[lk][4]
                    Lk_students = self.lp[lk][1][:]
                    # is sr in any of the ties from Lk's tail inwards? 
                    # If yes, set found to True
                    found = False                    
                    while lk_tail_pointer >= 0:
                        if sr in Lk_students[lk_tail_pointer]:
                            found = True
                            break
                        lk_tail_pointer -= 1
                        
                    if found:        
                        Lk_tail = Lk_students[-1] # copy the tail
                        self.lp[lk][1] = self.lp[lk][1][:-1] # delete the tail from Lk's preference list
                        self.lp[lk][4] -= 1 # decrement Lk's worst pointer 
                        p_k = set([i for i in self.lp[lk][2].keys()])  # all the projects that lk is offering
                        for st in Lk_tail:                            
                            a_t = set(self.sp_no_tie_deletions[st])  # the student's altered preference list without ties..
                            common_projects = p_k.intersection(a_t)
                            for pu in common_projects:
                                self.delete(st, pu, lk)
#                                print('line 29 delete', st, pu, lk)
                                    
#            update is useless here because no project or lecturer can be full after a deletion just took place
            self.update_project_tail()
            self.update_lecturer_tail()
                      
                                
    # =======================================================================
    # find set of projects in PR*
    # =======================================================================
    # unassigned, strict preference or indifference on student side
    def PR_star_condition_i(self, student, project, lecturer): # return True
        project_rank = self.sp[student][2][project][0]
        if self.G[student][0] == set():
            return True
        for p in self.G[student][0]:
            #p_lecturer = self.plc[p][0]
            p_rank = self.sp[student][2][p][0]
            # student prefers p_j to p_j'
            if project_rank <= p_rank:
                return True
        return False
    
    def PR_star_condition_ii(self, student, project, lecturer):
        if self.G[lecturer][2] > 0: # lk is undersubscribed
            return True
        if self.G[lecturer][2] <= 0: # lk is full and s_i \in G(l_k)
            if student in self.G[lecturer][0]:
                return True
            # l_k prefers s_i to some student assigned to her in G (i.e., tail)
            # since l_k is full, any tie at the tail must contain a student assigned to l_k in G
            # otherwise, this tie would have been removed !*******
            Lk_all_but_last = self.lp[lecturer][1][:-1]
            for tie in Lk_all_but_last:
                if student in tie:
                    return True
            
        return False

    # =======================================================================
    def find_PR_star(self):
        # PR contains projects that have been deleted from a student's list
        PR = [p for p in self.plc if self.plc[p][3] != []]
        PR_star = set()
        for p in PR:
            l = self.plc[p][0]
            for s in self.plc[p][3]:
#                if (self.PR_star_condition_ia(s,p,l) and self.PR_star_condition_iia(s,p,l)) or (self.PR_star_condition_ib(s,p,l) and self.PR_star_condition_iib(s,p,l)):
                if (self.PR_star_condition_i(s,p,l) and self.PR_star_condition_ii(s,p,l)):
                    PR_star.add(p)
        return PR_star
            

    # =======================================================================
    # form feasible matching self.M in self.G
    # ======================================================================= 
    def feasible_matching(self):
              
        students = set([s for s in self.sp if self.G[s][0] != set()])
        
        # deletions taking place in lines 31 - 36 0f algorithm
        # is s_i adjacent to p_j and p_j' offfered by different lecturers
        # via a bound and unbound edge?
#        for k,v in self.G.items():
#            print(k, '::>', v)
        restart = True
        while restart:
            restart = False
            self.update_bound_unbound()
#            print('***********************')
            for s in students:
                unbound_projects = self.G[s][2]
                bound_projects = list(self.G[s][1])[:]
#                print(s, unbound_projects, bound_projects)
                if bound_projects:
                    pj = bound_projects[0]
                    lk = self.plc[pj][0]
                    non_empty_projects = self.G[lk][1]            
                    # all unbound projects adjacent to s_i that are not offered by l_k
                    for p in unbound_projects.difference(non_empty_projects):
                        restart = True
                        l = self.plc[p][0]
                        self.delete(s,p,l) # line 36 deletion
    #                    print('line 36 delete', s, p, l)
                    
#        again, update is useless here after deletions had just taken place
        self.update_project_tail()
        self.update_lecturer_tail()     
        PR_star = self.find_PR_star()
        G_copy = deepcopy(self.G)
        G_copy = {k:v for k,v in G_copy.items() if v[0] != set()}
        projects = set([p for p in self.plc if p in G_copy])
        # make sure last elt for each project in G_copy is their true quota
        # no need to keep track of quota for lecturer here, cos 
        # critical set will be non-empty if lecturer cannot accommodate all students     
        for p in projects:
            cj = self.plc[p][1]
            G_copy[p].append(min(len(G_copy[p][0]), cj))

        clonedG_star = {}
        clonedG = {}
        for s in students:
            clonedG[s] = set([str(i)+':'+p for p in G_copy[s][0] for i in range(1,G_copy[p][2]+1)])
            if G_copy[s][0].intersection(PR_star) != set():
                clonedG_star[s] = set([str(i)+':'+p for p in G_copy[s][0].intersection(PR_star) for i in range(1,G_copy[p][2]+1)])                    

#        print('G_star:')        
#        for k,v in clonedG_star.items():
#            print('\t',k,v)
        clonedG_star_copy, maxmm2, us2, up2 = HopcroftKarp(clonedG_star).maximum_matching()
        clonedG_copy, self.M, us3, up3 = HopcroftKarp(clonedG, maxmm2).maximum_matching(keys_only=True)
        for k,v in self.M.items():
            project = v.split(':')[1]
            
            self.M[k] = project
            self.M_w_proj_lec[k] = project
            
            if project not in self.M_w_proj_lec:
                self.M_w_proj_lec[project] = [set([k]), 1] # set([k])
            else:
                self.M_w_proj_lec[project][0].add(k)
                self.M_w_proj_lec[project][1] += 1
                
            lecturer = self.plc[project][0]
            if lecturer not in self.M_w_proj_lec:
                self.M_w_proj_lec[lecturer] = [set([k]), 1]
            else:
                self.M_w_proj_lec[lecturer][0].add(k) # a student cannot be added twice, since M is a max matching
                self.M_w_proj_lec[lecturer][1] += 1 # thus, it is safe to increment for each k we encounter
    
    # =======================================================================    
    # blocking pair types
    # =======================================================================    
    def blockingpair_1bi(self, student, project, lecturer):
        #  project and lecturer capacity
        cj, dk = self.plc[project][1], self.lp[lecturer][0]
        # no of students assigned to project in M
        project_occupancy, lecturer_occupancy = 0, 0
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
        # no of students assigned to lecturer in M
        if lecturer in self.M_w_proj_lec: 
            lecturer_occupancy = self.M_w_proj_lec[lecturer][1]
        #  project and lecturer are both under-subscribed
        if project_occupancy < cj and lecturer_occupancy < dk:
            return True
        return False
    
    def blockingpair_1bii(self, student, project, lecturer):
        # p_j is undersubscribed, l_k is full and either s_i \in M(l_k)
        # or l_k prefers s_i to the worst student in M(l_k) or is indifferent between them
        cj, dk = self.plc[project][1], self.lp[lecturer][0]
        project_occupancy, lecturer_occupancy = 0, 0
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
            
        if lecturer in self.M_w_proj_lec: 
            lecturer_occupancy = self.M_w_proj_lec[lecturer][1]
            
        #  project is undersubscribed and lecturer is full
        if project_occupancy < cj and lecturer_occupancy == dk:
            Mlk_students = self.M_w_proj_lec[lecturer][0]
            if student in Mlk_students:
                return True
            remaining_Lk = self.lp[lecturer][1][:]
            for tie in remaining_Lk:
                if student in tie:
                    return True                
        return False
    
    def blockingpair_1biii(self, student, project, lecturer):
        # p_j is full and l_k prefers s_i to the worst student in M(p_j) or is indifferent between them
        cj, project_occupancy = self.plc[project][1], 0
        
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
            
        if project_occupancy == cj:
            remaining_Lkj = self.lp[lecturer][2][project][:]
            for tie in remaining_Lkj:
                if student in tie:
                    return True
        return False
    
    def blockingpair_2bi(self, student, project, lecturer):
        #  project and lecturer are both under-subscribed and s_i \notin M(l_k)
        cj, dk = self.plc[project][1], self.lp[lecturer][0]
        project_occupancy, lecturer_occupancy = 0, 0
        #Mlk_students = []
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
        if lecturer in self.M_w_proj_lec: 
            lecturer_occupancy = self.M_w_proj_lec[lecturer][1]
            #Mlk_students = self.M_w_proj_lec[lecturer][0]
#        if project_occupancy < cj and lecturer_occupancy < dk and student not in Mlk_students:
        if project_occupancy < cj and lecturer_occupancy < dk:
            return True
        return False
        
    
    def blockingpair_2bii(self, student, project, lecturer):
        # p_j is undersubscribed, l_k is full, s_i \notin M(l_k)
        # and l_k prefers s_i to the worst student in M(l_k)
        cj, dk = self.plc[project][1], self.lp[lecturer][0]
        project_occupancy, lecturer_occupancy = 0, 0
        Mlk_students = []
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
        if lecturer in self.M_w_proj_lec: 
            lecturer_occupancy = self.M_w_proj_lec[lecturer][1]
            Mlk_students = self.M_w_proj_lec[lecturer][0]
         
#        if project_occupancy < cj and lecturer_occupancy == dk and student not in Mlk_students:
        if project_occupancy < cj and lecturer_occupancy == dk:
            if student in Mlk_students:
                return True
            remaining_Lk = self.lp[lecturer][1][:-1] # excludes the tail
            for tie in remaining_Lk:
                if student in tie:
                    return True                
        return False
    
    def blockingpair_2biii(self, student, project, lecturer):
        # p_j is full, s_i \notin M(l_k) and l_k prefers s_i to the worst student in M(p_j)
        cj, project_occupancy = self.plc[project][1], 0
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
            
        if project_occupancy == cj:
            remaining_Lkj = self.lp[lecturer][2][project][:-1] #(excludes the tail)
            for tie in remaining_Lkj:
                if student in tie:
                    return True
        return False
    
    # =======================================================================    
    # Is M strongly stable? Check for blocking pair
    # self.blocking_pair is set to True if blocking pair exists
    # =======================================================================
    def check_stability(self):
        self.update_project_tail()
        self.update_lecturer_tail()
        for student in self.sp:
            preferred_projects, indifference = [], []
            if student not in self.M:
                preferred_projects = self.sp_no_tie[student]
            else:
                matched_project = self.M[student]
                rank_matched_project = self.sp_copy[student][2][matched_project][0]
                A_si = self.sp_copy[student][1]
                preferred_projects = [s for tie in A_si[:rank_matched_project] for s in tie]
                indifference = A_si[rank_matched_project][:] # deep copy the tie containing M(s_i)
                indifference.remove(matched_project)
        
            for project in preferred_projects:
                lecturer = self.plc[project][0]
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_1bi(student, project, lecturer)
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_1bii(student, project, lecturer)
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_1biii(student, project, lecturer)
                
                if self.blocking_pair:
#                    print(student, project, lecturer)
                    break
            
            for project in indifference:
                lecturer = self.plc[project][0]
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_2bi(student, project, lecturer)
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_2bii(student, project, lecturer)
                if self.blocking_pair is False:
                    self.blocking_pair = self.blockingpair_2biii(student, project, lecturer)
                
                if self.blocking_pair:                    
                    break
            
            if self.blocking_pair:
#                print(student, project, lecturer)
                
                break
    # =======================================================================
    # Three conditions in M/G to validate that no STSM exist student not in self.M or
    # =======================================================================
    # if a bound edge in G is not in M
    def student_checker(self):
        for student in self.sp:
            bound_projects = self.G[student][1]
            for project in bound_projects:
                if self.M[student] != project:
                    return True
#            unbound_projects = self.G[student][2]
        return False

    
    # replete lecturer is not full in M, or
    # non-replete lecturer has fewer assignee in M than in G
    def lecturer_checker(self):
        for lecturer,value in self.lp.items():
            # replete lecturer is not full in M
            dk, lecturer_occupancy_inM = self.lp[lecturer][0], 0
            if lecturer in self.M_w_proj_lec: 
                lecturer_occupancy_inM = self.M_w_proj_lec[lecturer][1]
                
            full_lk = value[3]
            if full_lk is True:
                if lecturer_occupancy_inM < dk :
#                    print('replete lecturer "{lk}" is not full in M'.format(lk=lecturer))
                    return True
            # non-replete lecturer has fewer assignee in M than in G
            else:
                lecturer_occupancy_inG = len(self.G[lecturer][0])                
                if lecturer_occupancy_inM < lecturer_occupancy_inG:
#                    print('non-replete lecturer "{lk}" has fewer assignees in M than G'.format(lk=lecturer))                    
                    return True
        return False
                        

    def project_checker(self):
        PR_star = self.find_PR_star()
#        print(PR_star)
        for p in PR_star:
            cj = self.plc[p][1]
            project_occupancy = 0
            if p in self.M_w_proj_lec: 
                project_occupancy = self.M_w_proj_lec[p][1]
            if project_occupancy < cj:
#                print(p)
                return True

        return False

        
    def run(self):
        
        self.outer_repeat()
        self.feasible_matching()
        self.check_stability()
    
#        print('student check :::>', self.student_checker())
#        print('project check :::>', self.project_checker())
#        print('lecturer check :::>', self.lecturer_checker())
#        print('blocking pair check :::>', self.blocking_pair)
#        print('===================================')
#        print('Feasible matching: ')
#        print(self.M)
#        for k,v in self.M_w_proj_lec.items():
#            print('\t', k, '::>', v)
#        print('===================================')
#        for k,v in s.G.items():
#            print(k, '::>', v)
#        print('***********************************')
# ###        
#         for k,v in s.sp.items():
#             print(k,v)
#         print('===================================')
#         for k,v in s.lp.items():
#             print(k, '::>', v)
#         print('===================================')
#         for k,v in s.plc.items():
#             print(k, '::>', v)
#         print('===================================')
        if self.student_checker() or self.lecturer_checker() or self.project_checker():
            
            self.found_stsm = 'N'
            
        elif self.blocking_pair:
            self.found_stsm = 'U'
            
        else:
            self.found_stsm ='Y'
        
        return self.found_stsm
"""
# count = 0
# for k in range(1, 10001):
#     filename = "CT/4/instance"+str(k)+".txt"
#     s = SPASTSTRONG(filename)
#     a = s.run()
#     if a == 'Y':
#         count +=1
#         print('instance'+str(k)+'.txt', a)
# print(count)
#filename = "../correctnessTesting/2/instance92144.txt"

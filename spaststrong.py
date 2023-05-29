#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last edited on Sat Feb  1 20:20:49 2020

@author: sofiat
"""


from readFile import SPASTFileReader
#from copy import deepcopy
import networkx as nx

class SPASTSTRONG:
    
    # initialise the instance variables
    def __init__(self, filename):
        self.filename = filename
        # read in the file
        r = SPASTFileReader(self.filename)
        r.read_file()
        # load in the number of students, projects, and lectures
        self.num_students = r.students
        self.num_projects = r.projects
        self.num_lecturers = r.lecturers
        # load in the data structures holding preference information
        self.sp = r.sp                
        self.sp_no_tie_deletions = r.sp_no_tie_deletions
        self.plc = r.plc
        self.lp = r.lp
        # ---------------------------------------------------------------------
        # setup data structures that need to be accessible across the program        
        # ---------------------------------------------------------------------
        self.unassigned_empty_list = list(self.sp.keys()) # keeps track of unassigned students who has an empty list
        self.G = {} #provisional assignment graph
        for student in self.sp.keys():            
            self.G[student] = {"projects": set(), "bound": set(), "unbound": set()} 
        for project in self.plc.keys():
            self.G[project] = {"students": set(), "replete": False, "rejected": [],"num_bound_edges": set(), "revised_quota": 0} 
        for lecturer in self.lp.keys():
            self.G[lecturer] = {"students": set(), "projects": set(), "replete": False, "num_bound_edges": set(), "revised_quota": 0} 
        # ---------------------------------------------------------------------
        # setup additional variables that need to be accessible across the program        
        # ---------------------------------------------------------------------
        self.build_Gr = False
        self.max_flow = {} # maximum matching in Gr to find Zp and Zs
        self.Zp = set()
        self.Zs = set()
        # ---------------------------------------------------------------------
        # self.M = {}
        # self.M_w_proj_lec = {}
        # self.not_assigned = 0
        # self.assigned = 0
        # self.student_check = False
        # self.project_check = False
        # self.lecturer_check = False
        # self.blocking_pair = False
        # self.full_projects = set()
        # self.found_stsm = 'N'

        
    def pquota(self, project): # q_j = min(c_j, d_G(p_j))
        return min(self.plc[project]["cap"], len(self.G[project]["students"]))
        
    
    def lquota(self, lecturer): # q_k = min(d_k, alpha_k)
        alpha_k = sum([self.pquota(project) for project in self.G[lecturer]["projects"]]) # can the if statement here be removed if self.G[lecturer]["projects"] is updated?
        return min(self.lp[lecturer]["cap"], alpha_k)
        
    # =======================================================================
    # add edge (s_i, p_j) to G
    # =======================================================================    
    def add_edge_to_G(self, student, project, lecturer):
        # set filters duplicates, so if a student is added to G(lk) 
        # multiple times, such student appears only once
        self.G[student]["projects"].add(project) 
        self.G[project]["students"].add(student)
        self.G[lecturer]["students"].add(student)        
        self.G[lecturer]["projects"].add(project)
        
    # =======================================================================
    # remove edge (s_i, p_j) from G
    # =======================================================================       
    def remove_edge_from_G(self, student, project, lecturer):
        if project in self.G[student]["projects"]:
            self.G[student]["projects"].remove(project)       
            self.G[project]["students"].remove(student) 
            self.G[project]["rejected"].append(student)  # keep track of students who were rejected from pj                      
            # if the project becomes an isolated vertex
            if len(self.G[project]["students"]) == 0:
                self.G[lecturer]["projects"].remove(project)
            # if in G, student no longer has any project in common w/ lecturer
            if len(self.G[student]["projects"].intersection(self.G[lecturer]["projects"])) == 0:
                self.G[lecturer]["students"].remove(student)        

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
            p_rank = self.sp[student]["list_rank"][project]
            self.sp[student]["list"][p_rank[0]][p_rank[1]] = 'dp'  # we replace the project in this position by a dummy project
            self.sp_no_tie_deletions[student].remove(project)
            
        
        # remove the edge from G
        self.remove_edge_from_G(student, project, lecturer)
        # if student is not adjacent to an edge in G and has a non-empty list
        # add her to the list of unassigned students
        if len(self.G[student]["projects"]) == 0 and student not in self.unassigned_empty_list and len(self.sp_no_tie_deletions[student]) > 0:
            self.unassigned_empty_list.append(student) 
       
    # =======================================================================
    # find dominated students in L_k^j
    # =======================================================================    
    def p_dominated_students(self, project):
        """
        param: p_j
        return: starting index of dominated students in Lkj as well as the students
        """        
        qpj = self.pquota(project)
        Lkj_students = self.plc[project]["list"]  # students who chose p_j according to Lk
        Lkj_tail_index = self.plc[project]["tail_idx"] # self.plc[project][4] is length of list while self.plc[project]["tail_idx"] is tail index after deletions
        Gpj = self.G[project]["students"] # students who are adjacent to pj in G
        dominated_index = None
        dominated_students = []
        count = 0 # this will increment the number of edges adjacent to students who are better than dominated students in Lkj
        for i in range(Lkj_tail_index+1):
            assigneees = Gpj.intersection(Lkj_students[i])
            count += len(assigneees)
            if count >= qpj:
                self.plc[project]["tail_idx"] = i # the updated tail of Lkj is now the tie in position i, again 0-based
                dominated_index = i+1
                dominated_students = self.plc[project]["list"][dominated_index:]
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
        Lk_students = self.lp[lecturer]["list"]  # students in L_k          
        qlk = self.lquota(lecturer)
        projects_inG = self.G[lecturer]["projects"] # non-isolated projects in P_k \cap \mathcal{P}
        p_reduced_quotas = {p:self.pquota(p) for p in projects_inG}
#        print(p_reduced_quotas)
        Lk_tail_index = self.lp[lecturer]["tail_idx"]
        dominated_index = None
        dominated_students = []
        count = 0
        for i in range(Lk_tail_index+1):
            for s in Lk_students[i]:
                for p in self.G[s]["projects"].intersection(projects_inG):
                    if p_reduced_quotas[p] > 0:
                        count += 1
                        p_reduced_quotas[p] -= 1
                    # print(lecturer, s, p, count, i)
                    if count >= qlk:
                        self.lp[lecturer]["tail_idx"] = i
                        dominated_index = i+1
                        dominated_students = Lk_students[dominated_index:]
                        # print('L_k dominated index is: ', dominated_index)
                        # print('L_k dominated student is: ', dominated_students)
                        return dominated_index, dominated_students
   
                
    # =======================================================================
    # get the next set of projects tied together at the head of a student's preference list
    # ======================================================================= 
    def next_tie_student_head(self, student):
        s_preference = self.sp[student]["list"]  # the projected preference list for this student.. this changes during the allocation process.
        # self.sp[student]["head_idx"] points to the tie at the head of s_i's list 
        # if tie pointer is not up to length of pref list --- length is not 0-based!
        head_index = self.sp[student]["head_idx"] # 0-based
        pref_list_length = self.sp[student]["list_len"] # !0-based
        if  head_index < pref_list_length:
            self.sp[student]["head_idx"] += 1  # we increment head_index pointer --> moves inward by 1
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
                    lecturer = self.plc[project]["lec"]
                    # add the edge (student, project) to G 
                    self.add_edge_to_G(student, project, lecturer)
                    # ----------- if the quota of pj in G  is equal to cj ----------
                    if self.pquota(project) == self.plc[project]["cap"]:  
                        self.G[project]["replete"] = True   # set replete(p_j) to True                        
                        dominated_index, dominated_students = self.p_dominated_students(project) # finds dominated students and the starting index on L_k^j
                        # we delete dominated students from L_k^j by replacing L_k_j with non-dominated students
                        #print(student, project, dominated_index, dominated_students)
                        if len(dominated_students) > 0:                            
                            self.plc[project]["list"] = self.plc[project]["list"][:dominated_index] 
                            #print('remaining L_k^j: ', self.plc[project]["list"])
                        #print()
                        # for each dominated student, delete (student, project)                          
                        for tie in dominated_students:
                            for st in tie:
                                self.delete(st, project, lecturer)
                                # print('line11 delete', st, project, lecturer)

                    # ----------- if the quota of lecturer is equal to dk -----------
                    if self.lquota(lecturer) == self.lp[lecturer]["cap"]:
                        # print(lecturer, " replete")
                        self.G[lecturer]["replete"] = True # set replete(l_k) to True 
                        dominated_index, dominated_students = self.l_dominated_students(lecturer) # finds dominated students and the starting index on l_k
                        
                        # print(lecturer, dominated_index, dominated_students, self.G[lecturer])
                        # print(self.lp[lecturer])
                        # print()
#                                for k,v in self.G.items():
#                                    print(k, '::>', v)
                        
                        # we delete dominated students from l_k by replacing l_k with non-dominated students
                        self.lp[lecturer]["list"] = self.lp[lecturer]["list"][:dominated_index]
                        p_k = self.lp[lecturer]["projects"] # all the projects that lecturer is offering
                        for tie in dominated_students:
                            for st in tie:
                                a_t = self.sp_no_tie_deletions[st]  # the student's altered preference list without ties..
                                common_projects = p_k.intersection(a_t)
                                for pu in common_projects:
                                    self.delete(st, pu, lecturer)
                                        # print('lne 16 delete', st,pu,lecturer)
    #                        print('***********************************')
                        # self.update_project_tail() # !important --- do not delete -- may not need afterall (verify!)
    #                        print('***********************************')
    #                        print(lecturer, self.lp[lecturer]["list"])
    #                        for k,v in s.G.items():
    #                            print('\t', k, '::>', v)
    #                        print('***********************************')

            #########################################################################################################################################
            if len(self.G[student]["projects"]) != 0 and student in self.unassigned_empty_list:
                self.unassigned_empty_list.remove(student)
                
            # !* if the current student is unassigned in the matching, with a non-empty preference list, we re-add the student to the unassigned list
            if len(self.G[student]["projects"]) == 0 and student not in self.unassigned_empty_list and len(self.sp_no_tie_deletions[student]) > 0:  
                self.unassigned_empty_list.append(student)
            #########################################################################################################################################
    

    
    
    # =========================================================================
    # is the edge (si, pj) bound? returns True or False
    # =========================================================================
    def is_bound(self, student, project, lecturer): 
        p_tail_idx = self.plc[project]["tail_idx"]
        l_tail_idx = self.lp[lecturer]["tail_idx"]
        # p_preflist = self.plc[project]["list"][:]
        # l_preflist = self.lp[lecturer]["list"][:]     
        qpj = self.pquota(project)
        qlk = self.lquota(lecturer)
        # (number of edges adjacent to p_j in G = qpj) or s_i is not in the tail of L_k_j (or both)
        p_boolean = (len(self.G[project]["students"]) == qpj) or (student not in self.plc[project]["list"][p_tail_idx])
#        if project == 'p3' and student == 's1':
#            print(project, p_preflist, p_boolean, self.G[project])
        # sum_{pj in Pk} (quota of pj) <= dk or (s_i is not in the tail of L_k) (or both)
        l_boolean = (sum([self.pquota(project) for project in self.G[lecturer]["projects"]]) == qlk) or student not in self.lp[lecturer]["list"][l_tail_idx]
        # if project == 'p2' and student == 's4':
        #     print(p_boolean, l_boolean, l_preflist[-1])
        return p_boolean and l_boolean
    
    # =========================================================================
    # reset the number of bound edges adjacent to each pj and lk to 0
    # as well as their revised quota in G from the previous iteration
    # ========================================================================= 
    def reset_bound(self):
        for project in self.plc:
            self.G[project]["num_bound_edges"], self.G[project]["revised_quota"] = 0, 0
        for lecturer in self.lp:
            self.G[lecturer]["num_bound_edges"], self.G[lecturer]["revised_quota"] = 0, 0
    
    # =========================================================================
    # update bound and unbound projects (edges) for each assigned student in self.G
    # ========================================================================= 
    def update_bound_unbound(self):
        self.build_Gr = False
        self.reset_bound()
        for s in self.sp:
            Gs = self.G[s]["projects"]
            bound_projects = set()
            for p in Gs:
                l = self.plc[p]["lec"]
                if self.is_bound(s, p, l):                    
                    bound_projects.add(p)
                    self.G[p]["num_bound_edges"] +=1
                    self.G[l]["num_bound_edges"] +=1
            unbound_projects = Gs.difference(bound_projects)
            self.G[s]["bound"] = bound_projects
            self.G[s]["unbound"] = unbound_projects
            if len(bound_projects) == 0 and len(unbound_projects) > 0:
                self.build_Gr = True
                    
    # =========================================================================
    # revised quota of project and lecturer in Gr
    # ========================================================================= 
    def revised_pquota(self, project):   
        qpj = self.pquota(project) # quota of pj in G
        # code below is dependent on self.update_bound_unbound() being called
        pbound_edges = self.G[project]["num_bound_edges"] # number of bound edges adjacent to pj in G 
        return qpj - pbound_edges
     
    def revised_lquota(self, lecturer):   
        qlk = self.lquota(lecturer) # quota of lk in G
        # code below is dependent on self.update_bound_unbound() being called
        lbound_edges = self.G[lecturer]["num_bound_edges"] # number of bound edges adjacent to projects offered by lk in G 
        return qlk - lbound_edges
    
    def update_revised_quota(self):
        #self.update_bound_unbound()
        for project in self.plc:
            self.G[project]["revised_quota"] = self.revised_pquota(project)
        for lecturer in self.lp:
            self.G[lecturer]["revised_quota"] = self.revised_lquota(lecturer)
            
      
    # =========================================================================
    # form reduced assignment graph Gr
    # ========================================================================= 
    def buildGr(self):
        '''
        This function forms the reduced assignment graph and uses networkx
        to find a maximum flow matching Mr

        Returns
        -------
        max_flow : dict of dict 
            The reduced assignment graph with their capacity in Mr.

        '''
        # the provisional assignment graph G passed into this function will have the
        # updated bound and unbound edges adjacent to each student, project, and lecturer
        # it will also have the information of each project and lecturer
        # whose revised quota is >= 1
        Gr = nx.DiGraph()
        for si in self.sp:
            # if the student is only adjacent to unbound edges
            if len(self.G[si]["bound"]) == 0 and len(self.G[si]["unbound"]) > 0:         
                Gr.add_edge('s', si, capacity=1) # source -> students
                # for each project that si is adjacent to via an unbound edge
                for pj in self.G[si]["unbound"]:
                    Gr.add_edge(si, pj, capacity=1) # students -> projects
        for lk in self.lp:
            # if the lecturer has positive revised quota 
            # then at least one of her projects will also have positive revised quota
            if self.G[lk]["revised_quota"] > 0:
                
                for pj in self.G[lk]["projects"]:
                    if self.G[pj]["revised_quota"] > 0: # if the project has a positive revised quota
                        Gr.add_edge(pj, lk, capacity=self.G[pj]["revised_quota"]) # lecturer (qlk*) -> sink
                Gr.add_edge(lk, 't', capacity = self.G[lk]["revised_quota"]) # lecturer (qlk*) -> sink
        max_flow = nx.max_flow_min_cost(Gr, 's', 't')
        #print(max_flow)
        # we can assert that the total lecturer revised quota is <= number of students
        # this is just for debugging, so that when the assert fails, something is massively wrong!
        return max_flow
    
    # =========================================================================
    # find unhappy projects in Gr relative to the matching Mr found by flow alg
    # =========================================================================
    def unhappy_projects(self):
        Gr_students = set(self.max_flow['s'].keys())        
        typeII_Us = set([])
        Up = set([]) # set of unhappy projects        
        for student in Gr_students:
            if self.max_flow['s'][student] == 1: # if the student is assigned a project in Mr
                #print(student, end=" ")
                # get the projects this student is adjacent to in Gr
                adjacent_projects = set(self.max_flow[student].keys())
                # out of these projects, which of them is the student assigned to by the flow algorithm            
                assigned_project = None
                for project in adjacent_projects:
                    if self.max_flow[student][project] == 1:
                        assigned_project = project
                        #print(assigned_project, end=" ")
                        break # stop the loop as soon as the assigned project is found
                lecturer = self.plc[assigned_project]["lec"] # lecturer who offers this project
                #print(lecturer, end=" ")
                # all the projects that the lecturer offers in G
                lecturer_projects = set(self.G[lecturer]["projects"])
                intersect = adjacent_projects.intersection(lecturer_projects)
                # first we need to remove the project that student is assigned in Mr
                intersect.remove(assigned_project)          
                #print(intersect)            
                # check that the number of assignees of the project in max_flow is
                # smaller than the revised quota of the project in Gr
                for project in intersect:
                    if self.max_flow[project][lecturer] < self.G[project]["revised_quota"]:
                        typeII_Us.add(student) # set filters duplicates
                        Up.add(project)
        print(f"typeII students in Gr {typeII_Us}\n")
        print(f"Unhappy projects in Gr {Up}\n")
        return Up, typeII_Us

    # =========================================================================
    # find critical set of projects in Gr (using unsaturated flow idea)
    # relative to the matching Mr found by flow alg, and unhappy projects
    # ========================================================================= 
    def criticalset_projects(self, Up, typeII_Us):
        ''' Every project that is unhappy will initiate the critical set        
        '''               
        unexplored_projects = set([p for p in Up])
        explored_projects, visited_students = set(), set()
        while unexplored_projects:
            project = unexplored_projects.pop() #pj
            lk = self.plc[project]["lec"] # lk
            explored_projects.add(project)
            unmatched_svertices = set([])
            # for each student who has not being visited before,
            # and who is adjacent to project via an unmatched edde
            for s in self.max_flow["s"]:                
                if s not in visited_students and project in self.max_flow[s] and self.max_flow[s][project] == 0:
                    unmatched_svertices.add(s)
            #print(f"project: {project} ---> unmatched_svertices: {unmatched_svertices}")
            # now follow all matched projects pt from each student in unmatched_svertices
            # and add them to unexplored project if they have not been explored
            # and if they are offered by lk 
            for s in unmatched_svertices:                
                visited_students.add(s)
                for p in self.max_flow[s]:
                    if self.max_flow[s][p] == 1 and self.plc[p]["lec"] == lk and p not in explored_projects:
                        unexplored_projects.add(p)
        print(f"critical set of projects: {explored_projects}")
        return explored_projects

    # =========================================================================
    # carry out the deletions in lines 19 - 21 of algorithm
    # by removing the edges from Mr, Gr, G, and the preference lists
    # =========================================================================   
    def Zp_deletions(self):        
        for project in self.Zp:
            lecturer = self.plc[project]["lec"]
            # all the flow going through project is removed
            # and we deduct this project's flow from that going into the lecturer
            project_flow = self.max_flow[project][lecturer]
            self.max_flow[lecturer]["t"] -= project_flow
            del self.max_flow[project]
            # for each student in Gr that is adjacent to a project in Zp
            # delete this edge in G. Also, if the student is matched to 
            # this project in Mr, then break this matching...            
            Gr_students = list(self.max_flow["s"].keys())
            for student in Gr_students:
                if project in self.max_flow[student]:
                    student_flow = self.max_flow[student][project]
                    self.max_flow["s"][student] -= student_flow
                    del self.max_flow[student][project]
                    # if student becomes an isolated vertex in Gr
                    if self.max_flow[student] == dict():
                        del self.max_flow["s"][student]
            # Now delete pj from the list of each student at the tail of Lkj, as well as from G
            Lkj_students = self.plc[project]["list"]  # reduced students who chose p_j according to Lk
            Lkj_tail_index = self.plc[project]["tail_idx"] # self.plc[project]["tail_idx"] is tail index after deletions
            Lkj_tail = Lkj_students[Lkj_tail_index]                
            for st in Lkj_tail:
                self.delete(st, project, lecturer)
                # print('line22  delete', st, project, lecturer)
                # if the current student is unassigned in G, with a non-empty preference list, 
                # we re-add the student to the unassigned list within the delete method
            # update the project tail as well as Lkj list
            self.plc[project]["list"] = self.plc[project]["list"][:Lkj_tail_index]
            self.plc[project]["tail_idx"] -= 1
      
    # =========================================================================
    # find unhappy students in Gr relative to the matching Mr found by flow alg
    # ========================================================================= 
    def unhappy_students(self):        
        Gr_students = set(self.max_flow['s'].keys())
        Us = set([student for student in Gr_students if self.max_flow['s'][student] == 0])
        print(f"typeI students in Gr {Us}\n")        
        return Us
    
    # =========================================================================
    # find critical set of students in Gr (using unsaturated flow idea)
    # relative to the matching Mr found by flow alg, and unhappy students
    # ========================================================================= 
    def criticalset_students(self, Us):
        ''' Every student who is unhappy will initiate the critical set        
        '''                       
        unexplored_students = set([s for s in Us])
        explored_students, visited_projects, explored_lecturers = set(), set(), set()
        while unexplored_students:
            student = unexplored_students.pop()
            explored_students.add(student)            
            # follow all projects adjacent to student via an unsaturated edge
            projects_to_explore = set([p for p in self.max_flow[student] if self.max_flow[student][p] == 0])
            #print(f"Initial projects to explore from {student}: {projects_to_explore}")
            # from above, we want to follow projects whose flow value is less 
            # than their revised quota in Gr to get to the lecturer
            new_projects = set()
            for p in projects_to_explore:
                lecturer = self.plc[p]["lec"]
                if lecturer not in explored_lecturers and self.max_flow[p][lecturer] < self.G[p]["revised_quota"]:
                    explored_lecturers.add(lecturer)
                    # find all projects that lecturer offers in Gr
                    PknPr = self.G[lecturer]["projects"].intersection(set(self.max_flow.keys()))
                    saturated_projects = set([pj for pj in PknPr if self.max_flow[pj][lecturer] == self.G[pj]["revised_quota"]])
                    new_projects.update(saturated_projects)
                    #print(f"{lecturer}:::::> PKnPr:::> {PknPr} ----- saturated_projects::::> {saturated_projects}")
            projects_to_explore.update(new_projects)
            #print(f"Updated projects to explore: {projects_to_explore}")
            for project in projects_to_explore:
                if project not in visited_projects:
                    visited_projects.add(project)
                    # follow all students adjacent to project via a matched edge
                    # there could be more than one of them since revised quota >= 1
                    all_matched_students = set()
                    # check all students adjacent to pj in G, if the student is in max_flow
                    # and max_flow[student][project] == 1, then student is matched to the project
                    # first get all students adjacent to project in Gr
                    Gr_students = self.G[project]["students"].intersection(set(self.max_flow.keys()))
                    for sr in Gr_students:
                        if self.max_flow[sr][project] == 1:
                            all_matched_students.add(sr)
                    for st in all_matched_students:
                        if st not in explored_students:
                            unexplored_students.add(st)
        print(f"Critical set Zs {explored_students}\n")
        return explored_students
    
    # =========================================================================
    # carry out the deletions in lines 22 - 25 of algorithm
    # by removing the edges from G and the preference lists
    # =========================================================================    
    def Zs_deletions(self):
        # find the projects adjacent in Gr to students in the critical set
        N_Zs = set() 
        for student in self.Zs:
            N_Zs.update(set(self.max_flow[student].keys()))                
        print(f"Neighbourhood of critical set N(Zs) {N_Zs}\n")
        for project in N_Zs:
            lecturer = self.plc[project]["lec"]                        
            Lkj_students = self.plc[project]["list"]  # reduced students who chose p_j according to Lk
            Lkj_tail_index = self.plc[project]["tail_idx"] # self.plc[project]["tail_idx"] is tail index after deletions
            Lkj_tail = Lkj_students[Lkj_tail_index]                
            for st in Lkj_tail:
                self.delete(st, project, lecturer)
                # print('line22  delete', st, project, lecturer)
                
                
            # update the project tail as well as Lkj list
            self.plc[project]["list"] = self.plc[project]["list"][:Lkj_tail_index]
            self.plc[project]["tail_idx"] -= 1
            
        
    # =========================================================================
    # inner repeat-until loop in line 3
    # keeps running until self.Zs union self.Zp is empty
    # ========================================================================= 
    def inner_repeat(self):
        self.Zs, self.Zp = set([None]), set([None]) # storing None allows the while loop to begin
        iteration = 0
        while self.Zs.union(self.Zp):
            self.Zs, self.Zp = set(), set()
            self.while_loop()
            iteration +=1
            print(f"Iteration {iteration} of while loop ...............................")
            self.update_bound_unbound() # classify edges in G as bound or unbound     
            
            # we form the reduced assignment graph if some student is adjacent to only unbound edges
            if self.build_Gr:                      
                self.update_revised_quota() # update revised quota based on # of bound edges
                self.max_flow = self.buildGr() # form reduced assigned graph and find Mr
                # -----------------------------------------------
                print("=========== Provisional assignment graph before Zp ============= \n")
                for s in self.sp:
                    print(f"{s} ---> {self.G[s]} \n --------> {self.sp[s]}")
                print()
                for p in self.plc:
                    print(f"{p} ---> {self.G[p]} \n --------> {self.plc[p]}")
                print()
                for l in self.lp:
                    print(f"{l} ---> {self.G[l]} \n --------> {self.lp[l]}")
                print()   
                # -----------------------------------------------  
                print("=========== Reduced assignment graph before Zp ============= \n")
                for k,v in self.max_flow.items():
                    print(f"{k} ---> {v}")
                # -------------------------------------------------------------
                # find the unhappy projects, critical set of projects, and delete
                Up, typeII_Us = self.unhappy_projects() 
                self.Zp = self.criticalset_projects(Up, typeII_Us)                 
                self.Zp_deletions()                              
                # -----------------------------------------------  
                print("=========== Reduced assignment graph after Zp ============= \n")
                for k,v in self.max_flow.items():
                    print(f"{k} ---> {v}")
                # -------------------------------------------------------------
                # find the unhappy students, critical set of students, and delete                
                Us = self.unhappy_students()
                self.Zs = self.criticalset_students(Us) 
                self.Zs_deletions()
                # -----------------------------------------------
                print("=========== Provisional assignment graph after Zs ============= \n")
                for s in self.sp:
                    print(f"{s} ---> {self.G[s]} \n --------> {self.sp[s]}")
                print()
                for p in self.plc:
                    print(f"{p} ---> {self.G[p]} \n --------> {self.plc[p]}")
                print()
                for l in self.lp:
                    print(f"{l} ---> {self.G[l]} \n --------> {self.lp[l]}")
                print()   
                
            
    

# filename = "caldam.txt"
#filename = "ex7.txt"
filename = "ex8.txt"
I = SPASTSTRONG(filename)
I.inner_repeat()
# I.update_bound_unbound()
# I.update_revised_quota()
# max_flow = I.buildGr()
# max_flow, US = I.unhappy_students(max_flow)
# Zs = I.criticalSet(max_flow, US)
# for s in I.sp:
#     print(f"{s} ;;;> {I.G[s]} ::: {I.sp[s][1]} ::: {I.sp_no_tie_deletions[s]}")
# print()
# for p in I.plc:
#     print(f"{p} ;;;> {I.G[p]}")
# print()
# for l in I.lp:
#     print(f"{l} ;;;> {I.G[l]}")
# print()    




"""
   

    
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
        if self.G[student]["projects"] == set():
            return True
        for p in self.G[student]["projects"]:
            #p_lecturer = self.plc[p]["lec"]
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
            Lk_all_but_last = self.lp[lecturer]["list"][:-1]
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
            l = self.plc[p]["lec"]
            for s in self.plc[p][3]:
#                if (self.PR_star_condition_ia(s,p,l) and self.PR_star_condition_iia(s,p,l)) or (self.PR_star_condition_ib(s,p,l) and self.PR_star_condition_iib(s,p,l)):
                if (self.PR_star_condition_i(s,p,l) and self.PR_star_condition_ii(s,p,l)):
                    PR_star.add(p)
        return PR_star
            

    # =======================================================================
    # form feasible matching self.M in self.G
    # ======================================================================= 
    def feasible_matching(self):
              
        students = set([s for s in self.sp if self.G[s]["projects"] != set()])
        
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
                unbound_projects = self.G[s]["unbound"]
                bound_projects = list(self.G[s]["bound"])[:]
#                print(s, unbound_projects, bound_projects)
                if bound_projects:
                    pj = bound_projects[0]
                    lk = self.plc[pj][0]
                    non_empty_projects = self.G[lk][1]            
                    # all unbound projects adjacent to s_i that are not offered by l_k
                    for p in unbound_projects.difference(non_empty_projects):
                        restart = True
                        l = self.plc[p]["lec"]
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
                
            lecturer = self.plc[project]["lec"]
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
        cj, dk = self.plc[project]["cap"], self.lp[lecturer]["cap"]
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
        cj, dk = self.plc[project]["cap"], self.lp[lecturer]["cap"]
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
            remaining_Lk = self.lp[lecturer]["list"][:]
            for tie in remaining_Lk:
                if student in tie:
                    return True                
        return False
    
    def blockingpair_1biii(self, student, project, lecturer):
        # p_j is full and l_k prefers s_i to the worst student in M(p_j) or is indifferent between them
        cj, project_occupancy = self.plc[project]["cap"], 0
        
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
            
        if project_occupancy == cj:
            remaining_Lkj = self.plc[project]["list"][:]
            for tie in remaining_Lkj:
                if student in tie:
                    return True
        return False
    
    def blockingpair_2bi(self, student, project, lecturer):
        #  project and lecturer are both under-subscribed and s_i \notin M(l_k)
        cj, dk = self.plc[project]["cap"], self.lp[lecturer]["cap"]
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
        cj, dk = self.plc[project]["cap"], self.lp[lecturer]["cap"]
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
            remaining_Lk = self.lp[lecturer]["list"][:-1] # excludes the tail
            for tie in remaining_Lk:
                if student in tie:
                    return True                
        return False
    
    def blockingpair_2biii(self, student, project, lecturer):
        # p_j is full, s_i \notin M(l_k) and l_k prefers s_i to the worst student in M(p_j)
        cj, project_occupancy = self.plc[project]["cap"], 0
        if project in self.M_w_proj_lec: 
            project_occupancy = self.M_w_proj_lec[project][1]
            
        if project_occupancy == cj:
            remaining_Lkj = self.plc[project]["list"][:-1] #(excludes the tail)
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
                lecturer = self.plc[project]["lec"]
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
                lecturer = self.plc[project]["lec"]
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
            dk, lecturer_occupancy_inM = self.lp[lecturer]["cap"], 0
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

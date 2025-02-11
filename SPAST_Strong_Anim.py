from copy import deepcopy
import networkx as nx
import matplotlib.pyplot as plt

from readFile import SPASTFileReader

class SPAST_STRONG_ANIM:
    
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

        self.unassigned_empty_list = list(self.sp.keys())
        self.G = {}
        for student in self.sp.keys():            
            self.G[student] = {"projects": set(), "bound": set(), "unbound": set()} 
        for project in self.plc.keys():
            self.G[project] = {"students": set(), "replete": False, "rejected": [],"num_bound_edges": set(), "revised_quota": 0} 
        for lecturer in self.lp.keys():
            self.G[lecturer] = {"students": set(), "projects": set(), "replete": False, "num_bound_edges": set(), "revised_quota": 0} 

        self.build_Gr = False
        self.max_flow = {}
        self.Zp = set()
        self.Zs = set()

        ### plotting ###
        self.figure, self.axes = plt.subplots(1, 2)
        self.dist = max(self.num_students, self.num_projects, self.num_lecturers)
        self.spacing = {"s": self.dist / (1 + self.num_students),
                   "p": self.dist / (1 + self.num_projects),
                   "l": self.dist / (1 + self.num_lecturers)}
        self.column = {"s": 1, "p": 2, "l": 3}

        ### plot control panel ###
        self.WAIT_PER_DRAW = 3
        self.style_info = {
            "with_labels" : True,
            "node_shape" : "s",
            "node_color" : "none",
            "node_size" : 600,
            "bbox" : {
                "facecolor" : (0.76, 0.69, 0.88),
                "edgecolor" : 'black',
                "boxstyle" : 'round,pad=0.1'
            }
        }

    def pquota(self, project):
        return min(self.plc[project]["cap"], len(self.G[project]["students"]))
    
    def lquota(self, lecturer):
        alpha_k = sum([self.pquota(project) for project in self.G[lecturer]["projects"]])
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

            if len(self.G[student]["projects"].intersection(self.G[lecturer]["projects"])) == 0:
                self.G[lecturer]["students"].remove(student)        

    def delete(self, student, project, lecturer):
        if project in self.sp_no_tie_deletions[student]:
            p_rank = self.sp[student]["list_rank"][project]
            self.sp[student]["list"][p_rank[0]][p_rank[1]] = 'dp'
            self.sp_no_tie_deletions[student].remove(project)

        self.remove_edge_from_G(student, project, lecturer)

        if len(self.G[student]["projects"]) == 0 and student not in self.unassigned_empty_list and len(self.sp_no_tie_deletions[student]) > 0:
            self.unassigned_empty_list.append(student) 
         
    def p_dominated_students(self, project):     
        qpj = self.pquota(project)
        Lkj_students = self.plc[project]["list"]
        Lkj_tail_index = self.plc[project]["tail_idx"]
        Gpj = self.G[project]["students"]
        dominated_index = None
        dominated_students = []
        count = 0
        for i in range(Lkj_tail_index+1):
            assigneees = Gpj.intersection(Lkj_students[i])
            count += len(assigneees)
            if count >= qpj:
                self.plc[project]["tail_idx"] = i
                dominated_index = i+1
                dominated_students = self.plc[project]["list"][dominated_index:]

                return dominated_index, dominated_students

    def l_dominated_students(self, lecturer):
        Lk_students = self.lp[lecturer]["list"]      
        qlk = self.lquota(lecturer)
        projects_inG = self.G[lecturer]["projects"]
        p_reduced_quotas = {p:self.pquota(p) for p in projects_inG}
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
                    if count >= qlk:
                        self.lp[lecturer]["tail_idx"] = i
                        dominated_index = i+1
                        dominated_students = Lk_students[dominated_index:]
                        return dominated_index, dominated_students

    def next_tie_student_head(self, student):
        s_preference = self.sp[student]["list"]
        head_index = self.sp[student]["head_idx"]
        pref_list_length = self.sp[student]["list_len"]
        if  head_index < pref_list_length:
            self.sp[student]["head_idx"] += 1
            return s_preference[head_index]
        return None

    def while_loop(self):
        while self.unassigned_empty_list:            
            student = self.unassigned_empty_list.pop(0)
            tie_at_head = self.next_tie_student_head(student)
            if tie_at_head is not None:
                for project in tie_at_head:
                    if project == 'dp':
                        continue                      
                    lecturer = self.plc[project]["lec"]
                    self.add_edge_to_G(student, project, lecturer)

                    if self.pquota(project) == self.plc[project]["cap"]:  
                        self.G[project]["replete"] = True                      
                        dominated_index, dominated_students = self.p_dominated_students(project)

                        if len(dominated_students) > 0:                            
                            self.plc[project]["list"] = self.plc[project]["list"][:dominated_index] 
                             
                        for tie in dominated_students:
                            for st in tie:
                                self.delete(st, project, lecturer)

                    if self.lquota(lecturer) == self.lp[lecturer]["cap"]:
                        self.G[lecturer]["replete"] = True
                        dominated_index, dominated_students = self.l_dominated_students(lecturer)
                        
                        self.lp[lecturer]["list"] = self.lp[lecturer]["list"][:dominated_index]
                        p_k = self.lp[lecturer]["projects"]
                        for tie in dominated_students:
                            for st in tie:
                                a_t = self.sp_no_tie_deletions[st]
                                common_projects = p_k.intersection(a_t)
                                for pu in common_projects:
                                    self.delete(st, pu, lecturer)

            if len(self.G[student]["projects"]) != 0 and student in self.unassigned_empty_list:
                self.unassigned_empty_list.remove(student)
           
            if len(self.G[student]["projects"]) == 0 and student not in self.unassigned_empty_list and len(self.sp_no_tie_deletions[student]) > 0:  
                self.unassigned_empty_list.append(student)


    def is_bound(self, student, project, lecturer): 
        p_tail_idx = self.plc[project]["tail_idx"]
        l_tail_idx = self.lp[lecturer]["tail_idx"]
  
        qpj = self.pquota(project)
        qlk = self.lquota(lecturer)

        p_boolean = (len(self.G[project]["students"]) == qpj) or (student not in self.plc[project]["list"][p_tail_idx])
        l_boolean = (sum([self.pquota(project) for project in self.G[lecturer]["projects"]]) == qlk) or student not in self.lp[lecturer]["list"][l_tail_idx]

        return p_boolean and l_boolean
    
    def reset_bound(self):
        for project in self.plc:
            self.G[project]["num_bound_edges"], self.G[project]["revised_quota"] = 0, 0
        for lecturer in self.lp:
            self.G[lecturer]["num_bound_edges"], self.G[lecturer]["revised_quota"] = 0, 0
    
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
                    self.G[p]["num_bound_edges"] +=1
                    self.G[L]["num_bound_edges"] +=1
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
                Gr.add_edge('s', si, capacity=1)
                for pj in self.G[si]["unbound"]:
                    Gr.add_edge(si, pj, capacity=1)

        for lk in self.lp:
            if self.G[lk]["revised_quota"] > 0:
                for pj in self.G[lk]["projects"]:
                    if self.G[pj]["revised_quota"] > 0:
                        Gr.add_edge(pj, lk, capacity=self.G[pj]["revised_quota"])
                Gr.add_edge(lk, 't', capacity = self.G[lk]["revised_quota"])

        max_flow = nx.max_flow_min_cost(Gr, 's', 't')
        return max_flow, Gr

    def unhappy_projects(self):
        Gr_students = set(self.max_flow['s'].keys())        
        typeII_Us = set([])
        Up = set([])    
        for student in Gr_students:
            if self.max_flow['s'][student] == 1:
                
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
                    if self.max_flow[project][lecturer] < self.G[project]["revised_quota"]:
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
                if s not in visited_students and project in self.max_flow[s] and self.max_flow[s][project] == 0:
                    unmatched_svertices.add(s)

            for s in unmatched_svertices:                
                visited_students.add(s)
                for p in self.max_flow[s]:
                    if self.max_flow[s][p] == 1 and self.plc[p]["lec"] == lk and p not in explored_projects:
                        unexplored_projects.add(p)

        return explored_projects

    def Zp_deletions(self):        
        for project in self.Zp:
            lecturer = self.plc[project]["lec"]

            project_flow = self.max_flow[project][lecturer]
            self.max_flow[lecturer]["t"] -= project_flow
            del self.max_flow[project]
           
            Gr_students = list(self.max_flow["s"].keys())
            for student in Gr_students:
                if project in self.max_flow[student]:
                    student_flow = self.max_flow[student][project]
                    self.max_flow["s"][student] -= student_flow
                    del self.max_flow[student][project]

                    if self.max_flow[student] == dict():
                        del self.max_flow["s"][student]

            Lkj_students = self.plc[project]["list"]
            Lkj_tail_index = self.plc[project]["tail_idx"]
            Lkj_tail = Lkj_students[Lkj_tail_index]                
            for st in Lkj_tail:
                self.delete(st, project, lecturer)

            self.plc[project]["list"] = self.plc[project]["list"][:Lkj_tail_index]
            self.plc[project]["tail_idx"] -= 1
      
    def unhappy_students(self):        
        Gr_students = set(self.max_flow['s'].keys())
        Us = set([student for student in Gr_students if self.max_flow['s'][student] == 0])
        return Us

    def criticalset_students(self, Us):               
        unexplored_students = set([s for s in Us])
        explored_students, visited_projects, explored_lecturers = set(), set(), set()
        while unexplored_students:
            student = unexplored_students.pop()
            explored_students.add(student)            
            projects_to_explore = set([p for p in self.max_flow[student] if self.max_flow[student][p] == 0])
            new_projects = set()

            for p in projects_to_explore:
                lecturer = self.plc[p]["lec"]
                if lecturer not in explored_lecturers and self.max_flow[p][lecturer] < self.G[p]["revised_quota"]:
                    explored_lecturers.add(lecturer)
                    PknPr = self.G[lecturer]["projects"].intersection(set(self.max_flow.keys()))
                    saturated_projects = set([pj for pj in PknPr if self.max_flow[pj][lecturer] == self.G[pj]["revised_quota"]])
                    new_projects.update(saturated_projects)

            projects_to_explore.update(new_projects)

            for project in projects_to_explore:
                if project not in visited_projects:
                    visited_projects.add(project)
                    all_matched_students = set()
                    Gr_students = self.G[project]["students"].intersection(set(self.max_flow.keys()))
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
        self.Zs, self.Zp = set([None]), set([None])
        while self.Zs.union(self.Zp):
            self.Zs, self.Zp = set(), set()
            self.while_loop()
            self.update_bound_unbound()

            self.figure.suptitle("Post-domination assignments.")
            self.draw_SPA_graph()
            
            if self.build_Gr:                      
                self.update_revised_quota()
                self.max_flow, G_reduced = self.buildGr()

                ### project ###
                Up, typeII_Us = self.unhappy_projects() 
                self.Zp = self.criticalset_projects(Up)
                self.Zp_deletions()

                self.figure.suptitle(f"Critical Projects Stage, Z_p = {self.Zp}.")
                self.draw_SPA_reduced(G_reduced)
                self.draw_SPA_graph()

                ### student ###
                Us = self.unhappy_students()
                self.Zs = self.criticalset_students(Us) 
                self.Zs_deletions()

                self.figure.suptitle(f"Critical Projects Stage, Z_s = {self.Zs}.")
                self.draw_SPA_reduced(G_reduced)
                self.draw_SPA_graph()

                self.axes[1].clear()
    
    def draw_graph_plot(self, G, labels, pos, ax):

        self.axes[0].set_title("G, the provisional assignment graph.")
        self.axes[1].set_title("G_r, the reduced assignment graph.")

        nx.draw(G,
                pos,
                with_labels = self.style_info["with_labels"],
                labels = labels,
                node_shape = self.style_info["node_shape"],
                node_color = self.style_info["node_color"],
                node_size = self.style_info["node_size"],
                bbox = self.style_info["bbox"],
                ax = ax)
        #nx.draw_networkx_edges(Gr, pos, ax=0, edgelist=M_r, edge_color='r')
        #nx.draw_networkx_edges(Gr, pos, ax=0, edgelist=G_r.edges-M_r)

        plt.show(block=False)
        plt.pause(self.WAIT_PER_DRAW)

    def draw_SPA_graph(self):
        G = nx.DiGraph()
        G.add_nodes_from(self.G.keys())

        for k in self.G.keys():
            if k[0] == 's':
                for pj in self.G[k]["projects"]:
                    G.add_edge(k, pj)
            elif k[0] == 'l':
                for pj in self.G[k]["projects"]:
                    G.add_edge(pj, k)

        pos = dict()
        labels = dict()

        for x in (self.G.keys()):
            letter = x[0]
            number = int(x[1:])
            pos[x] = (self.column[letter], number * self.spacing[letter])

            if letter == 'l':
                labels[x] = f"{x} : {self.lp[x]["cap"]}"
            elif letter == 'p':
                labels[x] = f"{x} : {self.plc[x]["cap"]}"
            else:
                labels[x] = x

        self.axes[0].clear()
        self.draw_graph_plot(G, labels, pos, self.axes[0])

    def draw_SPA_reduced(self, Gr):
        G_display = deepcopy(Gr)
        G_display.remove_node('s')
        G_display.remove_node('t')

        pos = dict()
        labels = dict()

        for x in (G_display.nodes()):
            letter = x[0]
            number = int(x[1:])
            pos[x] = (self.column[letter], number * self.spacing[letter])

            if letter in ('p','l'):
                labels[x] = f"{x} : {self.G[x]["revised_quota"]}"
            else:
                labels[x] = x

        self.axes[1].clear()
        self.draw_graph_plot(G_display, labels, pos, self.axes[1])

filename = "ex1.txt"
instance = SPAST_STRONG_ANIM(filename)
instance.inner_repeat()
print("Finished")
plt.pause(15)
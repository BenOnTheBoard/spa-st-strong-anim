import matplotlib.pyplot as plt
import networkx as nx

from bruteforce import STSMBruteForce
from spaststrong import SPAST_STRONG


class SPAST_STRONG_ANIM(SPAST_STRONG):
    def __init__(self, filename):
        super().__init__(filename=filename)

        self.ss_edges = []

        ### plotting ###
        self.figure, self.axes = plt.subplots(1, 3)
        self.dist = max(self.num_students, self.num_projects, self.num_lecturers)
        self.spacing = {
            "s": self.dist / (1 + self.num_students),
            "p": self.dist / (1 + self.num_projects),
            "l": self.dist / (1 + self.num_lecturers),
        }
        self.column = {"s": 1, "p": 2, "l": 3}

        ### plot control panel ###
        self.WAIT_PER_DRAW = 6
        self.END_WAIT = 400
        self.style_info = {
            "stable_edge_colour": "#2efd33",
            "with_labels": True,
            "node_shape": "s",
            "node_color": "none",
            "node_size": 600,
            "bbox": {
                "facecolor": "#a187cc",
                "edgecolor": "black",
                "boxstyle": "round,pad=0.1",
            },
        }

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
                self.max_flow, Gr_weights = self.buildGr()
                self.draw_matching_only()

                ### student ###
                Us = self.unhappy_students()
                self.Zs = self.criticalset_students(Us, Gr_weights)
                self.figure.suptitle(f"Critical Students Stage, Z_s = {self.Zs}.")
                self.draw_SPA_reduced()
                self.Zs_deletions()
                self.draw_SPA_graph()

                self.axes[1].clear()

    def draw_graph_plot(self, G, labels, pos, ax):
        self.axes[0].set_title("G, the provisional assignment graph.")
        self.axes[1].set_title("G_r, the reduced assignment graph.")
        self.axes[2].set_title("Maximum/feasible matching.")

        edge_color_list = []
        for e in G.edges:
            if e in self.ss_edges:
                edge_color_list.append(self.style_info["stable_edge_colour"])
            else:
                edge_color_list.append("#000000")

        nx.draw(
            G,
            pos,
            edge_color=edge_color_list,
            with_labels=self.style_info["with_labels"],
            labels=labels,
            node_shape=self.style_info["node_shape"],
            node_color=self.style_info["node_color"],
            node_size=self.style_info["node_size"],
            bbox=self.style_info["bbox"],
            ax=ax,
        )
        # nx.draw_networkx_edges(Gr, pos, ax=0, edgelist=M_r, edge_color='r')
        # nx.draw_networkx_edges(Gr, pos, ax=0, edgelist=G_r.edges-M_r)

        plt.show(block=False)
        plt.pause(self.WAIT_PER_DRAW)

    def draw_SPA_graph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(self.G.keys())

        for k in self.G.keys():
            if k[0] == "s":
                for pj in self.G[k]["projects"]:
                    graph.add_edge(k, pj)
            elif k[0] == "l":
                for pj in self.G[k]["projects"]:
                    graph.add_edge(pj, k)

        pos = dict()
        labels = dict()

        for x in self.G.keys():
            letter = x[0]
            number = int(x[1:])
            pos[x] = (self.column[letter], number * self.spacing[letter])

            if letter == "l":
                labels[x] = f"{x} : {self.lp[x]['cap']}"
            elif letter == "p":
                labels[x] = f"{x} : {self.plc[x]['cap']}"
            else:
                labels[x] = x

        self.axes[0].clear()
        self.draw_graph_plot(graph, labels, pos, self.axes[0])

    def max_flow_to_reduced_graph(self):
        flow_G = nx.DiGraph()

        for student in self.max_flow["s"].keys():
            for project in self.max_flow[student].keys():
                flow_G.add_edge(student, project)
                for lecturer in self.max_flow[project].keys():
                    flow_G.add_edge(project, lecturer)

        return flow_G

    def draw_SPA_reduced(self):
        G_display = self.max_flow_to_reduced_graph()

        pos = dict()
        labels = dict()

        for x in G_display.nodes():
            letter = x[0]
            number = int(x[1:])
            pos[x] = (self.column[letter], number * self.spacing[letter])

            labels[x] = x

        self.axes[1].clear()
        self.draw_graph_plot(G_display, labels, pos, self.axes[1])

    def max_flow_as_graph(self):
        flow_G = nx.DiGraph()

        if not self.max_flow:
            return flow_G

        for student in self.max_flow["s"].keys():
            for project, flow_val in self.max_flow[student].items():
                if flow_val == 1:
                    flow_G.add_edge(student, project)

        return flow_G

    def draw_matching_only(self):
        G_matching = self.max_flow_as_graph()

        pos = dict()
        labels = dict()

        for x in G_matching.nodes():
            letter = x[0]
            number = int(x[1:])
            pos[x] = (self.column[letter], number * self.spacing[letter])

            labels[x] = x

        self.axes[2].clear()
        self.draw_graph_plot(G_matching, labels, pos, self.axes[2])

    def find_ss_edges_dict(self):
        bruteforcer = STSMBruteForce(filename=filename)
        bruteforcer.choose()
        ssm_list = bruteforcer.get_ssm_list()

        self.ss_edges.clear()
        for matching in ssm_list:
            for k, v in matching.items():
                print(f"{k} : {v}")
                self.ss_edges.append((k, v))
            print("---" * 4)

    def run(self):
        self.find_ss_edges_dict()

        while not self.can_exit_outermost_loop():
            self.inner_repeat()  # lines 2 - 26
            self.repletion_deletions()  # lines 27 - 34

            self.figure.suptitle("Post-repletion deletions.")
            self.draw_SPA_graph()

        self.delete_lesser_unbound_edges()
        self.figure.suptitle("Post-unbound deletions.")
        self.draw_SPA_graph()

        double_bound = False
        for si in self.sp:
            if len(self.G[si]["bound"]) > 1:
                double_bound = True

        if double_bound:
            self.M = {si: "" for si in self.sp}
            self.figure.suptitle("End-state : No matching")
        else:
            self.get_feasible_matching()
            self.figure.suptitle("End-state : Feasible")
        self.draw_matching_only()
        plt.pause(self.END_WAIT)

        return self.M


# filename = "examples/misc/Zs_deletes_ssp.txt"
# filename = "examples/misc/too_few_del.txt"
# filename = "examples/small_breakers/no_del_2bii.txt"
filename = "alg_error.txt"
instance = SPAST_STRONG_ANIM(filename)
instance.run()
print("Finished")

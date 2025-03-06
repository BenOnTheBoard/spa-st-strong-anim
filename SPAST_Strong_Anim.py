import matplotlib.pyplot as plt
import networkx as nx

from spaststrong import SPAST_STRONG


class SPAST_STRONG_ANIM(SPAST_STRONG):
    def __init__(self, filename):
        super().__init__(filename)

        ### plotting ###
        self.figure, self.axes = plt.subplots(1, 2)
        self.dist = max(self.num_students, self.num_projects, self.num_lecturers)
        self.spacing = {
            "s": self.dist / (1 + self.num_students),
            "p": self.dist / (1 + self.num_projects),
            "l": self.dist / (1 + self.num_lecturers),
        }
        self.column = {"s": 1, "p": 2, "l": 3}

        ### plot control panel ###
        self.WAIT_PER_DRAW = 5
        self.style_info = {
            "with_labels": True,
            "node_shape": "s",
            "node_color": "none",
            "node_size": 600,
            "bbox": {
                "facecolor": (0.76, 0.69, 0.88),
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
                self.max_flow = self.buildGr()

                ### student ###
                Us = self.unhappy_students()
                self.Zs = self.criticalset_students(Us)
                self.figure.suptitle(f"Critical Students Stage, Z_s = {self.Zs}.")
                self.draw_SPA_reduced()
                self.Zs_deletions()
                self.draw_SPA_graph()

                if not self.Zs:
                    ### project ###
                    Up, typeII_Us = self.unhappy_projects()
                    self.Zp = self.criticalset_projects(Up)
                    self.figure.suptitle(f"Critical Projects Stage, Z_p = {self.Zp}.")
                    self.draw_SPA_reduced()
                    self.Zp_deletions()
                    self.draw_SPA_graph()

                self.axes[1].clear()

        self.figure.suptitle("End-state reached.")
        self.draw_SPA_graph()

    def draw_graph_plot(self, G, labels, pos, ax):
        self.axes[0].set_title("G, the provisional assignment graph.")
        self.axes[1].set_title("G_r, the reduced assignment graph.")

        nx.draw(
            G,
            pos,
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
        G = nx.DiGraph()
        G.add_nodes_from(self.G.keys())

        for k in self.G.keys():
            if k[0] == "s":
                for pj in self.G[k]["projects"]:
                    G.add_edge(k, pj)
            elif k[0] == "l":
                for pj in self.G[k]["projects"]:
                    G.add_edge(pj, k)

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
        self.draw_graph_plot(G, labels, pos, self.axes[0])

    def max_flow_as_graph(self):
        flow_G = nx.DiGraph()

        for student in self.max_flow["s"].keys():
            for project in self.max_flow[student].keys():
                flow_G.add_edge(student, project)
                for lecturer in self.max_flow[project].keys():
                    flow_G.add_edge(project, lecturer)

        return flow_G

    def draw_SPA_reduced(self):
        G_display = self.max_flow_as_graph()

        pos = dict()
        labels = dict()

        for x in G_display.nodes():
            letter = x[0]
            number = int(x[1:])
            pos[x] = (self.column[letter], number * self.spacing[letter])

            if letter in ("p", "l"):
                labels[x] = f"{x} : {self.G[x]['revised_quota']}"
            else:
                labels[x] = x

        self.axes[1].clear()
        self.draw_graph_plot(G_display, labels, pos, self.axes[1])


filename = "examples/K instances/K33.txt"
instance = SPAST_STRONG_ANIM(filename)
instance.inner_repeat()
print("Finished")
plt.pause(15)

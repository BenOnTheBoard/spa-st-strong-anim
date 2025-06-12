import matplotlib.pyplot as plt
import networkx as nx
from itertools import product

from readFile import SPASTFileReader


class ZPMaxMatchings:
    def __init__(self, filename):
        """
        Takes a file representing a reduced assignment graph and
        displays all the possible maximum matchings for that graph.
        """

        self.filename = filename
        r = SPASTFileReader(self.filename)
        r.read_file()

        self.num_students = r.students
        self.num_projects = r.projects
        self.num_lecturers = r.lecturers

        self.sp = r.sp
        self.plc = r.plc
        self.lp = r.lp

        for si, si_info in self.sp.items():
            if si_info["list_len"] > 1:
                return ValueError(f"{si} has more than one tie")

        for lk, lk_info in self.lp.items():
            if lk_info["list_len"] > 1:
                return ValueError(f"{lk} has more than one tie")

        self.construct_graph()
        self.set_style_info()

    def set_style_info(self):
        self.style_info = {
            "matched_edge_colour": "#00ff04",
            "with_labels": True,
            "node_shape": "s",
            "node_color": "none",
            "node_size": 600,
            "bbox": {
                "facecolor": "#a887dd",
                "edgecolor": "black",
                "boxstyle": "round,pad=0.1",
            },
        }

    def construct_graph(self):
        self.graph = nx.DiGraph()
        self.positions = dict()
        self.labels = dict()

        self.dist = max(self.num_students, self.num_projects, self.num_lecturers)
        self.spacing = {
            "s": self.dist / (1 + self.num_students),
            "p": self.dist / (1 + self.num_projects),
            "l": self.dist / (1 + self.num_lecturers),
        }
        self.column = {"s": 1, "p": 2, "l": 3}

        self.graph.add_nodes_from(self.sp.keys())
        self.graph.add_nodes_from(self.plc.keys())
        self.graph.add_nodes_from(self.lp.keys())

        for si, si_info in self.sp.items():
            for pj in si_info["list"][0]:
                self.graph.add_edge(si, pj)

        for lk, lk_info in self.lp.items():
            for pj in lk_info["projects"]:
                self.graph.add_edge(pj, lk)

        for x in self.graph.nodes:
            letter = x[0]
            number = int(x[1:])
            self.positions[x] = (self.column[letter], number * self.spacing[letter])

            if letter == "l":
                self.labels[x] = f"{x} : {self.lp[x]['cap']}"
            elif letter == "p":
                self.labels[x] = f"{x} : {self.plc[x]['cap']}"
            else:
                self.labels[x] = x

    def draw_graph_plot(self, match_edges, ax):
        projects_with_matches = set([pair[1] for pair in match_edges])
        edge_color_list = []
        for e in self.graph.edges:
            if e in match_edges:
                edge_color_list.append(self.style_info["matched_edge_colour"])
            elif e[0] in projects_with_matches:
                edge_color_list.append(self.style_info["matched_edge_colour"])
            else:
                edge_color_list.append("#000000")

        nx.draw(
            self.graph,
            self.positions,
            edge_color=edge_color_list,
            with_labels=self.style_info["with_labels"],
            labels=self.labels,
            node_shape=self.style_info["node_shape"],
            node_color=self.style_info["node_color"],
            node_size=self.style_info["node_size"],
            bbox=self.style_info["bbox"],
            ax=ax,
        )

    def brtueforce_max_matchings(self):
        """
        This proceeds by picking one project from every students' options which is
        allowed because every student must be saturated by the Z_p phase.
        """
        max_matchings = []
        si_lists = [si_info["list"][0] for si_info in self.sp.values()]

        for assignment in product(*si_lists):
            result = list(zip(self.sp.keys(), assignment))

            valid = True
            for pj in assignment:
                if assignment.count(pj) > self.plc[pj]["cap"]:
                    valid = False

            for lk, lk_info in self.lp.items():
                total_count = 0
                for pt in lk_info["projects"]:
                    total_count += assignment.count(pt)
                if total_count > lk_info["cap"]:
                    valid = False

            if valid:
                max_matchings.append(result)

        return max_matchings

    def run(self, row_length):
        all_max_matchings = self.brtueforce_max_matchings()
        for elt in all_max_matchings:
            print(elt)

        matching_rows = [
            all_max_matchings[i : i + 3] for i in range(0, len(all_max_matchings), 3)
        ]

        for row_contents in matching_rows:
            _, axes = plt.subplots(1, row_length)

            for i, matching in enumerate(row_contents):
                self.draw_graph_plot(matching, axes[i])

            plt.show()


if __name__ == "__main__":
    grmm = ZPMaxMatchings("examples/misc/5530simple.txt")
    grmm.run(row_length=3)

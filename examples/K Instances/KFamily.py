class KFamily:
    def __init__(self, n):
        self.n = n
        self.filename = f"examples/problematic/K{n}{n}.txt"

    def make_instance_file(self):
        header = [f"{self.n} {self.n * (self.n - 1)} {self.n}"]

        student_edges = {s + 1: [] for s in range(self.n)}
        project_lines = []
        for i in range(self.n):
            for j in range(self.n - 1):
                proj_num = i * (self.n - 1) + j + 1
                project_lines.append(f"{proj_num} 1 {i + 1}")
                # edges that create Z_p
                student_edges[i + 1].append(proj_num)

        # diagonal edges, create stable matching
        all_ports = [v.copy() for v in student_edges.values()]
        for i, lec_ports in enumerate(all_ports):
            lec_ports.insert(i, None)

        for student, prefs in student_edges.items():
            for lec_ports in all_ports:
                port_num = lec_ports[student - 1]
                if port_num is not None:
                    prefs.append(port_num)
            prefs.sort()
            student_edges[student] = list(map(str, prefs))

        student_lines = [f"{k} ({':'.join(v)})" for k, v in student_edges.items()]

        lecturer_lines = []
        preferences = f"({':'.join([str(s + 1) for s in range(self.n)])})"
        for i in range(self.n):
            lecturer_lines.append(f"{i + 1} 1 {preferences}")

        # combine and write
        all_lines = header + student_lines + project_lines + lecturer_lines
        contents = "\n".join(all_lines)

        with open(self.filename, "w") as inst_file:
            inst_file.writelines(contents)


k_gen = KFamily(6)
k_gen.make_instance_file()

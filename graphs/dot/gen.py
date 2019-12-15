from string import ascii_lowercase


def dense_static_wt():
    with open("g3_dense_w1.dot", "w") as f:
        f.write("digraph {\n")
        for c in ascii_lowercase:
            f.write(f"{c};\n")
        for c in ascii_lowercase:
            for d in ascii_lowercase:
                f.write(f"{c} -> {d} [weight=1];\n")

        f.write("}\n")


def dense_nosl_static_wt():
    with open("g3_dense_w1_nosl.dot", "w") as f:
        f.write("digraph {\n")
        for c in ascii_lowercase:
            f.write(f"{c};\n")
        for c in ascii_lowercase:
            for d in ascii_lowercase:
                if c != d:
                    f.write(f"{c} -> {d} [weight=1];\n")

        f.write("}\n")


def dense_dyn_wt():
    with open("g3_dense_wd.dot", "w") as f:
        f.write("digraph {\n")
        for c in ascii_lowercase:
            f.write(f"{c};\n")
        wt = 1
        for c in ascii_lowercase:
            for d in ascii_lowercase:
                f.write(f"{c} -> {d} [weight={wt}];\n")
                wt += 1

        f.write("}\n")


def dense_nosl_dyn_wt():
    with open("g3_dense_wd_nosl.dot", "w") as f:
        f.write("digraph {\n")
        for c in ascii_lowercase:
            f.write(f"{c};\n")
        wt = 1
        for c in ascii_lowercase:
            for d in ascii_lowercase:
                if c != d:
                    f.write(f"{c} -> {d} [weight={wt}];\n")
                    wt += 1

        f.write("}\n")


def g4():
    with open("g4.dot", "w") as f:
        f.write("digraph {\n")
        for c in ascii_lowercase:
            f.write(f"{c};\n")
        wt = 1
        for c in ascii_lowercase:
            for d in ascii_lowercase:
                if c != d:
                    f.write(f"{c} -> {d} [weight={wt}];\n")
            wt += 1

        f.write("}\n")


def g5():
    with open("g5.dot", "w") as f:
        f.write("digraph {\n")
        nodes=1024
        for i in range(nodes):
            f.write(f"node_{i};\n")
        wt = 1
        for i in range(nodes):
            for j in range(nodes):
                if i != j:
                    f.write(f"node_{i} -> node_{j} [weight={wt}];\n")
            wt += 1

        f.write("}\n")


def g6():
    with open("g6.dot", "w") as f:
        f.write("digraph {\n")
        nodes=256
        for i in range(nodes):
            f.write(f"{i};\n")
        for i in range(nodes):
            for j in range(nodes):
                if i != j:
                    f.write(f"{i} -> {j} [weight=1];\n")

        f.write("}\n")
g6()

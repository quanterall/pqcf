"""FORM scripts for the comparison with pq_symbolic: T1 power, T2 product, T3 permanent, T4 mesh unitary."""
import pathlib, sys
out = pathlib.Path(sys.argv[1])

def t1(n):
    xs = "+".join(f"x{i}" for i in range(1, 7))
    return f"Symbols x1,...,x6;\nLocal F = (1+{xs})^{n};\n.sort\nDrop;\n.end\n" if False else \
           f"Symbols x1,x2,x3,x4,x5,x6;\nLocal F = (1+{xs})^{n};\n.end\n"

def t2():
    return "Symbols x,y,z,w;\nLocal P = (1+x+y+z+w)^10;\n.sort\nLocal F = P*P;\n.sort\nDrop P;\n.end\n"

def t3(n):
    syms = ",".join(f"a{i}{j}" for i in range(1, n + 1) for j in range(1, n + 1))
    xs = ",".join(f"x{j}" for j in range(1, n + 1))
    s = f"Symbols {syms},{xs};\nLocal P = 1;\n.sort\n"
    for i in range(1, n + 1):
        lin = "+".join(f"a{i}{j}*x{j}" for j in range(1, n + 1))
        s += f"Multiply ({lin});\n" + "".join(f"id x{j}^2 = 0;\n" for j in range(1, n + 1)) + ".sort\n"
    s += "id " + "*".join(f"x{j}" for j in range(1, n + 1)) + " = 1;\n.end\n"
    return s

def mesh(k):
    """the Clements mesh's beamsplitters (m, m+1) in order, as in pq_nif's decomposition: k(k-1)/2 of them"""
    bs = []
    for layer in range(k):
        for m in range(layer % 2, k - 1, 2):
            bs.append(m)
    return bs[: k * (k - 1) // 2]

def t4(k):
    bs = mesh(k)
    n = len(bs)
    syms = ",".join([f"t{j}" for j in range(n)] + [f"p{j}" for j in range(n)])
    s = f"Symbols {syms};\n"
    # U starts as the identity; U{s}r{c}
    for r in range(k):
        for c in range(k):
            s += f"Local U0r{r}c{c} = {1 if r == c else 0};\n"
    s += ".sort\n"
    for j, m in enumerate(bs):
        a, b = m, m + 1
        cos = f"(t{j}+t{j}^-1)/2"
        sin = f"(-i_/2)*(t{j}-t{j}^-1)"
        gab = f"(-p{j}^-1)*{sin}"
        gba = f"p{j}*{sin}"
        for c in range(k):
            for r in range(k):
                if r == a:
                    e = f"{cos}*U{j}r{a}c{c} + {gab}*U{j}r{b}c{c}"
                elif r == b:
                    e = f"{gba}*U{j}r{a}c{c} + {cos}*U{j}r{b}c{c}"
                else:
                    e = f"U{j}r{r}c{c}"
                s += f"Local U{j+1}r{r}c{c} = {e};\n"
        s += ".sort\n" + "Drop " + ",".join(f"U{j}r{r}c{c}" for r in range(k) for c in range(k)) + ";\n.sort\n"
    s += ".end\n"
    return s

for n in (8, 12, 16, 20):
    (out / f"t1_{n}.frm").write_text(t1(n))
(out / "t2.frm").write_text(t2())
for n in (5, 6, 7, 8):
    (out / f"t3_{n}.frm").write_text(t3(n))
for k in (4, 6, 8):
    (out / f"t4_{k}.frm").write_text(t4(k))
print("written")

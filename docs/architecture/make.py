# The architecture diagram of the README, in GitHub's light and dark themes: three layers, L2 the applications, L1
# PQCF, L0 the FVLLMP federation of lab nodes.
#   python3 docs/architecture/make.py          writes docs/architecture/architecture-{light,dark}.svg
import os

W, H = 1060, 734
SANS = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
MONO = "ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, 'DejaVu Sans Mono', monospace"

THEMES = {
    "light": dict(bg="#ffffff", text="#1f2328", muted="#59636e", band="#f6f8fa", band_line="#d1d9e0",
                  tag="#1f2328", tag_text="#ffffff",
                  app="#ddf4ff", app_line="#54aeff", exe="#fff8c5", exe_line="#d4a72c",
                  sci="#dafbe1", sci_line="#4ac26b", node="#ffffff", node_line="#8250df",
                  inner="#fbefff", link="#8250df", arrow="#59636e"),
    "dark": dict(bg="#0d1117", text="#e6edf3", muted="#9198a1", band="#151b23", band_line="#3d444d",
                 tag="#e6edf3", tag_text="#0d1117",
                 app="#0c2d6b", app_line="#4493f8", exe="#3b2e0a", exe_line="#d29922",
                 sci="#0f3a1f", sci_line="#3fb950", node="#0d1117", node_line="#ab7df8",
                 inner="#271052", link="#ab7df8", arrow="#9198a1"),
}

X0 = 112                     # the content starts right of the layer tags
CW = W - X0 - 20             # content width


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=13, weight=400, fill="text", font=SANS, anchor="start"):
    return (f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" font-weight="{weight}" '
            f'fill="var(--{fill})" text-anchor="{anchor}">{esc(s)}</text>')


def rect(x, y, w, h, fill, stroke, r=10, dash=None, sw=1.5):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="var(--{fill})" '
            f'stroke="var(--{stroke})" stroke-width="{sw}"{d}/>')


def arrow(x1, y1, x2, y2):
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="var(--arrow)" stroke-width="1.5" '
            f'marker-end="url(#head)"/>')


def band(y, h, tag, name, sub):
    """A layer: a full-width band with its tag (L0, L1, L2) on the left."""
    return [rect(10, y, W - 20, h, "band", "band_line", r=14),
            rect(32, y + 14, 56, 30, "tag", "tag", r=15),
            text(60, y + 35, tag, 15, 700, fill="tag_text", anchor="middle"),
            text(60, y + 64, name, 12, 600, fill="text", anchor="middle"),
            *[text(60, y + 81 + 15 * i, s, 11, fill="muted", anchor="middle") for i, s in enumerate(sub)]]


def app(x, y, w, title, lines):
    out = [rect(x, y, w, 78, "app", "app_line"), text(x + w / 2, y + 26, title, 15, 600, anchor="middle")]
    for i, line in enumerate(lines):
        out.append(text(x + w / 2, y + 47 + 17 * i, line, 12, fill="muted", anchor="middle"))
    return out


def column(x, y, w, h, fill, stroke, title, subtitle, items):
    out = [rect(x, y, w, h, fill, stroke), text(x + 18, y + 28, title, 15, 600),
           text(x + 18, y + 46, subtitle, 12, fill="muted")]
    for i, (mods, what) in enumerate(items):
        yy = y + 74 + 33 * i
        out.append(text(x + 18, yy, mods, 12, 500, font=MONO))
        out.append(text(x + 18, yy + 16, what, 12, fill="muted"))
    return out


def svg(theme):
    t = THEMES[theme]
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="Three layers: L2 the applications (GenHub2, photonic simulation, scientific computing); L1 '
         f'PQCF, with distributed execution and the scientific engine; L0 the FVLLMP federation of NVIDIA DGX '
         f'Spark lab nodes, joined by mutual TLS with no central point">',
         '<defs><marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
         'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--arrow)"/></marker></defs>',
         f'<rect width="{W}" height="{H}" fill="var(--bg)"/>']

    # L2: the applications
    y2 = 10
    o += band(y2, 112, "L2", "Applications", ["built on", "pq_task"])
    aw = (CW - 2 * 20) / 3
    apps = [("GenHub2", ["virus classification for genetic labs", "every analysis a pq_task"]),
            ("Photonic simulation", ["Gaussian and Fock backends,", "boson and Gaussian boson sampling"]),
            ("Scientific computing", ["Feynman integrals, the QCD β-function,", "holographic QCD, quantum optics"])]
    for i, (title, lines) in enumerate(apps):
        o += app(X0 + i * (aw + 20), y2 + 17, aw, title, lines)

    # L1: PQCF, two capabilities
    y1, h1 = 146, 300
    o += band(y1, h1, "L1", "PQCF", ["Erlang/OTP", "and Rust"])
    colw = (CW - 30) / 2
    cy, ch = y1 + 16, h1 - 32
    o += column(X0, cy, colw, ch, "exe", "exe_line", "Distributed execution",
                "processes, cores and nodes: the same result on any number",
                [("pq_task", "the task behaviour: an application's work as plain data"),
                 ("pq_rexec, pq_rstore", "executors (this node, the owner's nodes, peers) and their store"),
                 ("pq_job, pq_parallel, pq_worker", "sampling jobs: chunks, retries, cancellation, timeouts"),
                 ("pq_node_registry, pq_spec_cache", "discovery, a lost node's work retried elsewhere"),
                 ("describe/1, map_traced/3", "a provenance record for every result"),
                 ("pq_service, python/", "a JSON service and a Python client")])
    sx = X0 + colw + 30
    o += column(sx, cy, colw, ch, "sci", "sci_line", "Scientific engine",
                "exact and numerical, checked against established tools",
                [("pq_gaussian, pq_fock, pq_sampling", "photonic states, permanents, hafnians, torontonians"),
                 ("pq_symbolic", "exact expressions, compiled programs, factor and gcd"),
                 ("pq_ff, pq_recon, pq_ibp", "reconstruction over finite fields; integral reduction"),
                 ("pq_ode, pq_jet, pq_series, pq_mp", "ODEs, series and limits, arbitrary precision"),
                 ("pq_geometry, pq_tensor", "curvature and tensors: holographic models"),
                 ("pq_qop, pq_lindblad, pq_qo_exact", "quantum optics: master equations, exact steady states")])
    my = cy + ch / 2
    o.append(arrow(sx - 2, my, X0 + colw + 2, my))
    o.append(text(X0 + colw + 15, my - 8, "tasks", 11, fill="muted", anchor="middle"))

    # L2 onto L1
    ab = y2 + 17 + 78
    c0, c1 = X0 + aw / 2, X0 + aw + 20 + aw / 2
    c2 = X0 + 2 * (aw + 20) + aw / 2
    o.append(arrow(c0, ab, c0, cy - 2))
    o.append(arrow(c1 - 40, ab, X0 + colw - 60, cy - 2))
    o.append(arrow(c1 + 40, ab, sx + 60, cy - 2))
    o.append(arrow(c2, ab, c2, cy - 2))

    # L0: the FVLLMP federation
    y0, h0 = y1 + h1 + 34, 248
    o += band(y0, h0, "L0", "FVLLMP", ["federation", "of lab nodes"])
    o.append(text(X0, y0 + 26, "Each lab's own NVIDIA DGX Spark nodes; mutual-TLS gates between them; "
                  "no central point", 13, 600))
    nw, nh, gap = 250, 132, (CW - 3 * 250) / 2
    ny = y0 + 64
    labs = ["lab A", "lab B", "lab C"]
    centres = []
    for i, lab in enumerate(labs):
        x = X0 + i * (nw + gap)
        centres.append((x, x + nw))
        o.append(rect(x, ny, nw, nh, "node", "node_line", sw=2))
        o.append(text(x + 16, ny + 26, f"NVIDIA DGX Spark · {lab}", 14, 600))
        o.append(rect(x + 14, ny + 40, nw - 28, 34, "inner", "node_line", r=6, sw=1))
        o.append(text(x + 26, ny + 62, "PQCF node (Erlang/OTP)", 12, 500))
        o.append(text(x + 16, ny + 96, "Rust kernels on CPU and GPU", 12, fill="muted"))
        o.append(text(x + 16, ny + 115, "the lab's data stays here", 12, fill="muted"))
    # mutual TLS between every pair: neighbours straight, A-C arched under the middle node
    ly = ny + nh / 2
    for (a, b) in [(0, 1), (1, 2)]:
        x1, x2 = centres[a][1], centres[b][0]
        o.append(f'<line x1="{x1}" y1="{ly}" x2="{x2}" y2="{ly}" stroke="var(--link)" stroke-width="2"/>')
        o.append(text((x1 + x2) / 2, ly - 8, "mTLS", 11, 600, fill="link", anchor="middle"))
    ax1, ax2 = (centres[0][0] + centres[0][1]) / 2, (centres[2][0] + centres[2][1]) / 2
    bot = ny + nh + 2
    o.append(f'<path d="M{ax1},{bot} C{ax1},{bot + 36} {ax2},{bot + 36} {ax2},{bot}" fill="none" '
             f'stroke="var(--link)" stroke-width="2"/>')
    o.append(text(X0 + CW / 2, bot + 17, "mTLS", 11, 600, fill="link", anchor="middle"))

    # L1 runs on L0
    o.append(arrow(X0 + CW / 2, y1 + h1 + 1, X0 + CW / 2, y0 - 1))
    o.append(text(X0 + CW / 2 + 10, y1 + h1 + 21, "PQCF runs on every node", 11, fill="muted"))

    o.append("</svg>")
    out = "\n".join(o) + "\n"
    for k in sorted(t, key=len, reverse=True):     # the colours inline: not every SVG renderer knows var()
        out = out.replace(f"var(--{k})", t[k])
    return out


here = os.path.dirname(os.path.abspath(__file__))
for theme in THEMES:
    with open(os.path.join(here, f"architecture-{theme}.svg"), "w") as f:
        f.write(svg(theme))

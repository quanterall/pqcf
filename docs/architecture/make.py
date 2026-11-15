# The architecture diagram of the README, in GitHub's light and dark themes.
#   python3 docs/architecture/make.py          writes docs/architecture/architecture-{light,dark}.svg
import os

W, H = 960, 618
SANS = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
MONO = "ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, 'DejaVu Sans Mono', monospace"

THEMES = {
    "light": dict(bg="#ffffff", text="#1f2328", muted="#59636e", line="#d1d9e0", box="#f6f8fa",
                  app="#ddf4ff", app_line="#54aeff", exe="#fff8c5", exe_line="#d4a72c",
                  sci="#dafbe1", sci_line="#4ac26b", node="#f6f8fa", node_line="#818b98",
                  core="#ffffff", core_line="#59636e", arrow="#59636e"),
    "dark": dict(bg="#0d1117", text="#e6edf3", muted="#9198a1", line="#3d444d", box="#151b23",
                 app="#0c2d6b", app_line="#4493f8", exe="#3b2e0a", exe_line="#d29922",
                 sci="#0f3a1f", sci_line="#3fb950", node="#151b23", node_line="#656c76",
                 core="#0d1117", core_line="#9198a1", arrow="#9198a1"),
}


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=13, weight=400, fill="text", font=SANS, anchor="start"):
    return (f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" font-weight="{weight}" '
            f'fill="var(--{fill})" text-anchor="{anchor}">{esc(s)}</text>')


def rect(x, y, w, h, fill, stroke, r=10, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="var(--{fill})" '
            f'stroke="var(--{stroke})" stroke-width="1.5"{d}/>')


def arrow(x1, y1, x2, y2):
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="var(--arrow)" stroke-width="1.5" '
            f'marker-end="url(#head)"/>')


def app(x, y, w, title, lines):
    out = [rect(x, y, w, 78, "app", "app_line"), text(x + w / 2, y + 26, title, 15, 600, anchor="middle")]
    for i, line in enumerate(lines):
        out.append(text(x + w / 2, y + 47 + 17 * i, line, 12, fill="muted", anchor="middle"))
    return out


def column(x, y, w, h, fill, stroke, title, subtitle, items):
    out = [rect(x, y, w, h, fill, stroke), text(x + 20, y + 30, title, 15, 600),
           text(x + 20, y + 49, subtitle, 12, fill="muted")]
    for i, (mods, what) in enumerate(items):
        yy = y + 80 + 34 * i
        out.append(text(x + 20, yy, mods, 12, 500, font=MONO))
        out.append(text(x + 20, yy + 16, what, 12, fill="muted"))
    return out


def svg(theme):
    t = THEMES[theme]
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="PQCF\'s architecture: applications on PQCF\'s two capabilities, distributed execution and '
         f'the scientific engine, on Erlang/OTP nodes with Rust kernels">',
         '<defs><marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
         'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--arrow)"/></marker></defs>',
         f'<rect width="{W}" height="{H}" fill="var(--bg)"/>']

    # applications
    o.append(text(24, 32, "APPLICATIONS", 11, 600, fill="muted"))
    o += app(24, 44, 288, "GenHub2", ["virus classification for genetic labs", "every analysis a pq_task"])
    o += app(336, 44, 288, "Photonic simulation", ["Gaussian and Fock backends,", "boson and Gaussian boson sampling"])
    o += app(648, 44, 288, "Scientific computing", ["Feynman integrals, the QCD β-function,", "holographic QCD, quantum optics"])

    # PQCF: two capabilities
    cy, ch = 180, 300
    o.append(rect(12, cy - 28, W - 24, ch + 46, "core", "core_line", r=14, dash="6 5"))
    o.append(text(W - 28, cy - 9, "PQCF", 13, 700, fill="muted", anchor="end"))
    o += column(24, cy + 8, 440, ch - 16, "exe", "exe_line", "Distributed execution",
                "processes, cores and nodes: the same result on any number",
                [("pq_task", "the task behaviour: an application's work as plain data"),
                 ("pq_rexec, pq_rstore", "executors (this node, the owner's nodes, federation) and their store"),
                 ("pq_job, pq_parallel, pq_worker", "sampling jobs: chunks, retries, cancellation, timeouts"),
                 ("pq_node_registry, pq_spec_cache", "discovery, a lost node's work retried elsewhere"),
                 ("describe/1", "a provenance record for every result"),
                 ("pq_service, python/", "a JSON service and a Python client")])
    o += column(496, cy + 8, 440, ch - 16, "sci", "sci_line", "Scientific engine",
                "exact and numerical, checked against established tools",
                [("pq_gaussian, pq_fock, pq_sampling", "photonic states, permanents, hafnians, torontonians"),
                 ("pq_symbolic", "exact expressions, compiled programs, factor and gcd"),
                 ("pq_ff, pq_recon, pq_ibp", "reconstruction over finite fields; integral reduction"),
                 ("pq_ode, pq_jet, pq_series, pq_mp", "ODEs, series and limits, arbitrary precision"),
                 ("pq_geometry, pq_tensor", "curvature and tensors: holographic models"),
                 ("pq_qop, pq_lindblad, pq_qo_exact", "quantum optics: master equations, exact steady states")])
    # the engine's work runs as tasks
    o.append(f'<line x1="496" y1="{cy + 160}" x2="466" y2="{cy + 160}" stroke="var(--arrow)" stroke-width="1.5" '
             f'marker-end="url(#head)"/>')
    o.append(text(480, cy + 150, "tasks", 11, fill="muted", anchor="middle"))

    # applications onto the capabilities
    o.append(arrow(168, 122, 168, cy + 6))
    o.append(arrow(430, 122, 330, cy + 6))
    o.append(arrow(530, 122, 630, cy + 6))
    o.append(arrow(792, 122, 792, cy + 6))

    # the runtime
    ny = cy + ch + 56
    o.append(text(24, ny - 8, "RUNTIME", 11, 600, fill="muted"))
    for i, (name, what) in enumerate([("node 1", "the coordinating node: tasks' data"),
                                      ("node 2", "workers on every core"),
                                      ("node n", "another machine, or a federation peer")]):
        x = 24 + i * 312
        o.append(rect(x, ny, 288, 74, "node", "node_line"))
        o.append(text(x + 16, ny + 26, f"Erlang/OTP {name}", 14, 600))
        o.append(text(x + 16, ny + 46, what, 12, fill="muted"))
        o.append(text(x + 16, ny + 63, "Rust kernels (pq_nif) on dirty schedulers", 11, fill="muted", font=MONO))
    o.append(arrow(W / 2, cy + ch + 18, W / 2, ny - 2))
    o.append("</svg>")
    out = "\n".join(o) + "\n"
    for k in sorted(t, key=len, reverse=True):     # the colours inline: not every SVG renderer knows var()
        out = out.replace(f"var(--{k})", t[k])
    return out


here = os.path.dirname(os.path.abspath(__file__))
for theme in THEMES:
    with open(os.path.join(here, f"architecture-{theme}.svg"), "w") as f:
        f.write(svg(theme))

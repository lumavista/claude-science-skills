#!/usr/bin/env python3
"""Standard figures for the report, drawn from the other scripts' outputs.

Every figure is driven by the files the pipeline already writes -- nothing
here knows about a particular genome, construct, sample or coordinate:

  coverage_<sample>_<construct>.svg
      Per-base depth over one vector contig (samtools depth -a), with the
      membership blocks of identity.json shaded underneath (construct-unique
      versus shared, labelled by which constructs carry each block), the
      regions the copy-number trim excludes, the single-genome baseline, and
      each block's mean retained depth in copies per haploid genome.
  genome_overview.svg
      Every host chromosome from a .fai as a length bar (any build, any
      organism -- no cytobands), and each called locus of called_sites.tsv as a
      tick, coloured by attributed construct and styled by confidence.
  copy_number.svg
      Copies per haploid genome per construct per sample, one bar per
      copy_number.py JSON supplied (e.g. two trim policies side by side), with
      the per-unique-block spread and the shared-block closure deficit.
  architecture.svg
      Per construct, the site-count bound from reconcile_sites.py (floor,
      attributed, ceiling, unresolved candidates) and the copies per site it
      implies.

A figure whose inputs are missing is skipped with a note; the run never fails
for that. `figures.json` in --out lists what was drawn and what was skipped.

`--fill REPORT.md` replaces every `[[FIG:<id>]]` placeholder of a report built
from assets/report_template.md with a Markdown image link to the figure (or a
visible "not generated" note), so the report carries the generated figures
rather than hand-drawn ones. `[[FIG:coverage_<construct>]]` expands to that
construct's track for every sample.

    python scripts/make_figures.py --out figures \\
        --identity identity.json --config config.json \\
        --copy-number accounting.json --copy-number legacy=accounting_ce.json \\
        --called called_sites.tsv --site-summary site_summary.json \\
        --fai genome.fa.fai --fill report.md --fill-out report_final.md

SVG only, Python 3.8+, standard library only.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

# --------------------------------------------------------------------------
# palette: categorical slots in fixed order, light and dark steps of the same
# hues. Colour follows the construct (identity.json order), never its rank.
# --------------------------------------------------------------------------

SERIES_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                "#e87ba4", "#008300", "#6250d6", "#e34948"]
SERIES_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500",
               "#d55181", "#0ca30c", "#9085e9", "#e66767"]
ROLES = {
    "surface": ("#fcfcfb", "#1a1a19"),
    "text": ("#0b0b0b", "#f0efec"),
    "text2": ("#52514e", "#c3c2b7"),
    "grid": ("#e4e3df", "#33332f"),
    "axis": ("#a3a29c", "#6b6a64"),
    "neutral": ("#8a8984", "#8f8e88"),
    "shared": ("#d9d8d3", "#3d3d39"),
    "trim": ("#52514e", "#c3c2b7"),
}
for _i, (_l, _d) in enumerate(zip(SERIES_LIGHT, SERIES_DARK)):
    ROLES["s%d" % _i] = (_l, _d)

FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


def fmt(x: Optional[float], nd: int = 1) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return "n/a"
    if abs(x) >= 1000:
        return "{:,.0f}".format(x)
    s = ("%." + str(nd) + "f") % x
    return s


def fmt_bp(n: float) -> str:
    n = float(n)
    for unit, div in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if n >= div:
            v = n / div
            return ("%.0f %s" if v >= 100 else "%.1f %s") % (v, unit)
    return "%d bp" % n


def nice_ticks(lo: float, hi: float, n: int = 5) -> List[float]:
    if hi <= lo:
        hi = lo + 1
    raw = (hi - lo) / max(1, n)
    mag = 10 ** math.floor(math.log10(raw))
    step = mag
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            step = m * mag
            break
    start = math.ceil(lo / step) * step
    out, v = [], start
    while v <= hi + step * 1e-9:
        out.append(round(v, 10))
        v += step
    if out and out[-1] < hi - step * 1e-9:
        out.append(round(out[-1] + step, 10))
    return out


def text_w(s: str, size: float = 11) -> float:
    """Rough rendered width of a string, for fitting labels."""
    return len(s) * size * 0.56


# --------------------------------------------------------------------------
# a minimal SVG builder
# --------------------------------------------------------------------------

class Svg:
    def __init__(self, w: int, h: int, fig_id: str, title: str, desc: str = ""):
        self.w, self.h, self.id = w, h, re.sub(r"[^A-Za-z0-9_-]", "_", fig_id)
        self.title, self.desc = title, desc
        self.body: List[str] = []
        self.defs: List[str] = []

    # paint by role: CSS class (so dark mode can swap it) plus the light value
    # as a presentation attribute, for renderers that ignore <style>
    @staticmethod
    def paint(fill: Optional[str] = None, stroke: Optional[str] = None,
              extra_class: Optional[str] = None) -> str:
        cls, attrs = ([extra_class] if extra_class else []), []
        if fill:
            if fill in ROLES:
                cls.append("f-" + fill)
                attrs.append('fill="%s"' % ROLES[fill][0])
            else:
                attrs.append('fill="%s"' % fill)
        if stroke:
            if stroke in ROLES:
                cls.append("k-" + stroke)
                attrs.append('stroke="%s"' % ROLES[stroke][0])
            else:
                attrs.append('stroke="%s"' % stroke)
        if cls:
            attrs.insert(0, 'class="%s"' % " ".join(cls))
        return " ".join(attrs)

    @staticmethod
    def _tip(title: Optional[str]) -> str:
        return "<title>%s</title>" % escape(title) if title else ""

    @staticmethod
    def _extra(kw: Dict) -> str:
        kw.pop("class", None)
        out = []
        for k, v in kw.items():
            if v is None:
                continue
            out.append('%s="%s"' % (k.replace("_", "-"), escape(str(v), {'"': "&quot;"})))
        return " ".join(out)

    def add(self, s: str) -> None:
        self.body.append(s)

    def rect(self, x, y, w, h, fill=None, stroke=None, title=None, **kw):
        self.add('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" %s %s>%s</rect>'
                 % (x, y, max(0.0, w), max(0.0, h), self.paint(fill, stroke, kw.get("class")),
                    self._extra(kw), self._tip(title)))

    def line(self, x1, y1, x2, y2, stroke="axis", title=None, **kw):
        kw.setdefault("stroke_width", 1)
        self.add('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" %s %s>%s</line>'
                 % (x1, y1, x2, y2, self.paint(None, stroke, kw.get("class")), self._extra(kw),
                    self._tip(title)))

    def poly(self, pts: Sequence[Tuple[float, float]], closed=False, fill=None,
             stroke=None, title=None, **kw):
        tag = "polygon" if closed else "polyline"
        if not closed:
            kw.setdefault("fill", "none")
        p = " ".join("%.1f,%.1f" % xy for xy in pts)
        self.add('<%s points="%s" %s %s>%s</%s>'
                 % (tag, p, self.paint(fill, stroke, kw.get("class")), self._extra(kw),
                    self._tip(title), tag))

    def circle(self, cx, cy, r, fill=None, stroke=None, title=None, **kw):
        self.add('<circle cx="%.2f" cy="%.2f" r="%.2f" %s %s>%s</circle>'
                 % (cx, cy, r, self.paint(fill, stroke, kw.get("class")), self._extra(kw),
                    self._tip(title)))

    def text(self, x, y, s, size=11, fill="text", anchor="start", weight=None,
             title=None, **kw):
        if weight:
            kw["font_weight"] = weight
        self.add('<text x="%.2f" y="%.2f" font-size="%s" text-anchor="%s" %s %s>%s%s</text>'
                 % (x, y, size, anchor, self.paint(fill, None, kw.get("class")), self._extra(kw),
                    escape(s), self._tip(title)))

    def open(self, **kw):
        cls = kw.get("class")
        self.add("<g%s %s>" % (' class="%s"' % escape(cls) if cls else "",
                               self._extra(kw)))

    def close(self):
        self.add("</g>")

    def hatch(self, name: str, role: str, angle: int = 45, gap: int = 6,
              width: float = 1.2) -> str:
        pid = "%s-%s" % (self.id, name)
        self.defs.append(
            '<pattern id="%s" patternUnits="userSpaceOnUse" width="%d" height="%d" '
            'patternTransform="rotate(%d)"><line x1="0" y1="0" x2="0" y2="%d" '
            'stroke-width="%.1f" %s/></pattern>'
            % (pid, gap, gap, angle, gap, width, self.paint(None, role)))
        return "url(#%s)" % pid

    def render(self) -> str:
        css = []
        for role, (lt, dk) in ROLES.items():
            css.append(".f-%s{fill:%s}.k-%s{stroke:%s}" % (role, lt, role, lt))
        dark = []
        for role, (lt, dk) in ROLES.items():
            dark.append(".f-%s{fill:%s}.k-%s{stroke:%s}" % (role, dk, role, dk))
        style = ("text{font-family:%s}%s"
                 "@media (prefers-color-scheme: dark){%s}" % (FONT, "".join(css),
                                                                "".join(dark)))
        head = ('<?xml version="1.0" encoding="UTF-8"?>\n'
                '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
                'viewBox="0 0 %d %d" role="img" aria-labelledby="%s-t %s-d" '
                'data-figure="%s">\n<title id="%s-t">%s</title>\n'
                '<desc id="%s-d">%s</desc>\n<style>%s</style>\n'
                % (self.w, self.h, self.w, self.h, self.id, self.id, self.id,
                   self.id, escape(self.title), self.id, escape(self.desc), style))
        defs = "<defs>%s</defs>\n" % "".join(self.defs) if self.defs else ""
        bg = '<rect x="0" y="0" width="%d" height="%d" %s/>\n' % (
            self.w, self.h, self.paint("surface"))
        return head + defs + bg + "\n".join(self.body) + "\n</svg>\n"

    def save(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.render())


def construct_roles(order: Sequence[str]) -> Dict[str, str]:
    """Colour role per construct, in identity order. Past eight constructs the
    rest fold to neutral rather than cycling hues (a cycled hue would claim two
    constructs are one)."""
    return {c: ("s%d" % i if i < len(SERIES_LIGHT) else "neutral")
            for i, c in enumerate(order)}


def legend_row(svg: Svg, x: float, y: float, items: List[Tuple[str, str]],
               max_x: float) -> float:
    """Draw swatch+label items left to right, wrapping. Each item is
    (label, svg-fragment drawing a 12x12 swatch at 0,0). Returns the y after."""
    cx = x
    for label, swatch in items:
        w = 18 + text_w(label, 11) + 16
        if cx + w > max_x and cx > x:
            cx, y = x, y + 18
        svg.add('<g transform="translate(%.1f,%.1f)">%s</g>' % (cx, y - 10, swatch))
        svg.text(cx + 18, y, label, size=11, fill="text2")
        cx += w
    return y + 18


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------

def load_json(path: str) -> Optional[Dict]:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def read_depth(path: str, want: Sequence[str]) -> Tuple[Optional[Dict[int, float]], str]:
    """Per-base depth of one contig from a `samtools depth` file.

    `want` lists acceptable contig names in order of preference. A file holding
    a single contig is accepted whatever it is called."""
    by: Dict[str, Dict[int, float]] = {}
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                by.setdefault(f[0], {})[int(f[1])] = float(f[2])
            except ValueError:
                continue
    for w in want:
        if w in by:
            return by[w], w
    if len(by) == 1:
        name = next(iter(by))
        return by[name], name
    return None, ("none of %s in %s (contigs: %s)"
                  % (", ".join(want), path, ", ".join(sorted(by)[:6]) or "none"))


def read_fai(path: str) -> List[Tuple[str, int]]:
    out = []
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                try:
                    out.append((f[0], int(f[1])))
                except ValueError:
                    continue
    return out


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def trim_regions(ident_c: Dict, trim: int, policy: str) -> List[Tuple[int, int]]:
    """The bases copy_number.py excludes before averaging (same rule)."""
    length = int(ident_c["length"])
    if trim <= 0:
        return []
    if policy == "contig-end":
        return [(1, min(trim, length)), (max(1, length - trim + 1), length)]
    blocks = ident_c["blocks"]
    circ = bool(ident_c.get("circular", False))
    out = []
    n = len(blocks)
    for i, b in enumerate(blocks):
        prev = blocks[i - 1] if i > 0 else (blocks[-1] if circ and n > 1 else None)
        nxt = blocks[i + 1] if i < n - 1 else (blocks[0] if circ and n > 1 else None)
        a, z = int(b["start"]), int(b["end"])
        if prev is None or prev["membership"] != b["membership"]:
            out.append((a, min(z, a + trim - 1)))
        if nxt is None or nxt["membership"] != b["membership"]:
            out.append((max(a, z - trim + 1), z))
    out.sort()
    merged: List[Tuple[int, int]] = []
    for a, z in out:
        if merged and a <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], z))
        else:
            merged.append((a, z))
    return merged


def inside(p: int, regions: Sequence[Tuple[int, int]]) -> bool:
    return any(a <= p <= z for a, z in regions)


# --------------------------------------------------------------------------
# figure 1: per-construct coverage track
# --------------------------------------------------------------------------

def fig_coverage(fig_id: str, sample: str, construct: str, ident: Dict,
                 roles: Dict[str, str], depth: Dict[int, float],
                 baseline: Optional[float], trim: int, policy: str) -> Tuple[Svg, str]:
    c = ident["constructs"][construct]
    length = int(c["length"]) or (max(depth) if depth else 1)
    blocks = c["blocks"]
    trims = trim_regions(c, trim, policy)

    W, L, R = 980, 64, 96
    PW = W - L - R
    top_title, strip_y, strip_h = 26, 46, 22
    py0, ph = strip_y + strip_h + 10, 230
    H = py0 + ph + 150

    nbin = min(length, PW)
    bins: List[Optional[Tuple[float, float, float]]] = []
    for i in range(nbin):
        a = 1 + int(i * length / nbin)
        z = max(a, int((i + 1) * length / nbin))
        vals = [depth[p] for p in range(a, z + 1) if p in depth]
        bins.append((sum(vals) / len(vals), min(vals), max(vals)) if vals else None)
    have = [b for b in bins if b]
    if not have:
        raise ValueError("no depth values inside the contig")
    ymax = max(b[2] for b in have)
    if baseline:
        ymax = max(ymax, baseline)
    ticks = nice_ticks(0, ymax * 1.08)
    ytop = ticks[-1]

    def X(pos: float) -> float:
        return L + (pos - 1) / max(1, length - 1) * PW if length > 1 else L

    def Y(v: float) -> float:
        return py0 + ph - v / ytop * ph

    own = roles.get(construct, "neutral")
    svg = Svg(W, H, fig_id, "Coverage of %s in %s" % (construct, sample),
              "Per-base depth over the %s vector contig (%s bp) in sample %s, "
              "with construct membership blocks shaded." % (construct, length, sample))
    hatch_shared = svg.hatch("shared", "neutral", 45, 6, 1)
    hatch_trim = svg.hatch("trim", "trim", 135, 5, 1)

    svg.text(L, top_title - 6, "%s · %s" % (construct, sample), size=15, weight="600")
    svg.text(L + text_w("%s · %s" % (construct, sample), 15) + 12, top_title - 6,
             "per-base depth over the vector contig, %s bp" % "{:,}".format(length),
             size=12, fill="text2")

    # membership strip + background shading
    block_means = []
    for i, b in enumerate(blocks):
        a, z = int(b["start"]), int(b["end"])
        mem = list(b["membership"])
        uniq = len(mem) == 1
        x0, x1 = X(a), X(z) if z < length else L + PW
        label = ("%s only" % mem[0]) if uniq else " + ".join(mem)
        tip = "%d-%d (%s bp): %s%s" % (a, z, "{:,}".format(z - a + 1),
                                       "unique to " if uniq else "shared by ", label)
        if b.get("features"):
            tip += " [%s]" % b["features"]
        svg.open(**{"class": "block", "data_start": a, "data_end": z,
                    "data_membership": ",".join(mem)})
        if uniq:
            svg.rect(x0, py0, x1 - x0, ph, fill=roles.get(mem[0], "neutral"),
                     fill_opacity=0.10)
            svg.rect(x0, strip_y, x1 - x0, strip_h, fill=roles.get(mem[0], "neutral"),
                     title=tip)
        else:
            svg.rect(x0, py0, x1 - x0, ph, fill="shared", fill_opacity=0.35)
            svg.rect(x0, strip_y, x1 - x0, strip_h, fill="shared", title=tip)
            svg.add('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"/>'
                    % (x0, strip_y, x1 - x0, strip_h, hatch_shared))
        # 2px surface gap between adjacent blocks
        if i:
            svg.line(x0, strip_y, x0, strip_y + strip_h, stroke="surface",
                     stroke_width=2)
        fitted = label
        if text_w(fitted, 10) > (x1 - x0) - 6:
            fitted = (mem[0] if uniq else "%d shared" % len(mem))
        if text_w(fitted, 10) <= (x1 - x0) - 6:
            svg.text((x0 + x1) / 2, strip_y + 15, fitted, size=10, anchor="middle",
                     fill="surface" if uniq else "text", weight="600")
        svg.close()

        # retained mean for this block
        lo = a if not inside(a, trims) else None
        vals = [depth[p] for p in range(a, z + 1) if p in depth and not inside(p, trims)]
        if vals:
            m = sum(vals) / len(vals)
            block_means.append((a, z, m, mem))

    # trimmed regions
    for a, z in trims:
        x0, x1 = X(a), X(z)
        svg.add('<rect class="trimmed" x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                'fill="%s" data-start="%d" data-end="%d"><title>%s</title></rect>'
                % (x0, py0, max(1.0, x1 - x0), ph, hatch_trim, a, z,
                   escape("%d-%d excluded from copy number (trim %d bp, %s)"
                          % (a, z, trim, policy))))

    # grid + y axis
    for t in ticks:
        y = Y(t)
        svg.line(L, y, L + PW, y, stroke="grid")
        svg.text(L - 6, y + 4, fmt(t, 0), size=10, fill="text2", anchor="end")
    svg.text(16, py0 + ph / 2, "depth (reads per base)", size=11, fill="text2",
             anchor="middle", transform="rotate(-90 16 %.1f)" % (py0 + ph / 2))

    # envelope (min-max per bin) and mean line
    env_hi, env_lo, mean_pts = [], [], []
    for i, b in enumerate(bins):
        if not b:
            continue
        x = L + (i + 0.5) / nbin * PW
        env_hi.append((x, Y(b[2])))
        env_lo.append((x, Y(b[1])))
        mean_pts.append((x, Y(b[0])))
    if env_hi:
        svg.poly(env_hi + env_lo[::-1], closed=True, fill=own, fill_opacity=0.22,
                 **{"class": "envelope"})
        svg.poly(mean_pts, stroke=own, stroke_width=1.6, stroke_linejoin="round",
                 **{"class": "depth"})

    # per-block retained mean, in copies per haploid genome
    for a, z, m, mem in block_means:
        x0, x1 = X(a), X(z)
        y = Y(m)
        cp = (m / baseline) if baseline else None
        tip = "%s %d-%d: mean retained depth %sx%s" % (
            "+".join(mem), a, z, fmt(m), (" = %s copies/haploid (%s)"
                                          % (fmt(cp, 2), " + ".join(mem)))
            if cp is not None else "")
        svg.line(x0, y, x1, y, stroke="text", stroke_width=2, title=tip,
                 **{"class": "block-mean"})
        if cp is not None:
            lab = "%s copies" % fmt(cp, 1 if cp >= 10 else 2)
            if text_w(lab, 10) < (x1 - x0) - 4:
                svg.text((x0 + x1) / 2, y - 5, lab, size=10, anchor="middle",
                         fill="text", weight="600")

    # baseline
    if baseline:
        y = Y(baseline)
        svg.line(L, y, L + PW, y, stroke="text2", stroke_width=1.2,
                 stroke_dasharray="5 4", **{"class": "baseline"},
                 title="single-genome baseline %sx = 1 copy per haploid genome"
                 % fmt(baseline, 2))
        svg.text(L + PW + 6, y - 2, "baseline", size=10, fill="text2")
        svg.text(L + PW + 6, y + 10, "%sx = 1 copy" % fmt(baseline), size=10,
                 fill="text2")

    # x axis
    svg.line(L, py0 + ph, L + PW, py0 + ph, stroke="axis")
    for t in nice_ticks(1, length, 8):
        if t < 1 or t > length:
            continue
        x = X(t)
        svg.line(x, py0 + ph, x, py0 + ph + 4, stroke="axis")
        svg.text(x, py0 + ph + 16, "{:,}".format(int(t)), size=10, fill="text2",
                 anchor="middle")
    svg.text(L + PW / 2, py0 + ph + 32, "position on %s (bp)" % construct, size=11,
             fill="text2", anchor="middle")

    # legend
    items = []
    present_unique = []
    for b in blocks:
        if len(b["membership"]) == 1 and b["membership"][0] not in present_unique:
            present_unique.append(b["membership"][0])
    for u in present_unique:
        items.append(("unique to %s" % u, '<rect width="12" height="12" rx="2" %s/>'
                      % Svg.paint(roles.get(u, "neutral"))))
    if any(len(b["membership"]) > 1 for b in blocks):
        items.append(("shared (strip names the carriers)",
                      '<rect width="12" height="12" rx="2" %s/><rect width="12" '
                      'height="12" fill="%s"/>' % (Svg.paint("shared"), hatch_shared)))
    if trims:
        items.append(("excluded by trim (%d bp, %s)" % (trim, policy),
                      '<rect width="12" height="12" fill="%s" %s stroke-width="0.5"/>'
                      % (hatch_trim, Svg.paint(None, "axis"))))
    items.append(("depth, mean per bin", '<line x1="0" y1="6" x2="12" y2="6" '
                  'stroke-width="2" %s/>' % Svg.paint(None, own)))
    items.append(("min-max per bin", '<rect width="12" height="12" fill-opacity="0.22" %s/>'
                  % Svg.paint(own)))
    items.append(("block mean (retained bases)", '<line x1="0" y1="6" x2="12" y2="6" '
                  'stroke-width="2" %s/>' % Svg.paint(None, "text")))
    if baseline:
        items.append(("baseline", '<line x1="0" y1="6" x2="12" y2="6" stroke-width="1.2" '
                      'stroke-dasharray="3 2" %s/>' % Svg.paint(None, "text2")))
    y = legend_row(svg, L, py0 + ph + 56, items, L + PW)
    notes = []
    if not baseline:
        notes.append("No baseline depth available: block means are shown in reads, "
                     "not copies. Supply --copy-number or --baseline.")
    notes.append("Each construct's copy number is the mean depth over its own "
                 "unique blocks; shared blocks carry the sum of their carriers.")
    for n in notes:
        svg.text(L, y + 4, n, size=10, fill="text2")
        y += 15
    svg.h = int(y + 10)
    caption = ("Per-base depth over %s in %s. Strip: membership blocks from the "
               "identity map (coloured = unique to one construct, hatched grey = "
               "shared, labelled by carriers). Hatched overlay: bases excluded by "
               "the %d bp %s trim. Black bars: mean retained depth per block%s."
               % (construct, sample, trim, policy,
                  ", in copies per haploid genome against the %sx baseline (dashed)"
                  % fmt(baseline) if baseline else ""))
    return svg, caption


# --------------------------------------------------------------------------
# figure 2: genome integration-site overview
# --------------------------------------------------------------------------

CONF_ORDER = ["certain", "inferred", "unresolved"]


def locus_marker(cx: float, cy: float, conf: str, role: str) -> str:
    if conf == "certain":
        return '<circle cx="%.1f" cy="%.1f" r="4.5" %s stroke-width="1.5"/>' % (
            cx, cy, Svg.paint(role, "surface"))
    if conf == "inferred":
        return '<circle cx="%.1f" cy="%.1f" r="4" %s stroke-width="2"/>' % (
            cx, cy, Svg.paint("surface", role))
    r = 4.5
    return ('<polygon points="%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f" %s '
            'stroke-width="1.5"/>' % (cx, cy - r, cx + r, cy, cx, cy + r, cx - r, cy,
                                      Svg.paint("surface", "neutral")))


def fig_genome(fai: List[Tuple[str, int]], called: List[Dict[str, str]],
               order: Sequence[str], roles: Dict[str, str], exclude: Sequence[str],
               min_frac: float, max_contigs: int,
               summary: Optional[Dict]) -> Tuple[Svg, str, List[str]]:
    notes: List[str] = []
    excl = set(exclude)
    host = [(n, l) for n, l in fai if n not in excl]
    if len(host) < len(fai):
        notes.append("%d vector contig(s) in the .fai left out of the host genome."
                     % (len(fai) - len(host)))
    loci_contigs = {r.get("contig", "") for r in called}
    longest = max((l for _, l in host), default=1)
    keep_names = {n for n, l in host if l >= min_frac * longest or n in loci_contigs}
    if len(keep_names) > max_contigs:
        by_len = sorted(((l, n) for n, l in host if n in keep_names
                         and n not in loci_contigs), reverse=True)
        room = max(0, max_contigs - len(loci_contigs & keep_names))
        keep_names = (loci_contigs & keep_names) | {n for _, n in by_len[:room]}
    shown = [(n, l) for n, l in host if n in keep_names]
    omitted = len(host) - len(shown)
    if omitted:
        notes.append("%d short contig(s) without called loci not drawn (below %.1f%% "
                     "of the longest, or past %d contigs)."
                     % (omitted, min_frac * 100, max_contigs))
    lengths = dict(shown)
    off = [r for r in called if r.get("contig") not in lengths]
    if off:
        notes.append("%d called locus/loci on contigs absent from the .fai: %s."
                     % (len(off), ", ".join(sorted({r["contig"] for r in off}))))

    name_w = max([text_w(n, 11) for n, _ in shown] + [40])
    L = 20 + name_w + 10
    W, R = 980, 90
    PW = W - L - R
    row = 22
    y0 = 64
    H = y0 + row * max(1, len(shown)) + 140

    svg = Svg(W, H, "genome_overview", "Integration sites across the host genome",
              "Host chromosomes drawn to length from the .fai, with each called "
              "locus marked by attributed construct and confidence.")
    svg.text(20, 22, "Called integration sites", size=15, weight="600")
    n_by = {}
    for r in called:
        n_by[r.get("confidence", "")] = n_by.get(r.get("confidence", ""), 0) + 1
    sub = "%d called loci on %d contigs drawn to scale" % (len(called), len(shown))
    if n_by:
        sub += " · " + ", ".join("%d %s" % (n_by[k], k) for k in CONF_ORDER if k in n_by)
    svg.text(20, 40, sub, size=12, fill="text2")

    by_contig: Dict[str, List[Dict[str, str]]] = {}
    for r in called:
        by_contig.setdefault(r.get("contig", ""), []).append(r)

    for i, (name, ln) in enumerate(shown):
        y = y0 + i * row
        svg.open(**{"class": "chrom", "data_name": name, "data_length": ln})
        svg.text(L - 10, y + 4, name, size=11, fill="text2", anchor="end")
        bw = ln / longest * PW
        svg.rect(L, y - 4, bw, 8, fill="shared", rx=4,
                 title="%s: %s" % (name, fmt_bp(ln)))
        svg.close()
        for r in by_contig.get(name, []):
            try:
                pos = (int(r["start"]) + int(r.get("end") or r["start"])) / 2.0
            except (KeyError, ValueError):
                continue
            x = L + pos / longest * PW
            conf = r.get("confidence", "") or "unresolved"
            cons = r.get("construct", "")
            role = roles.get(cons, "neutral") if conf != "unresolved" else "neutral"
            dash = ' stroke-dasharray="2 2"' if conf == "unresolved" else ""
            tip = "%s:%s  %s (%s%s)  support %s%s" % (
                name, "{:,}".format(int(pos)), cons or "unattributed", conf,
                (", " + r["evidence"]) if r.get("evidence") else "",
                r.get("total_support", "?"),
                ("  candidates: " + r["candidates"]) if conf == "unresolved"
                and r.get("candidates") else "")
            svg.add('<g class="locus" data-contig="%s" data-pos="%d" data-construct="%s" '
                    'data-confidence="%s"><title>%s</title>'
                    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke-width="2"%s %s/>%s</g>'
                    % (escape(name), pos, escape(cons), escape(conf), escape(tip),
                       x, y - 9, x, y + 7, dash, Svg.paint(None, role),
                       locus_marker(x, y - 10, conf, role)))

    # scale bar
    yb = y0 + len(shown) * row + 6
    for t in nice_ticks(0, longest, 6):
        if t > longest:
            continue
        x = L + t / longest * PW
        svg.line(x, yb, x, yb + 4, stroke="axis")
        svg.text(x, yb + 16, fmt_bp(t) if t else "0", size=10, fill="text2",
                 anchor="middle")
    svg.line(L, yb, L + PW, yb, stroke="axis")

    # legend, from the data
    items = []
    for c in order:
        n = sum(1 for r in called if r.get("construct") == c
                and r.get("confidence") != "unresolved")
        if n:
            items.append(("%s (%d)" % (c, n), '<rect width="12" height="12" rx="2" %s/>'
                          % Svg.paint(roles.get(c, "neutral"))))
    for c in sorted({r.get("construct") for r in called
                     if r.get("construct") and r.get("construct") not in order}):
        items.append(("%s (not in identity map)" % c,
                      '<rect width="12" height="12" rx="2" %s/>' % Svg.paint("neutral")))
    for conf in CONF_ORDER:
        if n_by.get(conf):
            label = {"certain": "certain: construct-unique sequence",
                     "inferred": "inferred: vector end or copy ceiling",
                     "unresolved": "unresolved: no construct attributable"}[conf]
            items.append(("%s (%d)" % (label, n_by[conf]),
                          locus_marker(6, 6, conf, "text2" if conf != "unresolved"
                                       else "neutral")))
    y = legend_row(svg, 20, yb + 40, items, W - 20)
    if summary:
        extra = []
        if summary.get("n_sub") is not None:
            extra.append("%s sub-threshold locus/loci not drawn (the noise tail)"
                         % summary["n_sub"])
        if summary.get("n_entries") is not None and summary.get("n_unified") is not None:
            extra.append("%s per-mapping entries collapsed to %s loci"
                         % (summary["n_entries"], summary["n_unified"]))
        if extra:
            notes.insert(0, "; ".join(extra) + ".")
    if not called:
        notes.insert(0, "No called loci in called_sites.tsv.")
    for n in notes:
        svg.text(20, y + 4, n, size=10, fill="text2")
        y += 15
    svg.h = int(y + 10)
    caption = ("Host genome overview: %d contigs drawn to scale from the .fai; "
               "%d called loci from reconcile_sites.py, coloured by attributed "
               "construct (filled = certain, open = inferred, grey diamond = "
               "unresolved)." % (len(shown), len(called)))
    return svg, caption, notes


# --------------------------------------------------------------------------
# figure 3: copy number
# --------------------------------------------------------------------------

def cn_label(path: str, data: Dict, used: Sequence[str]) -> str:
    t = data.get("trim") or {}
    lab = t.get("policy") or "copy number"
    if t.get("trim_boundaries") is not None:
        lab = "%s %s bp" % (lab, t["trim_boundaries"])
    if lab in used:
        lab = os.path.splitext(os.path.basename(path))[0]
    return lab


def fig_copy_number(cns: List[Tuple[str, Dict]], order: Sequence[str],
                    roles: Dict[str, str]) -> Tuple[Svg, str]:
    samples: List[str] = []
    constructs: List[str] = []
    for _, d in cns:
        for s, sd in (d.get("samples") or {}).items():
            if s not in samples:
                samples.append(s)
            for c in (sd.get("copies_per_haploid") or {}):
                if c not in constructs:
                    constructs.append(c)
    constructs.sort(key=lambda c: order.index(c) if c in order else len(order))
    bars = []  # (construct, sample, mode_index, label, value, lo, hi, source)
    vmax = 0.0
    for c in constructs:
        for s in samples:
            for mi, (lab, d) in enumerate(cns):
                sd = (d.get("samples") or {}).get(s)
                if sd is None:
                    continue
                v = (sd.get("copies_per_haploid") or {}).get(c)
                solve = sd.get("solve") or {}
                ub = (solve.get("unique_blocks") or {}).get(c) or {}
                lo = ub.get("min") if (ub.get("n_blocks") or 0) > 1 else None
                hi = ub.get("max") if (ub.get("n_blocks") or 0) > 1 else None
                src = (solve.get("source") or {}).get(c)
                bars.append((c, s, mi, lab, v, lo, hi, src))
                for x in (v, hi):
                    if isinstance(x, (int, float)):
                        vmax = max(vmax, x)
    if not bars:
        raise ValueError("no copies_per_haploid in the copy-number JSON")
    ticks = nice_ticks(0, max(vmax, 1) * 1.12)
    ytop = ticks[-1]

    bw, gap, ggap = 30, 2, 40
    per_group = {c: sum(1 for b in bars if b[0] == c) for c in constructs}
    L, R = 64, 24
    PW = sum(per_group[c] * (bw + gap) - gap for c in constructs) + ggap * (len(constructs) + 1)
    W = max(720, L + PW + R)
    if L + PW + R < W:
        bars_w = PW - ggap * (len(constructs) + 1)
        ggap = min(160, (W - L - R - bars_w) / (len(constructs) + 1))
    py0, ph = 64, 260
    svg = Svg(W, 900, "copy_number", "Copy number per construct",
              "Copies per haploid genome for each construct, sample and "
              "copy-number run.")
    textures = [None, svg.hatch("m1", "surface", 45, 5, 1.6),
                svg.hatch("m2", "surface", 135, 5, 1.6),
                svg.hatch("m3", "surface", 90, 4, 1.4)]
    svg.text(L, 22, "Copy number", size=15, weight="600")
    svg.text(L, 40, "copies per haploid genome, from unique-block depth"
             + (" · %d runs side by side" % len(cns) if len(cns) > 1 else ""),
             size=12, fill="text2")

    def Y(v):
        return py0 + ph - v / ytop * ph

    for t in ticks:
        svg.line(L, Y(t), W - R, Y(t), stroke="grid")
        svg.text(L - 6, Y(t) + 4, fmt(t, 0 if ytop >= 10 else 1), size=10,
                 fill="text2", anchor="end")
    svg.text(16, py0 + ph / 2, "copies per haploid genome", size=11, fill="text2",
             anchor="middle", transform="rotate(-90 16 %.1f)" % (py0 + ph / 2))

    x = L + ggap
    flagged = False
    for c in constructs:
        role = roles.get(c, "neutral")
        gx0 = x
        for b in [b for b in bars if b[0] == c]:
            _, s, mi, lab, v, lo, hi, src = b
            tip = "%s · %s · %s: %s copies/haploid%s%s" % (
                c, s, lab, fmt(v, 2) if isinstance(v, (int, float)) else "n/a",
                (" (unique blocks %s-%s)" % (fmt(lo, 2), fmt(hi, 2))) if lo is not None else "",
                ("; source: %s" % src) if src else "")
            svg.open(**{"class": "cn-bar", "data_construct": c, "data_sample": s,
                        "data_run": lab, "data_value": v})
            if isinstance(v, (int, float)) and v > 0:
                h = v / ytop * ph
                y = py0 + ph - h
                r = min(4.0, bw / 2, h / 2)
                # rounded data end, square at the baseline
                svg.add('<path d="M%.1f,%.1f v%.1f a%.1f,%.1f 0 0 1 %.1f,%.1f h%.1f '
                        'a%.1f,%.1f 0 0 1 %.1f,%.1f v%.1f z" %s><title>%s</title></path>'
                        % (x, py0 + ph, -(h - r), r, r, r, -r, bw - 2 * r, r, r, r, r,
                           h - r, Svg.paint(role), escape(tip)))
                tex = textures[mi % len(textures)]
                if tex:
                    svg.add('<path d="M%.1f,%.1f v%.1f a%.1f,%.1f 0 0 1 %.1f,%.1f h%.1f '
                            'a%.1f,%.1f 0 0 1 %.1f,%.1f v%.1f z" fill="%s"/>'
                            % (x, py0 + ph, -(h - r), r, r, r, -r, bw - 2 * r, r, r, r,
                               r, h - r, tex))
                svg.text(x + bw / 2, y - 5 if lo is None else Y(hi) - 5,
                         fmt(v, 1 if v >= 10 else 2), size=10, anchor="middle",
                         weight="600")
                if lo is not None and hi is not None:
                    xm = x + bw / 2
                    svg.line(xm, Y(lo), xm, Y(hi), stroke="text", stroke_width=1.5)
                    svg.line(xm - 5, Y(lo), xm + 5, Y(lo), stroke="text", stroke_width=1.5)
                    svg.line(xm - 5, Y(hi), xm + 5, Y(hi), stroke="text", stroke_width=1.5)
                if src and src != "unique block":
                    flagged = True
                    svg.text(x + bw / 2, py0 + ph - 6, "†", size=12, anchor="middle",
                             fill="surface", weight="700")
            else:
                svg.text(x + bw / 2, py0 + ph - 6, "n/a", size=10, anchor="middle",
                         fill="text2")
            svg.close()
            if len(samples) > 1:
                svg.text(x + bw / 2, py0 + ph + 14, s, size=9, fill="text2",
                         anchor="middle")
            x += bw + gap
        x -= gap
        svg.text((gx0 + x) / 2, py0 + ph + (32 if len(samples) > 1 else 18), c,
                 size=12, anchor="middle", weight="600")
        x += ggap
    svg.line(L, py0 + ph, W - R, py0 + ph, stroke="axis")

    y = py0 + ph + 58
    items = []
    if len(cns) > 1:
        for mi, (lab, _) in enumerate(cns):
            tex = textures[mi % len(textures)]
            sw = '<rect width="12" height="12" rx="2" %s/>' % Svg.paint("text2")
            if tex:
                sw += '<rect width="12" height="12" fill="%s"/>' % tex
            items.append((lab, sw))
    if any(b[5] is not None for b in bars):
        items.append(("range of per-unique-block estimates",
                      '<line x1="6" y1="0" x2="6" y2="12" stroke-width="1.5" %s/>'
                      '<line x1="2" y1="0" x2="10" y2="0" stroke-width="1.5" %s/>'
                      '<line x1="2" y1="12" x2="10" y2="12" stroke-width="1.5" %s/>'
                      % ((Svg.paint(None, "text"),) * 3)))
    if flagged:
        items.append(("† not from a unique block (see source in the report)", ""))
    if items:
        y = legend_row(svg, L, y, items, W - R)

    # shared-block closure: the additivity check
    svg.text(L, y + 6, "Shared-block closure (observed ÷ sum of carriers' unique-block "
             "depth − 1)", size=12, weight="600")
    y += 24
    for s in samples:
        for lab, d in cns:
            sd = (d.get("samples") or {}).get(s)
            if sd is None:
                continue
            solve = sd.get("solve") or {}
            cl = solve.get("closure") or {}
            head = "%s%s: " % (s, (" · " + lab) if len(cns) > 1 else "")
            if cl.get("deficit_pct") is not None:
                msg = "%s%% (observed %sx vs %sx predicted)" % (
                    ("%+.1f" % cl["deficit_pct"]), fmt(cl.get("observed")),
                    fmt(cl.get("predicted")))
                bad = abs(cl["deficit_pct"]) >= 10
            elif solve:
                msg, bad = "n/a (no shared class with every carrier measured)", False
            else:
                msg, bad = "n/a (copy number not measured from depth: %s)" % (
                    solve.get("method") or "fallback route"), False
            svg.add('<g class="closure" data-sample="%s" data-run="%s" data-deficit="%s">'
                    % (escape(s), escape(lab), cl.get("deficit_pct")))
            svg.add('<text x="%.1f" y="%.1f" font-size="11" %s>%s<tspan %s%s>%s'
                    '</tspan></text>' % (L, y, Svg.paint("text2"), escape(head),
                                         Svg.paint("text"),
                                         ' font-weight="600"' if bad else "",
                                         escape(msg)))
            svg.close()
            y += 16
            for cav in (solve.get("caveats") or [])[:2]:
                first = cav.split(". ")[0].strip().rstrip(".") + "."
                svg.text(L + 16, y, "! " + first[:150], size=10, fill="text2")
                y += 14
    svg.text(L, y + 8, "A closure well below 0 means shared sequence carries less depth "
             "than additivity predicts (e.g. backbone loss in concatemers).",
             size=10, fill="text2")
    svg.h = int(y + 22)
    caption = ("Copies per haploid genome per construct and sample, from each "
               "construct's unique-block depth%s. Whiskers: range of individual "
               "unique-block estimates where a construct has more than one. Below: "
               "shared-block closure from the additivity cross-check."
               % ("; runs compared side by side: " + ", ".join(l for l, _ in cns)
                  if len(cns) > 1 else ""))
    return svg, caption


# --------------------------------------------------------------------------
# figure 4: architecture / reconciliation
# --------------------------------------------------------------------------

def fig_architecture(summary: Dict, order: Sequence[str], roles: Dict[str, str],
                     cn: Optional[Dict]) -> Tuple[Svg, str]:
    pc = summary.get("per_construct") or {}
    if not pc:
        raise ValueError("site summary has no per_construct block")
    cons = sorted(pc, key=lambda c: order.index(c) if c in order else len(order))
    params = summary.get("parameters") or {}
    ploidy = params.get("ploidy")
    mcps = params.get("min_copies_per_site")

    measured: Dict[str, float] = {}
    if cn:
        prim = cn.get("primary_mapping")
        for s, sd in (cn.get("samples") or {}).items():
            cps = (sd.get("architecture") or {}).get("copies_per_site")
            if prim and isinstance(cps, (int, float)) and prim not in measured:
                measured[prim] = float(cps)

    def g(c, k):
        v = pc[c].get(k)
        return v if isinstance(v, (int, float)) else 0

    smax = max([max(g(c, "ceiling"), g(c, "attributed") + g(c, "unresolved_candidate_for"),
                    g(c, "site_high")) for c in cons] + [1])
    sticks = nice_ticks(0, smax * 1.1, min(8, int(smax * 1.1) + 1))
    stop = sticks[-1] if sticks[-1] > 0 else 1

    rows = []
    cmax = 0.0
    for c in cons:
        cph = pc[c].get("copies_per_haploid")
        cpc = cph * ploidy if isinstance(cph, (int, float)) and ploidy else None
        floor, high, att = g(c, "site_floor"), g(c, "site_high"), g(c, "attributed")
        lo = cpc / high if cpc is not None and high else None
        hi = cpc / floor if cpc is not None and floor else None
        pt = cpc / att if cpc is not None and att else None
        rows.append((c, cpc, lo, hi, pt))
        for v in (lo, hi, pt, measured.get(c), mcps):
            if isinstance(v, (int, float)):
                cmax = max(cmax, v)
    cticks = nice_ticks(0, max(cmax, 1) * 1.1)
    ctop = cticks[-1]

    name_w = max([text_w(c, 12) for c in cons] + [50])
    L = 20 + name_w + 12
    W = 980
    gapc = 70
    P1 = (W - L - 24 - gapc) * 0.55
    P2 = (W - L - 24 - gapc) - P1
    X2 = L + P1 + gapc
    row = 58
    y0 = 112
    svg = Svg(W, 900, "architecture", "Site-count bound and copies per site",
              "Per construct: site floor, attributed sites, copy-number ceiling and "
              "unresolved candidates from reconcile_sites.py, and implied copies per site.")
    hatch_ve = svg.hatch("ve", "surface", 45, 5, 1.6)
    svg.text(20, 22, "Site-count bound and copies per site", size=15, weight="600")
    svg.text(20, 40, "from reconcile_sites.py%s" % (
        " · ploidy %s, min copies per site %s" % (fmt(ploidy, 1), fmt(mcps, 0))
        if ploidy else ""), size=12, fill="text2")
    svg.text(L, 66, "Integration sites", size=12, weight="600")
    svg.text(X2, 66, "Copies per site (copies per cell ÷ sites)", size=12, weight="600")

    def SX(v):
        return L + v / stop * P1

    def CX(v):
        return X2 + v / ctop * P2

    yb = y0 + row * len(cons) - 14
    for t in sticks:
        svg.line(SX(t), y0 - 16, SX(t), yb, stroke="grid")
        svg.text(SX(t), yb + 14, fmt(t, 0), size=10, fill="text2", anchor="middle")
    for t in cticks:
        svg.line(CX(t), y0 - 16, CX(t), yb, stroke="grid")
        svg.text(CX(t), yb + 14, fmt(t, 0 if ctop >= 5 else 1), size=10,
                 fill="text2", anchor="middle")
    if isinstance(mcps, (int, float)):
        svg.line(CX(mcps), y0 - 16, CX(mcps), yb, stroke="text2", stroke_dasharray="4 3",
                 title="min copies per site assumed by the ceiling: %s" % fmt(mcps, 1))
        svg.text(CX(mcps), yb + 28, "assumed min %s" % fmt(mcps, 0), size=9,
                 fill="text2", anchor="middle")

    for i, (c, cpc, lo, hi, pt) in enumerate(rows):
        y = y0 + i * row
        role = roles.get(c, "neutral")
        p = pc[c]
        svg.text(L - 12, y + 4, c, size=12, anchor="end", weight="600")
        svg.open(**{"class": "arch-row", "data_construct": c,
                    "data_floor": p.get("site_floor"), "data_high": p.get("site_high"),
                    "data_ceiling": p.get("ceiling"),
                    "data_unresolved": p.get("unresolved_candidate_for")})
        x = L
        for k, lab, fill, tex in (("certain", "certain", role, None),
                                  ("vector_end", "vector end", role, hatch_ve),
                                  ("copy_ceiling", "copy ceiling", role, "light")):
            n = g(c, k)
            if not n:
                continue
            w = n / stop * P1
            if tex == "light":
                svg.rect(x, y - 8, w - 2, 16, fill=fill, fill_opacity=0.35, rx=2,
                         title="%s: %d %s" % (c, n, lab))
            else:
                svg.rect(x, y - 8, w - 2, 16, fill=fill, rx=2,
                         title="%s: %d %s" % (c, n, lab))
                if tex:
                    svg.add('<rect x="%.1f" y="%.1f" width="%.1f" height="16" fill="%s"/>'
                            % (x, y - 8, w - 2, tex))
            x += w
        un = g(c, "unresolved_candidate_for")
        if un:
            w = un / stop * P1
            svg.rect(x, y - 8, w - 2, 16, fill="surface", stroke="neutral", rx=2,
                     stroke_dasharray="3 2", stroke_width=1.2,
                     title="%s is a candidate for %d unresolved locus/loci" % (c, un))
            x += w
        ceil = p.get("ceiling")
        if isinstance(ceil, (int, float)):
            svg.line(SX(ceil), y - 13, SX(ceil), y + 13, stroke="text", stroke_width=2.5,
                     title="copy-number ceiling: %s site(s)" % ceil)
            svg.text(SX(ceil), y - 16, "ceiling %s" % ceil, size=9, fill="text2",
                     anchor="middle")
        fl = p.get("site_floor")
        if isinstance(fl, (int, float)) and fl > 0:
            svg.line(SX(fl), y + 9, SX(fl), y + 14, stroke="text", stroke_width=1.5)
        rng = "sites %s" % (p.get("site_range") or "?")
        if un:
            rng += " · +%d unresolved" % un
        if ceil is None:
            rng += " · no ceiling (no copy number)"
        svg.text(L, y + 24, rng, size=10, fill="text2")
        # copies per site
        if lo is not None or hi is not None or pt is not None:
            a = lo if lo is not None else 0
            b = hi if hi is not None else ctop
            svg.line(CX(a), y, CX(b), y, stroke=role, stroke_width=3,
                     stroke_linecap="round",
                     title="%s: %s copies per cell over %s site(s) = %s-%s copies/site"
                     % (c, fmt(cpc, 1), p.get("site_range"), fmt(lo, 1),
                        fmt(hi, 1) if hi is not None else "unbounded"))
            if hi is None:
                svg.poly([(CX(b) - 6, y - 5), (CX(b), y), (CX(b) - 6, y + 5)],
                         stroke=role, stroke_width=2)
            if pt is not None:
                svg.circle(CX(pt), y, 5, fill=role, stroke="surface", stroke_width=2,
                           title="%s copies per site at %d attributed site(s)"
                           % (fmt(pt, 1), g(c, "attributed")))
            lab = "%s-%s" % (fmt(lo, 1) if lo is not None else "?",
                             fmt(hi, 1) if hi is not None else "∞")
            svg.text(CX(b) + 8 if hi is not None else CX(a), y + 4 if hi is not None
                     else y + 20, lab, size=10, fill="text2")
        else:
            svg.text(X2, y + 4, "no copy number", size=10, fill="text2")
        if c in measured:
            m = measured[c]
            svg.poly([(CX(m), y - 7), (CX(m) + 7, y), (CX(m), y + 7), (CX(m) - 7, y)],
                     closed=True, fill="surface", stroke="text", stroke_width=1.5,
                     title="copies per site measured by copy_number.py from "
                     "flank/internal junction reads: %s" % fmt(m, 2))
        svg.close()

    y = yb + 52
    items = [("certain", '<rect width="12" height="12" rx="2" %s/>' % Svg.paint("text2")),
             ("vector end", '<rect width="12" height="12" rx="2" %s/><rect width="12" '
              'height="12" fill="%s"/>' % (Svg.paint("text2"), hatch_ve)),
             ("copy ceiling", '<rect width="12" height="12" rx="2" fill-opacity="0.35" %s/>'
              % Svg.paint("text2")),
             ("unresolved candidate", '<rect x="0.5" y="0.5" width="11" height="11" rx="2" '
              'stroke-dasharray="3 2" %s/>' % Svg.paint("surface", "neutral")),
             ("ceiling", '<line x1="6" y1="0" x2="6" y2="12" stroke-width="2.5" %s/>'
              % Svg.paint(None, "text")),
             ("copies/site range (dot: at attributed sites)",
              '<line x1="0" y1="6" x2="12" y2="6" stroke-width="3" %s/>'
              % Svg.paint(None, "text2"))]
    if measured:
        items.append(("measured from junctions", '<polygon points="6,0 12,6 6,12 0,6" '
                      'stroke-width="1.5" %s/>' % Svg.paint("surface", "text")))
    y = legend_row(svg, 20, y, items, W - 20)
    notes = []
    if summary.get("n_entries") is not None:
        notes.append("%s per-mapping entries collapsed to %s loci; %s called, %s "
                     "sub-threshold; %d unresolved locus/loci attributed to no construct."
                     % (summary.get("n_entries"), summary.get("n_unified"),
                        summary.get("n_called"), summary.get("n_sub"),
                        len(summary.get("unresolved") or [])))
    notes.append("Floor = certain + vector-end loci; ceiling = ceil(copies x ploidy ÷ "
                 "min copies per site). The site count is the range, never one number.")
    for n in notes:
        svg.text(20, y + 4, n, size=10, fill="text2")
        y += 15
    svg.h = int(y + 10)
    caption = ("Per construct: called sites by evidence tier with the copy-number "
               "ceiling and unresolved candidates (left), and the copies per site "
               "the bounded site count implies (right).")
    return svg, caption


# --------------------------------------------------------------------------
# filling the report template
# --------------------------------------------------------------------------

PLACEHOLDER = re.compile(r"\[\[FIG:([^\]]+)\]\]")


def fill_report(src: str, dst: str, manifest: Dict, fig_dir: str) -> Tuple[int, List[str]]:
    with open(src, encoding="utf-8") as fh:
        text = fh.read()
    figs = {f["id"]: f for f in manifest.get("figures", [])}
    skipped = {s["id"]: s["reason"] for s in manifest.get("skipped", []) if s.get("id")}
    rel_base = os.path.dirname(os.path.abspath(dst))
    missing: List[str] = []
    count = [0]

    def link(f):
        rel = os.path.relpath(os.path.join(fig_dir, f["file"]), rel_base)
        count[0] += 1
        return "![%s](%s)\n\n*Figure: %s*" % (f["title"], rel.replace(os.sep, "/"),
                                               f["caption"])

    def sub(m):
        key = m.group(1).strip()
        if key in figs:
            return link(figs[key])
        if key.startswith("coverage_"):
            cons = key[len("coverage_"):]
            hits = [f for f in manifest.get("figures", [])
                    if f.get("kind") == "coverage" and f.get("construct") == cons]
            if hits:
                return "\n\n".join(link(f) for f in hits)
        missing.append(key)
        if "<" in key:
            why = "the placeholder still holds a template name; replace it " \
                  "with the construct name from identity.json and rerun"
        else:
            why = skipped.get(key) or "no such figure was generated"
        return "> **Figure `%s` not generated:** %s" % (key, why)

    # placeholders inside HTML comments (the template's own instructions) are
    # documentation, not figure slots
    parts = re.split(r"(<!--.*?-->)", text, flags=re.S)
    out = "".join(p if p.startswith("<!--") else PLACEHOLDER.sub(sub, p)
                  for p in parts)
    with open(dst, "w", encoding="utf-8") as fh:
        fh.write(out)
    return count[0], missing


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def parse_kv(items: Sequence[str], what: str) -> List[Tuple[str, str]]:
    out = []
    for it in items:
        if "=" in it:
            k, v = it.split("=", 1)
            out.append((k.strip(), v.strip()))
        else:
            out.append(("", it.strip()))
    return out


def discover_depth(d: str, constructs: Sequence[str]) -> Dict[str, Dict[str, str]]:
    """depth_<sample>_<construct>.txt, depth_<construct>_<sample>.txt or
    depth_<construct>.txt (single sample) in a directory."""
    out: Dict[str, Dict[str, str]] = {}
    for p in sorted(glob.glob(os.path.join(d, "depth_*.txt"))):
        stem = os.path.basename(p)[len("depth_"):-len(".txt")]
        for c in sorted(constructs, key=len, reverse=True):
            if stem == c:
                out.setdefault("sample", {})[c] = p
            elif stem.endswith("_" + c):
                out.setdefault(stem[:-len(c) - 1], {})[c] = p
            elif stem.startswith(c + "_"):
                out.setdefault(stem[len(c) + 1:], {})[c] = p
            else:
                continue
            break
    return out


def apply_run(a: argparse.Namespace) -> Optional[Tuple[Dict, str]]:
    """Fill unset inputs from a run manifest; return (config, its base dir)
    built from the manifest when no --config was given."""
    if not a.run:
        return None
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from via_run import RunError, copy_number_config, load_run
    try:
        run = load_run(a.run)
    except RunError as e:
        print("  warning: %s" % e)
        return None
    st = lambda name, key: (run.resolve(run.stage(name)[key])  # noqa: E731
                            if run.stage(name).get(key) else None)
    a.out = a.out or os.path.join(run.base, "figures")
    a.identity = a.identity or run.identity_path()
    if not a.copy_number and st("copy_number", "json"):
        a.copy_number = [st("copy_number", "json")]
    a.called = a.called or st("reconcile", "called")
    a.site_summary = a.site_summary or st("reconcile", "summary")
    a.fai = a.fai or st("host", "fai")
    if a.config:
        return None
    try:
        return copy_number_config(run)[0], run.base
    except RunError as e:
        print("  warning: no depth configuration from the manifest (%s)" % e)
        return None


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Draw the report's standard figures (SVG) from the pipeline "
                    "outputs, and fill [[FIG:...]] placeholders in a report.")
    ap.add_argument("--run", default="", metavar="RUN.json",
                    help="run manifest: take every input below that is not "
                         "given explicitly from what the stages recorded "
                         "(identity, depth files, copy-number JSON, called "
                         "sites, site summary, host .fai)")
    ap.add_argument("--out", default="", help="output directory (default: "
                    "figures, beside the manifest with --run)")
    ap.add_argument("--identity", help="identity.json from construct_identity.py")
    ap.add_argument("--config", help="the copy_number.py config: supplies each "
                    "sample's depth_files, vector contig names and baseline_depth")
    ap.add_argument("--depth", action="append", default=[], metavar="SAMPLE:CONSTRUCT=PATH",
                    help="one samtools depth file; repeatable. Overrides --config.")
    ap.add_argument("--depth-dir", help="discover depth_<sample>_<construct>.txt "
                    "(or depth_<construct>_<sample>.txt, depth_<construct>.txt) here")
    ap.add_argument("--copy-number", action="append", default=[], metavar="[LABEL=]JSON",
                    help="copy_number.py --json output; repeatable, one bar per run "
                         "(e.g. two trim policies). The first also supplies baselines "
                         "and the trim used for the coverage tracks.")
    ap.add_argument("--called", help="called_sites.tsv from reconcile_sites.py")
    ap.add_argument("--site-summary", help="site_summary.json from reconcile_sites.py")
    ap.add_argument("--fai", help="host genome .fai (samtools faidx) for the overview")
    ap.add_argument("--baseline", action="append", default=[], metavar="SAMPLE=DEPTH",
                    help="single-genome baseline depth per sample (overrides JSON/config)")
    ap.add_argument("--trim", type=int, help="bases trimmed at block boundaries "
                    "(default: from the copy-number JSON, else config, else 150)")
    ap.add_argument("--min-contig-frac", type=float, default=0.01,
                    help="overview: draw contigs at least this fraction of the longest, "
                         "plus every contig with a called locus (default 0.01)")
    ap.add_argument("--max-contigs", type=int, default=60,
                    help="overview: at most this many contigs (default 60)")
    ap.add_argument("--fill", metavar="REPORT.md",
                    help="replace [[FIG:<id>]] placeholders in this report")
    ap.add_argument("--fill-out", metavar="OUT.md",
                    help="where to write the filled report (default: overwrite --fill)")
    a = ap.parse_args(argv)
    run_cfg = apply_run(a)
    a.out = a.out or "figures"

    os.makedirs(a.out, exist_ok=True)
    figures: List[Dict] = []
    skipped: List[Dict] = []

    def skip(fid: str, reason: str) -> None:
        skipped.append({"id": fid, "reason": reason})
        print("  skipped %s: %s" % (fid, reason))

    def emit(fid: str, kind: str, svg: Svg, caption: str, **extra) -> None:
        fn = fid + ".svg"
        svg.save(os.path.join(a.out, fn))
        rec = {"id": fid, "kind": kind, "file": fn, "title": svg.title,
               "caption": caption, "placeholder": "[[FIG:%s]]" % fid}
        rec.update(extra)
        figures.append(rec)
        print("  wrote %s" % os.path.join(a.out, fn))

    ident = None
    if a.identity:
        try:
            ident = load_json(a.identity)
        except (OSError, ValueError) as e:
            print("  warning: identity map unreadable (%s)" % e)
    order: List[str] = list((ident or {}).get("construct_order")
                            or (ident or {}).get("constructs", {}).keys())

    cfg, cfg_dir = run_cfg if run_cfg else (None, "")
    if a.config:
        try:
            cfg = load_json(a.config)
            cfg_dir = os.path.dirname(os.path.abspath(a.config))
        except (OSError, ValueError) as e:
            print("  warning: config unreadable (%s)" % e)

    cns: List[Tuple[str, Dict]] = []
    for lab, path in parse_kv(a.copy_number, "copy-number"):
        try:
            d = load_json(path)
        except (OSError, ValueError) as e:
            skip("copy_number", "copy-number JSON %s unreadable (%s)" % (path, e))
            continue
        cns.append((lab or cn_label(path, d, [l for l, _ in cns]), d))
    for _, d in cns:
        for s in (d.get("samples") or {}).values():
            for c in (s.get("copies_per_haploid") or {}):
                if c not in order:
                    order.append(c)

    summary = None
    if a.site_summary:
        try:
            summary = load_json(a.site_summary)
        except (OSError, ValueError) as e:
            print("  warning: site summary unreadable (%s)" % e)
    called: Optional[List[Dict[str, str]]] = None
    if a.called:
        try:
            called = read_tsv(a.called)
        except OSError as e:
            print("  warning: called-site table unreadable (%s)" % e)
    if summary:
        for c in summary.get("mappings") or []:
            if c not in order:
                order.append(c)
    for r in called or []:
        c = r.get("construct")
        if c and c not in order:
            order.append(c)
    roles = construct_roles(order)

    # ---- 1. coverage tracks
    depth_paths: Dict[str, Dict[str, str]] = {}
    contig_names: Dict[str, str] = {}
    if cfg:
        for s, sd in (cfg.get("samples") or {}).items():
            for c, p in (sd.get("depth_files") or {}).items():
                depth_paths.setdefault(s, {})[c] = p if os.path.isabs(p) \
                    else os.path.join(cfg_dir, p)
        for c, m in (cfg.get("mappings") or {}).items():
            if isinstance(m, dict) and m.get("contig"):
                contig_names[c] = m["contig"]
    if a.depth_dir and order:
        for s, m in discover_depth(a.depth_dir, order).items():
            for c, p in m.items():
                depth_paths.setdefault(s, {}).setdefault(c, p)
    single = None
    if cns and len((cns[0][1].get("samples") or {})) == 1:
        single = next(iter(cns[0][1]["samples"]))
    for k, p in parse_kv(a.depth, "depth"):
        if ":" in k:
            s, c = k.split(":", 1)
        else:
            s, c = (single or "sample"), k
        depth_paths.setdefault(s, {})[c] = p
    # a discovered single-sample "sample" key adopts the copy-number sample name
    if "sample" in depth_paths and single and single not in depth_paths:
        depth_paths[single] = depth_paths.pop("sample")

    baselines: Dict[str, float] = {}
    if cfg:
        for s, sd in (cfg.get("samples") or {}).items():
            if sd.get("baseline_depth"):
                baselines[s] = float(sd["baseline_depth"])
    if cns:
        for s, sd in (cns[0][1].get("samples") or {}).items():
            b = sd.get("baseline_depth") or sd.get("genome_coverage")
            if b:
                baselines[s] = float(b)
    for s, v in parse_kv(a.baseline, "baseline"):
        try:
            baselines[s or (single or "sample")] = float(v)
        except ValueError:
            print("  warning: --baseline %s=%s is not a number" % (s, v))

    trim, policy = 150, "trim_boundaries"
    if cfg:
        trim = int(cfg.get("trim_boundaries") or cfg.get("read_length") or trim)
    if cns:
        t = cns[0][1].get("trim") or {}
        if t.get("trim_boundaries") is not None:
            trim = int(t["trim_boundaries"])
        policy = t.get("policy") or policy
    if a.trim is not None:
        trim = a.trim

    if not ident:
        skip("coverage", "no identity map (--identity), so no membership blocks")
    elif not depth_paths:
        skip("coverage", "no depth files (give --config, --depth or --depth-dir)")
    else:
        for s, m in depth_paths.items():
            for c in [c for c in order if c in m] + [c for c in m if c not in order]:
                fid = "coverage_%s_%s" % (s, c)
                if c not in ident.get("constructs", {}):
                    skip(fid, "%s is not in the identity map" % c)
                    continue
                if not os.path.exists(m[c]):
                    skip(fid, "depth file %s not found" % m[c])
                    continue
                want = [x for x in (contig_names.get(c), c) if x]
                depth, got = read_depth(m[c], want)
                if depth is None:
                    skip(fid, got)
                    continue
                try:
                    svg, cap = fig_coverage(fid, s, c, ident, roles, depth,
                                            baselines.get(s), trim, policy)
                except ValueError as e:
                    skip(fid, str(e))
                    continue
                emit(fid, "coverage", svg, cap, sample=s, construct=c)

    # ---- 2. genome overview
    if not a.fai:
        skip("genome_overview", "no host .fai (--fai, or run host_reference.py)")
    elif not os.path.exists(a.fai):
        skip("genome_overview", ".fai %s not found" % a.fai)
    else:
        fai = read_fai(a.fai)
        if not fai:
            skip("genome_overview", ".fai %s has no contigs" % a.fai)
        else:
            if called is None:
                print("  note: no called_sites.tsv; drawing chromosomes without loci")
            exclude = list(order) + list(contig_names.values())
            svg, cap, notes = fig_genome(fai, called or [], order, roles, exclude,
                                         a.min_contig_frac, a.max_contigs, summary)
            for n in notes:
                print("  note: %s" % n)
            emit("genome_overview", "genome_overview", svg, cap)

    # ---- 3. copy number
    if not cns:
        if not any(s["id"] == "copy_number" for s in skipped):
            skip("copy_number", "no copy-number JSON (--copy-number)")
    else:
        try:
            svg, cap = fig_copy_number(cns, order, roles)
            emit("copy_number", "copy_number", svg, cap,
                 runs=[l for l, _ in cns])
        except ValueError as e:
            skip("copy_number", str(e))

    # ---- 4. architecture
    if not summary:
        skip("architecture", "no reconcile_sites.py summary (--site-summary)")
    else:
        try:
            svg, cap = fig_architecture(summary, order, roles, cns[0][1] if cns else None)
            emit("architecture", "architecture", svg, cap)
        except ValueError as e:
            skip("architecture", str(e))

    manifest = {"figures": figures, "skipped": skipped,
                "constructs": order, "trim": {"bp": trim, "policy": policy}}
    with open(os.path.join(a.out, "figures.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")
    print("figures: %d written, %d skipped -> %s"
          % (len(figures), len(skipped), os.path.join(a.out, "figures.json")))

    if a.fill:
        dst = a.fill_out or a.fill
        try:
            n, missing = fill_report(a.fill, dst, manifest, a.out)
        except OSError as e:
            print("error: cannot fill %s (%s)" % (a.fill, e), file=sys.stderr)
            return 2
        print("report: %d figure link(s) placed in %s%s"
              % (n, dst, ("; not generated: " + ", ".join(missing)) if missing else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

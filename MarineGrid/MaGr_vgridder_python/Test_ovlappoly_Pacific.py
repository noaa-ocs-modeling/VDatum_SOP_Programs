#!/usr/bin/env python
# coding: utf-8

import math
from dataclasses import dataclass
from pathlib import Path
from itertools import combinations

import shapefile
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon


# ============================================================
# USER SETTINGS
# ============================================================

INPUT_SHP = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/MarinGrid/MarinGrid_polygons_3832_v6.shp"
OUTDIR = Path("maringrid_overlap_report_v6")
ALLOWED_PARENT_CHILD = set()

# Meter-based CRS usually requires WRAP_TO_360 = False
WRAP_TO_360 = False
MAKE_PLOTS = True


# ============================================================
# CONSTANTS (Optimized for Meters)
# ============================================================

SMALL_COMMON = 0.1   # 10cm tolerance
DEL_PERT = 0.05      
RAD = 0.01745329
SEG_TOL = 0.1        # 10cm boundary tolerance


# ============================================================
# DATA STRUCTURES
# ============================================================

@dataclass
class PolyData:
    name: str
    x: list
    y: list
    record_index: int
    part_index: int

    @classmethod
    def from_points(cls, name, points, record_index, part_index, wrap360=False):
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        if wrap360:
            xs = [x if x >= 0 else x + 360.0 for x in xs]
        return cls(name, xs, ys, record_index, part_index)

    def validate_like_fortran(self):
        if len(self.x) < 4: return
        if not (self.x[0] == self.x[-1] and self.y[0] == self.y[-1]): return


# ============================================================
# GEOMETRY HELPERS
# ============================================================

def part_points(shape):
    pts = shape.points
    idx = list(shape.parts) + [len(pts)]
    return [pts[a:b] for a, b in zip(idx[:-1], idx[1:])]

def point_common_dist(xi, yi, xt, yt):
    return math.sqrt((xi - xt) ** 2 + (yi - yt) ** 2)

def inside_pert_fortran(xi, yi, x, y):
    ys = [val + DEL_PERT if yi == val else val for val in y]
    nc = 0
    for n in range(len(x) - 1):
        if min(ys[n], ys[n + 1]) < yi < max(ys[n], ys[n + 1]):
            xc = x[n] if x[n] == x[n + 1] else x[n] + (x[n + 1] - x[n]) * (yi - ys[n]) / (ys[n + 1] - ys[n])
            if xc > xi:
                nc += 1
    return nc % 2

def point_on_segment(px, py, x1, y1, x2, y2, tol=SEG_TOL):
    if px < min(x1, x2) - tol or px > max(x1, x2) + tol: return False
    if py < min(y1, y2) - tol or py > max(y1, y2) + tol: return False
    dx, dy = x2 - x1, y2 - y1
    if abs(dx) < 1e-12 and abs(dy) < 1e-12: return math.sqrt((px-x1)**2 + (py-y1)**2) < tol
    cross = abs((px - x1) * dy - (py - y1) * dx)
    dist = math.sqrt(dx**2 + dy**2)
    return (cross / dist) < tol

def mark_points_on_any_border(skip, pxs, pys, oxs, oys, tol=SEG_TOL):
    for i in range(len(pxs) - 1):
        if skip[i]: continue
        for j in range(len(oxs) - 1):
            if point_on_segment(pxs[i], pys[i], oxs[j], oys[j], oxs[j+1], oys[j+1], tol=tol):
                skip[i] = 1
                break


# ============================================================
# ANALYSIS & REPORTING
# ============================================================

def analyze_overlaps(record_infos, part_polys):
    results = []
    by_record = {info['record_index']: [] for info in record_infos}
    for poly in part_polys:
        by_record[poly.record_index].append(poly)

    for i, j in combinations(range(len(record_infos)), 2):
        name_i, name_j = record_infos[i]["MATNAME"], record_infos[j]["MATNAME"]
        parts_i, parts_j = by_record[i], by_record[j]
        pair_hits = []

        for p1 in parts_i:
            for p2 in parts_j:
                ni, nt = len(p1.x), len(p2.x)
                skip1, skip2 = [0] * ni, [0] * nt

                # Boundary snapping across all parts
                for op2 in parts_j: mark_points_on_any_border(skip1, p1.x, p1.y, op2.x, op2.y)
                for op1 in parts_i: mark_points_on_any_border(skip2, p2.x, p2.y, op1.x, op1.y)

                nt1 = sum(1 for idx in range(ni-1) if not skip1[idx] and inside_pert_fortran(p1.x[idx], p1.y[idx], p2.x, p2.y))
                nt2 = sum(1 for idx in range(nt-1) if not skip2[idx] and inside_pert_fortran(p2.x[idx], p2.y[idx], p1.x, p1.y))

                s1, s2 = sum(skip1), sum(skip2)
                
                # Only record if there is an actual interaction
                if (nt1 + nt2) > 0 or (s1 + s2) > 0:
                    cls = "real_overlap" if (nt1 + nt2) > 0 else "touch_only"
                    pair_hits.append({
                        "part_i": p1.part_index, "part_j": p2.part_index, 
                        "class": cls, "ntot1": nt1, "ntot2": nt2, "skip1": s1, "skip2": s2
                    })

        pair_class = "no_overlap"
        if any(h["class"] == "real_overlap" for h in pair_hits): pair_class = "real_overlap"
        elif pair_hits: pair_class = "touch_only"

        results.append({"record_i": i, "name_i": name_i, "record_j": j, "name_j": name_j, "pair_class": pair_class, "details": pair_hits})
    return results


def write_summary_txt(results, out_txt):
    with open(out_txt, "w", encoding="utf-8") as f:
        for cls in ["real_overlap", "touch_only", "no_overlap"]:
            f.write(f"{cls.upper()}\n" + "-" * len(cls) + "\n")
            rows = [r for r in results if r["pair_class"] == cls]
            if not rows: f.write("None\n\n"); continue
            for r in rows:
                f.write(f'{r["name_i"]} vs {r["name_j"]} -> {r["pair_class"]}\n')
                for d in r["details"]:
                    f.write(f'    part {d["part_i"]} vs part {d["part_j"]}: ntot1={d["ntot1"]}, ntot2={d["ntot2"]}, skip1={d["skip1"]}, skip2={d["skip2"]}\n')
            f.write("\n")

# ============================================================
# PLOTTING UTILS
# ============================================================

def plot_overview(record_infos, results, out_png):
    fig, ax = plt.subplots(figsize=(16, 9))
    allx, ally = [], []
    for info in record_infos:
        parts = [p for p in part_points(info["shape"])]
        for pts in parts:
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            allx.extend(xs); ally.extend(ys)
            ax.add_patch(MplPolygon(pts, closed=True, fill=True, alpha=0.15, edgecolor="black"))
        ax.text(info["center"][0], info["center"][1], info["MATNAME"], fontsize=8, ha='center')
    
    ax.set_xlim(min(allx), max(allx)); ax.set_ylim(min(ally), max(ally))
    ax.set_aspect("equal"); fig.savefig(out_png, dpi=150); plt.close(fig)

def main():
    OUTDIR.mkdir(exist_ok=True, parents=True)
    sf = shapefile.Reader(INPUT_SHP)
    
    # Load data
    record_infos, part_polys = [], []
    for idx, (rec, shp) in enumerate(zip(sf.records(), sf.shapes())):
        name = rec["MATNAME"]
        pts = [p for part in part_points(shp) for p in part]
        record_infos.append({"record_index": idx, "MATNAME": name, "shape": shp, "center": (sum(p[0] for p in pts)/len(pts), sum(p[1] for p in pts)/len(pts))})
        for p_idx, p_pts in enumerate(part_points(shp), 1):
            part_polys.append(PolyData.from_points(name, p_pts, idx, p_idx))

    results = analyze_overlaps(record_infos, part_polys)
    write_summary_txt(results, OUTDIR / "pair_summary.txt")
    if MAKE_PLOTS: plot_overview(record_infos, results, OUTDIR / "overview.png")
    print(f"Complete. Reports in {OUTDIR.resolve()}")

if __name__ == "__main__":
    main()






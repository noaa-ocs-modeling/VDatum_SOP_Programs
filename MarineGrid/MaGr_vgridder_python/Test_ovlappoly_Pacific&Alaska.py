#!/usr/bin/env python
# coding: utf-8

# "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/MarinGrid/MarinGrid_polygons_4326_v3.shp"
# "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/MarinGrid/Alaska/AK_BPs_V16.shp"

import math
from dataclasses import dataclass
from pathlib import Path

import shapefile
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon


# ============================================================
# USER SETTINGS
# ============================================================

MARINGRID_SHP = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/MarinGrid/MarinGrid_polygons_3832_v3.shp"
ALASKA_SHP = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/MarinGrid/Alaska/AK_BPs_V16.shp"

OUTDIR = Path("/mnt/c/Users/mojgan.rostaminia/Documents/bin/maringrid_vs_alaska_report")

WRAP_TO_360 = True
MAKE_PLOTS = True


# ============================================================
# CONSTANTS
# ============================================================

SMALL_COMMON = 5e-6
DEL_PERT = 2e-6
RAD = 0.01745329
SEG_TOL = 1e-8


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
        if len(self.x) < 4:
            raise ValueError(f"{self.name}: polygon has too few points")
        if not (self.x[0] == self.x[-1] and self.y[0] == self.y[-1]):
            raise ValueError(f"{self.name}: polygon is not closed")
        for i in range(len(self.x) - 1):
            if self.x[i] == self.x[i + 1] and self.y[i] == self.y[i + 1]:
                raise ValueError(f"{self.name}: repeated consecutive point at index {i}")


# ============================================================
# GEOMETRY HELPERS
# ============================================================

def wrap360_pts(pts):
    if not WRAP_TO_360:
        return pts
    return [(x if x >= 0 else x + 360.0, y) for x, y in pts]


def part_points(shape):
    pts = shape.points
    idx = list(shape.parts) + [len(pts)]
    return [pts[a:b] for a, b in zip(idx[:-1], idx[1:])]


def point_common_dist(xi, yi, xt, yt):
    c = math.cos(RAD * yi) ** 2
    return math.sqrt(c * (xi - xt) ** 2 + (yi - yt) ** 2)


def inside_pert_fortran(xi, yi, x, y):
    ys = list(y)

    for n in range(len(ys)):
        if yi == ys[n]:
            ys[n] = ys[n] + DEL_PERT

    nc = 0
    for n in range(len(x) - 1):
        if yi > min(ys[n], ys[n + 1]) and yi < max(ys[n], ys[n + 1]):
            if x[n] == x[n + 1]:
                xc = x[n]
            else:
                xc = x[n] + (x[n + 1] - x[n]) * (yi - ys[n]) / (ys[n + 1] - ys[n])
            if xc > xi:
                nc += 1

    return nc % 2


def point_on_segment(px, py, x1, y1, x2, y2, tol=SEG_TOL):
    if px < min(x1, x2) - tol or px > max(x1, x2) + tol:
        return False
    if py < min(y1, y2) - tol or py > max(y1, y2) + tol:
        return False

    dx = x2 - x1
    dy = y2 - y1

    if abs(dx) < tol and abs(dy) < tol:
        return abs(px - x1) < tol and abs(py - y1) < tol

    cross = (px - x1) * dy - (py - y1) * dx
    if abs(cross) > tol:
        return False

    return True


def mark_points_on_any_border(skip, pxs, pys, oxs, oys, tol=SEG_TOL):
    npts = len(pxs)
    nother = len(oxs)

    for i in range(npts - 1):
        if skip[i]:
            continue
        for j in range(nother - 1):
            if point_on_segment(pxs[i], pys[i], oxs[j], oys[j], oxs[j + 1], oys[j + 1], tol=tol):
                skip[i] = 1
                break


# ============================================================
# OVERLAP TEST
# ============================================================

def ov_test_fortran_equiv(p1: PolyData, p2: PolyData):
    p1.validate_like_fortran()
    p2.validate_like_fortran()

    ni = len(p1.x)
    nt = len(p2.x)

    skip1 = [0] * ni
    skip2 = [0] * nt

    # Nearly common points
    for n in range(ni - 1):
        for i in range(1, nt):
            d = point_common_dist(p1.x[n], p1.y[n], p2.x[i], p2.y[i])
            if d <= SMALL_COMMON:
                skip1[n] = 1
                skip2[i] = 1

    # Horizontal/vertical border checks from original Fortran logic
    for n in range(ni - 1):
        if skip1[n]:
            continue
        for m in range(nt - 1):
            if p1.y[n] == p2.y[m]:
                mp = m + 1
                if p1.y[n] == p2.y[mp]:
                    ibet = (
                        (p1.x[n] >= p2.x[m] and p1.x[n] < p2.x[mp]) or
                        (p1.x[n] < p2.x[m] and p1.x[n] >= p2.x[mp])
                    )
                    if ibet:
                        skip1[n] = 1

    for n in range(ni - 1):
        if skip1[n]:
            continue
        for m in range(nt - 1):
            if p1.x[n] == p2.x[m]:
                mp = m + 1
                if p1.x[n] == p2.x[mp]:
                    ibet = (
                        (p1.y[n] >= p2.y[m] and p1.y[n] < p2.y[mp]) or
                        (p1.y[n] < p2.y[m] and p1.y[n] >= p2.y[mp])
                    )
                    if ibet:
                        skip1[n] = 1

    for m in range(nt - 1):
        if skip2[m]:
            continue
        for n in range(ni - 1):
            if p1.y[n] == p2.y[m]:
                np1 = n + 1
                if p1.y[np1] == p2.y[m]:
                    ibet = (
                        (p2.x[m] >= p1.x[n] and p2.x[m] < p1.x[np1]) or
                        (p2.x[m] < p1.x[n] and p2.x[m] >= p1.x[np1])
                    )
                    if ibet:
                        skip2[m] = 1

    for m in range(nt - 1):
        if skip2[m]:
            continue
        for n in range(ni - 1):
            if p1.x[n] == p2.x[m]:
                np1 = n + 1
                if p1.x[np1] == p2.x[m]:
                    ibet = (
                        (p2.y[m] >= p1.y[n] and p2.y[m] < p1.y[np1]) or
                        (p2.y[m] < p1.y[n] and p2.y[m] >= p1.y[np1])
                    )
                    if ibet:
                        skip2[m] = 1

    # Extra general border check for slanted shared boundaries
    mark_points_on_any_border(skip1, p1.x, p1.y, p2.x, p2.y, tol=SEG_TOL)
    mark_points_on_any_border(skip2, p2.x, p2.y, p1.x, p1.y, tol=SEG_TOL)

    ntot1 = 0
    ntot2 = 0
    inside_pts1 = []
    inside_pts2 = []

    for n in range(ni - 1):
        if not skip1[n]:
            inside = inside_pert_fortran(p1.x[n], p1.y[n], p2.x, p2.y)
            ntot1 += inside
            if inside:
                inside_pts1.append((n, p1.x[n], p1.y[n]))

    for n in range(nt - 1):
        if not skip2[n]:
            inside = inside_pert_fortran(p2.x[n], p2.y[n], p1.x, p1.y)
            ntot2 += inside
            if inside:
                inside_pts2.append((n, p2.x[n], p2.y[n]))

    overlap = (ntot1 + ntot2) > 0

    return {
        "overlap": overlap,
        "ntot1": ntot1,
        "ntot2": ntot2,
        "skip1": sum(skip1),
        "skip2": sum(skip2),
        "inside_pts1": inside_pts1,
        "inside_pts2": inside_pts2,
    }


# ============================================================
# LOAD SHAPEFILES
# ============================================================

def load_parts(sf, source_name):
    records = sf.records()
    shapes = sf.shapes()

    record_infos = []
    part_polys = []

    for rec_idx, (rec, shp) in enumerate(zip(records, shapes)):
        rec_dict = rec.as_dict()

        if "MATNAME" in rec_dict and rec_dict["MATNAME"] not in (None, ""):
            name = rec_dict["MATNAME"]
        else:
            name = f"{source_name}_{rec_idx}"

        parts = part_points(shp)
        valid_parts = []
        valid_part_count = 0

        for part_idx, pts in enumerate(parts, start=1):
            if pts is None or len(pts) < 4:
                print(f"skipping short part: {name}, part {part_idx}, npts={0 if pts is None else len(pts)}")
                continue

            if pts[0] != pts[-1]:
                pts = list(pts) + [pts[0]]

            if len(pts) < 4:
                print(f"skipping invalid closed part: {name}, part {part_idx}, npts={len(pts)}")
                continue

            valid_parts.append(pts)

            poly = PolyData.from_points(
                name=name,
                points=pts,
                record_index=rec_idx,
                part_index=part_idx,
                wrap360=WRAP_TO_360,
            )
            part_polys.append(poly)
            valid_part_count += 1

        if valid_part_count == 0:
            print(f"skipping record with no valid polygon parts: {name}")
            continue

        wrapped_parts = [wrap360_pts(p) for p in valid_parts]
        xs = [x for part in wrapped_parts for x, _ in part]
        ys = [y for part in wrapped_parts for _, y in part]
        center = (sum(xs) / len(xs), sum(ys) / len(ys))

        record_infos.append({
            "record_index": rec_idx,
            "MATNAME": name,
            "shape": shp,
            "center": center,
            "source_name": source_name,
        })

    return record_infos, part_polys


# ============================================================
# CROSS-SHAPEFILE ANALYSIS
# ============================================================

def analyze_cross_overlaps(record_infos_a, part_polys_a, record_infos_b, part_polys_b):
    results = []

    by_record_a = {}
    for poly in part_polys_a:
        by_record_a.setdefault(poly.record_index, []).append(poly)

    by_record_b = {}
    for poly in part_polys_b:
        by_record_b.setdefault(poly.record_index, []).append(poly)

    for info_a in record_infos_a:
        for info_b in record_infos_b:
            rec_a = info_a["record_index"]
            rec_b = info_b["record_index"]
            name_a = info_a["MATNAME"]
            name_b = info_b["MATNAME"]

            pair_hits = []

            for p1 in by_record_a.get(rec_a, []):
                for p2 in by_record_b.get(rec_b, []):
                    ov = ov_test_fortran_equiv(p1, p2)

                    if ov["overlap"]:
                        cls = "real_overlap"
                    elif ov["skip1"] > 0 or ov["skip2"] > 0:
                        cls = "touch_only"
                    else:
                        cls = "no_overlap"

                    if cls != "no_overlap":
                        pair_hits.append({
                            "part_a": p1.part_index,
                            "part_b": p2.part_index,
                            "class": cls,
                            "ntot1": ov["ntot1"],
                            "ntot2": ov["ntot2"],
                            "skip1": ov["skip1"],
                            "skip2": ov["skip2"],
                            "inside_pts1": ov.get("inside_pts1", []),
                            "inside_pts2": ov.get("inside_pts2", []),
                        })

            if pair_hits:
                pair_class = "touch_only"
                if any(h["class"] == "real_overlap" for h in pair_hits):
                    pair_class = "real_overlap"

                results.append({
                    "record_a": rec_a,
                    "name_a": name_a,
                    "record_b": rec_b,
                    "name_b": name_b,
                    "pair_class": pair_class,
                    "details": pair_hits,
                })
            else:
                results.append({
                    "record_a": rec_a,
                    "name_a": name_a,
                    "record_b": rec_b,
                    "name_b": name_b,
                    "pair_class": "no_overlap",
                    "details": [],
                })

    return results


# ============================================================
# WRITERS
# ============================================================

def write_cross_csv(results, out_csv):
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("record_a,name_a,record_b,name_b,pair_class,details\n")
        for r in results:
            safe_details = str(r["details"]).replace('"', "'")
            f.write(
                f'{r["record_a"]},{r["name_a"]},'
                f'{r["record_b"]},{r["name_b"]},'
                f'{r["pair_class"]},"{safe_details}"\n'
            )


def write_cross_summary_txt(results, out_txt, title="MARINGRID vs ALASKA CHECK"):
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(title + "\n")
        f.write("=" * len(title) + "\n\n")

        for cls in ["real_overlap", "touch_only", "no_overlap"]:
            f.write(f"{cls.upper()}\n")
            f.write("-" * len(cls) + "\n")
            rows = [r for r in results if r["pair_class"] == cls]
            if not rows:
                f.write("None\n\n")
                continue

            for r in rows:
                f.write(f'{r["name_a"]} vs {r["name_b"]} -> {r["pair_class"]}\n')
                for d in r["details"]:
                    f.write(
                        f'    part {d["part_a"]} vs part {d["part_b"]}: '
                        f'class={d["class"]}, '
                        f'ntot1={d["ntot1"]}, ntot2={d["ntot2"]}, '
                        f'skip1={d["skip1"]}, skip2={d["skip2"]}, '
                        f'inside_pts1={d.get("inside_pts1", [])}, '
                        f'inside_pts2={d.get("inside_pts2", [])}\n'
                    )
            f.write("\n")


# ============================================================
# PLOTTING
# ============================================================

def get_record_parts(record_infos):
    out = {}
    for info in record_infos:
        shp = info["shape"]
        parts = [wrap360_pts(p) for p in part_points(shp)]
        out[info["record_index"]] = parts
    return out


def plot_cross_overview(record_infos_a, record_infos_b, results, out_png):
    parts_a = get_record_parts(record_infos_a)
    parts_b = get_record_parts(record_infos_b)

    real_pairs = [f'{r["name_a"]}-{r["name_b"]}' for r in results if r["pair_class"] == "real_overlap"]
    touch_pairs = [f'{r["name_a"]}-{r["name_b"]}' for r in results if r["pair_class"] == "touch_only"]

    fig, ax = plt.subplots(figsize=(16, 9))
    allx = []
    ally = []

    for info in record_infos_a:
        for pts in parts_a[info["record_index"]]:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            allx.extend(xs)
            ally.extend(ys)
            ax.add_patch(MplPolygon(pts, closed=True, fill=True, alpha=0.22, edgecolor="black", linewidth=1.0))
        cx, cy = info["center"]
        ax.text(cx, cy, info["MATNAME"], fontsize=8,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.0))

    for info in record_infos_b:
        for pts in parts_b[info["record_index"]]:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            allx.extend(xs)
            ally.extend(ys)
            ax.add_patch(MplPolygon(pts, closed=True, fill=True, alpha=0.18, edgecolor="red", linewidth=1.0))
        cx, cy = info["center"]
        ax.text(cx, cy, info["MATNAME"], fontsize=8,
                bbox=dict(facecolor="mistyrose", alpha=0.85, edgecolor="none", pad=1.0))

    txt = []
    txt.append("Real overlaps: " + (", ".join(real_pairs) if real_pairs else "None"))
    txt.append("Touch only: " + (", ".join(touch_pairs) if touch_pairs else "None"))

    ax.text(
        0.01, 0.99,
        "\n".join(txt),
        transform=ax.transAxes,
        va="top",
        bbox=dict(facecolor="white", alpha=0.92),
    )

    xmin, xmax = min(allx), max(allx)
    ymin, ymax = min(ally), max(ally)
    padx = 0.03 * (xmax - xmin)
    pady = 0.03 * (ymax - ymin)

    ax.set_xlim(xmin - padx, xmax + padx)
    ax.set_ylim(ymin - pady, ymax + pady)

    ax.set_title("MarinGrid vs Alaska overview")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, alpha=0.25)
    ax.set_aspect("equal", adjustable="box")

    fig.savefig(str(out_png), dpi=100)
    print("wrote:", Path(out_png).resolve())
    plt.close(fig)


def plot_cross_real_overlap_pairs(record_infos_a, record_infos_b, results, outdir):
    parts_a = get_record_parts(record_infos_a)
    parts_b = get_record_parts(record_infos_b)

    for r in results:
        if r["pair_class"] != "real_overlap":
            continue

        info_a = next(x for x in record_infos_a if x["record_index"] == r["record_a"])
        info_b = next(x for x in record_infos_b if x["record_index"] == r["record_b"])

        fig, ax = plt.subplots(figsize=(10, 8))

        for pts in parts_a[info_a["record_index"]]:
            ax.add_patch(MplPolygon(pts, closed=True, fill=True, alpha=0.30, edgecolor="black", linewidth=1.2))
            ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", markersize=2)

        for pts in parts_b[info_b["record_index"]]:
            ax.add_patch(MplPolygon(pts, closed=True, fill=True, alpha=0.25, edgecolor="red", linewidth=1.2))
            ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="s", markersize=2)

        allx = [x for pts in parts_a[info_a["record_index"]] + parts_b[info_b["record_index"]] for x, _ in pts]
        ally = [y for pts in parts_a[info_a["record_index"]] + parts_b[info_b["record_index"]] for _, y in pts]
        padx = max(1.0, 0.08 * (max(allx) - min(allx)))
        pady = max(1.0, 0.08 * (max(ally) - min(ally)))

        ax.set_xlim(min(allx) - padx, max(allx) + padx)
        ax.set_ylim(min(ally) - pady, max(ally) + pady)

        lines = [f'{r["name_a"]} vs {r["name_b"]}']
        for d in r["details"]:
            lines.append(
                f'part {d["part_a"]} vs part {d["part_b"]}: '
                f'ntot1={d["ntot1"]}, ntot2={d["ntot2"]}, '
                f'skip1={d["skip1"]}, skip2={d["skip2"]}'
            )
            if d.get("inside_pts1"):
                lines.append(f'inside_pts1={d["inside_pts1"]}')
            if d.get("inside_pts2"):
                lines.append(f'inside_pts2={d["inside_pts2"]}')

        ax.text(
            0.01, 0.99,
            "\n".join(lines),
            transform=ax.transAxes,
            va="top",
            bbox=dict(facecolor="white", alpha=0.92),
        )

        ax.set_title(f'Real overlap: {r["name_a"]} vs {r["name_b"]}')
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.25)
        ax.set_aspect("equal", adjustable="box")

        out_png = outdir / f'maringrid_vs_alaska_real_overlap_{r["name_a"]}__{r["name_b"]}.png'
        fig.savefig(str(out_png), dpi=100)
        print("wrote:", Path(out_png).resolve())
        plt.close(fig)


# ============================================================
# MAIN
# ============================================================

def main():
    OUTDIR.mkdir(exist_ok=True, parents=True)

    sf_mg = shapefile.Reader(MARINGRID_SHP)
    mg_infos, mg_parts = load_parts(sf_mg, "MarinGrid")

    sf_ak = shapefile.Reader(ALASKA_SHP)
    ak_infos, ak_parts = load_parts(sf_ak, "Alaska")

    results = analyze_cross_overlaps(mg_infos, mg_parts, ak_infos, ak_parts)

    print("\n=== MARINGRID vs ALASKA REAL OVERLAPS ===")
    real_rows = [r for r in results if r["pair_class"] == "real_overlap"]
    if real_rows:
        for r in real_rows:
            print(f'{r["name_a"]} vs {r["name_b"]}')
    else:
        print("None")

    print("\n=== MARINGRID vs ALASKA TOUCH ONLY ===")
    touch_rows = [r for r in results if r["pair_class"] == "touch_only"]
    if touch_rows:
        for r in touch_rows:
            print(f'{r["name_a"]} vs {r["name_b"]}')
    else:
        print("None")

    write_cross_csv(results, OUTDIR / "maringrid_vs_alaska.csv")
    write_cross_summary_txt(results, OUTDIR / "maringrid_vs_alaska_summary.txt")

    if MAKE_PLOTS:
        plot_cross_overview(mg_infos, ak_infos, results, OUTDIR / "maringrid_vs_alaska_overview.png")
        plot_cross_real_overlap_pairs(mg_infos, ak_infos, results, OUTDIR)

    print(f"\nDone. Outputs saved in: {OUTDIR.resolve()}")


if __name__ == "__main__":
    main()





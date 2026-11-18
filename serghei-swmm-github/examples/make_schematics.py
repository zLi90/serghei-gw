"""Generate a side-view schematic for each drainage demo case."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrow, Polygon, Rectangle
import os


GROUND = 9.0


def new_ax(title):
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_xlim(0, 10)
    ax.set_ylim(3.5, 10.6)
    ax.axis("off")
    return fig, ax


def ground(ax):
    ax.plot([0, 10], [GROUND, GROUND], color="saddlebrown", lw=2)
    ax.add_patch(Rectangle((0, 3.5), 10, GROUND - 3.5, color="tan", alpha=0.35))
    ax.text(9.6, GROUND + 0.08, "ground z=9.0", ha="right", fontsize=8, color="saddlebrown")


def rain(ax, xmin=0.5, xmax=9.5):
    for x in (xmin, (xmin + xmax) / 2, xmax):
        ax.add_patch(FancyArrow(x, 10.35, 0, -0.5, width=0.04,
                                head_width=0.16, head_length=0.16,
                                length_includes_head=True, color="steelblue"))
    ax.text((xmin + xmax) / 2, 10.52, "rainfall -> surface runoff",
            ha="center", fontsize=8, color="steelblue")


def manhole(ax, x, z_bottom, z_top, label, water=0.0, open_lid=True):
    w = 0.55
    ax.add_patch(Rectangle((x - w / 2, z_bottom), w, z_top - z_bottom,
                           fill=False, edgecolor="black", lw=1.5))
    if water > 0:
        ax.add_patch(Rectangle((x - w / 2, z_bottom), w, min(water, z_top - z_bottom),
                               color="deepskyblue", alpha=0.65))
    if open_lid:
        ax.plot([x - w / 2 - 0.08, x + w / 2 + 0.08], [z_top, z_top],
                color="black", lw=2)
    ax.text(x, z_bottom - 0.25, label, ha="center", fontsize=8)
    ax.text(x, z_bottom + 0.12, f"z={z_bottom}", ha="center", fontsize=6.5, color="dimgray")


def storage_tank(ax, x, z_bottom, z_top, label, water=0.0):
    w = 1.4
    ax.add_patch(Rectangle((x - w / 2, z_bottom), w, z_top - z_bottom,
                           fill=False, edgecolor="black", lw=1.5))
    if water > 0:
        ax.add_patch(Rectangle((x - w / 2, z_bottom), w, min(water, z_top - z_bottom),
                               color="deepskyblue", alpha=0.65))
    ax.plot([x - w / 2 - 0.08, x + w / 2 + 0.08], [z_top, z_top], color="black", lw=2)
    ax.text(x, z_bottom - 0.25, label, ha="center", fontsize=8)
    ax.text(x, z_bottom + 0.12, f"z={z_bottom}", ha="center", fontsize=6.5, color="dimgray")


def pipe(ax, x1, z1, x2, z2, label=None, mid=None):
    ax.plot([x1, x2], [z1, z2], color="dimgray", lw=3, solid_capstyle="round")
    if label:
        mx, my = mid or ((x1 + x2) / 2, (z1 + z2) / 2)
        ax.text(mx, my - 0.25, label, ha="center", fontsize=8, color="dimgray")


def pump(ax, x, z, label="pump", direction=1):
    ax.add_patch(Circle((x, z), 0.28, fill=True, facecolor="gold", edgecolor="black", lw=1.5))
    ax.add_patch(Polygon([(x - 0.13, z - 0.12), (x + 0.16, z), (x - 0.13, z + 0.12)],
                         color="black"))
    ax.text(x, z - 0.5, label, ha="center", fontsize=8)


def river(ax, x, z_bottom, stage, label="river"):
    ax.add_patch(Rectangle((x, z_bottom), 10 - x, stage - z_bottom,
                           color="deepskyblue", alpha=0.6))
    ax.plot([x, 10], [stage, stage], color="steelblue", lw=1.5, ls="--")
    ax.text(x + 0.4, stage + 0.08, f"{label} stage={stage:.1f}", fontsize=7.5, color="steelblue")


def flow_arrow(ax, x1, z1, x2, z2, color="crimson"):
    ax.annotate("", xy=(x2, z2), xytext=(x1, z1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8))


def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def case01(out):
    fig, ax = new_ax("Case 1 - Gravity drainage, no pump (FREE outfall)")
    ground(ax)
    rain(ax)
    manhole(ax, 1.5, 7.0, 9.0, "J0 (sealed=0)", water=0.65)
    manhole(ax, 4.0, 7.0, 9.0, "J1 (sealed=1)", water=0.3)
    pipe(ax, 1.5, 7.0, 4.0, 7.0, "C0")
    pipe(ax, 4.0, 7.0, 6.5, 6.0, "C1")
    ax.plot([6.5, 10], [6.0, 6.0], color="black", lw=2)
    river(ax, 6.5, 5.7, 6.0, "FREE")
    ax.text(6.5, 5.35, "O2 (outlet, z=6.0)", ha="left", fontsize=8)
    flow_arrow(ax, 0.8, 8.5, 1.5, 8.5)
    flow_arrow(ax, 2.8, 6.6, 2.8, 7.0, color="royalblue")
    flow_arrow(ax, 5.2, 6.1, 5.2, 6.5, color="royalblue")
    ax.text(5.0, 10.0, "surface(2D) -> manhole -> pipes -> free outfall",
            ha="center", fontsize=9, style="italic")
    save(fig, out)


def case02(out):
    fig, ax = new_ax("Case 2 - Wet well + TYPE3 pump, FIXED outfall stage")
    ground(ax)
    rain(ax)
    manhole(ax, 1.2, 7.0, 9.0, "J1 (sealed=0)", water=0.5)
    pipe(ax, 1.2, 7.0, 3.4, 5.0, "C0")
    storage_tank(ax, 3.4, 5.0, 9.0, "Wet well S0\n25 m2, z=5.0", water=1.3)
    pump(ax, 5.3, 7.2, "P0 (TYPE3)")
    pipe(ax, 5.3, 5.0, 5.3, 7.2, label=None)
    pipe(ax, 5.3, 7.2, 6.8, 6.5, label=None)
    ax.plot([6.8, 10], [6.5, 6.5], color="black", lw=2)
    river(ax, 6.8, 6.2, 6.5, "FIXED")
    ax.text(6.8, 5.85, "O2 (outlet, z=6.5)", ha="left", fontsize=8)
    flow_arrow(ax, 0.6, 8.4, 1.2, 8.4)
    flow_arrow(ax, 2.3, 6.4, 2.3, 6.9, color="royalblue")
    flow_arrow(ax, 3.4, 4.3, 5.3, 4.3, color="crimson")
    ax.text(5.0, 10.0, "runoff -> manhole -> wet well -> pump -> fixed-stage outfall",
            ha="center", fontsize=9, style="italic")
    save(fig, out)


def case03(out):
    fig, ax = new_ax("Case 3 - Wet well + TYPE3 pump, time-varying river stage")
    ground(ax)
    rain(ax)
    manhole(ax, 1.2, 7.0, 9.0, "J1 (sealed=0)", water=0.5)
    pipe(ax, 1.2, 7.0, 3.4, 5.0, "C0")
    storage_tank(ax, 3.4, 5.0, 9.0, "Wet well S0", water=1.1)
    pump(ax, 5.3, 7.2, "P0 (TYPE3)")
    pipe(ax, 5.3, 5.0, 5.3, 7.2)
    pipe(ax, 5.3, 7.2, 6.8, 6.5)
    river(ax, 6.8, 6.2, 8.5, "high tide")
    river(ax, 7.6, 6.2, 6.5, "low tide")
    ax.plot([6.8, 10], [6.5, 6.5], color="black", lw=2)
    ax.text(6.8, 5.85, "O2 (outlet, z=6.5)", ha="left", fontsize=8)
    flow_arrow(ax, 0.6, 8.4, 1.2, 8.4)
    flow_arrow(ax, 3.4, 4.3, 5.3, 4.3, color="crimson")
    ax.text(5.0, 10.0, "rising river raises head -> pump flow drops",
            ha="center", fontsize=9, style="italic")
    save(fig, out)


def case04(out):
    fig, ax = new_ax("Case 4 - IDEAL pump + time-varying river + flap gate")
    ground(ax)
    rain(ax)
    manhole(ax, 1.3, 7.0, 9.0, "J0 (sealed=0)", water=0.8)
    pump(ax, 2.6, 7.6, "P0 (IDEAL)")
    pipe(ax, 1.3, 7.0, 2.6, 7.6)
    pipe(ax, 2.6, 7.6, 3.7, 7.0)
    manhole(ax, 3.7, 7.0, 9.0, "J1 (sealed=1)", water=0.45)
    pipe(ax, 3.7, 7.0, 6.3, 6.0, "C0")
    river(ax, 6.3, 5.7, 8.5, "high tide")
    ax.plot([6.3, 10], [6.0, 6.0], color="black", lw=2)
    ax.text(6.3, 5.3, "O2 (outlet z=6.0, gated=yes)", ha="left", fontsize=8)
    ax.text(6.4, 8.8, "flap gate blocks backflow", fontsize=8, color="purple")
    flow_arrow(ax, 0.7, 8.4, 1.3, 8.4)
    flow_arrow(ax, 1.9, 7.6, 2.6, 7.6, color="crimson")
    flow_arrow(ax, 5.0, 6.5, 5.0, 7.0, color="royalblue")
    save(fig, out)


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    case01(f"{base}/case01_no_pump/schematic.png")
    case02(f"{base}/case02_type3_fixed/schematic.png")
    case03(f"{base}/case03_type3_timeseries/schematic.png")
    case04(f"{base}/case04_ideal_timeseries/schematic.png")
    print("schematics written")


if __name__ == "__main__":
    main()

"""Write the ARC geometry and optics block of an ARCTracing card from the k4geo sources (no external file needed at run time).

  python3 arc_card_block.py <RadiatorCell_FinalOptimisation xml> <optical tables> > block.tcl

<RadiatorCell_FinalOptimisation xml>: per-cell mirror and sensor parameters of the ARC (k4geo
  FCCee/CLD/compact/CLD_o3_v01/RadiatorCell_FinalOptimisation_o1_v01.xml), in m and rad.
<optical tables>: "TABLE <name> <n>" followed by n lines "energy[eV] value", from the Geant4 materials of the geometry.
The block is pasted inside the ARCTracing module of a card. The cell layout (which cell types sit where, 27 barrel phi
copies, 6 endcap sectors with mirrored copies) is that of k4geo ARC_geo_o1_v01; the module rebuilds every cell from it."""
import re
import sys

xml, tables = sys.argv[1], sys.argv[2]
P = {}
for m in re.finditer(r'name="(\w+?)_(Curvature|XPosition|ZPosition|DetPosition|DetTilt)"\s+value="([-\d.eE+]+)\s*\*\s*(m|rad)"', open(xml).read()):
    P.setdefault(m.group(1), {})[m.group(2)] = float(m.group(3)) * (1000.0 if m.group(4) == "m" else 1.0)
KEYS = ("Curvature", "XPosition", "ZPosition", "DetPosition", "DetTilt")
fmt = lambda v: f"{v:.6g}"

print("  ## ---- ARC cell optics (k4geo ARC_geo_o1_v01 construction; mm, rad) ----")
print("  ## barrel cell types: column row mirror-radius mirror-x mirror-z sensor-offset sensor-tilt")
for name in sorted((n for n in P if n.startswith("Radiator_")), key=lambda n: tuple(int(x) for x in re.findall(r"\d+", n))):
    col, row = (int(x) for x in re.findall(r"\d+", name))
    print(f"  add BarrelCellTypes {{{col} {row} " + " ".join(fmt(P[name][k]) for k in KEYS) + "}")
# endcap: the cells of one 60-degree sector as placed by the driver (row, column, x in hexagon sides, y in apothems, mirrored copy)
EC = [(1, 2, 0, 4, 0), (1, 3, 0, 6, 0), (1, 4, 0, 8, 0), (1, 5, 0, 10, 0), (1, 6, 0, 12, 0), (2, 2, -1.5, 3, 0), (2, 3, -1.5, 5, 1),
      (2, 4, -1.5, 7, 1), (2, 5, -1.5, 9, 1), (2, 6, -1.5, 11, 1), (2, 7, -1.5, 13, 1), (3, 3, -3, 6, 0), (3, 4, -3, 8, 1),
      (3, 5, -3, 10, 1), (3, 6, -3, 12, 1), (4, 5, -4.5, 9, 0), (4, 6, -4.5, 11, 1)]
print("  ## endcap cells of one sector: row column x/side y/apothem mirrored-copy mirror-radius mirror-x mirror-z sensor-offset sensor-tilt")
for row, col, x, y, ref in EC:
    q = P[f"EndCapRadiator_c{col}_r{row}"]
    print(f"  add EndcapCells {{{row} {col} {x:g} {y:g} {ref} " + " ".join(fmt(q[k]) for k in KEYS) + "}")

print("  ## ---- optical tables (Geant4 materials of the geometry): name, then energy [eV] and value pairs ----")
# the tables an ARCTracing card can use (indices, attenuation, SiPM detection efficiency, mirror reflectivity)
USED = {"ARC_SiPM_QuantumEfficiency", "REFLECTIVITY__ARC_Mirror", "RINDEX__Xenon_ARC", "ABSLENGTH__Xenon_ARC", "RINDEX__C4F10_PFRICH",
        "ABSLENGTH__C4F10_PFRICH", "RINDEX__Aerogel_PFRICH", "ABSLENGTH__Aerogel_PFRICH", "RAYLEIGH__Aerogel_BelleII", "RAYLEIGH__Aerogel_PFRICH"}
lines = open(tables).read().splitlines()
i = 0
while i < len(lines):
    f = lines[i].split()
    if f and f[0] == "TABLE":
        name, n = f[1], int(f[2])
        if name not in USED:
            i += n + 1
            continue
        pairs = [lines[i + 1 + k].split()[:2] for k in range(n)]
        print(f"  add OpticsTable {{{name} " + " ".join(f"{e} {v}" for e, v in pairs) + "}")
        i += n + 1
    else:
        i += 1

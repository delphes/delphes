"""Dump the ARC cell optics of a k4geo compact file for the Delphes ARCTracing module (ray-tracing fast simulation).

One line per placed cell, world and cell-frame quantities in mm / deg. The cell frame has its origin at the TGeo cell node
origin, z along the cell axis (radial in the barrel, beam axis in the endcaps) and x along the mirror-centre offset, so that
the mirror centre, the sensor and the aerogel all lie in the x-z plane (the optics of Martin Tat's stand-alone optimisation).
Usage: python3 dump_arc_cells.py <compact.xml> <out.txt>
"""
import math
import sys

import numpy as np
import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gSystem.Load("libDDCore")
compact, out = sys.argv[1], sys.argv[2]
desc = ROOT.dd4hep.Detector.getInstance()
desc.fromXML(compact)
mgr = desc.manager()


def point(m, v):
    loc, glob = np.array(v, float), np.zeros(3)
    m.LocalToMaster(loc, glob)
    return glob * 10.0  # cm -> mm


def vect(m, v):
    loc, glob = np.array(v, float), np.zeros(3)
    m.LocalToMasterVect(loc, glob)
    return glob


cells = []


def walk(node, mat, depth):
    name = node.GetName()
    if depth > 0 and "arc" not in name.lower():
        return
    if "_cell" in name:
        cells.append((name, ROOT.TGeoHMatrix(mat), node))
        return
    for i in range(node.GetNdaughters()):
        d = node.GetDaughter(i)
        m = ROOT.TGeoHMatrix(mat)
        m.Multiply(d.GetMatrix())
        walk(d, m, depth + 1)


walk(mgr.GetTopNode(), ROOT.TGeoHMatrix(), 0)
print("cells found:", len(cells))

f = open(out, "w")
f.write("# ARC cells from %s\n" % compact)
f.write("# region name O(3) ex(3) ey(3) ez(3) mirror(3) Rmin Rmax sensor(3) sex(3) sey(3) sez(3) hsx hsy hsz "
        "aerogel(3) aero_normal(3) aero_half aero_rmax hex_phi0_deg cell_zmin cell_zmax hex_apothem  [mm, deg; cell frame]\n")
count = {"barrel": 0, "endcap": 0}
shown = set()
for name, M, node in cells:
    region = "barrel" if "barrel" in name.lower() else "endcap"
    O, ez = point(M, (0, 0, 0)), vect(M, (0, 0, 1))
    mir = sen = aer = None
    for k in range(node.GetNdaughters()):
        d = node.GetDaughter(k)
        nm, sh = d.GetName(), d.GetVolume().GetShape()
        Md = ROOT.TGeoHMatrix(M)
        Md.Multiply(d.GetMatrix())
        if "mirror" in nm:
            bn = sh.GetBoolNode()
            t = bn.GetRightMatrix().GetTranslation()
            sph = bn.GetRightShape()
            mir = dict(C=point(Md, (t[0], t[1], t[2])), Rmin=sph.GetRmin() * 10, Rmax=sph.GetRmax() * 10)
        elif "sensor" in nm:
            sen = dict(S=point(Md, (0, 0, 0)), ex=vect(Md, (1, 0, 0)), ey=vect(Md, (0, 1, 0)), ez=vect(Md, (0, 0, 1)),
                       h=(sh.GetDX() * 10, sh.GetDY() * 10, sh.GetDZ() * 10))
        elif "aerogel" in nm:
            bn = sh.GetBoolNode()
            Ma = ROOT.TGeoHMatrix(Md)
            Ma.Multiply(bn.GetRightMatrix())
            tube = bn.GetRightShape()
            aer = dict(A=point(Ma, (0, 0, 0)), n=vect(Ma, (0, 0, 1)), half=tube.GetDz() * 10, rmax=tube.GetRmax() * 10)
    if mir is None or sen is None or aer is None:
        print("incomplete cell", name)
        continue
    # cell frame: z = cell axis, x = mirror-centre offset (sensor offset if the mirror is on the axis)
    v = mir["C"] - O
    v -= (v @ ez) * ez
    if np.linalg.norm(v) < 0.5:
        v = sen["S"] - O
        v -= (v @ ez) * ez
    ex = v / np.linalg.norm(v)
    ey = np.cross(ez, ex)
    loc = lambda P: np.array([(P - O) @ ex, (P - O) @ ey, (P - O) @ ez])
    lv = lambda w: np.array([w @ ex, w @ ey, w @ ez])
    # cell solid: hexagon orientation and z range
    sh = node.GetVolume().GetShape()
    cls = sh.ClassName()
    if cls == "TGeoTessellated":
        vs = np.array([[sh.GetVertex(i).x(), sh.GetVertex(i).y(), sh.GetVertex(i).z()] for i in range(sh.GetNvertices())]) * 10
        zmin, zmax = vs[:, 2].min(), vs[:, 2].max()
        inner = vs[np.abs(vs[:, 2] - zmin) < 1e-3]
        apothem = np.abs(inner[:, 0]).max()  # flat sides face local x; vertex along local y
        hexv = vect(M, (0, 1, 0))
    else:  # endcap: polygon prism, vertex along local x
        zmin, zmax = -sh.GetDZ() * 10, sh.GetDZ() * 10
        apothem = sh.GetRmax(0) * 10 if hasattr(sh, "GetRmax") and cls == "TGeoPgon" else float("nan")
        hexv = vect(M, (1, 0, 0))
    hv = lv(hexv)
    phi0 = math.degrees(math.atan2(hv[1], hv[0]))
    sl, al = loc(sen["S"]), loc(aer["A"])
    row = ([region, name] + list(O) + list(ex) + list(ey) + list(ez) + list(loc(mir["C"])) + [mir["Rmin"], mir["Rmax"]]
           + list(sl) + list(lv(sen["ex"])) + list(lv(sen["ey"])) + list(lv(sen["ez"])) + list(sen["h"])
           + list(al) + list(lv(aer["n"])) + [aer["half"], aer["rmax"], phi0, zmin, zmax, apothem])
    f.write(" ".join(f"{x:.7f}" if isinstance(x, float) else str(x) for x in row) + "\n")
    count[region] += 1
    key = (region, cls)
    if key not in shown:
        shown.add(key)
        print(f"\n{name}: shape {cls}, O world {np.round(O, 1)}, ez {np.round(ez, 3)}, ex {np.round(ex, 3)}")
        print(f"  mirror centre local {np.round(loc(mir['C']), 1)}  Rmin {mir['Rmin']:.1f}  apex z {loc(mir['C'])[2] + mir['Rmin']:.1f}")
        print(f"  sensor centre local {np.round(sl, 1)}  normal local {np.round(lv(sen['ez']), 3)}  half {np.round(sen['h'], 1)}")
        print(f"  aerogel centre local {np.round(al, 1)} normal {np.round(lv(aer['n']), 3)} half {aer['half']:.1f} rmax {aer['rmax']:.1f}")
        print(f"  hex phi0 {phi0:.1f} deg, cell z {zmin:.1f}..{zmax:.1f}, apothem {apothem:.1f}")
f.close()
print("\nwritten", out, count)

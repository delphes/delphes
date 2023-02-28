#!/usr/bin/env python
#################################################################################################
#
# Requires python v2.7
#
# Simple script for generation of a LHE file containing events distributed with a flat (log) e, pt and eta spectrum.
#
# For help type:
# python flatGunLHEventProducer.py -h
#
#################################################################################################
import os, random, math, argparse, sys

allowed_pdgCodes = [
    1,
    2,
    3,
    4,
    5,
    6,
    11,
    12,
    13,
    14,
    15,
    16,
    21,
    22,
    23,
    24,
    25,
    211,
    130,
]

# define particle masses
mass = dict()
mass[1] = 0.0
mass[2] = 0.0
mass[3] = 0.095
mass[4] = 1.27e00
mass[5] = 4.7
mass[6] = 173.0
mass[11] = 5.11e-04
mass[12] = 0.0
mass[13] = 1.056600e-01
mass[14] = 0.0
mass[15] = 1.777000e00
mass[16] = 0.0
mass[21] = 0.0
mass[22] = 0.0
mass[23] = 9.118800e01
mass[24] = 80.419002
mass[25] = 1.250000e02
mass[211] = 0.13957
mass[130] = 0.497611

charge = dict()
charge[1] = 1.0
charge[2] = 1.0
charge[3] = 1.0
charge[4] = 1.0
charge[5] = 1.0
charge[6] = 1.0
charge[11] = 1.0
charge[12] = 0.0
charge[13] = 1.0
charge[14] = 0.0
charge[15] = 1.0
charge[16] = 0.0
charge[21] = 0.0
charge[22] = 0.0
charge[23] = 1.0
charge[24] = 0.0
charge[25] = 0.0
charge[211] = 1.0
charge[130] = 0.0

color = dict()
color[1] = (501, 0)
color[2] = (501, 0)
color[3] = (501, 0)
color[4] = (501, 0)
color[5] = (501, 0)
color[6] = (501, 0)
color[11] = (0, 0)
color[12] = (0, 0)
color[13] = (0, 0)
color[14] = (0, 0)
color[15] = (0, 0)
color[16] = (0, 0)
color[21] = (501, 502)
color[22] = (0, 0)
color[23] = (0, 0)
color[24] = (0, 0)
color[25] = (0, 0)
color[211] = (0, 0)
color[130] = (0, 0)


# __________________________________________________________
def print_params(args):

    print("pdg codes  : ", args.pdg)
    print("pmax       : ", args.pmax)
    print("pmin       : ", args.pmin)
    print("etamin     : ", args.etamin)
    print("etamax     : ", args.etamax)
    print("nevts      : ", args.nevts)
    print("Seed       : ", args.seed)
    print("flat log   : ", args.log)
    print("output     : ", args.output)


# __________________________________________________________
def write_init(args):

    out = open(args.output, "a")

    ebeam = args.pmax
    out.write('<LesHouchesEvents version="3.0">\n')
    out.write("<init>\n")
    out.write("-11 11 {:+.8e} {:+.8e} 0 0 247000 247000 -4 1\n".format(ebeam, ebeam))
    # out.write('1000. 1. 1000. 1\n')
    out.write("4.032802e+01 5.428322e-02 4.032802e+01 1\n")
    out.write("</init>\n")
    out.close()


# __________________________________________________________
def write_event(args, px, py, pz, e, m):
    def alpha_s(q):
        return 12 * math.pi / ((33.0 - 6.0) * math.log(q ** 2 / 0.3 ** 2))

    pdg = args.pdg

    p1 = [0.0, 0.0, e, e, mass[11]]
    p2 = [0.0, 0.0, -e, e, mass[11]]
    p3 = [px, py, pz, e, m]
    p4 = [-px, -py, -pz, e, m]

    out = open(args.output, "a")
    out.write("<event>\n")
    out.write(
        " 4      1 +1000. {:.8e} 7.54677100e-03 {:.8e}\n".format(2 * e, alpha_s(e))
    )

    prefix_antip = dict()
    prefix_antip[0] = ""
    prefix_antip[1] = "-"

    out.write(
        "       -11 -1    0    0  0    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n".format(
            p1[0], p1[1], p1[2], p1[3], p1[4]
        )
    )
    out.write(
        "        11 -1    0    0  0    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n".format(
            p2[0], p2[1], p2[2], p2[3], p2[4]
        )
    )
    out.write(
        "         {}  1    1    2  {}    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n".format(
            pdg, color[pdg][0], color[pdg][1], p3[0], p3[1], p3[2], p3[3], p3[4], 0.0
        )
    )
    out.write(
        "        {}{}  1    1    2  {}    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n".format(
            prefix_antip[charge[pdg]],
            pdg,
            color[pdg][1],
            color[pdg][0],
            p4[0],
            p4[1],
            p4[2],
            p4[3],
            p4[4],
            0.0,
        )
    )

    out.write("</event>\n")


# __________________________________________________________

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--pdg", type=int)
    parser.add_argument(
        "--pmin", type=float, help="minimum pt/e [GeV] (default: 1.)", default=1.0
    )
    parser.add_argument(
        "--pmax",
        type=float,
        help="maximum pt/e [GeV] (default: 50000.)",
        default=50000.0,
    )
    parser.add_argument(
        "--etamin", type=float, help="minimum eta (default: -2.5)", default=-6.0
    )
    parser.add_argument(
        "--etamax", type=float, help="maximum eta (default: 2.5)", default=6.0
    )
    parser.add_argument(
        "--nevts",
        type=int,
        help="number of events to generate (default: 1000)",
        default=1000,
    )
    parser.add_argument(
        "--seed", type=int, help="random seed (default uses cpu time)", default=None
    )
    parser.add_argument(
        "--log", dest="log", help="flat in log pt (default: no)", action="store_true"
    )
    parser.set_defaults(log=False)
    parser.add_argument(
        "--output", help="output LHE file (default: events.lhe)", default="events.lhe"
    )

    args = parser.parse_args()

    # check if provided pdgCode is allowed
    if args.pdg not in allowed_pdgCodes:
        print(
            args.pdg,
            "Please provide a list of supported pdgCodes : {}".format(allowed_pdgCodes),
        )
        sys.exit(0)

    # print user-defined parameters
    print_params(args)
    print("")

    # intialize file and write LHE file header
    out = open(args.output, "w+")
    out.close()
    write_init(args)

    # initialize random seed
    random.seed(args.seed)

    print("Start event generation ...")

    count = 0

    # start event loop
    nfail = 0
    while count < args.nevts:

        phi = random.uniform(0.0, math.pi)
        eta = random.uniform(args.etamin, args.etamax)
        theta = 2 * math.atan(math.exp(-eta))
        m = mass[args.pdg]

        # flat in p or in log(p)
        if args.log:
            p = math.pow(
                10, random.uniform(math.log10(args.pmin), math.log10(args.pmax))
            )
        else:
            p = random.uniform(args.pmin, args.pmax)

        pt = p * math.sin(theta)

        px = p * math.sin(theta) * math.cos(phi)
        py = p * math.sin(theta) * math.sin(phi)
        pz = p * math.cos(theta)
        e = math.sqrt(m ** 2 + p ** 2)

        write_event(args, px, py, pz, e, m)
        count += 1
        if (count + 1) % 500 == 0:
            print(" ... processed {} events ...".format(count + 1))

    out = open(args.output, "a")
    out.write("</LesHouchesEvents>\n")
    out.close()

    # os.system("gzip {}".format(args.output))

    print("")
    print("Event generation completed.")
    print("Output file:")
    # print("{}".format(os.path.abspath(args.output + ".gz")))
    print("{}".format(os.path.abspath(args.output)))

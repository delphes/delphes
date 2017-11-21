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

#__________________________________________________________
def print_params(args):

    print 'sqrt(s)    : ', args.ecm
    print 'pdg codes  : ', args.pdg
    print 'gun mode   : ', args.guntype
    print 'pmax       : ', args.pmax
    print 'pmin       : ', args.pmin
    print 'etamin     : ', args.etamin
    print 'etamax     : ', args.etamax
    print 'nevts      : ', args.nevts
    print 'Seed       : ', args.seed
    print 'flat log   : ', args.log
    print 'output     : ', args.output
#__________________________________________________________
def write_init(args):

    ebeam = args.ecm/2.

    out = open(args.output, "a")
    out.write('<LesHouchesEvents version="3.0">\n')
    out.write('<init>\n')
    out.write('2212 2212 {0} {1} 0 0 1 1 -4 1\n'.format(ebeam, ebeam))
    out.write('1000. 1. 1000. 1\n')
    out.write('</init>\n')
    out.close()
#__________________________________________________________
def write_event(args, pt, eta, phi):

    # define particle masses
    mass = dict()
    mass[1] = 0.
    mass[2] = 0.
    mass[3] = 0.
    mass[4] = 1.29
    mass[5] = 4.7
    mass[6] = 173.
    mass[11] = 0.
    mass[12] = 0.
    mass[13] = 0.
    mass[14] = 0.
    mass[15] = 0.
    mass[16] = 0.
    mass[21] = 0.
    mass[22] = 0.
    mass[23] = 9.118800e+01
    mass[24] = 80.419002
    mass[25] = 1.250000e+02

    def alpha_s(q):
       return 12*math.pi/((33.-6.)*math.log(q**2/0.3**2))

    # randomly pick a pdg code from supplied list
    pdg = random.choice(args.pdg)
    
    # compute particles 4-vectors: px, py, pz, e

    px = pt*math.cos(phi)
    py = pt*math.sin(phi)
    pz = pt*math.sinh(eta)
    e = math.sqrt(mass[pdg]**2 + px**2 + py**2 + pz**2)
    
    p1 = [0., 0., e, e]
    p2 = [0., 0., -e, e]
    p3 = [px, py, pz, e]
    p4 = [-px, -py, -pz, e]
    
    cf_list = []

    # initialize stuff
    color = [501, 501, 0, 0]

    out = open(args.output, "a")
    out.write('<event>\n')
    out.write(' 4      1 +1000. {:.8e} 7.54677100e-03 {:.8e}\n'.format(2*e, alpha_s(e)))


    ########## g g -> g g ###################

    if pdg == 21:
    
       cf_list.append([503, 501, 504, 502, 503, 502, 504, 501])
       cf_list.append([504, 501, 503, 502, 503, 501, 504, 502])

       color = cf_list[random.randint(0, 1)]

       out.write('       21 -1    0    0  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[0], color[1], p1[0], p1[1], p1[2], p1[3], 0.))
       out.write('       21 -1    0    0  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[2], color[3], p2[0], p2[1], p2[2], p2[3], 0.))
       out.write('       21  1    1    2  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[4], color[5], p3[0], p3[1], p3[2], p3[3], 0.))
       out.write('       21  1    1    2  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[6], color[7], p4[0], p4[1], p4[2], p4[3], 0.))

    ########## q q -> q' q' ###############

    elif pdg > 0 and pdg < 7:

       cf_list.append([501, 501, 502, 502])
       cf_list.append([501, 502, 501, 502])

       color = cf_list[random.randint(0, 1)]
       
       out.write('        1 -1    0    0  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[0], p1[0], p1[1], p1[2], p1[3], 0.))
       out.write('       -1 -1    0    0  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n'.format(color[1], p2[0], p2[1], p2[2], p2[3], 0.))
       out.write('        {}  1    1    2  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} 1.0000e+00\n'.format(pdg, color[2], p3[0], p3[1], p3[2], p3[3], mass[pdg], 0.))
       out.write('       -{}  1    1    2  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} -1.0000e+00\n'.format(pdg, color[3], p4[0], p4[1], p4[2], p4[3], mass[pdg], 0.))

    ########## q q -> l l (l = e, mu, tau, ve, vm , vtau)  ###############


    elif pdg > 10 and pdg < 17:

        out.write('        1 -1    0    0  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[0], p1[0], p1[1], p1[2], p1[3], 0.))
        out.write('       -1 -1    0    0  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n'.format(color[1], p2[0], p2[1], p2[2], p2[3], 0.))
        out.write('        {}  1    1    2  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} -1.0000e+00\n'.format(pdg, color[2], p3[0], p3[1], p3[2], p3[3], mass[pdg], 0.))
        out.write('       -{}  1    1    2  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} 1.0000e+00\n'.format(pdg, color[3], p4[0], p4[1], p4[2], p4[3], mass[pdg], 0.))

    ########## q q -> B B  (B = gamma, Z, W, H)  ###############

    elif pdg > 21 and pdg < 26:

        out.write('        1 -1    0    0  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[0], p1[0], p1[1], p1[2], p1[3], 0.))
        out.write('       -1 -1    0    0  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n'.format(color[1], p2[0], p2[1], p2[2], p2[3], 0.))
        out.write('        {}  1    1    2  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} 0.0000e+00\n'.format(pdg, color[2], p3[0], p3[1], p3[2], p3[3], mass[pdg], 0.))
        if pdg == 24:
           out.write('       -{}  1    1    2  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} 0.0000e+00\n'.format(pdg, color[3], p4[0], p4[1], p4[2], p4[3], mass[pdg], 0.))
        else:
           out.write('       {}  1    1    2  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} {:.8e} 0.0000e+00\n'.format(pdg, color[3], p4[0], p4[1], p4[2], p4[3], mass[pdg], 0.))

    out.write('</event>\n')




#__________________________________________________________

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdg", nargs='+', type=int)
    parser.add_argument('--guntype', dest='guntype', help='pt or e gun. The parameters pmin(max) are then interpreted as ptmin(max) or emin(max) depending on the specified gunmode', default='pt')
    parser.add_argument("--pmin", type=float, help="minimum pt/e [GeV] (default: 1.)", default=1.)
    parser.add_argument("--pmax", type=float, help="maximum pt/e [GeV] (default: 50000.)", default=50000.)
    parser.add_argument("--etamin", type=float, help="minimum eta (default: -2.5)", default=-6.)
    parser.add_argument("--etamax", type=float, help="maximum eta (default: 2.5)", default=6.)
    parser.add_argument("--ecm,", dest='ecm', type=float, help="center of mass energy (default: 13000)", default=13000)
    parser.add_argument("--nevts", type=int, help="number of events to generate (default: 1000)", default=1000)
    parser.add_argument("--seed", type=int, help="random seed (default uses cpu time)", default=None)
    parser.add_argument('--log', dest='log', help="flat in log pt (default: yes)", action='store_true')
    parser.add_argument('--nolog', dest='log', help="flat in pt (default: false)", action='store_false')
    parser.set_defaults(log=True)
    parser.add_argument("--output", help="output LHE file (default: events.lhe)", default='events.lhe')

    args = parser.parse_args()
    
    # check if provided pdgCode is allowed
    allowed_pdgCodes = [1,2,3,4,5,6,11,12,13,14,15,16,21,22,23,24,25]
    for p1 in args.pdg:
        if p1 not in allowed_pdgCodes:
           print args.pdg ,'Please provide a list of supported pdgCodes : 1, 2, 3, 4, 5 or 21'
           sys.exit(0)
       
    # print user-defined parameters
    print_params(args)
    print ''
   
    # intialize file and write LHE file header
    out = open(args.output, "w+")
    out.close()
    write_init(args)
    
    # initialize random seed
    random.seed(args.seed)

    print 'Start event generation ...'
    
    ebeam = args.ecm/2.
    count = 0

    if args.guntype not in ['pt', 'e']:
       print 'Please specify either gun mode "pt" or "e".'
       sys.exit(0)

    # start event loop
    nfail = 0
    while count < args.nevts:
       
       phi = random.uniform(0., math.pi)
       eta = random.uniform(args.etamin, args.etamax)
          
       # flat in e/pt or in log e/pt
       if args.log:
          gun = math.pow(10, random.uniform(math.log10(args.pmin), math.log10(args.pmax)))
       else:
          gun = random.uniform(args.pmin, args.pmax)

       # generating "balanced" collision, i.e x1 = x2 = 2*energy/sqrt(s)
       if args.guntype == 'pt':
          pt = gun
          e = pt*math.cosh(eta)
       elif args.guntype == 'e':
          e = gun
          pt = e/math.cosh(eta)

       # write event corresponding to required process
       if e <= ebeam:
          write_event(args, pt, eta, phi)
          count += 1
          if (count+1)%500 == 0:
             print ' ... processed {} events ...'.format(count+1)
       else:
          nfail += 1

       if nfail > 10*args.nevts:
          print 'Too many events fail phase space requirements. Usually means ptmax is too high given eta cuts.'
          print 'To fix this either increase center of mass/ change cuts.'
          out.close()
          os.system('rm {}'.format(args.output))
          sys.exit(0)

    out = open(args.output, "a")
    out.write('</LesHouchesEvents>\n')
    out.close()

    os.system('gzip {}'.format(args.output))

    print ''
    print 'Event generation completed.'
    print 'Output file:'
    print '{}'.format(os.path.abspath(args.output +'.gz'))
   

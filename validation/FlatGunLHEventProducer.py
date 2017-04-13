#!/usr/bin/env python
################################################################################
#
# Requires python v2.7
#
# Simple script for generation of a LHE file containing events distributed with a flat (log)pt and eta spectrum.
#
# For help type:
# python FlatSpectrumProducer.py -h
#
################################################################################

import os, random, math, argparse, sys

def print_params(args):

    print 'sqrt(s)    : ', args.ecm
    print 'ptmin      : ', args.ptmin
    print 'ptmax      : ', args.ptmax
    print 'etamin     : ', args.etamin
    print 'etamax     : ', args.etamax
    print 'size       : ', args.size
    print 'seed       : ', args.seed
    print 'flat log   : ', args.log
    print 'output     : ', args.output

def write_init(args):

    ebeam = args.ecm/2.

    out = open(args.output, "a")
    out.write('<LesHouchesEvents version="3.0">\n')
    out.write('<init>\n')
    out.write('2212 2212 {0} {1} 0 0 1 1 -4 1\n'.format(ebeam, ebeam))
    out.write('1000. 1. 1000. 1\n')
    out.write('</init>\n')
    out.close()

def write_event(args, pt, eta, phi):

    e = pt*math.cosh(eta)

    # compute particles 4-vectors (for massless particles): px, py, pz, e
    
    p1 = [0., 0., e, e]
    p2 = [0., 0., -e, e]
    p3 = [pt*math.cos(phi), pt*math.sin(phi), pt*math.sinh(eta), e]
    p4 = [- pt*math.cos(phi), - pt*math.sin(phi), - pt*math.sinh(eta), e]

    if args.pdg == 21:
       pdg = 0
    else: 
       pdg = args.pdg

    cf_list = []

    ########## quarks
    if pdg > 0:

       cf_list.append([501, 501, 502, 502])
       cf_list.append([501, 502, 501, 502])

       color = cf_list[random.randint(0, 1)]
       
       # assume massless quarks
       out = open(args.output, "a")
       out.write('<event>\n')
       out.write(' 4      1 +1000. {:.8e} 0.78186083E-02 0.11800000E+00\n'.format(2*e))
       out.write('        1 -1    0    0  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[0], p1[0], p1[1], p1[2], p1[3], 0.))
       out.write('       -1 -1    0    0  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n'.format(color[1], p2[0], p2[1], p2[2], p2[3], 0.))
       out.write('        {}  1    1    2  {}    0 {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(pdg, color[2], p3[0], p3[1], p3[2], p3[3], 0.))
       out.write('       -{}  1    1    2  0    {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 -1.0000e+00\n'.format(pdg, color[3], p4[0], p4[1], p4[2], p4[3], 0.))
       out.write('</event>\n')
    
    ########## gluons
    else:
    
       cf_list.append([503, 501, 504, 502, 503, 502, 504, 501])
       cf_list.append([504, 501, 503, 502, 503, 501, 504, 502])

       color = cf_list[random.randint(0, 1)]

       out = open(args.output, "a")
       out.write('<event>\n')
       out.write(' 4      1 +1000. {:.8e} 0.78186083E-02 0.11800000E+00\n'.format(2*e))
       out.write('       21 -1    0    0  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[0], color[1], p1[0], p1[1], p1[2], p1[3], 0.))
       out.write('       21 -1    0    0  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[2], color[3], p2[0], p2[1], p2[2], p2[3], 0.))
       out.write('       21  1    1    2  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[4], color[5], p3[0], p3[1], p3[2], p3[3], 0.))
       out.write('       21  1    1    2  {}   {} {:+.8e} {:+.8e} {:+.8e} {:.8e} {:.8e} 0.0000e+00 1.0000e+00\n'.format(color[6], color[7], p4[0], p4[1], p4[2], p4[3], 0.))
       out.write('</event>\n')



#__________________________________________________________

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdg", type=int, help='available: 1/2/3/4/5/21 i.e: [uubar/ddbar/ssbar/ccbar/bbbar/gg] (default: 1)', default='1')
    parser.add_argument("--ptmin", type=float, help="minimum pt [GeV] (default: 1.)", default=1.)
    parser.add_argument("--ptmax", type=float, help="maximum pt [GeV] (default: 50000.)", default=50000.)
    parser.add_argument("--etamin", type=float, help="minimum eta (default: 6.)", default=-6.)
    parser.add_argument("--etamax", type=float, help="maximum eta (default: 6.)", default=6.)
    parser.add_argument("--ecm,", dest='ecm', type=float, help="center of mass energy (default: 13000)", default=13000)
    parser.add_argument("--size", type=int, help="number of events to generate (default: 1000)", default=1000)
    parser.add_argument("--seed", type=int, help="random seed (default: 0)", default=0)
    parser.add_argument('--log', dest='log', help="flat in log pt (default: yes)", action='store_true')
    parser.add_argument('--nolog', dest='log', help="flat in pt (default: false)", action='store_false')
    parser.set_defaults(log=True)
    parser.add_argument("--output", help="output LHE file (default: events.lhe)", default='events.lhe')

    args = parser.parse_args()
    
    # check if provided pdgCode is allowed
    allowed_pdgCodes = [1, 2, 3, 4, 5, 21]
    if args.pdg not in allowed_pdgCodes:
       print args.pdg ,'Please provide a supported pdgCode : 1, 2, 3, 4, 5 or 21'
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

    # start event loop
    while count < args.size:
       
       phi = random.uniform(0., math.pi)
       eta = random.uniform(args.etamin, args.etamax)
          
       # flat in pt or in logpt
       if args.log:
          pt = math.pow(10, random.uniform(math.log10(args.ptmin), math.log10(args.ptmax)))
       else:
          pt = random.uniform(args.ptmin, args.ptmax)

       # generating "balanced" collision, i.e x1 = x2 = 2*energy/sqrt(s)
       e = pt*math.cosh(eta)

       # write event corresponding to required process
       if e < ebeam/2.:
          write_event(args, pt, eta, phi)
          count += 1
          if (count+1)%1000 == 0:
             print ' ... processed {} events ...'.format(count+1)

    print ''
    print 'Event generation completed.'
    print 'Output file:'
    print '{}'.format(os.path.abspath(args.output))
   

    out = open(args.output, "a")
    out.write('</LesHouchesEvents>\n')
    out.close()

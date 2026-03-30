import DelphesPython as delphes
import matplotlib.pyplot as plt
import ROOT


delphes.load("libDelphesPythia.so")

processor = delphes.Delphes()
#processor.reader = delphes.Reader("LHEF", inputFiles = ["/home/laurent/work/dev/cepgen/build/output.lhe"])
processor.reader = delphes.Reader("Pythia8",
    pythiaConfig = [
        'Beams:idA = 2212               ! proton beam 1',
        'Beams:idB = 2212               ! proton beam 2',
        'Beams:eCM = 13600.             ! 13.6 TeV center-of-mass energy',
        'HardQCD:all = on               ! Enable all 2-to-2 QCD processes (qg, gg, qq, etc.)',
        'PhaseSpace:pThatMin = 50.      ! Minimum pT of hard scattering (GeV)',
    ],
)
processor.loadTCL("../cards/delphes_card_CMS.tcl")

num_events = 5000

jet_algos = [1, 2, 3, 4, 5, 6]
labels = ["CDFJetClu", "MidPoint", "SIScone", "kt", "Cambridge/Aachen", "antikt"]
colours = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue-2, ROOT.kGreen+1, ROOT.kOrange+1, ROOT.kMagenta+2]

jet_pts_vs_jetalgos = []
jet_etas_vs_jetalgos = []

algo_id = 0
for jet_algo in jet_algos:
    processor.modules['GenJetFinder']['JetAlgorithm'] = jet_algo
    print(f"Jet algo:={processor.modules['GenJetFinder']['JetAlgorithm']}")

    jet_pts = ROOT.TH1D(f"jet_pts{algo_id}", labels[algo_id]+";Jet p_{T} (GeV);Frequency", 50, 0., 250.)
    jet_etas = ROOT.TH1D(f"jet_etas{algo_id}", f"{labels[algo_id]};Jet #eta;Frequency", 50, -5., 5.)
    i = 0
    try:
        while i < num_events:
            event = processor.next()
            for jet in event["UniqueObjectFinder/jets"]:
                jet_pts.Fill(jet.pt)
                jet_etas.Fill(jet.eta)
            i += 1
    except KeyboardInterrupt:
        pass
    jet_pts_vs_jetalgos.append(jet_pts)
    jet_etas_vs_jetalgos.append(jet_etas)
    algo_id += 1

c_eta = ROOT.TCanvas()
hs_eta = ROOT.THStack()
for i in range(len(jet_algos)):
    jet_etas_vs_jetalgos[i].SetLineColor(colours[i])
    jet_etas_vs_jetalgos[i].SetFillColorAlpha(colours[i], 0.25)
    hs_eta.Add(jet_etas_vs_jetalgos[i])
hs_eta.Draw("nostack")
hs_eta.GetHistogram().SetTitle(";Jet #eta;Frequency")
l = c_eta.BuildLegend(.75, .55, .95, .9, "", "f")
l.Draw()
c_eta.SaveAs("jet_eta_vs_algo.pdf")

c_pt = ROOT.TCanvas()
hs_pt = ROOT.THStack()
for i in range(len(jet_algos)):
    jet_pts_vs_jetalgos[i].SetLineColor(colours[i])
    jet_pts_vs_jetalgos[i].SetFillColorAlpha(colours[i], 0.25)
    hs_pt.Add(jet_pts_vs_jetalgos[i])
hs_pt.Draw("nostack")
hs_pt.GetHistogram().SetTitle(";Jet p_{T} (GeV);Frequency")
c_pt.SetLogy()
l = c_pt.BuildLegend(.75, .55, .95, .9, "", "f")
l.Draw()
c_pt.SaveAs("jet_pt_vs_algo.pdf")


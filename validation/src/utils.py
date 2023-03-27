import sys, os
import numpy as np
from array import array
import ROOT
import math
from collections import OrderedDict
from ROOT import TVector3, TLorentzVector, TVectorD, TMath
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import matplotlib

# matplotlib.use("tkagg")

plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.labelweight"] = "bold"

from ROOT import gROOT

gROOT.SetBatch(True)

# plt.gcf().subplots_adjust(bottom=0.15)

# _______________________________________________________________________________
""" Particle class """


class Particle:
    def __init__(self, name, label, pid):
        self.name = name
        self.pid = pid
        self.label = label


# _______________________________________________________________________________
""" Observable class """


class Observable:
    def __init__(self, name, varname, label, unit, nbins, xmin, xmax, opt):
        self.name = name
        self.varname = varname
        self.label = label
        self.unit = unit
        self.nbins = nbins
        self.xmin = xmin
        self.xmax = xmax
        self.opt = opt  # abs (reco/gen) or rel (reco-gen)
        self.reso_label = "$\sigma({})$ {}".format(self.label, self.unit)
        if opt == "mean":
            self.reso_label = "$< {} >$ {}".format(self.label, self.unit)

        if opt == "rel":
            self.reso_label = "$\sigma({})/{}$".format(self.label, self.label)
        self.eff_label = "${}$ {}".format(self.label, self.unit)
        self.slice_label = "${}$".format(self.label)


# _______________________________________________________________________________
""" Slice 1D class """


class Slices1D:
    def __init__(self, obs, bins, scale):
        self.obs = obs
        self.bins = bins
        self.scale = scale  # lin/log


# _______________________________________________________________________________
""" Slice 2D class """


class Slices2D:
    def __init__(self, name, obsX, obsY, binsX, binsY):
        self.name = name
        self.obsX = obsX
        self.obsY = obsY
        self.binsX = binsX
        self.binsY = binsY
        bins = []
        for i in binsX:
            for j in binsY:
                bins.append((i, j))
        # print(bins)
        self.bins = bins


# _______________________________________________________________________________
""" Resolution histo class """


class ResolutionHisto:
    def __init__(self, particle, collection, observable, binning, slice, label):
        self.particle = particle
        self.observable = observable
        self.collection = collection
        self.binning = binning
        self.slice = slice
        self.label = label
        self.histograms = OrderedDict()
        self.histogram_names = OrderedDict()
        self.finalhistname = "res_{}_{}_{}_{}_{}_".format(
            collection,
            particle.name,
            observable.name,
            self.observable.opt,
            binning.obs.name,
        )

        for obs, bin2 in slice.items():
            self.finalhistname += "{}{}_{}".format(obs.name, bin2[0], bin2[1])

        for bin in binning.bins:
            histname = "res_{}_{}_{}_{}_{}{}_{}_".format(
                collection,
                particle.name,
                observable.name,
                self.observable.opt,
                binning.obs.name,
                bin[0],
                bin[1],
            )

            for obs, bin2 in slice.items():
                histname += "{}{}_{}".format(obs.name, bin2[0], bin2[1])

            # print(histname)
            self.histogram_names[bin] = histname

    def construct(self):
        for bin in self.binning.bins:
            # print(bin)
            self.histograms[bin] = ROOT.TH1F(
                self.histogram_names[bin],
                self.histogram_names[bin],
                self.observable.nbins,
                self.observable.xmin,
                self.observable.xmax,
            )
            # print("constructing {}".format(self.histogram_names[bin]))

    def fill(self, gen_coll, reco_coll):
        for gen in gen_coll:
            if self.particle.pid not in [1, 2, 3, 4, 5, 21]:
                if abs(gen.PID) != self.particle.pid:
                    continue

            funcname_gen = "get_gen{}".format(self.observable.varname)
            funcname_sliceX = "get_gen{}".format(self.binning.obs.varname)
            funcname_reco = "get_reco{}".format(self.observable.varname)
            val_gen = globals()[funcname_gen](gen)

            val_sliceX = globals()[funcname_sliceX](gen)
            for binX in self.binning.bins:
                if val_sliceX > binX[0] and val_sliceX <= binX[1]:
                    for obs, binY in self.slice.items():
                        funcname_sliceY = "get_gen{}".format(obs.varname)
                        val_sliceY = globals()[funcname_sliceY](gen)
                        if val_sliceY > binY[0] and val_sliceY <= binY[1]:
                            reco = get_reco(gen, reco_coll)
                            if reco:
                                val_reco = globals()[funcname_reco](reco)
                                # print(val_reco, val_gen)
                                if self.observable.opt == "rel":
                                    self.histograms[binX].Fill(val_reco / val_gen)
                                else:
                                    self.histograms[binX].Fill(val_reco - val_gen)

    def write(self):
        for bin in self.binning.bins:
            self.histograms[bin].Write()


# _______________________________________________________________________________
""" ResolutionPlot class """


class ResolutionPlot:
    def __init__(
        self,
        name,
        res_histos,
        text,
    ):
        self.name = name
        self.res_histos = res_histos
        self.text = text
        self.datasets_reso = []

    def construct(self, pid):
        for hist in self.res_histos:
            if hist.particle.pid == pid:
                hist.construct()

    def write_histos(self, pid):
        for hist in self.res_histos:
            if hist.particle.pid == pid:
                hist.write()

    def plot(self, validation_files, outdir):

        ## compute ratios and store them in a new file
        reso_filename = "{}/reso_{}.root".format(
            outdir,
            self.name,
        )

        reso_file = ROOT.TFile(reso_filename, "RECREATE")
        self.datasets = []

        for h in self.res_histos:
            input_file = validation_files[h.particle]

            file = ROOT.TFile(input_file)
            file.cd()

            nbins = len(h.binning.bins)
            bins = array("d", [x[0] for x in h.binning.bins])
            bins.append(h.binning.bins[-1][1])

            final_histogram = ROOT.TH1F(h.finalhistname, h.finalhistname, nbins, bins)

            i = 1
            for bin in h.binning.bins:
                hist = file.Get(h.histogram_names[bin])
                x = (bin[0] + bin[1]) * 0.5

                ## extract place where maximum is
                mode = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())

                ## extract resolution
                # sigma = getEffSigma(hist, wmin=0.0, wmax=2.0, epsilon=0.01)
                # sigma = getEffSigma2(hist)
                # sigma = hist.GetRMS()

                sigma = getSigmaGaus(hist, 1.0)
                sigma_err = hist.GetRMSError()
                mean = hist.GetMean()
                mean_err = hist.GetMeanError()

                value = 0.0
                value_err = 0.0

                ### plot raw sigma from residual plot
                if h.observable.opt == "abs":
                    value = sigma
                    value_err = sigma_err

                ### plot sigma corrected by the mode (relative resolution)
                elif h.observable.opt == "rel":
                    if mode > 0:
                        value = sigma
                        value_err = sigma_err
                    else:
                        # print("did not find histogram maximum in: {}".format(h.histogram_names[bin]))
                        value = 1
                        value_err = 0.0

                ### plot absolute scale
                elif h.observable.opt == "mean":
                    value = mean
                    value_err = mean_err

                debug_str = "{} {}: x={:.2f}, mode={:.2f}, sigma={:.2f}".format(
                    h.histogram_names[bin], h.observable.opt, x, value, value_err
                )

                final_histogram.SetBinContent(i, value)
                final_histogram.SetBinError(i, value_err)
                i += 1

            reso_file.cd()
            final_histogram.Write()
            file.Close()

            self.datasets_reso.append(
                Dataset(
                    reso_filename,
                    h.finalhistname,
                    h.label,
                )
            )

        reso_file.Close()

        ## produce actual plots

        x_label = self.res_histos[0].binning.obs.eff_label
        y_label = self.res_histos[0].observable.reso_label

        self.plot_path = "{}/{}.pdf".format(outdir, self.name)
        plot_reso = PlotHisto1D(self.datasets_reso, self.plot_path, x_label, y_label)

        if "log" in self.res_histos[0].binning.scale:
            plot_reso.set_xscale("log")

        plot_reso.set_yscale("log")

        plot_reso.add_text(self.text)
        plot_reso.set_xmin(self.res_histos[0].binning.bins[0][0])
        plot_reso.plot()

        print("plotted {}".format(self.plot_path))
        return self.plot_path


# _______________________________________________________________________________
""" ResolutionBlock function """


def ResolutionBlock(particle_list, branch, observables, bins_mom, bins_eta, reso_plots):
    reso_dict = dict()
    for obs in observables:
        for p in particle_list:
            reso_pt = []
            for bin in bins_eta[1].bins:
                res_histo_pt = ResolutionHisto(
                    p,
                    branch,
                    obs,
                    bins_mom[0],
                    {bins_eta[1].obs: bin},
                    "{} < ${}$ < {}".format(bin[0], bins_eta[1].obs.label, bin[1]),
                )
                reso_pt.append(res_histo_pt)
            reso_pt_plot = ResolutionPlot(
                "reso_{}_{}_{}_pt".format(branch.lower(), p.name, obs.name),
                reso_pt,
                Text("", (0.5, 0.5)),
            )
            reso_dict[(p, obs, "pt")] = reso_pt_plot
            reso_plots.append(reso_pt_plot)

            reso_eta = []
            for bin in bins_mom[1].bins:
                res_histo_eta = ResolutionHisto(
                    p,
                    branch,
                    obs,
                    bins_eta[0],
                    {bins_mom[1].obs: bin},
                    "{} < ${}$ < {}".format(bin[0], bins_mom[1].obs.label, bin[1]),
                )
                reso_eta.append(res_histo_eta)
            reso_eta_plot = ResolutionPlot(
                "reso_{}_{}_{}_eta".format(branch.lower(), p.name, obs.name),
                reso_eta,
                Text("", (0.5, 0.5)),
            )
            reso_dict[(p, obs, "eta")] = reso_eta_plot
            reso_plots.append(reso_eta_plot)

    return reso_dict


# _______________________________________________________________________________
""" Efficiency histo class """


class EfficiencyHisto:
    def __init__(self, particle, collection, binning, slice, label):
        self.particle = particle
        self.collection = collection
        self.binning = binning
        self.slice = slice
        self.label = label

        self.histname = "{}_{}_{}_".format(
            particle.name,
            collection,
            binning.obs.name,
        )

        for obs, bin in slice.items():
            self.histname += "{}{}_{}".format(obs.name, bin[0], bin[1])

        self.num_histname = "num_{}".format(self.histname)
        self.den_histname = "den_{}".format(self.histname)
        self.eff_histname = "eff_{}".format(self.histname)

        # takes 2D slices as input, plots as function of X and slicing in Y
        self.obs = binning.obs
        self.nbins = len(binning.bins)
        self.bins = binning.bins
        self.bin_array = array("d", [x[0] for x in self.bins])
        self.bin_array.append(self.bins[-1][1])

    def construct(self):
        self.num_hist = ROOT.TH1F(self.num_histname, self.num_histname, self.nbins, self.bin_array)
        self.den_hist = ROOT.TH1F(self.den_histname, self.den_histname, self.nbins, self.bin_array)
        self.num_hist.Sumw2()
        self.den_hist.Sumw2()

    def fill(self, gen_coll, reco_coll):
        funcname_gen_obs = "get_gen{}".format(self.obs.varname)
        funcname_gen_cuts = dict()
        for obs, bin in self.slice.items():
            funcname_gen_cuts[obs] = "get_gen{}".format(obs.varname)
        for gen in gen_coll:
            if abs(gen.PID) != self.particle.pid:
                continue
            val_x = globals()[funcname_gen_obs](gen)
            val_y = dict()
            ok = True

            for obs, bin in self.slice.items():
                val_y = globals()[funcname_gen_cuts[obs]](gen)
                # print(val_x, obs.name, val_y, bin)
                if val_y < bin[0] or val_y > bin[1]:
                    ok = False
            if ok:
                self.den_hist.Fill(val_x)
                reco = get_reco(gen, reco_coll)
                # print("find gen")
                if reco:
                    self.num_hist.Fill(val_x)
                    # print("find reco")

    def write(self):
        self.num_hist.Write()
        self.den_hist.Write()


# _______________________________________________________________________________
""" EfficiencyPlot1D class """


class EfficiencyPlot1D:
    def __init__(self, name, eff_histos, text):
        self.name = name
        self.eff_histos = eff_histos
        self.text = text

    def construct(self, pid):
        for hist in self.eff_histos:
            if hist.particle.pid == pid:
                hist.construct()

    def write_histos(self, pid):
        for hist in self.eff_histos:
            if hist.particle.pid == pid:
                hist.write()

    def plot(self, validation_files, outdir):

        ## compute ratios and store them in a new file
        eff_filename = "{}/{}.root".format(
            outdir,
            self.name,
        )

        eff_file = ROOT.TFile(eff_filename, "RECREATE")

        histnames_labels = OrderedDict()
        self.datasets = []

        for hist in self.eff_histos:
            input_file = validation_files[hist.particle]
            # input_file = "/Users/michele/Programs/delphes/new_validation/python_bis/validation_delphes_card_FCChh_I/eff_track_mu_eta.root"
            file = ROOT.TFile(input_file)
            file.cd()
            hist_num = file.Get(hist.num_histname)
            hist_den = file.Get(hist.den_histname)
            hist_eff = hist_num.Clone()
            # hist_num.Sumw2()
            # hist_den.Sumw2()

            hist_eff.Sumw2()
            hist_eff.Divide(hist_den)
            # hist_eff.Sumw2()
            hist_eff.SetNameTitle(hist.eff_histname, hist.eff_histname)
            eff_file.cd()
            hist_eff.Write()
            file.Close()
            self.datasets.append(
                Dataset(
                    eff_filename,
                    hist.eff_histname,
                    hist.label,
                )
            )
        eff_file.Close()

        self.x_label = self.eff_histos[0].binning.obs.eff_label
        self.y_label = "efficiency"
        self.plot_path = "{}/{}.pdf".format(outdir, self.name)

        plot_eff = PlotHisto1D(self.datasets, self.plot_path, self.x_label, self.y_label)
        plot_eff.set_xscale(self.eff_histos[0].binning.scale)

        # text.set_weight("bold")
        plot_eff.add_text(self.text)
        plot_eff.plot()

        print("plotted {}".format(self.plot_path))
        return self.plot_path


# _______________________________________________________________________________
""" Efficiency histo class """


class EfficiencyTaggingHisto:
    def __init__(self, particle, collection, gen_cond, reco_cond, binning, slice, label):
        self.particle = particle
        self.collection = collection
        self.binning = binning
        self.slice = slice
        self.gen_cond = gen_cond
        self.reco_cond = reco_cond
        self.label = label

        self.histname = "{}_{}_{}_{}_{}_".format(
            particle.name,
            collection,
            gen_cond.replace(" ", ""),
            reco_cond.replace(" ", ""),
            binning.obs.name,
        )

        for obs, bin in slice.items():
            self.histname += "{}{}_{}".format(obs.name, bin[0], bin[1])

        self.num_histname = "num_{}".format(self.histname)
        self.den_histname = "den_{}".format(self.histname)
        self.eff_histname = "eff_{}".format(self.histname)

        # takes 2D slices as input, plots as function of X and slicing in Y
        self.obs = binning.obs
        self.nbins = len(binning.bins)
        self.bins = binning.bins
        self.bin_array = array("d", [x[0] for x in self.bins])
        self.bin_array.append(self.bins[-1][1])

    def construct(self):
        self.num_hist = ROOT.TH1F(self.num_histname, self.num_histname, self.nbins, self.bin_array)
        self.den_hist = ROOT.TH1F(self.den_histname, self.den_histname, self.nbins, self.bin_array)

    def fill(self, jet_coll):
        # need only one jet collection for tagging, since gen info is extracted from them

        gen_coll = [p for p in jet_coll if eval(self.gen_cond)]
        reco_coll = [p for p in jet_coll if eval(self.reco_cond)]

        funcname_gen_obs = "get_gen{}".format(self.obs.varname)
        funcname_gen_cuts = dict()
        for obs, bin in self.slice.items():
            funcname_gen_cuts[obs] = "get_gen{}".format(obs.varname)
        for gen in gen_coll:
            val_x = globals()[funcname_gen_obs](gen)
            val_y = dict()
            ok = True

            for obs, bin in self.slice.items():
                val_y = globals()[funcname_gen_cuts[obs]](gen)
                # print(val_x, obs.name, val_y, bin)
                if val_y < bin[0] or val_y > bin[1]:
                    ok = False
            if ok:
                self.den_hist.Fill(val_x)
                reco = get_reco(gen, reco_coll)
                # print("find gen")
                if reco:
                    self.num_hist.Fill(val_x)
                    # print("find reco")

    def write(self):
        self.num_hist.Write()
        self.den_hist.Write()


# _______________________________________________________________________________
def EfficiencyParticleBlock(particle_list, branch, mom, eta, eff_plots):
    eff_dict = OrderedDict()
    for p in particle_list:

        eff_pt_histos = []
        for bin in eta[1].bins:
            eff_histo = EfficiencyHisto(
                p,
                branch,
                mom[0],
                {eta[1].obs: bin},
                "{} < ${}$ < {}".format(bin[0], eta[1].obs.label, bin[1]),
            )
            eff_pt_histos.append(eff_histo)

        eff_pt_plot = EfficiencyPlot1D(
            "eff_{}_{}_{}".format(branch.lower(), p.name, mom[0].obs.name),
            eff_pt_histos,
            Text("", (0.5, 0.5)),
        )

        eff_dict[(p, "pt")] = eff_pt_plot
        eff_plots.append(eff_pt_plot)

        eff_eta_histos = []
        for bin in mom[1].bins:
            eff_histo = EfficiencyHisto(
                p,
                branch,
                eta[0],
                {mom[1].obs: bin},
                "{} < ${}$ < {}".format(bin[0], mom[1].obs.label, bin[1]),
            )
            eff_eta_histos.append(eff_histo)

        eff_eta_plot = EfficiencyPlot1D(
            "eff_{}_{}_{}".format(branch.lower(), p.name, eta[0].obs.name),
            eff_eta_histos,
            Text("", (0.5, 0.5)),
        )
        eff_dict[(p, "eta")] = eff_eta_plot
        eff_plots.append(eff_eta_plot)

    return eff_dict


# _______________________________________________________________________________
def EfficiencyTaggingBlock(jet_coll, flavor_tag, mom, eta, gen_cond, wps, tag_flag, eff_tag_plots):

    eff_tag_dict = dict()
    for part, cond in gen_cond.items():
        for wp in wps:
            eff_tag_pt_histos = []
            for etabin in eta[1].bins:
                eff_tag_pt = EfficiencyTaggingHisto(
                    part,
                    jet_coll,
                    cond,  # gen condition
                    "p.{} & (1 << {})".format(tag_flag, wp),  # reco condition
                    mom[0],
                    {eta[1].obs: etabin},
                    "{} < ${}$ < {}".format(etabin[0], eta[1].obs.label, etabin[1]),
                )
                eff_tag_pt_histos.append(eff_tag_pt)
            eff_tag_pt_plot = EfficiencyPlot1D(
                "eff_{}_{}_{}_pt".format(flavor_tag, part.name, wp),
                eff_tag_pt_histos,
                Text("", (0.5, 0.5)),
            )
            eff_tag_dict[(part, wp, "pt")] = eff_tag_pt_plot
            eff_tag_plots.append(eff_tag_pt_plot)

            eff_tag_eta_histos = []
            for ptbin in mom[1].bins:
                eff_tag_eta = EfficiencyTaggingHisto(
                    part,
                    jet_coll,
                    cond,  # gen condition
                    "p.{} & (1 << {})".format(tag_flag, wp),  # reco condition
                    eta[0],
                    {mom[1].obs: ptbin},
                    "{} < ${}$ < {}".format(ptbin[0], mom[1].obs.label, ptbin[1]),
                )
                eff_tag_eta_histos.append(eff_tag_eta)

            eff_tag_eta_plot = EfficiencyPlot1D(
                "eff_{}_{}_{}_eta".format(flavor_tag, part.name, wp),
                eff_tag_eta_histos,
                Text("", (0.5, 0.5)),
            )
            eff_tag_dict[(part, wp, "eta")] = eff_tag_eta_plot
            eff_tag_plots.append(eff_tag_eta_plot)

    return eff_tag_dict


# _______________________________________________________________________________
""" 1D plot class """


class PlotHisto1D:
    def __init__(self, data, name, title_x, title_y):
        self.data = data  ## list of Datasets
        self.name = name
        self.title_x = title_x
        self.title_y = title_y
        self.xscale = "linear"
        self.yscale = "linear"
        self.leg_loc = "best"
        self.size = 16
        self.normalize = False
        self.text = []

    def set_title_x(self, title_x):
        self.title_x = title_x

    def set_title_y(self, title_y):
        self.title_y = title_y

    def set_xscale(self, xscale):
        self.xscale = xscale

    def set_yscale(self, yscale):
        self.yscale = yscale

    def set_leg_loc(self, leg_loc):
        self.leg_loc = leg_loc

    def set_normalize(self):
        self.normalize = True

    def set_ymin(self, ymin):
        self.ymin = ymin

    def set_ymax(self, ymax):
        self.ymax = ymax

    def set_xmin(self, xmin):
        self.xmin = xmin

    def set_xmax(self, xmax):
        self.xmax = xmax

    def add_text(self, text):
        self.text.append(text)

    def plot(self):

        fig, ax = plt.subplots()
        samples = self.data

        j = 0

        for sample in samples:

            filename = sample.filename
            file = ROOT.TFile(filename)
            hname = sample.histname
            hist = file.Get(hname)

            # print(hist, hname)
            if not hist:
                continue
            # print(hname)
            x, y = [], []
            yerr = []
            integral = 1.0
            if self.normalize:
                integral = hist.Integral(0, hist.GetNbinsX() + 1)

            nbins = hist.GetNbinsX() + 1

            for i in range(1, nbins):
                # print(hist.GetBinCenter(i), hist.GetBinContent(i) / integral)
                x.append(hist.GetBinCenter(i))
                y.append(hist.GetBinContent(i) / integral)
                yerr.append(hist.GetBinError(i))

            file.Close()
            """
            ax.hist(
                x,
                len(x),
                weights=y,
                label="{}".format(sample.label),
                histtype="step",
                linewidth=2,
            )
            """
            """
            ax.step(
                x,
                y,
                label="{}".format(sample.label),
                # histtype="step",
                linewidth=2,
            )
            """
            ax.errorbar(x, y, yerr=yerr, fmt="o", label="{}".format(sample.label), linewidth=2)

        # Create new legend handles but use the colors from the existing ones

        handles, labels = ax.get_legend_handles_labels()
        # new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
        fontsize = self.size
        ax.legend(
            # handles=new_handles,
            labels=labels,
            frameon=False,
            loc=self.leg_loc,
            fontsize=fontsize,
        )

        # add text to plot
        for text in self.text:
            # weight = text.weight
            size = self.size
            ax.text(
                text.location[0],
                text.location[1],
                text.content,
                verticalalignment="center",
                horizontalalignment="center",
                transform=ax.transAxes,
                weight="bold",
                fontsize=size,
            )

        ax.set_xlabel(self.title_x, fontsize=self.size)
        ax.set_ylabel(self.title_y, fontsize=self.size)

        ax.tick_params(axis="both", labelsize=14)
        ax.grid(linestyle="dashed")
        fig.tight_layout()
        fig_file = "{}".format(self.name)
        ax.set_xscale(self.xscale)
        ax.set_yscale(self.yscale)
        fig.savefig(fig_file)
        fig.clf()
        ax.cla()
        plt.close("all")


# _______________________________________________________________________________
""" dataset class """


class Dataset:
    def __init__(self, filename, histname, label):
        self.filename = filename  # root file name
        self.histname = histname
        self.label = label


# _______________________________________________________________________________
""" text plot class """


class Text:
    def __init__(self, content, location):
        self.content = content
        # self.weight = "normal"
        # self.size = 12
        self.location = location

    def set_size(self, size):
        self.size = size

    def set_weight(self, weight):
        self.weight = weight


# _______________________________________________________________________________
""" Latex Report Class """


class LatexReport:
    def __init__(self, name, outdir):

        self.name = name
        self.outdir = outdir

        ## replace underscores and remove extension
        self.title = name.replace("_", " ")
        self.tex_lines = "\n".join(
            "{}".format(ln)
            for ln in r"""\documentclass[8pt]{beamer}
        \setbeamertemplate{frametitle}{
        \insertframetitle}
        \setbeamertemplate{footline}{%
          \raisebox{5pt}{\makebox[\paperwidth]{\hfill\makebox[10pt]{\scriptsize\insertframenumber}}}}
        \setbeamertemplate{navigation symbols}{}
        \setbeamersize{text margin left=0mm,text margin right=0mm}
        \usepackage{graphicx}
        \usepackage{bookmark}
        \usepackage{caption,subcaption}
        \graphicspath{ {.} }
        \usepackage[utf8]{inputenc}

        \title{DUMMYNAME}

        %\author{RTB}
        %\institute{CMS}
        \date{\today}

        \begin{document}

        \section{Cover}

        \frame{\titlepage}

        """.replace(
                "DUMMYNAME", self.title
            ).split(
                "\n"
            )
        )

        # self.tex_lines += "\n" + r"\section{Jetpuppi}" + "\n" + r"\subsection{Response}"

    def section(self, name):
        self.tex_lines += "\n" + r"\section{" + name + "}"

    def subsection(self, name):
        self.tex_lines += "\n" + r"\subsection{" + name + "}"

    def subsubsection(self, name):
        self.tex_lines += "\n" + r"\subsubsection{" + name + "}"

    def begin_frame(self, name):
        """Return the beginning lines for a new page."""
        tex_line = "\n" r"\begin{frame}" + "\n" + r"\frametitle{" + name
        tex_line += "\n" + r"}" + r"\begin{figure}"
        # remove the numbering under the subfigures
        tex_line += "\n" + r"\captionsetup[subfigure]{labelformat=empty}"
        self.tex_lines += tex_line

    def subfigure(self, figure, caption):
        """Return tex lines for subfigures."""
        tex_line = r"\begin{subfigure}{0.45\textwidth}" + "\n" + r"\includegraphics[width=\linewidth]{"
        tex_line += figure
        tex_line += r"}" + "\n" + r"\vspace*{-0.15cm}" + "\n" + r"\caption{"

        tex_line += r"\text{{\tiny " + caption + r"}}}"
        tex_line += "\n" + r"\end{subfigure}" + "\n" + r"\hfil"
        return tex_line

    def add_figures(self, figures):
        """Add figure structure to the script."""
        tex_line = "\n"
        plots_per_page = False
        for figure in figures:
            if plt == "resolution":
                new_batch_size = batch_size
                if batch_size > 9:  # at max. 9 figures can be displayed per page
                    new_batch_size = round(batch_size / 2)
                if i != 0 and i % new_batch_size == 0:
                    tex_line += "\n" + r"\end{figure}" + "\n" + r"\end{frame}"
                    tex_line += "\n" + beginFrame(" cont'd")
            elif len(figures) > 9:
                plots_per_page = round(float(len(figure_dict)) / 2)
                if i != 0 and plots_per_page and i % plots_per_page == 0:
                    tex_line += "\n" + r"\end{figure}" + "\n" + r"\end{frame}"
                    tex_line += "\n" + self.begin_frame(" cont'd")
            tex_line += self.subfigure(figure, "") + "\n"
        tex_line += r"\end{figure}" + "\n" + r"\end{frame}" + "\n" + r"\newpage" + 2 * "\n"
        self.tex_lines += tex_line

    def compile(self):
        self.tex_lines += "\n" + r"\end{document}"

        texfile_path = "{}/{}.tex".format(self.outdir, self.name)
        pdffile_path = texfile_path.replace("tex", "pdf")

        with open(texfile_path, "w") as tex_output:
            tex_output.write(self.tex_lines)

        print("\n {} file is created!\n".format(texfile_path))

        print("Creating pdf file... ".format(pdffile_path))

        os.system("pdflatex {}".format(texfile_path))
        os.system("rm *.aux *.log *.toc *.snm *.nav ")
        # print("\nCreated {}".format(pdffile_path))


# _______________________________________________________________________________
""" compute minimal interval that contains 68% area under curve """


def getEffSigma(theHist, wmin=0.2, wmax=1.8, epsilon=0.01):

    point = wmin
    weight = 0.0
    points = []
    thesum = theHist.Integral(0, theHist.GetNbinsX() + 1)

    if thesum == 0:
        return 0.0
    # fill list of bin centers and the integral up to those point
    for i in range(theHist.GetNbinsX()):
        weight += theHist.GetBinContent(i)
        points.append([theHist.GetBinCenter(i), weight / thesum])

    low = wmin
    high = wmax

    width = wmax - wmin
    for i in range(len(points)):
        for j in range(i, len(points)):
            wy = points[j][1] - points[i][1]
            # print(wy)
            if abs(wy - 0.683) < epsilon:
                # print("here")
                wx = points[j][0] - points[i][0]
                if wx < width:
                    low = points[i][0]
                    high = points[j][0]
                    # print(points[j][0], points[i][0], wy, wx)
                    width = wx
    # print(low, high)
    return 0.5 * (high - low)


# _______________________________________________________________________________
""" compute minimal interval that contains 68% area under curve """


def getEffSigma2(theHist):

    weight = 0.0
    points = []
    thesum = theHist.Integral(0, theHist.GetNbinsX() + 1)

    if thesum == 0:
        return 0.0
    # fill list of bin centers and the integral up to those point
    integrals = dict()
    for i in range(theHist.GetNbinsX() + 1):
        for j in range(i, theHist.GetNbinsX() + 1):
            delta_x = 0.5 * (theHist.GetBinCenter(j) - theHist.GetBinCenter(i))
            integrals[delta_x] = abs(0.683 - (theHist.Integral(i, j) / thesum))

    return min(integrals, key=integrals.get)


# _______________________________________________________________________________
""" compute minimal interval that contains 68% area under curve """


def getSigmaGaus(histo, sigma_range=2):

    ## extract place where maximum is
    x0 = histo.GetXaxis().GetBinCenter(histo.GetMaximumBin())
    d = histo.GetRMS()

    # now perform gaussian fit in [x_max_sigm, x_max_sigp]
    f = ROOT.TF1("gaus", "gaus", 0.0, 2.0)

    s = sigma_range
    histo.Fit("gaus", "Q", "", x0 - s * d, x0 + s * d)

    return f.GetParameter(2)


# _______________________________________________________________________________
""" computes bins for resolution slicing """


def compute_bins(xmin, xmax, nbins, precision, opt="lin"):
    list = []
    if opt == "lin":
        # print(xmin, xmax, nbins + 1, precision)
        list = np.round(np.linspace(xmin, xmax, nbins + 1), precision)
    if opt == "log":
        list = np.round(np.exp(np.linspace(np.log(xmin), np.log(xmax), nbins + 1)), precision)
    bins = []
    # print(list)
    for i in range(nbins):
        bin = (list[i], list[i + 1])
        # print(bin)
        bins.append(bin)
    return bins


# _______________________________________________________________________________
""" returns gen track reference """


def get_gentrack(part):
    gen = None
    if hasattr(part, "Particle"):
        gen = part.Particle.GetObject()
    return gen


# _______________________________________________________________________________
""" returns gen tower reference """


def get_gentower(part):
    gen = None
    if hasattr(part, "Particles"):
        if len(part.Particles) > 0:
            gen = part.Particles.At(0)
    return gen


# _______________________________________________________________________________
""" returns gen jet reference """


def get_genjet(part, genjets):
    gen = None
    dr = 999.0
    for genjet in genjets:
        if genjet.P4().DeltaR(part.P4()) < dr:
            dr = genjet.P4().DeltaR(part.P4())
            gen = genjet
    if dr > 0.4:
        gen = None
    return gen


# _______________________________________________________________________________
""" returns closest reco to gen in a given collection"""


def get_reco(genpart, coll):
    best_reco = None
    dr = 999.0
    for reco in coll:
        if genpart.P4().DeltaR(reco.P4()) < dr:
            dr = genpart.P4().DeltaR(reco.P4())
            best_reco = reco
    if dr > 0.4:
        best_reco = None
    return best_reco


# _______________________________________________________________________________
""" returns gen track parameters (D0, phi, C, z, cotan_th) """


def get_genTrackParam(part):

    x = TVector3(part.X, part.Y, part.Z) * 1.0e-03  # in meters
    p = TVector3(part.Px, part.Py, part.Pz)
    Q = part.Charge

    Par = TVectorD(5)

    from src.init import Bz

    a = -Q * Bz * 0.2998  # Units are Tesla, GeV and meters
    pt = p.Pt()
    if a > 1.0e-06:
        # print(Q, Bz, a)
        # Half curvature
        C = a / (2 * pt)
        r2 = x.Perp2()
        cross = x(0) * p(1) - x(1) * p(0)
        T = TMath.Sqrt(pt * pt - 2 * a * cross + a * a * r2)
        phi0 = TMath.ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T)  # Phi0
        D = (T - pt) / a  # Impact parameter D
        if pt > 10.0:
            D = (-2 * cross + a * r2) / (T + pt)

        Par[0] = D  # Store D
        Par[1] = phi0  # Store phi0
        Par[2] = C  # Store C
        # Longitudinal parameters
        B = C * TMath.Sqrt(TMath.Max(r2 - D * D, 0.0) / (1 + 2 * C * D))
        st = TMath.ASin(B) / C
        ct = p(2) / pt
        z0 = x(2) - ct * st
        #
        Par[3] = z0  # Store z0
        Par[4] = ct  # Store cot(theta)
        #

    else:
        Par[0] = -999.0
        Par[1] = p.Phi()
        Par[2] = -999.0
        Par[3] = -999.0
        if p.Pt() < 1.0e-06:
            Par[4] = TMath.Sign(part.Pz, 1.0) * 999999
        else:
            Par[4] = 1 / math.tan(p.Theta())

    return Par


# ________________________________________________________________________________
def getC(Q, pt):
    from src.init import Bz

    a = -Q * Bz * 0.2998  # Units are Tesla, GeV and meters
    C = 0
    if a > 1.0e-06:
        # print(Q, Bz, a)
        # Half curvature
        C = a / (2 * pt) * 1e03  # in mm^-1
    return C


# _______________________________________________________________________________
""" returns arccotan for x in [0;PI] """


def arccotan(x):
    return math.pi / 2 - math.atan2(x, 1)


# _______________________________________________________________________________
""" returns gen D0 """


def get_genD0(part):
    # return get_genTrackParam(part)[0] * 1.0e03
    return 0.0


# _______________________________________________________________________________
""" returns reco D0 """


def get_recoD0(part):
    return part.D0


# _______________________________________________________________________________
""" returns gen DZ """


def get_genDZ(part):
    # return get_genTrackParam(part)[3] * 1.0e03
    return 0.0


# _______________________________________________________________________________
""" returns reco DZ """


def get_recoDZ(part):
    return part.DZ


# _______________________________________________________________________________
""" returns gen Phi """


def get_genPhi(part):
    # return get_genTrackParam(part)[1]
    return part.P4().Phi()


# _______________________________________________________________________________
""" returns reco Phi """


def get_recoPhi(part):
    return part.Phi


# _______________________________________________________________________________
""" returns gen CtgTheta """


def get_genCtgTheta(part):
    # return get_genTrackParam(part)[4]
    return math.cos(get_genTheta(part)) / math.sin(get_genTheta(part))


# _______________________________________________________________________________
""" returns reco CtgTheta """


def get_recoCtgTheta(part):
    return part.CtgTheta


# _______________________________________________________________________________
""" returns gen Theta """


def get_genTheta(part):
    # return arccotan(get_genTrackParam(part)[4])
    # return part.P4().Theta()
    if abs(part.Eta) > 10:
        return math.copysign(1, part.Eta) * 3.14
    else:
        return 2 * math.atan(math.exp(-part.Eta))


# _______________________________________________________________________________
""" returns reco Theta """


def get_recoTheta(part):
    if hasattr(part, "CtgTheta"):
        # return arccotan(part.CtgTheta)
        # return part.P4().Theta()
        return 2 * math.atan(math.exp(-part.Eta))
    else:
        # return part.P4().Theta()
        return 2 * math.atan(math.exp(-part.Eta))


# _______________________________________________________________________________
""" returns gen P """


def get_genP(part):
    return part.P4().P()


# _______________________________________________________________________________
""" returns reco P """


def get_recoP(part):
    return part.P


# _______________________________________________________________________________
""" returns gen C """


def get_genC(part):
    return getC(part.Charge, part.PT)


# _______________________________________________________________________________
""" returns reco C """


def get_recoC(part):
    return returnC(part.Charge, part.PT)


# _______________________________________________________________________________
""" returns gen PT """


def get_genPT(part):
    return part.PT


# _______________________________________________________________________________
""" returns reco PT """


def get_recoPT(part):
    return part.PT


# _______________________________________________________________________________
""" returns gen E """


def get_genE(part):
    if hasattr(part, "E"):
        return part.E
    else:
        return part.P4().E()


# _______________________________________________________________________________
""" returns reco E """


def get_recoE(part):
    if hasattr(part, "E"):
        return part.E
    else:
        return part.P4().E()


# _______________________________________________________________________________
""" returns gen Eta """


def get_genEta(part):
    return part.Eta


# _______________________________________________________________________________
""" returns reco Eta """


def get_recoEta(part):
    return part.Eta


# ________________________________________________________________________________
""" compute correlation element i, j  """


def get_corr(part, i, j):

    cov_xy = part.CovarianceMatrix()(i, j)
    sigma_x = math.sqrt(part.CovarianceMatrix()(i, i))
    sigma_y = math.sqrt(part.CovarianceMatrix()(j, j))
    return cov_xy / (sigma_x * sigma_y)


# ________________________________________________________________________________
""" returns corr matrix:  D0, Pḧi, C, DZ, CtgTheta """


def get_recoCorrD0Phi(part):
    return get_corr(part, 0, 1)


def get_genCorrD0Phi(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrD0C(part):
    return get_corr(part, 0, 2)


def get_genCorrD0C(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrD0DZ(part):
    return get_corr(part, 0, 3)


def get_genCorrD0DZ(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrD0CtgTheta(part):
    return get_corr(part, 0, 4)


def get_genCorrD0CtgTheta(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrPhiC(part):
    return get_corr(part, 1, 2)


def get_genCorrPhiC(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrPhiDZ(part):
    return get_corr(part, 1, 3)


def get_genCorrPhiDZ(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrPhiCtgTheta(part):
    return get_corr(part, 1, 4)


def get_genCorrPhiCtgTheta(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrCDZ(part):
    return get_corr(part, 2, 3)


def get_genCorrCDZ(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrCCtgTheta(part):
    return get_corr(part, 2, 4)


def get_genCorrCCtgTheta(part):
    return 0.0


# ________________________________________________________________________________


def get_recoCorrDZCtgTheta(part):
    return get_corr(part, 3, 4)


def get_genCorrDZCtgTheta(part):
    return 0.0


# ________________________________________________________________________________

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
plt.gcf().subplots_adjust(bottom=0.15)

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
        self.opt = opt  # absolute (reco/gen) or relative (reco-gen)
        self.reso_label = "$\sigma$({}) {}".format(self.label, self.unit)
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
    def __init__(self, particle, observable, slices, label):
        self.particle = particle
        self.observable = observable
        self.slices = slices
        self.label = label
        self.histograms = OrderedDict()
        self.histogram_names = OrderedDict()

        for bin in slices.bins:
            # print(particle.name, observable.name, slices.obsX.name,slices.obsY.name, bin)
            binX = bin[0]
            binY = bin[1]
            histname = "res_{}_{}_{}_{}_{}{}_{}_{}{}_{}_{}".format(
                label,
                particle.name,
                observable.name,
                self.observable.opt,
                slices.obsX.name,
                binX[0],
                binX[1],
                slices.obsY.name,
                binY[0],
                binY[1],
                slices.name,
            )
            # print(histname)
            self.histogram_names[bin] = histname

            self.histograms[bin] = ROOT.TH1F(
                histname, histname, observable.nbins, observable.xmin, observable.xmax
            )

    def fill(self, reco_particle, gen_particle):
        funcname_gen = "get_gen{}".format(self.observable.varname)
        funcname_sliceX = "get_gen{}".format(self.slices.obsX.varname)
        funcname_sliceY = "get_gen{}".format(self.slices.obsY.varname)
        funcname_reco = "get_reco{}".format(self.observable.varname)
        val_gen = globals()[funcname_gen](gen_particle)
        if val_gen > 0:
            val_sliceX = globals()[funcname_sliceX](gen_particle)
            val_sliceY = globals()[funcname_sliceY](gen_particle)
            val_reco = globals()[funcname_reco](reco_particle)
            # print(self.observable.name, val_gen, val_reco)
            for bin in self.slices.bins:
                binX = bin[0]
                binY = bin[1]
                if (
                    val_sliceX > binX[0]
                    and val_sliceX <= binX[1]
                    and val_sliceY > binY[0]
                    and val_sliceY <= binY[1]
                ):
                    if self.observable.opt == "rel":
                        self.histograms[bin].Fill(val_reco / val_gen)
                    elif self.observable.opt == "abs":
                        self.histograms[bin].Fill(val_reco - val_gen)

    def write(self):
        for bin in self.slices.bins:
            self.histograms[bin].Write()


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

    def construct(self):
        for hist in self.eff_histos:
            hist.construct()

    def write_histos(self):
        for hist in self.eff_histos:
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
            hist_eff.Divide(hist_den)
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

        self.x_label = "${}$".format(self.eff_histos[0].binning.obs.label)
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
def EfficiencyBlock(particle_list, branch, bins_mom, bins_eta, eff_plots):
    eff_histos_mom = OrderedDict()
    eff_histos_eta = OrderedDict()

    eff_plot_mom = OrderedDict()
    eff_plot_eta = OrderedDict()

    for p in particle_list:
        # for p in [pion]:
        eff_histos_mom[p] = []
        eff_histos_eta[p] = []
        for bin in bins_eta[1].bins:
            eff_histos_mom[p].append(
                EfficiencyHisto(
                    p,
                    branch,
                    bins_mom[0],
                    {bins_eta[1].obs: bin},
                    "{} < ${}$ < {}".format(bin[0], bins_eta[1].obs.label, bin[1]),
                )
            )

        for bin in bins_mom[1].bins:
            eff_histos_eta[p].append(
                EfficiencyHisto(
                    p,
                    branch,
                    bins_eta[0],
                    {bins_mom[1].obs: bin},
                    "{} < ${}$ < {}".format(bin[0], bins_mom[1].obs.label, bin[1]),
                )
            )

        eff_plot_mom[p] = EfficiencyPlot1D(
            "eff_{}_{}_{}".format(branch.lower(), p.name, bins_mom[0].obs.name),
            eff_histos_mom[p],
            Text("{} {}".format(p.label, branch.lower()), (0.5, 0.5)),
        )
        eff_plot_eta[p] = EfficiencyPlot1D(
            "eff_{}_{}_{}".format(branch.lower(), p.name, bins_eta[0].obs.name),
            eff_histos_eta[p],
            Text("{} {}".format(p.label, branch.lower()), (0.5, 0.5)),
        )

        eff_plots.append(eff_plot_mom[p])
        eff_plots.append(eff_plot_eta[p])

    return eff_plot_mom, eff_plot_eta


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
            integral = 1.0
            if self.normalize:
                integral = hist.Integral(0, hist.GetNbinsX() + 1)

            for i in range(1, hist.GetNbinsX() + 1):
                x.append(hist.GetBinCenter(i))
                y.append(hist.GetBinContent(i) / integral)

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

            ax.step(
                x,
                y,
                label="{}".format(sample.label),
                # histtype="step",
                linewidth=2,
            )

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

        if hasattr(self, "ymin") and hasattr(self, "ymax"):
            ax.set_ylim(self.ymin, self.ymax)
        if hasattr(self, "xmin") and hasattr(self, "xmax"):
            ax.set_xlim(self.xmin, self.xmax)

        ax.set_xscale(self.xscale)
        ax.set_yscale(self.yscale)
        fig.tight_layout()
        fig_file = "{}".format(self.name)
        fig.savefig(fig_file)
        return fig, ax


# _______________________________________________________________________________
""" ResolutionPlot class """


class ResolutionPlot:
    def __init__(self, particle, colstring, observable, slices, input_file, outdir):
        self.particle = particle
        self.colstring = colstring
        self.observable = observable
        self.slices = slices
        self.input_file = input_file
        self.reso_histos = ResolutionHisto(particle, observable, slices, colstring)
        ## slice resolution plots according to Y dimension

        self.x_label = self.reso_histos.slices.obsX.eff_label
        self.y_label = self.reso_histos.observable.reso_label
        self.reso_plot_name = "{}/reso_{}_{}_{}_{}_{}".format(
            outdir,
            self.reso_histos.label,
            self.reso_histos.particle.name,
            self.reso_histos.observable.name,
            self.reso_histos.observable.opt,
            self.reso_histos.slices.name,
        )
        self.plot_path = "{}.pdf".format(self.reso_plot_name)
        # write tree
        resolution_file = "{}.root".format(self.reso_plot_name)
        root_reso = ROOT.TFile(resolution_file, "RECREATE")

        slices = self.reso_histos.slices.binsY
        bins = self.reso_histos.slices.binsX

        histnames_labels = OrderedDict()
        for slice in slices:
            nbinsX = len(self.reso_histos.slices.binsX)
            binsX = [x[0] for x in self.reso_histos.slices.binsX]
            binsX.append(self.reso_histos.slices.binsX[-1][1])
            histname = "{}_{}_{}_{}_{}{}_{}_{}".format(
                self.reso_histos.label,
                self.reso_histos.particle.name,
                self.reso_histos.slices.obsX.name,
                self.reso_histos.observable.opt,
                self.reso_histos.slices.obsY.name,
                slice[0],
                slice[1],
                self.reso_histos.slices.name,
            )
            histname = "reso_{}".format(histname)
            histogram = ROOT.TH1F(histname, histname, nbinsX, array("d", binsX))
            label = "{} < {} < {}".format(
                slice[0], self.reso_histos.slices.obsY.slice_label, slice[1]
            )
            histnames_labels[label] = histname

            for bin in bins:
                file = ROOT.TFile(input_file)
                hist = file.Get(self.reso_histos.histogram_names[(bin, slice)])
                print(input_file, self.reso_histos.histogram_names[(bin, slice)])
                # hist.Rebin(2)
                x = (bin[0] + bin[1]) * 0.5

                ## extract place where maximum is
                mode = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())

                ## extract resolution
                # sigma = getEffSigma(hist, wmin=0.0, wmax=2.0, epsilon=0.01)
                # sigma = getFWHM(hist) / 2.35
                sigma = hist.GetRMS()

                debug_str = "{}: x={:.2f}, mode={:.2f}, sigma={:.2f}".format(
                    histname, x, mode, sigma
                )
                # print(debug_str)
                if mode > 0:
                    sigma = sigma / mode
                else:
                    print("Did not find histogram maximum ...")

                histogram.Fill(x, sigma)

            root_reso.cd()
            histogram.Write()
        root_reso.Close()

        self.datasets_reso = []
        ## prduce resolution plots

        for label, histname in histnames_labels.items():
            self.datasets_reso.append(
                Dataset(
                    resolution_file,
                    histname,
                    label,
                )
            )

    def plot(self):
        # print(datasets_reso, self.reso_plot_name, self.x_label, self.y_label)
        plot_reso = PlotHisto1D(self.datasets_reso, self.reso_plot_name, self.x_label, self.y_label)

        xobs_str = self.reso_histos.slices.name.partition("_")[0]
        yobs_str = self.reso_histos.slices.name.partition("_")[2]

        if "log" in xobs_str:
            plot_reso.set_xscale("log")

        text_str = "{} ({})".format(self.colstring, self.particle.label)
        text = Text(text_str, (0.5, 0.5))
        text.set_weight("bold")
        plot_reso.text.append(text)
        plot_reso.plot()
        print("plotted {} ... ".format(self.plot_path))
        return self.plot_path


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
        tex_line = (
            r"\begin{subfigure}{0.45\textwidth}" + "\n" + r"\includegraphics[width=\linewidth]{"
        )
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

    x = TVector3(part.X, part.Y, part.Y) * 1.0e-03  # in meters
    p = TVector3(part.Px, part.Py, part.Pz)
    Q = part.Charge

    Par = TVectorD(5)

    from DelphesValidationInit import Bz

    a = -Q * Bz * 0.2998  # Units are Tesla, GeV and meters
    pt = p.Pt()
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
    return Par


# _______________________________________________________________________________
""" returns arccotan for x in [0;PI] """


def arccotan(x):
    return math.pi / 2 - math.atan(x)


# _______________________________________________________________________________
""" returns gen D0 """


def get_genD0(part):
    return get_genTrackParam(part)[0] * 1.0e03


# _______________________________________________________________________________
""" returns reco D0 """


def get_recoD0(part):
    return part.D0


# _______________________________________________________________________________
""" returns gen DZ """


def get_genDZ(part):
    return get_genTrackParam(part)[3] * 1.0e03


# _______________________________________________________________________________
""" returns reco DZ """


def get_recoDZ(part):
    return part.DZ


# _______________________________________________________________________________
""" returns gen Phi """


def get_genPhi(part):
    return get_genTrackParam(part)[1]


# _______________________________________________________________________________
""" returns reco Phi """


def get_recoPhi(part):
    return part.Phi


# _______________________________________________________________________________
""" returns gen CtgTheta """


def get_genCtgTheta(part):
    return get_genTrackParam(part)[4]


# _______________________________________________________________________________
""" returns reco CtgTheta """


def get_recoCtgTheta(part):
    return part.CtgTheta


# _______________________________________________________________________________
""" returns gen Theta """


def get_genTheta(part):
    return arccotan(get_genTrackParam(part)[4])


# _______________________________________________________________________________
""" returns reco Theta """


def get_recoTheta(part):
    return arccotan(part.CtgTheta)


# _______________________________________________________________________________
""" returns gen P """


def get_genP(part):
    return part.P


# _______________________________________________________________________________
""" returns reco P """


def get_recoP(part):
    return part.P


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

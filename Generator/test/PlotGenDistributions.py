"""Draw distributions of generated events for H -> Tau Tau analysis"""
"""usage example: python PlotGenDistributions.py -d RootFiles/ higgs_pt"""

import ROOT as rt
from optparse import OptionParser
from os import listdir


def get_binning(variable):
    """Return a tuple (nbins, min_bin, max_bin) for a given variable"""
    bins = {}
    bins['higgs_pt'] = (100, 0, 500)
    bins['higgs_eta'] = (80, -10, 10)
    bins['higgs_phi'] = (70, -3.5, 3.5)
    bins['higgs_mass'] = (20, 124., 126.)
    bins['tau_pt'] = (100, 0., 500.)
    bins['tau_eta'] = (80, -10., 10.)
    bins['tau_phi'] = (70, -3.5, 3.5)
    bins['tau_nprongsK'] = (8, 0, 8)
    bins['taudecay_vr'] = (40, 0., 4.)
    bins['taudecay_vz'] = (40, -2., 2.)
    bins['tau_nvisstabledaug'] = (15, 0, 15)
    bins['tau_vispt'] = (50, 0., 200.)
    bins['stabledaug_pt'] = (25, 0., 125.)
    bins['stabledaug_eta'] = (80, -10, 10)
    bins['stabledaug_phi'] = (70, -3.5, 3.5)
    bins['stabledaug_mass'] = (60, -0.05, 0.55)
    bins['recohiggs_mass'] = (80, 0., 160.)
    bins['stabledaug_id'] = (140, 210., 350.)
        
    return bins[variable]


def make_histogram(variable, tree, color):
    """Create histogram from TTree 'tree' and TLeaf 'variable'"""

    nbins, min_bin, max_bin = get_binning(variable)
    hist = rt.TH1D(variable, variable, nbins, min_bin, max_bin)
    tree.Project(variable, variable)
    hist.SetLineColor(int(color+1))

    if int(color+1)==10:
        hist.SetLineColor(28)

    return hist


def find_max_histogram(hists):
    """Returns histogram index with highest maximum"""
    n_max = hists[0].GetMaximum()
    max_hist = 0
    for i in range(1, len(hists)):
        if hists[i].GetMaximum() > n_max:
            n_max = hists[i].GetMaximum()
            max_hist = i

    return max_hist


def plot_histograms(hists):
    """Plot histograms into a new canvas; takes a list of histograms"""
    max_h = find_max_histogram(hists)

    c1 = rt.TCanvas()
    hists[max_h].Draw()
    for i in range(len(hists)):
        if i == max_h:
            continue
        hists[i].Draw("same")
    return c1


def create_legend(histos):
    """Returns a TLegend from a list of histograms"""
    leg = rt.TLegend(0.7, 0.6, 0.89, 0.89)
    leg.SetHeader("Higgs productions")
    for legend in histos:
        leg.AddEntry(legend, legend.GetName(), "l")
    return leg

def get_norm(hist):
    """return normalization factor according to difference cross section; assume 100 fb^-1"""
    norm = 0
    totentry = hist.Integral()
    if (hist.GetName()=='ffbar2H'):
        norm = 6.365e-10*1e14/totentry
    if (hist.GetName()=='gg2H'):
        norm = 1.09e-8*1e14/totentry
    if (hist.GetName()=='ffbar2HZ'):
        norm = 2.878e-10*1e14/totentry
    if (hist.GetName()=='ffbar2HW'):
        norm = 5.215e-10*1e14/totentry
    if (hist.GetName()=='ff2Hff_ZZ'):
        norm = 1.313e-09*1e14/totentry
    if (hist.GetName()=='ff2Hff_WW'):
        norm = 3.455e-09*1e14/totentry
    if (hist.GetName()=='gg2Httbar'):
        norm = 1.919e-09*1e14/totentry
    if (hist.GetName()=='qqbar2Httbar'):
        norm =  1.717e-09*1e14/totentry
    return norm


if __name__ == '__main__':

    usage = ("usage: %prog [-d <input-dir> -f <output-file>] <variable>\n\n"
             "You can also specify a single file with the [-i <input>] option"
             " instead of\n[-d <directory>]")
    PARSER = OptionParser(usage)
    PARSER.add_option('-f', '--file', dest="filename", metavar="FILE",
                      default=None, help="output file name")
    PARSER.add_option('-i', '--input', dest="open_file", type="string",
                      help="input file name, if use only one file")
    PARSER.add_option('-d', '--dir', dest="directory", type="string",
                      help="directory with all input files")

    (OPTIONS, ARGS) = PARSER.parse_args()

    # If you want to run multiple files
    if OPTIONS.directory:
        rt.gStyle.SetOptStat(0)
        print ARGS

        my_dir = listdir(OPTIONS.directory)
        histos = []
        for i in range(len(my_dir)):
            if my_dir[i].find('All') is -1:
                print "opening", my_dir[i]
                process_name = my_dir[i][15:-5]

                ROOT_FILE = rt.TFile.Open(OPTIONS.directory + '/' +
                                          str(my_dir[i]))
                TREE = ROOT_FILE.Get("GeneratorNtupleMaker/H2TauTauTree")
                my_hist = make_histogram(ARGS[0], TREE, i)
                my_hist.SetName(process_name)
                my_norm = get_norm(my_hist)
                my_hist.Scale(my_norm)
                histos.append(my_hist)

        my_canvas = plot_histograms(histos)
        leg = create_legend(histos)
        leg.Draw()
        if OPTIONS.filename:
            my_canvas.SaveAs(OPTIONS.filename)
        else:
            my_canvas.SaveAs(ARGS[0] + ".png")
            my_canvas.SaveAs(ARGS[0] + ".C")

    # If you want to run a single root file
    elif OPTIONS.open_file:
        ROOT_FILE = rt.TFile.Open(OPTIONS.open_file)
        TREE = ROOT_FILE.Get("GeneratorNtupleMaker/H2TauTauTree")
        CANVAS = rt.TCanvas()
        my_hist = make_histogram(OPTIONS.variable, TREE, 1)
        my_hist.Draw()
        CANVAS.Print('Dunno.png')

import ROOT
import math
from DataFormats.FWLite import Events, Handle
from array import array
import csv
import numpy
import argparse
import sys


''' example of how to run this: python fitPhiK.py testing_phi_k_RMS.root rms > output_rms.txt '''



parser = argparse.ArgumentParser()
parser.add_argument('file_name', type=str,help='file_name.root')
parser.add_argument('type_plot', choices=['phik','rms'], help='either phik or rms plot to fit parameters')
args = parser.parse_args()

#turn the command line arguments into variables
file_name = args.file_name
type_plot = args.type_plot


def getAll(d,basepath=''):
    "Go into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            for i in getAll(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)

#fit each histogram in the root file
#ROOT.gROOT.SetBatch(True)
f = ROOT.TFile(file_name)
for name, obj in getAll(f):
    if type_plot == 'phik':
	print
	print
	print
	print 'fitting parameters for ',name
	#file_phik = open("fit_params_phik","w")
	canvas = ROOT.TCanvas(name,name)
	obj_pfx = obj.ProfileX()
	fit_func = ROOT.TF1('fit_func','x*[0]/(1+[1]*abs(x)) + [2]',obj_pfx.GetXaxis().GetXmin(),obj_pfx.GetXaxis().GetXmax())
	fitted = obj_pfx.Fit('fit_func','S',)
	canvas.Print(("{}_phik.png").format(name))
	canvas.SaveAs(("{}_phik.root").format(name))
    elif type_plot == 'rms':
	print 
	print
	print
	print 'fitting parameters for ',name
	canvas = ROOT.TCanvas(name,name) 
	fit_func = ROOT.TF1('fit_func','sqrt([0]*[0] + ([1]*[1])*(x*x))',obj.GetXaxis().GetXmin(),obj.GetXaxis().GetXmax())
	fitted = obj.Fit('fit_func','S')
	#fitted.Draw("ALP")
	#canvas.cd()
	obj.Draw("AP")
	fit_func.Draw("same")
	canvas.Print(("{}.png").format(name))
	canvas.SaveAs(("{}.root").format(name))



#!/usr/bin/python
import sys,os
import array
import time
import math
from optparse import OptionParser
from subprocess import call

DEBUG=1

if(DEBUG>0):print "----- BEGIN -----"

if(DEBUG>0):print "-PARSING OPTIONS-"

if __name__=="__main__":
	usage = "usage: %prog [options] arg1 arg2"
	parser=OptionParser(usage=usage)
	parser.add_option("-f","--fileName" ,dest='fileName',type='string',help="File Name with response matrix",default="UnfoldResults.root")
	(options,args)=parser.parse_args()


import ROOT
ROOT.gROOT.SetBatch()

if __name__=="__main__":
	variables=["nCharged","nNeutral","ptD","axis1","axis2"]
	etaRegions=["centr","fwd"]
	if(DEBUG>0): print "--> Load FitTools"
	ROOT.gSystem.Load("fitTools.so")

	ROOT.gROOT.ProcessLine("\
		const int nPtBins = 10;\
		Double_t ptBins[nPtBins+1];\
	")
	from ROOT import nPtBins,ptBins
	ROOT.fitTools.getBins_int( nPtBins, ptBins, 20., 1000.)
	ptBins[nPtBins] = 3500.;
	fRoot=ROOT.TFile.Open(options.fileName)
	for v in variables:
	   for e in etaRegions:	
	      for ptBin in range(0,nPtBins):
			baseName=v+"_"+e+"_pt%d_%d"%(ptBins[ptBin],ptBins[ptBin+1])
			directory=v+"/"+e+"/"
			unfold=fRoot.Get(directory+"unfold"+baseName)
			gen=fRoot.Get(directory+baseName+"_h1g")		
			reco=fRoot.Get(directory+baseName+"_h1r")		
			matrix=fRoot.Get(directory+baseName)	

			C1=ROOT.TCanvas("C1","C1",800,800)
			matrix.GetYaxis().SetTitle("Gen " + v)
			matrix.GetXaxis().SetTitle("Reco " + v)
			matrix.Draw("Box")

			l=ROOT.TLatex()
			l.SetNDC()
			l.SetTextSize(0.04)
			l.SetTextAlign(22)
			l.DrawLatex(0.40,.88,"CMS Simulation, Work in Progress")
			l.DrawLatex(0.40,.82,"#sqrt{s}=8TeV")

			C2=ROOT.TCanvas("C2","C2",800,800)
			gen.SetLineColor(ROOT.kRed)
			reco.SetLineColor(ROOT.kBlue)
			reco.SetLineStyle(ROOT.kDashed)
			reco.SetLineWidth(ROOT.kDashed)
			unfold.SetMarkerStyle(20)
			gen.GetXaxis().SetTitle(v)
			gen.SetMaximum(reco.GetMaximum()*1.2 )
			gen.Draw("HIST")
			reco.Draw("HIST SAME")
			unfold.Draw("P SAME")

			L=ROOT.TLegend(0.65,.65,.89,.89,e+" "+ "%d #leq p_{T}[GeV] #leq%d"%(ptBins[ptBin],ptBins[ptBin+1]) )
			L.AddEntry(gen,"gen","L")
			L.AddEntry(reco,"reco+bkg","L")
			L.AddEntry(unfold,"unfold-closure","P")
			L.Draw()
			l.DrawLatex(0.40,.88,"CMS Simulation, Work in Progress")
			l.DrawLatex(0.40,.82,"#sqrt{s}=8TeV")

			call(["mkdir","plots"])			
			C1.SaveAs("plots/"+baseName+".pdf")
			C2.SaveAs("plots/gen_"+baseName+".pdf")



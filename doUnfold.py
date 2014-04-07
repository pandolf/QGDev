#!/usr/bin/python
import sys,os
import array
import time
import math
from optparse import OptionParser

DEBUG=1

if(DEBUG>0):print "----- BEGIN -----"

if(DEBUG>0):print "-PARSING OPTIONS-"

if __name__=="__main__":
	usage = "usage: %prog [options] arg1 arg2"
	parser=OptionParser(usage=usage)
	parser.add_option("-f","--fileName" ,dest='fileName',type='string',help="File Name with response matrix",default="")
	parser.add_option("-o","--outName" ,dest='outName',type='string',help="OutputFileName. \n\t default=%default",default="UnfoldResults.root")
	#parser.add_option("-l","--libRooUnfold" ,dest='libRooUnfold',type='string',help="Shared RooUnfoldLibrary",default="/afs/cern.ch/user/a/amarini/work/RooUnfold-1.1.1/libRooUnfold.so")
	parser.add_option("-l","--libRooUnfold" ,dest='libRooUnfold',type='string',help="Shared RooUnfoldLibrary. \n\r\t default=%default \n\r\t slc5=/afs/cern.ch/user/a/amarini/work/RooUnfold-1.1.1/libRooUnfold.so",default="/afs/cern.ch/user/a/amarini/work/RooUnfoldSLC6/libRooUnfold.so")
	(options,args)=parser.parse_args()


import ROOT
ROOT.gROOT.SetBatch()

if(DEBUG>0): print "--> Load RooUnfold Library"
ROOT.gSystem.Load(options.libRooUnfold)

def Unfold(R,G,M,H,par,typ="BAYES"):
	''' Macro that performs unfolding:
	\n\t R=Reco hist (MC)
	\n\t G= Gen hist (MC)
	\n\t M=Matrix (MC)
	\n\t H=TargetDistribution Reco (Data/MC)
	\n\t par= reg parameter
	\n\t typ= reg typ'''
	if R==None : print "Ignoring bkg subtraction"
	if G==None : print "Ignoring eff."
	Response= ROOT.RooUnfoldResponse(R,G,M,"Response","Response")

	if typ.lower() == "svd":
                U=ROOT.RooUnfoldSvd(Response,H,par,1000)
        elif typ.lower() == "invert" or typ.lower() == "inversion":
                U=ROOT.RooUnfoldInvert(Response,H)
        elif typ.lower() == "bayes":
                U=ROOT.RooUnfoldBayes(Response,H,par,False) # last= smooth

	#use 1000
	U.SetNToys(10)
        u=U.Hreco(ROOT.RooUnfold.kCovToy)
        c=U.Ereco(ROOT.RooUnfold.kCovToy)
        return (u,c)
	

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
	fOut=ROOT.TFile.Open(options.outName,"RECREATE")
	for v in variables:
	   for e in etaRegions:	
	      for ptBin in range(0,nPtBins):
			#gen=fRoot.Get("nNeutral_fwd_pt1000_3500_h1g")
			baseName=v+"_"+e+"_pt%d_%d"%(ptBins[ptBin],ptBins[ptBin+1])
			gen=fRoot.Get(baseName+"_h1g")		
			reco=fRoot.Get(baseName+"_h1r")
			matrix=fRoot.Get(baseName)
			###def Unfold(R,G,M,H,par,typ="BAYES"):
			fOut.mkdir(v)
			fOut.mkdir(v+"/"+e)
			fOut.cd(v+"/"+e)
			h=Unfold(reco,gen,matrix,reco,3,"BAYES")[0]
			h.SetName("unfold"+baseName)
			h.Write()
			gen.Clone().Write()
			reco.Clone().Write()
			matrix.Clone().Write()


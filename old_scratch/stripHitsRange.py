import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle
from array import array
import csv
import numpy






def median(lst):
    return numpy.median(numpy.array(lst))

def fetchDTSegmentsPhi(event,ontime=True,twinMux=True):
    phiSeg    = Handle  ('L1MuDTChambPhContainer')
    if twinMux:
        event.getByLabel('simTwinMuxDigis',phiSeg)
    else:
        event.getByLabel('simDtTriggerPrimitiveDigis',phiSeg)

    #take only hits from wheel=2
    prefiltered = filter(lambda x: abs(x.whNum())==2 and x.stNum()!=4,phiSeg.product().getContainer())
    if ontime:
        filtered=filter(lambda x: x.bxNum()==0, prefiltered)
        return filtered
    else:
        return prefiltered


def cscSegmentsPhi(event,ontime=True,twinMux=True):

    phiSeg = Handle('std::vector<std::pair<CSCDetId,CSCCorrelatedLCTDigi> >')
    event.getByLabel('cscTriggerHits',phiSeg)
    digis = phiSeg.product()
    return digis


def fetchDTGEANT(event):
    geantH  = Handle  ('vector<PSimHit>')
    event.getByLabel('g4SimHits:MuonDTHits',geantH)
    geant=filter(lambda x: x.pabs()>0.5 and abs(x.particleType())==13,geantH.product())
    return geant

def fetchCSCGEANT(event):
    geantH  = Handle  ('vector<PSimHit>')
    event.getByLabel('g4SimHits:MuonCSCHits',geantH)
    geant=filter(lambda x: x.pabs()>0.5 and abs(x.particleType())==13,geantH.product())
    return geant


class CSCSeg(object):
    def __init__(self,s):
        self.seg=s
        self.truePhiB=-1
    def __getattr__(self,name):
        return getattr(self.seg,name)

def matchedCSCStubs(muon,segments,geant):
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
    chambers=[]
    for p in thisMuonGEANT:
        detid=ROOT.CSCDetId(p.detUnitId())
        chambers.append(p.detUnitId())

    assocSeg=[]
    for s in segments:
        minLayer=8;
        maxLayer=0;
        phiMin=0;
        phiMax=0;
        associated=False
        for c in thisMuonGEANT:
            detid=ROOT.CSCDetId(c.detUnitId())
            if detid.layer()<minLayer:
                minLayer=detid.layer()
                phiMin=c.entryPoint().phi().value()
            if detid.layer()>maxLayer:
                maxLayer=detid.layer()
                phiMax=c.entryPoint().phi().value()
            if detid.endcap()==s.first.endcap() and detid.station()==s.first.station() and detid.ring()==s.first.ring() and detid.chamber()==s.first.chamber():
                associated=True
        if associated:
            seg = CSCSeg(s)
            seg.truePhiB = phiMax-phiMin
            assocSeg.append(seg)
    return assocSeg



def fetchGEN(event,etaMax=2.4):
    genH  = Handle  ('vector<reco::GenParticle>')
    event.getByLabel('genParticles',genH)
    genMuons=filter(lambda x: abs(x.pdgId())==13 and x.status()==1 and abs(x.eta())<etaMax,genH.product())
    return genMuons

def segINT(seg,f1=1,f2=1):
    return seg.phi()*f1,seg.phiB()*f2


def qPTInt(qPT,bits=14):
    lsb = lsBIT(bits)
    floatbinary = int(math.floor(abs(qPT)/lsb))
    return int((qPT/abs(qPT))*floatbinary)

def lsBIT(bits=14):
    maximum=1.25
    lsb = 1.25/pow(2,bits-1)
    return lsb


def matchedDTStubs(muon,segments,geant):
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
    chambers=[]
    for p in thisMuonGEANT:
        detid=ROOT.DTChamberId(p.detUnitId())
        chambers.append(p.detUnitId())

    chambers=list(set(chambers))


    assocSeg=[]
    for s in segments:
        for c in chambers:
            detid=ROOT.DTChamberId(c)
            if s.whNum()==detid.wheel() and s.stNum()==detid.station() and s.scNum()==detid.sector()-1:
                if not (s in assocSeg):
                    assocSeg.append(s)

    return assocSeg



#dictionary that links the (encap,station,ring) to the correct getPhiGlobal inputs
#ec_s_r_dict = {(1,1,3):(0,36,64), (1,1,2):(0,36,80),(1,2,2):(0,36,80),(1,3,2):(0,36,80),(1,4,2):(0,36,80)}


#events = Events(['file:/scratch3/MTF/data/singleMuonWithCSC.root'])
events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

striphits_1_1_3 = ROOT.TH1D("Striphits_1_1_3","Striphits_1_1_3",100,0,200)
striphits_1_1_2 = ROOT.TH1D("Striphits_1_1_2","Striphits_1_1_2",100,0,200)
striphits_1_2_1 = ROOT.TH1D("Striphits_1_2_1","Striphits_1_2_1",100,0,200)
striphits_1_2_2 = ROOT.TH1D("Striphits_1_2_2","Striphits_1_2_2",100,0,200)
striphits_1_3_1 = ROOT.TH1D("Striphits_1_3_1","Striphits_1_3_1",100,0,200)
striphits_1_3_2 = ROOT.TH1D("Striphits_1_3_2","Striphits_1_3_2",100,0,200)
striphits_1_4_1 = ROOT.TH1D("Striphits_1_4_1","Striphits_1_4_1",100,0,200)
striphits_1_4_2 = ROOT.TH1D("Striphits_1_4_2","Striphits_1_4_2",100,0,200)



N=0

for event in events:
    N=N+1
    if N==100000:
        break;
    #real muons
    genMuons=fetchGEN(event)
    segmentsPhiDT=fetchDTSegmentsPhi(event)
    geantDT=fetchDTGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)
    for s in segmentsCSC:
	ec_s_r = [s.first.endcap(),s.first.station(),s.first.ring()]
	num_chamber = float(s.first.chamber())
	num_half_strip = float(s.second.getStrip())
	if ec_s_r == [1,1,3]:
	    striphits_1_1_3.Fill(num_half_strip)
	elif ec_s_r == [1,1,2]:
	    striphits_1_1_2.Fill(num_half_strip)
	elif ec_s_r == [1,2,1]:
	    striphits_1_2_1.Fill(num_half_strip)
	elif ec_s_r == [1,2,2]:
	    striphits_1_2_2.Fill(num_half_strip)
	elif ec_s_r == [1,3,1]:
	    striphits_1_3_1.Fill(num_half_strip)
	elif ec_s_r == [1,3,2]:
	    striphits_1_3_2.Fill(num_half_strip)
	elif ec_s_r == [1,4,1]:
	    striphits_1_4_1.Fill(num_half_strip)
	elif ec_s_r == [1,4,2]:
	    striphits_1_4_2.Fill(num_half_strip)




#RESULTS: saved on desktop/research/1/23/19... ALL except 1_1_3 have 160 half strips.

f=ROOT.TFile("plots_v1.root","RECREATE")
f.cd()
striphits_1_1_3.Write()
striphits_1_1_2.Write()
striphits_1_2_1.Write()
striphits_1_2_2.Write()
striphits_1_3_1.Write()
striphits_1_3_2.Write()
striphits_1_4_1.Write()
striphits_1_4_2.Write() 
f.Close()
    
#    for g in genMuons:
#        # drop real muons that do not land in our detectors
#        if abs(g.eta())<0.9 or abs(g.eta())>1.2:
#            continue



#        dtMatchedSegments=matchedDTStubs(g,segmentsPhIDT,geantDT)

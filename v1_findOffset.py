import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle
from array import array

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

def matchedCSCStubs(muon,segments,geant):
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
    chambers=[]
    for p in thisMuonGEANT:        
        detid=ROOT.CSCDetId(p.detUnitId())
        chambers.append(p.detUnitId())
        
    chambers=list(set(chambers))
    assocSeg=[]   
    for s in segments:
        for c in chambers:
            detid=ROOT.CSCDetId(c)
            if detid.endcap()==s.first.endcap() and detid.station()==s.first.station() and detid.ring()==s.first.ring() and detid.chamber()==s.first.chamber():
                if not (s in assocSeg):
                    assocSeg.append(s)

    return assocSeg



def fetchGEN(event,etaMax=1.2):
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


#dictionary that links the (encap,station,ring) to the corresponding [offset, chambers, strips_per_chamber] 
ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,2):[0,36,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80]}


def getPhi(s):
    '''takes in hit and returns global phi between -pi and pi radians'''
    if (s.first.endcap(),s.first.station(),s.first.ring()) in ec_s_r_dict:
        offset, chambers, strips_per_chamber = ec_s_r_dict[s.first.endcap(), s.first.station(),s.first.ring()]
        num_chamber = float(s.first.chamber())
        num_half_strip = float(s.second.getStrip())
        delta_phi_chamber = 360/chambers
        half_strips_per_chamber = strips_per_chamber*2
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
        phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
        return math.radians(phi) - math.pi



hist_1_1_3 = ROOT.TH2D("offset_1_1_3","Offset_1_1_3; Curvature; Difference in Phi",100,-2,2,100,-2,2)
hist_1_1_2 = ROOT.TH2D("offset_1_1_2","Offset_1_1_2; Curvature; Difference in Phi",100,-2,2,100,-2,2)
hist_1_2_2 = ROOT.TH2D("offset_1_2_2","Offset_1_2_2; Curvature; Difference in Phi",100,-2,2,100,-2,2)
hist_1_3_2 = ROOT.TH2D("offset_1_3_2","Offset_1_3_2; Curvature; Difference in Phi",100,-2,2,100,-2,2)
hist_1_4_2 = ROOT.TH2D("offset_1_4_2","Offset_1_4_2; Curvature; Difference in Phi",100,-2,2,100,-2,2)


events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])



N=0
#make histograms for each (encap,station,ring) to find individual offsets
#x-axis is curvature (charge/pT); y-axis is phi(calculated by me)-phi(system) 

for event in events:
    N=N+1
    #if N==10000000:
    #    break;


    #real muons
    genMuons=fetchGEN(event)
    segmentsPhiDT=fetchDTSegmentsPhi(event)
    geantDT=fetchDTGEANT(event)
    geantCSC=fetchCSCGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)

    for g in genMuons:
        # drop real muons that do not land in our detectors
        if abs(g.eta())>0.9 and abs(g.eta())<1.2:
	#find the matched CSC Segments on this muon
	    cscMatchedSegments = matchedCSCStubs(g,segmentsCSC,geantCSC) 
            #loop on matched segments
            #print 'New generated muon',g.pt(),g.eta(),g.phi()
            for s in cscMatchedSegments:
                phitrack = g.phi()
                phi = getPhi(s)
                if phi:
                    if [s.first.endcap(),s.first.station(),s.first.ring()] == [1,1,3]:
                    	hist_1_1_3.Fill((g.charge()/g.pt()),phi-phitrack)
                    elif [s.first.endcap(),s.first.station(),s.first.ring()] == [1,1,2]:
                    	hist_1_1_2.Fill((g.charge()/g.pt()),phi-phitrack)
                    elif [s.first.endcap(),s.first.station(),s.first.ring()] == [1,2,2]:
                    	hist_1_2_2.Fill((g.charge()/g.pt()),phi-phitrack)
                    elif [s.first.endcap(),s.first.station(),s.first.ring()] == [1,3,2]:
                    	hist_1_3_2.Fill((g.charge()/g.pt()),phi-phitrack)
                    elif [s.first.endcap(),s.first.station(),s.first.ring()] == [1,4,2]:
                    	hist_1_4_2.Fill((g.charge()/g.pt()),phi-phitrack)

f=ROOT.TFile("plots.root","RECREATE")
f.cd()
hist_1_1_3.Write()
hist_1_1_2.Write()
hist_1_2_2.Write()
hist_1_3_2.Write()
hist_1_4_2.Write()
#hist_1_1_3.FitSlicesY()
#hist_1_1_3_1.Fit("pol1")
#hist_1_1_3.Draw()
#hist_1_1_3_1.Draw("same")
#hist_1_1_3.Draw()
f.Close()






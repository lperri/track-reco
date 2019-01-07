import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle
from array import array
ROOT.FWLiteEnabler.enable()
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
ec_s_r_dict = {(1,1,3):[0,36,64], (1,1,2):[0,36,80],(1,2,2):[0,36,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80]}

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

events = Events(['file:/scratch3/MTF/data/singleMuonWithCSC.root'])
N=0


phi_global_vals = []
phi_hist = ROOT.TH1D("GlobalPhi","GlobalPhi",100,-2*math.pi,2*math.pi)

for event in events:
    N=N+1
    if N==10000:
        break;
    #real muons
    genMuons=fetchGEN(event)
    segmentsPhiDT=fetchDTSegmentsPhi(event)
    geantDT=fetchDTGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)
    for s in segmentsCSC:
	phi = getPhi(s)
        if phi:
            phi_global_vals.append(phi)
            phi_hist.Fill(phi)

f=ROOT.TFile("plots.root","RECREATE")
f.cd()
phi_hist.Write()
f.Close()
# note for tues-- getting (1,1,4) as a key its looking for i.e. (1,1,4) get hit -- endcap1,station1,ring4?? DNE?

    
   
#    for g in genMuons:
#        # drop real muons that do not land in our detectors
#        if abs(g.eta())<0.9 or abs(g.eta())>1.2:
#            continue



#        dtMatchedSegments=matchedDTStubs(g,segmentsPhIDT,geantDT)

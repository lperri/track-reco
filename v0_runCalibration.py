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


def getPhiGlobal(offset, chambers, delta_phi_chamber, strips_per_chamber, num_chamber, num_half_strip):
    ''' offset: offset angle, 
	chambers: number of chambers in ring,
	delta_phi_chamber: the angle subtended by a chamber (10 for rings with 36 chambers and 20 for rings with 18 chambers),
	strips_per_chamber: the number of cathode strips per chamber
    	num_chamber: chamber number of hit
	num_half_strip: half strip number of hit
    '''
    half_strips_per_chamber = strips_per_chamber*2
    delta_phi_half_strip = float(delta_phi_chamber/num_half_strip)
    return offset + (chambers - 1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip)


events = Events(['file:/scratch3/MTF/data/singleMuonWithCSC.root'])

#hist_chambers_1_1_3 = ROOT.TH1D("Chambers_1_1_3","Chambers_1_1_3",100,1,37)
#hist_chambers_1_1_2 = ROOT.TH1D("Chambers_1_1_2","Chambers_1_1_2",100,1,37)
#hist_chambers_1_2_2 = ROOT.TH1D("Chambers_1_2_2","Chambers_1_2_2",100,1,37)
#hist_chambers_1_3_2 = ROOT.TH1D("Chambers_1_3_2","Chambers_1_3_2",100,1,37)
#hist_chambers_1_4_2 = ROOT.TH1D("Chambers_1_4_2","Chambers_1_4_2",100,1,37)

#hist_chambers_2_1_3 = ROOT.TH1D("Chambers_2_1_3","Chambers_2_1_3",100,1,37)
#hist_chambers_2_1_2 = ROOT.TH1D("Chambers_2_1_2","Chambers_2_1_2",100,1,37)
#hist_chambers_2_2_2 = ROOT.TH1D("Chambers_2_2_2","Chambers_2_2_2",100,1,37)
#hist_chambers_2_3_2 = ROOT.TH1D("Chambers_2_3_2","Chambers_2_3_2",100,1,37)
#hist_chambers_2_4_2 = ROOT.TH1D("Chambers_2_4_2","Chambers_2_4_2",100,1,37)












N=0

stripnum = []
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
	ec_s_r = [s.first.endcap(),s.first.station(),s.first.ring()]
	num_chamber = float(s.first.chamber())
	num_half_strip = float(s.second.getStrip())
	
	if ec_s_r == [1,1,2]:
	    #hist_chambers_1_1_2.Fill(chamber)
	if ec_s_r == [1,1,3]:
	    #hist_chambers_1_1_3.Fill(chamber)
	if ec_s_r == [1,2,2]:
	    #hist_chambers_1_2_2.Fill(chamber)
	if ec_s_r == [1,3,2]:
	    #hist_chambers_1_3_2.Fill(chamber)
	if ec_s_r == [1,4,2]:
	    #hist_chambers_1_4_2.Fill(chamber)
	if ec_s_r == [2,1,3]:
	    #hist_chambers_2_1_3.Fill(chamber)
	if ec_s_r == [2,1,2]:
	    #hist_chambers_2_1_2.Fill(chamber)
	if ec_s_r == [2,2,2]:
	    #hist_chambers_2_2_2.Fill(chamber)
	if ec_s_r == [2,3,2]:
	    #hist_chambers_2_3_2.Fill(chamber)
	if ec_s_r == [2,4,2]:
	    #hist_chambers_2_4_2.Fill(chamber)
    

f=ROOT.TFile("plots.root","RECREATE")
f.cd()
hist_chambers_1_1_2.Write()
hist_chambers_1_1_3.Write()
hist_chambers_1_2_2.Write()
hist_chambers_1_3_2.Write()
hist_chambers_1_4_2.Write()
hist_chambers_2_1_3.Write()
hist_chambers_2_1_2.Write()
hist_chambers_2_2_2.Write()
hist_chambers_2_3_2.Write()
hist_chambers_2_4_2.Write()
f.Close()
    
#    for g in genMuons:
#        # drop real muons that do not land in our detectors
#        if abs(g.eta())<0.9 or abs(g.eta())>1.2:
#            continue



#        dtMatchedSegments=matchedDTStubs(g,segmentsPhIDT,geantDT)

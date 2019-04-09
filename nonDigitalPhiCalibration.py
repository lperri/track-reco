import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle
from array import array
import csv
import numpy
import subprocess



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
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==-13) or ((muon.charge()<0) and x.particleType()==13),geant)
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


#dictionary that links the (encap,station,ring) to the corresponding [offset, chambers, strips_per_chamber] 
ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,1):[0,18,80],(1,2,2):[0,36,80],(1,3,1):[0,18,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80],(2,1,3):[0,36,64],(2,1,2):[0,36,80],(2,2,1):[0,18,80],(2,2,2):[0,36,80],(2,3,1):[0,18,80],(2,3,2):[0,36,80],(2,4,2):[0,36,80],(1,4,1):[0,18,80],(2,4,1):[0,18,80]}




def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
    if (ec,station,ring) in ec_s_r_dict:
        offset, chambers, strips_per_chamber = ec_s_r_dict[ec,station,ring]
        num_chamber = float(stub.first.chamber())
        num_half_strip = float(stub.second.getStrip())
        delta_phi_chamber = float((2*math.pi)/chambers)
        half_strips_per_chamber = float(strips_per_chamber)*2
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
        if (ec,station) == (2,1) or (ec,station) == (2,2) or (ec,station) == (1,3) or (ec,station) == (1,4):
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + ((half_strips_per_chamber+1-num_half_strip)*delta_phi_half_strip))
        elif (ec,station) == (2,3) or (ec,station) == (2,4) or (ec,station) == (1,1) or (ec,station) == (1,2):
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
        while phi<(-1*math.pi):
            phi += 2*math.pi
        while phi>math.pi:
            phi -= 2*math.pi
            
        return phi

def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def cosTheta(theta):
    return math.cos(theta)


events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0

hists = {}
#ec_s_r_eta gives cutoffs at each five regions in eta-space that the chambers occupy -- starts at [minimum_region1,minimum_region2...]
ec_s_r_eta = {(1,1,2):[1.2,1.3,1.4,1.5,1.6],(1,1,3):[0.9,0.95,1.0,1.05,1.10],(1,2,1):[1.6,1.78,1.96,2.14,2.32],(1,2,2):[1.0,1.12,1.24,1.36,1.48],(1,3,1):[1.7,1,86,2.02,2.18,2.34],(1,3,2):[1.1,1.22,1.34,1.46,1.58],(1,4,2):[1.15,1.28,1.41,1.54,1.67],(2,1,3):[-1.1,-1.05,-1.0,-0.95,-0.9],(2,1,2):[-1.6,-1.5,-1.4,-1.3,-1.2],(2,2,1):[-2.32,-2.14,-1.96,-1.78,-1.6],(2,2,2):[-1.48,-1.36,-1.24,-1.12,-1.0],(2,3,1):[-2.34,-2.18,-2.02,-1.86,-1.7],(2,3,2):[-1.58,-1.46,-1.34,-1.22,-1.1],(2,4,2):[-1.67,-1.54,-1.41,-1.28,-1.15],(1,4,1):[1.8,1.94,2.08,2.22,2.36],(2,4,1):[-2.36,-2.22,-2.08,-1.94,-1.80]}
ec_s = [(1,1),(1,2),(1,3),(1,4),(-1,1),(-1,2),(-1,3),(-1,4)]

for combo in ec_s_r_dict:
    for eta_region in [1,2,3,4,5]:
        name = 'phi_'+str(combo)[1]+str(combo)[4]+str(combo)[7]+'etaRegion'+str(eta_region)
        hists[(combo,eta_region)] = ROOT.TH2D(name,name+';curvature; phi_fp - gen_phi',100,-0.4,0.4,100,-1.4,1.4)

#note upper eta means actully lower on the chamber... so the lower half of the chamber is upper eta

for event in events:
    N=N+1
    if N in [50000,100000,150000,200000]:
	print ('analyzing event {}').format(N)
    if N > 1000:
        break;
    #real muons
    genMuons=fetchGEN(event)
    segmentsPhiDT=fetchDTSegmentsPhi(event)
    geantDT=fetchDTGEANT(event)
    geantCSC=fetchCSCGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)

    for g in genMuons:
        # drop real muons that do not land in our detectors
        if abs(g.eta())>0.9 and abs(g.eta())<2.4:
	#find the matched CSC Segments on this muon
	    cscMatchedSegments = matchedCSCStubs(g,segmentsCSC,geantCSC)
            #loop on matched segments
	    for s in cscMatchedSegments:
               	ec, station, ring = (s.first.endcap(),s.first.station(),s.first.ring())
      		if (ec,station,ring) in ec_s_r_dict:
#		if (ec,station,ring) == (2,2,2):
       		    k = float(g.charge()/g.pt())
		    cos_theta = cosTheta(getTheta(g.eta()))
		    diff_phi = (getPhi(s)) - g.phi()
		    kcos = k/(cos_theta)
    	            while diff_phi < (-1*math.pi):
   	    	        diff_phi += 2*math.pi
    	    	    while diff_phi > math.pi:
    	    	        diff_phi -= 2*math.pi	
                    if g.eta()>= ec_s_r_eta[(ec,station,ring)][0] and g.eta() <= ec_s_r_eta[(ec,station,ring)][1]:
                        hists[((ec,station,ring),1)].Fill(kcos,diff_phi)
                    elif g.eta() >= ec_s_r_eta[(ec,station,ring)][1] and g.eta() <= ec_s_r_eta[(ec,station,ring)][2]:
                        hists[((ec,station,ring),2)].Fill(kcos,diff_phi)
                    elif g.eta() >= ec_s_r_eta[(ec,station,ring)][2] and g.eta() <= ec_s_r_eta[(ec,station,ring)][3]:
                        hists[((ec,station,ring),3)].Fill(kcos,diff_phi)
                    elif g.eta() >= ec_s_r_eta[(ec,station,ring)][3] and g.eta() <= ec_s_r_eta[(ec,station,ring)][4]:
                        hists[((ec,station,ring),4)].Fill(kcos,diff_phi)
                    elif g.eta() >= ec_s_r_eta[(ec,station,ring)][4]:
                        hists[((ec,station,ring),5)].Fill(kcos,diff_phi)
                     


g=ROOT.TFile("offset_phi_split_eta_in5___fixed_fitting.root","RECREATE")
g.cd()

for (combo,eta_region) in hists:
    hists[(combo,eta_region)].Write()
    obj_pfx = hists[(combo,eta_region)].ProfileX()
    fit_func = ROOT.TF1('fit_func','x*[0]/(1+[1]*abs(x)) + [2]',-0.4,0.4)
    print combo,eta_region
    obj_pfx.Fit('fit_func','M')
    # check for non-convering fits
#    while "CALL LIMIT" in ROOT.gMinuit.fCstatu.Data():
#        print 'refitting'
    import pdb;pdb.set_trace()
    #    obj_pfx.Fit('fit_func','M')
    #obj_pfx.Write()
g.Close()


#
#

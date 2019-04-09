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
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==-13) or ((muon.charge()<0) and x.particleType()==13),geant)
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


ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}


def wireCalib(x,y,z,eta):
    return x+ y*eta + z*(eta**2)

def resolution(u,v,k): 
    return math.sqrt(u**2+(v**2)*k**2) 

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0

wire_diff_212 = ROOT.TH1D("wire_212_without_res","wire_212_without_res",8,-4,4)
wire_diff_113 = ROOT.TH1D("wire_113_without_res","wire_113_without_res",8,-4,4)
wire_diff_121 = ROOT.TH1D("wire_121_without_res","wire_121_without_res",8,-4,4)
wire_diff_122 = ROOT.TH1D("wire_122_without_res","wire_122_without_res",8,-4,4)
wire_diff_131 = ROOT.TH1D("wire_131_without_res","wire_131_without_res",8,-4,4)
wire_diff_132 = ROOT.TH1D("wire_132_without_res","wire_132_without_res",8,-4,4)
wire_diff_142 = ROOT.TH1D("wire_142_without_res","wire_142_without_res",8,-4,4)
wire_diff_213 = ROOT.TH1D("wire_213_without_res","wire_213_without_res",8,-4,4)
wire_diff_221 = ROOT.TH1D("wire_221_without_res","wire_221_without_res",8,-4,4)
wire_diff_222 = ROOT.TH1D("wire_222_without_res","wire_222_without_res",8,-4,4)
wire_diff_231 = ROOT.TH1D("wire_231_without_res","wire_231_without_res",8,-4,4)
wire_diff_232 = ROOT.TH1D("wire_232_without_res","wire_232_without_res",8,-4,4)
wire_diff_242 = ROOT.TH1D("wire_242_without_res","wire_242_without_res",8,-4,4)
wire_diff_141 = ROOT.TH1D("wire_141_without_res","wire_141_without_res",8,-4,4)
wire_diff_112 = ROOT.TH1D("wire_112_without_res","wire_112_without_res",8,-4,4)
wire_diff_241 = ROOT.TH1D("wire_241_without_res","wire_241_without_res",8,-4,4)
 
dict_of_hist = {(1,1,2):wire_diff_112,(1,1,3):wire_diff_113,(1,2,2):wire_diff_122,(1,3,2):wire_diff_132,(1,4,2):wire_diff_142,(2,1,2):wire_diff_212,(2,1,3):wire_diff_213,(2,2,2):wire_diff_222,(2,3,2):wire_diff_232,(2,4,2):wire_diff_242,(2,2,1):wire_diff_221,(1,2,1):wire_diff_121,(1,3,1):wire_diff_131,(1,4,1):wire_diff_141,(2,4,1):wire_diff_241,(2,3,1):wire_diff_231}

dict_wire_calib = {(1, 2, 1): [499.9, -303.3, 39.59], (2, 1, 3): [273.5, 364.9, 106.9], (2, 3, 2): [336.7, 334.3, 79.78], (1, 3, 2): [332.3, -327.8, 77.39], (2, 2, 2): [312, 330.2, 83.64], (1, 1, 3): [285, -388, 118], (1, 4, 2): [339.5, -315.5, 69.73], (1, 4, 1): [812.3, -564.8, 94.43], (1, 3, 1): [524.1, -324.8, 44.81], (1, 1, 2): [197.2, -86.25, -21.13], (1, 2, 2): [310, -327, 82.4], (2, 1, 2): [360.4, 316.2, 59.59], (2, 4, 1): [803, 555.9, 92.31], (2, 4, 2): [349, 328.5, 74.13], (2, 3, 1): [473.1, 276.8, 33.59], (2, 2, 1): [470, 270.9, 30.84]}

dict_wire_rms = {(2, 1, 2): [3.78391, 14.105], (2, 1, 3): [1.3635, 16.8427], (2, 2, 2): [1.85818, 16.8152], (2, 3, 2): [1.7054, 16.2489], (1, 3, 2): [1.69797, 17.1057], (1, 2, 1): [4.2155, 18.7755], (1, 2, 2): [1.85454, 16.7872], (1, 4, 1): [3.21621, 21.5728], (1, 3, 1): [3.32301, 17.47], (2, 4, 2): [1.65894, 16.0229], (2, 2, 1): [4.41787, 14.9874], (1, 4, 2): [1.61513, 17.1709], (1, 1, 3): [1.40545, 16.968], (2, 4, 1): [3.25701, 22.0568], (1, 1, 2): [3.7663, 12.4721], (2, 3, 1): [3.48749, 15.1532]}

for event in events:
    N=N+1
#    if N > 50000:
#	break
    if N in [50000,100000,150000,200000]:
	print ('analyzing event {}').format(N)
   # if N > 100000:
   # 	break
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
                    k = float(g.charge()/g.pt())
                    eta_tracker = g.eta()
                    wire_num = s.second.getKeyWG()
                    
                    a,b,c = dict_wire_calib[(ec,station,ring)] 
                    u,v = dict_wire_rms[(ec,station,ring)]
                    res = resolution(u,v,k) 
                    #eta_calib = etaCalib(a,b,wire_num)
                    #eta_diff = eta_calib-eta_tracker

                    wire_calib = wireCalib(a,b,c,eta_tracker)
                    wire_diff = int(wire_calib) -  wire_num

                    wire_diff_res = wire_diff/res
#                    print 'eta_tracker',eta_tracker
#                    print 'wire_num',wire_num
#                    print 'wire_calib',wire_calib
#                    print 'wire_diff',wire_diff
#                    import pdb;pdb.set_trace()
                    (dict_of_hist[(ec,station,ring)]).Fill(wire_diff_res) 
                

f=ROOT.TFile("pull_wire.root","RECREATE")
f.cd()
for combo in dict_of_hist:
    (dict_of_hist[combo]).Write()
f.Close()
#
#

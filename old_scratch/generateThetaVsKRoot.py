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

tan_ab_dict = {(1,1,3):[0.74,0.0072],(1,1,2):[0.40,0.0038],(1,2,1):[0.18,0.0020],(1,2,2):[0.44,0.0061],(1,3,1):[0.18,0.0018],(1,3,2):[0.39,0.0054],(1,4,2):[0.36,0.0049],(2,1,3):[-0.74,-0.0072],(2,1,2):[-0.40,-0.0038],(2,2,1):[-0.18,-0.0020],(2,2,2):[-0.44,-0.0061],(2,3,1):[-0.18,-0.0018],(2,3,2):[-0.39,-0.0054],(2,4,2):[-0.36,-0.0050],(1,4,1):[0.19,0.0015],(2,4,1):[-0.19,-0.0015]}

theta_rms_dict = {(1,1,2):[0.0185,0.06],(1,1,3):[0.0108,0.09],(1,2,1):[0.0244,0.04],(1,2,2):[0.0143,0.08],(1,3,1):[0.0171,0.0],(1,3,2):[0.0117,0.08],(1,4,1):[0.0118,0.0],(1,4,2):[0.0109,0.08],(2,1,2):[0.0181,0.07],(2,1,3):[0.0104,0.11],(2,2,1):[0.0249,0.04],(2,2,2):[0.0133,0.09],(2,3,1):[0.0171,0],(2,3,2):[0.0114,0.09],(2,4,1):[0.0118,0.02],(2,4,2):[0.0113,0.07]}

def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def tanTheta(theta):
    return math.tan(theta)

def cosTheta(theta):
    return math.cos(theta)

def tanThetaCalib(x,y,wire_num):
    return x+y*wire_num

def graphRMSofProfileY(hist_2d):
    ''' must be 2d hist '''
    graph_error = ROOT.TGraphErrors()
    for i in range(1,hist_2d.GetNbinsX()+1):
        proj_y = hist_2d.ProjectionY("",i,i)
	k = (hist_2d.GetXaxis()).GetBinCenter(i)
	rms = proj_y.GetRMS()
	rms_error = proj_y.GetRMSError()
	graph_error.SetPoint(i,k,rms)
	graph_error.SetPointError(i,0,rms_error)	
    return graph_error



def resolution(u,v,k):
    return math.sqrt(u**2+(v**2)*k**2)

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0
theta_diff_112 = ROOT.TH2D("theta_diff_112","theta_diff_112",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_113 = ROOT.TH2D("theta_diff_113","theta_diff_113",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_122 = ROOT.TH2D("theta_diff_122","theta_diff_122",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_132 = ROOT.TH2D("theta_diff_132","theta_diff_132",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_142 = ROOT.TH2D("theta_diff_142","theta_diff_142",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_212 = ROOT.TH2D("theta_diff_212","theta_diff_212",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_213 = ROOT.TH2D("theta_diff_213","theta_diff_213",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_222 = ROOT.TH2D("theta_diff_222","theta_diff_222",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_232 = ROOT.TH2D("theta_diff_232","theta_diff_232",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_242 = ROOT.TH2D("theta_diff_242","theta_diff_242",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_221 = ROOT.TH2D("theta_diff_221","theta_diff_221",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_121 = ROOT.TH2D("theta_diff_121","theta_diff_121",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_131 = ROOT.TH2D("theta_diff_131","theta_diff_131",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_141 = ROOT.TH2D("theta_diff_141","theta_diff_141",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_231 = ROOT.TH2D("theta_diff_231","theta_diff_231",25,-0.4,0.4,25,-0.1,0.1)
theta_diff_241 = ROOT.TH2D("theta_diff_241","theta_diff_241",25,-0.4,0.4,25,-0.1,0.1)

dict_of_hist = {(1,1,2):theta_diff_112,(1,1,3):theta_diff_113,(1,2,2):theta_diff_122,(1,3,2):theta_diff_132,(1,4,2):theta_diff_142,(2,1,2):theta_diff_212,(2,1,3):theta_diff_213,(2,2,2):theta_diff_222,(2,3,2):theta_diff_232,(2,4,2):theta_diff_242,(2,2,1):theta_diff_221,(1,2,1):theta_diff_121,(1,3,1):theta_diff_131,(1,4,1):theta_diff_141,(2,4,1):theta_diff_241,(2,3,1):theta_diff_231}


for event in events:
    N=N+1
#    if N > 50:
#	break
    if N%5==0:
        print ('analyzing event {}').format(N)
#    if N > 10000:
#	break
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
                    tan_theta_tracker = tanTheta(getTheta(g.eta()))
                    wire_num = s.second.getKeyWG()
                    a,b = tan_ab_dict[(ec,station,ring)]
                    
                    tan_theta_calib = tanThetaCalib(a,b,wire_num)
        
                    (dict_of_hist[(ec,station,ring)]).Fill(k,(tan_theta_calib-tan_theta_tracker))



g=ROOT.TFile("test__theta_diff_k.root","RECREATE")
g.cd()
for combo in dict_of_hist:
    dict_of_hist[combo].Write()
g.Close()


#
#

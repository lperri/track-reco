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

#tan_ab_dict = {(1,1,3):[0.74,0.0072],(1,1,2):[0.40,0.0038],(1,2,1):[0.18,0.0020],(1,2,2):[0.44,0.0061],(1,3,1):[0.18,0.0018],(1,3,2):[0.39,0.0054],(1,4,2):[0.36,0.0049],(2,1,3):[-0.74,-0.0072],(2,1,2):[-0.40,-0.0038],(2,2,1):[-0.18,-0.0020],(2,2,2):[-0.44,-0.0061],(2,3,1):[-0.18,-0.0018],(2,3,2):[-0.39,-0.0054],(2,4,2):[-0.36,-0.0050],(1,4,1):[0.19,0.0015],(2,4,1):[-0.19,-0.0015]}

#trying to include more digits
#tan_ab_dict = {(1,1,3):[0.742239,0.00722698],(1,1,2):[0.399156,0.00378522],(1,2,1):[0.17798,0.00195661],(1,2,2):[0.442669,0.00610899],(1,3,1):[0.180084,0.00182782],(1,3,2):[0.391205,0.00539775],(1,4,2):[0.357142,0.00492794],(2,1,3):[-0.741665,-0.00724789],(2,1,2):[-0.398662,-0.00379369],(2,2,1):[-0.178575,-0.00194616],(2,2,2):[-0.44244,-0.00611859],(2,3,1):[-0.180227,-0.00182506],(2,3,2):[-0.391448,-0.00539299],(2,4,2):[-0.357112,-0.00492924],(1,4,1):[0.186378,0.0014822],(2,4,1):[-0.18632,-0.00148197]}

tan_ab_dict_quad = {(1,1,3):[0,0,0],(1,1,2):[0.406769,0.00275293,1.79252e-05],(1,2,1):[0,0,0],(1,2,2):[0,0,0],(1,3,1):[0,0,0],(1,3,2):[0,0,0],(1,4,2):[0,0,0],(2,1,3):[0,0,0],(2,1,2):[-0.405766,-0.00280487,-0.0000173429],(2,2,1):[0,0,0],(2,2,2):[0,0,0],(2,3,1):[0,0,0],(2,3,2):[0,0,0],(2,4,2):[0,0,0],(1,4,1):[0,0,0],(2,4,1):[0,0,0]}


# this one caused issues
#theta_rms_dict = {(1,1,2):[0.0185,0.06],(1,1,3):[0.0108,0.09],(1,2,1):[0.0244,0.04],(1,2,2):[0.0143,0.08],(1,3,1):[0.0171,0.0],(1,3,2):[0.0117,0.08],(1,4,1):[0.0118,0.0],(1,4,2):[0.0109,0.08],(2,1,2):[0.0181,0.07],(2,1,3):[0.0104,0.11],(2,2,1):[0.0249,0.04],(2,2,2):[0.0133,0.09],(2,3,1):[0.0171,0],(2,3,2):[0.0114,0.09],(2,4,1):[0.0118,0.02],(2,4,2):[0.0113,0.07]}

# this one is with the added k^4 term
theta_rms_dict = {  (1,1,2):[float("1.15504e-02"),float("7.56813e-02"),float("-3.00578e-07")],
                    (1,1,3):[float("6.60352e-03"),float("1.07894e-01"),float("8.30937e-01")],
                    (1,2,1):[float("6.94269e-03"),float("2.47499e-02"),float("1.17663e-01")],
                    (1,2,2):[float("1.05986e-02"),float("7.92184e-02"),float("2.49972e-01")],
                    (1,3,1):[float("5.82480e-03"),float("1.38042e-02"),float("1.35796e-01")],
                    (1,3,2):[float("8.38573e-03"),float("6.72039e-02"),float("2.49933e-01")],
                    (1,4,2):[float("7.24900e-03"),float("6.37295e-02"),float("2.50591e-01")],
                    (2,1,3):[float("6.70201e-03"),float("1.08154e-01"),float("7.48731e-01")],
                    (2,1,2):[float("1.25786e-02"),float("4.73764e-02"),float("1.67428e-01")],
                    (2,2,1):[float("6.88419e-03"),float("2.96877e-02"),float("1.05547e-01")],
                    (2,2,2):[float("1.05741e-02"),float("7.77776e-02"),float("2.59509e-01")],
                    (2,3,1):[float("5.76880e-03"),float("2.23054e-02"),float("1.27920e-01")],
                    (2,3,2):[float("8.46642e-03"),float("7.01040e-02"),float("8.59636e-02")],
                    (2,4,2):[float("7.06954e-03"),float("6.18149e-02"),float("2.34171e-01")],
                    (1,4,1):[float("4.07216e-03"),float("2.96945e-02"),float("1.02614e-01")],
                    (2,4,1):[float("4.02844e-03"),float("3.28505e-02"),float("6.93297e-02")]  }


def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def tanTheta(theta):
    return math.tan(theta)

def cosTheta(theta):
    return math.cos(theta)

def tanThetaCalib(x,y,z,wire_num):
    return x+y*wire_num + z*(wire_num**2)

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

tantheta_k_112 = ROOT.TH2D("tantheta_k_121","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
tantheta_k_212 = ROOT.TH2D("tantheta_k_221","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)


for event in events:
    N=N+1
#    if N > 50000:
#	break
    if N in [50000,100000,150000,200000]:
	print ('analyzing event {}').format(N)
#    if N > 100000:
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
#  	        if (ec,station,ring) in ec_s_r_dict:
                if (ec,station,ring) == (2,1,2) or (ec,station,ring) == (1,1,2):	            
                    k = float(g.charge()/g.pt())
                    tan_theta_tracker = tanTheta(getTheta(g.eta()))
                    wire_num = s.second.getKeyWG()
                    a,b,c = tan_ab_dict_quad[(ec,station,ring)]
                    
                    tan_theta_calib = tanThetaCalib(a,b,c,wire_num)
        
                    diff_tan_theta = tan_theta_calib - tan_theta_tracker
                    u,v,x = theta_rms_dict[(ec,station,ring)]
                    res = resolution(u,v,k)
  	            if (ec,station,ring) == (2,1,2):
                        tantheta_k_212.Fill(k,diff_tan_theta)
  	            elif (ec,station,ring) == (1,1,2):
                        tantheta_k_112.Fill(k,diff_tan_theta)

tantheta_k_error_112 = graphRMSofProfileY(tantheta_k_112)
tantheta_k_error_112.SetNameTitle("rms_112","rms_112")

tantheta_k_error_212 = graphRMSofProfileY(tantheta_k_212)
tantheta_k_error_212.SetNameTitle("rms_212","rms_212")

f=ROOT.TFile("testing_212_and_112_tantheta_k_quad.root","RECREATE")
f.cd()
tantheta_k_112.Write()
tantheta_k_212.Write()
f.Close()

g=ROOT.TFile("testing_212_and_112_tantheta_k_quad_RMS.root","RECREATE")
tantheta_k_error_112.Write()
tantheta_k_error_212.Write()
g.Close()

#
#

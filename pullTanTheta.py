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

tan_ab_dict = {(1,1,3):[0.74,0.0072],(1,1,2):[0.40,0.0038],(1,2,1):[0.18,0.0200],(1,2,2):[0.44,0.0061],(1,3,1):[0.18,0.0018],(1,3,2):[0.39,0.0054],(1,4,2):[0.36,0.0049],(2,1,3):[-0.74,-0.0072],(2,1,2):[-0.40,-0.0038],(2,2,1):[-0.18,-0.0020],(2,2,2):[-0.44,-0.0061],(2,3,1):[-0.18,-0.0018],(2,3,2):[-0.39,-0.0054],(2,4,2):[-0.36,-0.0050],(1,4,1):[0.19,0.0015],(2,4,1):[-0.19,-0.0015]}

theta_rms_dict = {(1,1,2):[0.0185,0.06],(1,1,3):[0.0108,0.09],(1,2,1):[0.0244,0.04],(1,2,2):[0.0143,0.08],(1,3,1):[0.0171,0.0],(1,3,2):[0.0117,0.08],(1,4,1):[0.0118,0.0],(1,4,2):[0.0109,0.08],(2,1,2):[0.181,0.07],(2,1,3):[0.0104,0.11],(2,2,1):[0.0249,0.04],(2,2,2):[0.0133,0.09],(2,3,1):[0.0171,0],(2,3,2):[0.0114,0.09],(2,4,1):[0.0118,0.02],(2,4,2):[0.0113,0.07]}

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
#possible_paths = [((2,4,2),(2,3,2)),((2,4,2),(2,2,2)),((2,4,2),(2,1,2)),((2,4,2),(2,2,1)),((2,3,2),(2,2,2)),((2,3,2),(2,1,3)),((2,3,2),(2,1,2)),((2,3,2),(2,2,1)),((2,2,2),(2,1,3)),((2,2,2),(2,1,2)),((1,4,2),(1,3,2)),((1,4,2),(1,2,2)),((1,4,2),(1,1,2)),((1,4,2),(1,2,1)),((1,3,2),(1,2,2)),((1,3,2),(1,1,3)),((1,3,2),(1,1,2)),((1,3,2),(1,2,1)),((1,2,2),(1,1,3)),((1,2,2),(1,1,2)),((1,4,1),(1,3,1)),((1,4,1),(1,2,1)),((2,4,1),(2,3,1)),((2,4,1),(2,2,1))]

pull_theta_112 = ROOT.TH1D("pull_phi_112","pull_phi_112",100,-math.pi,math.pi)
pull_theta_113 = ROOT.TH1D("pull_phi_113","pull_phi_113",100,-math.pi,math.pi)
pull_theta_122 = ROOT.TH1D("pull_phi_122","pull_phi_122",100,-math.pi,math.pi)
pull_theta_132 = ROOT.TH1D("pull_phi_132","pull_phi_132",100,-math.pi,math.pi)
pull_theta_142 = ROOT.TH1D("pull_phi_142","pull_phi_142",100,-math.pi,math.pi)
pull_theta_212 = ROOT.TH1D("pull_phi_212","pull_phi_212",100,-math.pi,math.pi)
pull_theta_213 = ROOT.TH1D("pull_phi_213","pull_phi_213",100,-math.pi,math.pi)
pull_theta_222 = ROOT.TH1D("pull_phi_222","pull_phi_222",100,-math.pi,math.pi)
pull_theta_232 = ROOT.TH1D("pull_phi_232","pull_phi_232",100,-math.pi,math.pi)
pull_theta_242 = ROOT.TH1D("pull_phi_242","pull_phi_242",100,-math.pi,math.pi)
pull_theta_221 = ROOT.TH1D("pull_phi_221","pull_phi_221",100,-math.pi,math.pi)
pull_theta_121 = ROOT.TH1D("pull_phi_121","pull_phi_121",100,-math.pi,math.pi)
pull_theta_131 = ROOT.TH1D("pull_phi_131","pull_phi_131",100,-math.pi,math.pi)
pull_theta_141 = ROOT.TH1D("pull_phi_141","pull_phi_141",100,-math.pi,math.pi)
pull_theta_231 = ROOT.TH1D("pull_phi_231","pull_phi_231",100,-math.pi,math.pi)
pull_theta_241 = ROOT.TH1D("pull_phi_241","pull_phi_241",100,-math.pi,math.pi)


for event in events:
    N=N+1
#    if N > 10000:
#	break
    if N in [50000,100000,150000,200000]:
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
                    tantheta_diff = tan_theta_calib - tan_theta_tracker

  	            u,v = theta_rms_dict[(ec,station,ring)]
  	            
  	            res = resolution(u,v,k)
                    tantheta_diff_res = (tantheta_diff)/res

                    if (ec,station,ring) == (2,4,2):
  	                pull_theta_242.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,2,2):
  	                pull_theta_222.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,1,2):
  	                pull_theta_212.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,2,1):
  	                pull_theta_221.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,3,2):
  	                pull_theta_232.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,1,3):
  	                pull_theta_213.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,4,1):
  	                pull_theta_241.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (2,3,1):
  	                pull_theta_231.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,4,2):
  	                pull_theta_142.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,2,2):
  	                pull_theta_122.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,1,2):
  	                pull_theta_112.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,2,1):
  	                pull_theta_121.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,3,2):
  	                pull_theta_132.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,1,3):
  	                pull_theta_113.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,4,1):
  	                pull_theta_141.Fill(tantheta_diff_res)
  	            elif (ec,station,ring) == (1,3,1):
  	                pull_theta_131.Fill(tantheta_diff_res)


g=ROOT.TFile("theta_pull.root","RECREATE")
g.cd()
pull_theta_221.Write()
pull_theta_231.Write()
pull_theta_241.Write()
pull_theta_212.Write()
pull_theta_213.Write()
pull_theta_222.Write()
pull_theta_232.Write()
pull_theta_242.Write()
pull_theta_121.Write()
pull_theta_131.Write()
pull_theta_141.Write()
pull_theta_112.Write()
pull_theta_113.Write()
pull_theta_122.Write()
pull_theta_132.Write()
pull_theta_142.Write()
g.Close()


#
#

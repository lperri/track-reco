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


#post-endcap switch!!! looks correct:
ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}

tan_ab_dict = {(1,1,3):[0.74,0.0072],(1,1,2):[0.40,0.0038],(1,2,1):[0.18,0.0020],(1,2,2):[0.44,0.0061],(1,3,1):[0.18,0.0018],(1,3,2):[0.39,0.0054],(1,4,2):[0.36,0.0049],(2,1,3):[-0.74,-0.0072],(2,1,2):[-0.40,-0.0038],(2,2,1):[-0.18,-0.0020],(2,2,2):[-0.44,-0.0061],(2,3,1):[-0.18,-0.0018],(2,3,2):[-0.39,-0.0054],(2,4,2):[-0.36,-0.0050],(1,4,1):[0.19,0.0015],(2,4,1):[-0.19,-0.0015]}

def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def tanTheta(theta):
    return math.tan(theta)

def tanThetaCalib(a,b,wire_num):
    return a+b*wire_num

def cosTheta(theta):
    return math.cos(theta)


def graphRMSofProfileY(hist_2d):
    ''' must be 2d hist '''
    graph_error = ROOT.TGraphErrors()
    graph_error.SetTitle(("RMS for {}").format(hist_2d))
    for i in range(1,hist_2d.GetNbinsX()+1):
        proj_y = hist_2d.ProjectionY("",i,i)
	k = (hist_2d.GetXaxis()).GetBinCenter(i)
	rms = proj_y.GetRMS()
	rms_error = proj_y.GetRMSError()
	graph_error.SetPoint(i,k,rms)
	graph_error.SetPointError(i,0,rms_error)	
    return graph_error


events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0

#tantheta_k_112 = ROOT.TH2D("tantheta_k_112","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_113 = ROOT.TH2D("tantheta_k_113","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
tantheta_k_121 = ROOT.TH2D("tantheta_k_121","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_122 = ROOT.TH2D("tantheta_k_122","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_131 = ROOT.TH2D("tantheta_k_131","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_132 = ROOT.TH2D("tantheta_k_132","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_141 = ROOT.TH2D("tantheta_k_141","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_142 = ROOT.TH2D("tantheta_k_142","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)


#tantheta_k_212 = ROOT.TH2D("tantheta_k_212","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_213 = ROOT.TH2D("tantheta_k_213","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
tantheta_k_221 = ROOT.TH2D("tantheta_k_221","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_222 = ROOT.TH2D("tantheta_k_222","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_231 = ROOT.TH2D("tantheta_k_231","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_232 = ROOT.TH2D("tantheta_k_232","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_241 = ROOT.TH2D("tantheta_k_241","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)
#tantheta_k_242 = ROOT.TH2D("tantheta_k_242","$\Delta$$\Tantheta$ vs  k; k; $\Delta\Theta$",100,-0.4,0.4,100,-0.5,0.5)


for event in events:
    N=N+1
    if N in [100,1000,50000,100000,150000,200000]:
	print ('analyzing event {}').format(N)
#    if N>10000:
#    	break    
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
      		#if (ec,station,ring) in ec_s_r_dict:
		if (ec,station,ring) == (1,2,1) or (ec,station,ring) == (2,2,1):
       		    k = float(g.charge()/g.pt())
                    tan_theta_tracker = tanTheta(getTheta(g.eta()))
                    wire_num = s.second.getKeyWG()
                    a,b = tan_ab_dict[(ec,station,ring)]
                    tan_theta_calib = tanThetaCalib(a,b,wire_num)
                    diff_tan_theta = tan_theta_calib - tan_theta_tracker
#                    if (ec,station,ring) == (1,1,2):
#			tantheta_k_112.Fill(k,diff_tan_theta)       			
#       	    	    if (ec,station,ring) == (1,1,3):
#			tantheta_k_113.Fill(k,diff_tan_theta)       			
       	    	    if (ec,station,ring) == (1,2,1):
			tantheta_k_121.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (1,2,2):
#			tantheta_k_122.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (1,3,1):
#			tantheta_k_131.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (1,3,2):
#			tantheta_k_132.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (1,4,1):
#			tantheta_k_141.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (1,4,2):
#			tantheta_k_142.Fill(k,diff_tan_theta)

#       	    	    if (ec,station,ring) == (2,1,2):
#			tantheta_k_212.Fill(k,diff_tan_theta)       			
#       	    	    if (ec,station,ring) == (2,1,3):
#			tantheta_k_213.Fill(k,diff_tan_theta)       			
       	    	    if (ec,station,ring) == (2,2,1):
			tantheta_k_221.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (2,2,2):
#			tantheta_k_222.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (2,3,1):
#			tantheta_k_231.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (2,3,2):
#			tantheta_k_232.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (2,4,1):
#			tantheta_k_241.Fill(k,diff_tan_theta)
#       	    	    if (ec,station,ring) == (2,4,2):
#			tantheta_k_242.Fill(k,diff_tan_theta)


print 'computing RMS...'
#tantheta_k_error_112 = graphRMSofProfileY(tantheta_k_112)
#tantheta_k_error_112.SetNameTitle("rms_112","rms_112")
#tantheta_k_error_113 = graphRMSofProfileY(tantheta_k_113)
#tantheta_k_error_113.SetNameTitle("rms_113","rms_113")
tantheta_k_error_121 = graphRMSofProfileY(tantheta_k_121)
tantheta_k_error_121.SetNameTitle("rms_121","rms_121")
#tantheta_k_error_122 = graphRMSofProfileY(tantheta_k_122)
#tantheta_k_error_122.SetNameTitle("rms_122","rms_122")
#tantheta_k_error_131 = graphRMSofProfileY(tantheta_k_131)
#tantheta_k_error_131.SetNameTitle("rms_131","rms_131")
#tantheta_k_error_132 = graphRMSofProfileY(tantheta_k_132)
#tantheta_k_error_132.SetNameTitle("rms_132","rms_132")
#tantheta_k_error_141 = graphRMSofProfileY(tantheta_k_141)
#tantheta_k_error_141.SetNameTitle("rms_141","rms_141")
#tantheta_k_error_142 = graphRMSofProfileY(tantheta_k_142)
#tantheta_k_error_142.SetNameTitle("rms_142","rms_142")
#tantheta_k_error_212 = graphRMSofProfileY(tantheta_k_212)
#tantheta_k_error_212.SetNameTitle("rms_212","rms_212")
#tantheta_k_error_213 = graphRMSofProfileY(tantheta_k_213)
#tantheta_k_error_213.SetNameTitle("rms_213","rms_213")
tantheta_k_error_221 = graphRMSofProfileY(tantheta_k_221)
tantheta_k_error_221.SetNameTitle("rms_221","rms_221")
#tantheta_k_error_222 = graphRMSofProfileY(tantheta_k_222)
#tantheta_k_error_222.SetNameTitle("rms_222","rms_222")
#tantheta_k_error_231 = graphRMSofProfileY(tantheta_k_231)
#tantheta_k_error_231.SetNameTitle("rms_231","rms_231")
#tantheta_k_error_232 = graphRMSofProfileY(tantheta_k_232)
#tantheta_k_error_232.SetNameTitle("rms_232","rms_232")
#tantheta_k_error_241 = graphRMSofProfileY(tantheta_k_241)
#tantheta_k_error_241.SetNameTitle("rms_241","rms_241")
#tantheta_k_error_242 = graphRMSofProfileY(tantheta_k_242)
#tantheta_k_error_242.SetNameTitle("rms_242","rms_242")

print 'writing files...'
f=ROOT.TFile("141_221_tantheta_k_plots_ecSwitch.root","RECREATE")
f.cd()
#tantheta_k_112.Write()
#tantheta_k_113.Write()
tantheta_k_121.Write()
#tantheta_k_122.Write()
#tantheta_k_131.Write()
#tantheta_k_132.Write()
#tantheta_k_141.Write()
#tantheta_k_142.Write()
#tantheta_k_212.Write()
#tantheta_k_213.Write()
tantheta_k_221.Write()
#tantheta_k_222.Write()
#tantheta_k_231.Write()
#tantheta_k_232.Write()
#tantheta_k_241.Write()
#tantheta_k_242.Write()
f.Close()

g=ROOT.TFile("121_221_tantheta_rms_plots_ecSwitch.root","RECREATE")
g.cd()
#tantheta_k_error_112.Write()
#tantheta_k_error_113.Write()
tantheta_k_error_121.Write()
#tantheta_k_error_122.Write()
#tantheta_k_error_131.Write()
#tantheta_k_error_132.Write()
#tantheta_k_error_141.Write()
#tantheta_k_error_142.Write()
#tantheta_k_error_212.Write()
#tantheta_k_error_213.Write()
tantheta_k_error_221.Write()
#tantheta_k_error_222.Write()
#tantheta_k_error_231.Write()
#tantheta_k_error_232.Write()
#tantheta_k_error_241.Write()
#tantheta_k_error_242.Write()
g.Close()


#
#

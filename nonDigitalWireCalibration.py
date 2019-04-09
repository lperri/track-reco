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


def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def tanTheta(theta):
    return math.tan(theta)

def cosTheta(theta):
    return math.cos(theta)

def wireNumCalib(stub):
    dict_wire_calib = {(1, 2, 1): [499.9, -303.3, 39.59], (2, 1, 3): [273.5, 364.9, 106.9], (2, 3, 2): [336.7, 334.3, 79.78], (1, 3, 2): [332.3, -327.8, 77.39], (2, 2, 2): [312, 330.2, 83.64], (1, 1, 3): [285, -388, 118], (1, 4, 2): [339.5, -315.5, 69.73], (1, 4, 1): [812.3, -564.8, 94.43], (1, 3, 1): [524.1, -324.8, 44.81], (1, 1, 2): [197.2, -86.25, -21.13], (1, 2, 2): [310, -327, 82.4], (2, 1, 2): [360.4, 316.2, 59.59], (2, 4, 1): [803, 555.9, 92.31], (2, 4, 2): [349, 328.5, 74.13], (2, 3, 1): [473.1, 276.8, 33.59], (2, 2, 1): [460, 270.9, 30.84]}
    ec,station,ring = stub.first.endcap(),stub.first.station(),stub.first.ring()
    x,y,w = dict_wire_calib[(ec,station,ring)]
    eta = stub.eta()
    return x+(y*eta)+(w*(eta**2))

def graphRMSofProfileY(hist_2d):
    graph_error = dict_of_hist_rms[hist_2d]
    for i in range(1,hist_2d.GetNbinsX()+1):
        proj_y = hist_2d.ProjectionY("",i,i)
        k = (hist_2d.GetXaxis()).GetBinCenter(i)
        rms = proj_y.GetRMS()
        rms_error = proj_y.GetRMSError()
        graph_error.SetPoint(i,k,rms)
        graph_error.SetPointError(i,0,rms_error)
#        print 'str(hist_2d)',str(hist_2d)
        graph_error.SetNameTitle("rms_"+str(hist_2d)[-29:-16],"rms_"+str(hist_2d)[-29:-16])
    return graph_error
    
events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0

#theta_vs_theta_all = ROOT.TH2D("theta_vs_theta_all","theta_vs_theta_all",80,-1,1,80,-1,1)
wire_diff_112 = ROOT.TH2D("wire_diff_112","wire_diff_112",25,-0.4,0.4,60,-30,30)
wire_diff_113 = ROOT.TH2D("wire_diff_113","wire_diff_113",25,-0.4,0.4,60,-30,30)
wire_diff_122 = ROOT.TH2D("wire_diff_122","wire_diff_122",25,-0.4,0.4,60,-30,30)
wire_diff_132 = ROOT.TH2D("wire_diff_132","wire_diff_132",25,-0.4,0.4,60,-30,30)
wire_diff_142 = ROOT.TH2D("wire_diff_142","wire_diff_142",25,-0.4,0.4,60,-30,30)
wire_diff_212 = ROOT.TH2D("wire_diff_212","wire_diff_212",25,-0.4,0.4,60,-30,30)
wire_diff_213 = ROOT.TH2D("wire_diff_213","wire_diff_213",25,-0.4,0.4,60,-30,30)
wire_diff_222 = ROOT.TH2D("wire_diff_222","wire_diff_222",25,-0.4,0.4,60,-30,30)
wire_diff_232 = ROOT.TH2D("wire_diff_232","wire_diff_232",25,-0.4,0.4,60,-30,30)
wire_diff_242 = ROOT.TH2D("wire_diff_242","wire_diff_242",25,-0.4,0.4,60,-30,30)
wire_diff_221 = ROOT.TH2D("wire_diff_221","wire_diff_221",25,-0.4,0.4,60,-30,30)
wire_diff_121 = ROOT.TH2D("wire_diff_121","wire_diff_121",25,-0.4,0.4,60,-30,30)
wire_diff_131 = ROOT.TH2D("wire_diff_131","wire_diff_131",25,-0.4,0.4,60,-30,30)
wire_diff_141 = ROOT.TH2D("wire_diff_141","wire_diff_141",25,-0.4,0.4,60,-30,30)
wire_diff_231 = ROOT.TH2D("wire_diff_231","wire_diff_231",25,-0.4,0.4,60,-30,30)
wire_diff_241 = ROOT.TH2D("wire_diff_241","wire_diff_241",25,-0.4,0.4,60,-30,30)

rms_wire_112 = ROOT.TGraphErrors()
rms_wire_113 = ROOT.TGraphErrors()
rms_wire_122 = ROOT.TGraphErrors()
rms_wire_132 = ROOT.TGraphErrors()
rms_wire_142 = ROOT.TGraphErrors()
rms_wire_212 = ROOT.TGraphErrors()
rms_wire_213 = ROOT.TGraphErrors()
rms_wire_222 = ROOT.TGraphErrors()
rms_wire_232 = ROOT.TGraphErrors()
rms_wire_242 = ROOT.TGraphErrors()
rms_wire_221 = ROOT.TGraphErrors()
rms_wire_121 = ROOT.TGraphErrors()
rms_wire_131 = ROOT.TGraphErrors()
rms_wire_141 = ROOT.TGraphErrors()
rms_wire_231 = ROOT.TGraphErrors()
rms_wire_241 = ROOT.TGraphErrors()



dict_of_hist = {(1,1,2):wire_diff_112,(1,1,3):wire_diff_113,(1,2,2):wire_diff_122,(1,3,2):wire_diff_132,(1,4,2):wire_diff_142,(2,1,2):wire_diff_212,(2,1,3):wire_diff_213,(2,2,2):wire_diff_222,(2,3,2):wire_diff_232,(2,4,2):wire_diff_242,(2,2,1):wire_diff_221,(1,2,1):wire_diff_121,(1,3,1):wire_diff_131,(1,4,1):wire_diff_141,(2,4,1):wire_diff_241,(2,3,1):wire_diff_231}

#associate TgraphError to each 2d histogram so that can calculate rms using 2d hist
dict_of_hist_rms = {wire_diff_112:rms_wire_112,wire_diff_113:rms_wire_113,wire_diff_122:rms_wire_122,wire_diff_132:rms_wire_132,wire_diff_142:rms_wire_142,wire_diff_212:rms_wire_212,wire_diff_213:rms_wire_213,wire_diff_222:rms_wire_222,wire_diff_232:rms_wire_232,wire_diff_242:rms_wire_242,wire_diff_221:rms_wire_221,wire_diff_121:rms_wire_121,wire_diff_131:rms_wire_131,wire_diff_141:rms_wire_141,wire_diff_241:rms_wire_241,wire_diff_231:rms_wire_231}



for event in events:
    N=N+1
    #if N > 50:
#	break
#    if N%500==0:
    print ('analyzing event {}').format(N)
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
                    gen_eta = g.eta()
                    wire_num = s.second.getKeyWG()
                    a,b,c = dict_wire_calib[(ec,station,ring)]
                    
                    wire_calib = wireNumCalib(a,b,c,gen_eta)
                    wire_diff = int(wire_calib) - wire_num
                    (dict_of_hist[(ec,station,ring)]).Fill(k,(int(wire_calib)-wire_num))

for hist_2d in dict_of_hist_rms:
    dict_of_hist_rms[hist_2d] = graphRMSofProfileY(hist_2d)

#g=ROOT.TFile("wire_diff_k.root","RECREATE")
#g.cd()
#for combo in dict_of_hist:
#    dict_of_hist[combo].Write()
#g.Close()

f = ROOT.TFile("wire_diff_k_rms.root","RECREATE")
f.cd()
for hist_2d in dict_of_hist_rms:
    dict_of_hist_rms[hist_2d].Write()
f.Close()
#
#

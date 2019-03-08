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



def getAll(d,basepath=''):
    "Go into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            for i in getAll(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)


#dictionary that links the (encap,station,ring) to the corresponding [offset, chambers, strips_per_chamber] 
#ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,1):[0,18,80],(1,2,2):[0,36,80],(1,3,1):[0,18,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80],(2,1,3):[0,36,64],(2,1,2):[0,36,80],(2,2,1):[0,18,80],(2,2,2):[0,36,80],(2,3,1):[0,18,80],(2,3,2):[0,36,80],(2,4,2):[0,36,80],(1,4,1):[0,18,80],(2,4,1):[0,18,80]}
#ec_s_r_dict = {(1,1,3):[0.086,36,64],(1,1,2):[0.082,36,80],(1,2,2):[0.086,36,80],(1,3,2):[0.069,36,80],(1,4,2):[0.073,36,80],(2,1,3):[0.095,36,64],(2,1,2):[0.079,36,80],(2,2,2):[0.074,36,80],(2,3,2):[0.087,36,80],(2,4,2):[0.086,36,80]}


#ec_s_r_dict [CORRESPONDING TO 'NEWEST_PLOTS',PRE-ECSWITCH] = {(1,1,3):[0.096,36,64],(1,1,2):[0.094,36,80],(1,2,1):[0.12,18,80],(1,2,2):[0.098,36,80],(1,3,1):[0.091,18,80],(1,3,2):[0.084,36,80],(1,4,2):[0.084,36,80],(2,1,3):[0.089,36,64],(2,1,2):[0.086,36,80],(2,2,1):[0.093,18,80],(2,2,2):[0.074,36,80],(2,3,1):[0.11,18,80],(2,3,2):[0.081,36,80],(2,4,2):[0.081,36,80],(1,4,1):[0.080,18,80],(2,4,1):[0.12,18,80]}

#post-endcap switch!!! looks correct:
ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}



#need to update detector_z_dict to include new stations
#detector_z_dict = {(1,1,3):-7.0,(1,1,2):-7.0,(1,2,2):-8.3,(1,3,2):-9.3,(1,4,2):-10.3,(2,1,3):7.0,(2,1,2):7.0,(2,2,2):8.3,(2,3,2):9.3,(2,4,2):10.3}


def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def tanTheta(theta):
    return math.tan(theta)

def cosTheta(theta):
    return math.cos(theta)


events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0


tan_wire_112 = ROOT.TH2D("theta_wire_112","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_113 = ROOT.TH2D("theta_wire_113","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_121 = ROOT.TH2D("theta_wire_121","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_122 = ROOT.TH2D("theta_wire_122","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_131 = ROOT.TH2D("theta_wire_131","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_132 = ROOT.TH2D("theta_wire_132","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_141 = ROOT.TH2D("theta_wire_141","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)
tan_wire_142 = ROOT.TH2D("theta_wire_142","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,0,1)


tan_wire_212 = ROOT.TH2D("theta_wire_212","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_213 = ROOT.TH2D("theta_wire_213","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_221 = ROOT.TH2D("theta_wire_221","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_222 = ROOT.TH2D("theta_wire_222","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_231 = ROOT.TH2D("theta_wire_231","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_232 = ROOT.TH2D("theta_wire_232","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_241 = ROOT.TH2D("theta_wire_241","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)
tan_wire_242 = ROOT.TH2D("theta_wire_242","Tan$\Theta$ vs wire_num; wire_num; Tan$\Theta$",116,0,116,100,-1,0)


for event in events:
    N=N+1
    if N in [100,1000,50000,100000,150000,200000]:
	print ('analyzing event {}').format(N)
  #  if N>10000:
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
#		if (ec,station,ring) == (2,2,2):
       		    k = float(g.charge()/g.pt())
		    tan_theta = tanTheta(getTheta(g.eta()))
                    wire_num = s.second.getKeyWG()
                    if (ec,station,ring) == (1,1,2):
			tan_wire_112.Fill(wire_num,tan_theta)       			
       	    	    if (ec,station,ring) == (1,1,3):
			tan_wire_113.Fill(wire_num,tan_theta)       			
       	    	    if (ec,station,ring) == (1,2,1):
			tan_wire_121.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (1,2,2):
			tan_wire_122.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (1,3,1):
			tan_wire_131.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (1,3,2):
			tan_wire_132.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (1,4,1):
			tan_wire_141.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (1,4,2):
			tan_wire_142.Fill(wire_num,tan_theta)

       	    	    if (ec,station,ring) == (2,1,2):
			tan_wire_212.Fill(wire_num,tan_theta)       			
       	    	    if (ec,station,ring) == (2,1,3):
			tan_wire_213.Fill(wire_num,tan_theta)       			
       	    	    if (ec,station,ring) == (2,2,1):
			tan_wire_221.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (2,2,2):
			tan_wire_222.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (2,3,1):
			tan_wire_231.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (2,3,2):
			tan_wire_232.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (2,4,1):
			tan_wire_241.Fill(wire_num,tan_theta)
       	    	    if (ec,station,ring) == (2,4,2):
			tan_wire_242.Fill(wire_num,tan_theta)



print 'writing files...'
f=ROOT.TFile("tan_wire_plots_ecSwitch.root","RECREATE")
f.cd()
tan_wire_112.Write()
tan_wire_113.Write()
tan_wire_121.Write()
tan_wire_122.Write()
tan_wire_131.Write()
tan_wire_132.Write()
tan_wire_141.Write()
tan_wire_142.Write()
tan_wire_212.Write()
tan_wire_213.Write()
tan_wire_221.Write()
tan_wire_222.Write()
tan_wire_231.Write()
tan_wire_232.Write()
tan_wire_241.Write()
tan_wire_242.Write()


for name, obj in getAll(f):
    print 'fitting parameters for ',name
    canvas = ROOT.TCanvas(name,name)
    obj.FitSlicesY()
    #with open('output_tan_wire.txt','w') as output_file:
    #    print >> output_file, name
    print 
    print name
    obj.Fit('pol1','S')

        

f.Close()


        

#g=ROOT.TFile("fitted_tan_wire_plots_ecSwitch.root","RECREATE")
#g.cd()
#tan_wire_error_112.Write()
#tan_wire_error_113.Write()
#tan_wire_error_121.Write()
#tan_wire_error_122.Write()
#tan_wire_error_131.Write()
#tan_wire_error_132.Write()
#tan_wire_error_141.Write()
#tan_wire_error_142.Write()
#tan_wire_error_212.Write()
#tan_wire_error_213.Write()
#tan_wire_error_221.Write()
#tan_wire_error_222.Write()
#tan_wire_error_231.Write()
#tan_wire_error_232.Write()
#tan_wire_error_241.Write()
#tan_wire_error_242.Write()
#g.Close()


#
#

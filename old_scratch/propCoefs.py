import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle
from array import array

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

##def matchedCSCStubs(muon,segments,geant):
#    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
#    chambers=[]
#    for p in thisMuonGEANT:        
#        detid=ROOT.CSCDetId(p.detUnitId())
#        chambers.append(p.detUnitId())
#        
#    chambers=list(set(chambers))
#    assocSeg=[]   
#    for s in segments:
#        for c in chambers:
#            detid=ROOT.CSCDetId(c)
#            if detid.endcap()==s.first.endcap() and detid.station()==s.first.station() and detid.ring()==s.first.ring() and detid.chamber()==s.first.chamber():
#                if not (s in assocSeg):
#                    assocSeg.append(s)
#
#    return assocSeg
class CSCSeg(object):
    def __init__(self,s):
        self.seg=s
        self.truePhiB=-1
    def __getattr__(self,name):
        return getattr(self.seg,name)



def matchedCSCStubs(muon,segments,geant):
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
    chambers=[]
    for p in thisMuonGEANT:
        detid=ROOT.CSCDetId(p.detUnitId())
        chambers.append(p.detUnitId())

    assocSeg=[]
    for s in segments:
        minLayer=8;
        maxLayer=0;
        xin=0;

        yin=0;
        xout=0;
        yout=0;

        associated=False
        for c in thisMuonGEANT:
            detid=ROOT.CSCDetId(c.detUnitId())
#            print 'Det=',detid.station(),detid.ring(),detid.chamber(),detid.layer(),'Angle=',c.entryPoint().phi().value(),'entry=',c.entryPoint().x(),c.entryPoint().y(),c.entryPoint().z()
            if detid.layer()<minLayer:
                minLayer=detid.layer()
                xin =c.entryPoint().x()
                yin =c.entryPoint().y()
            if detid.layer()>maxLayer:
                maxLayer=detid.layer()
                xout =c.entryPoint().x()
                yout =c.entryPoint().y()
            if detid.endcap()==s.first.endcap() and detid.station()==s.first.station() and detid.ring()==s.first.ring() and detid.chamber()==s.first.chamber():
                associated=True
        if associated:
            seg = CSCSeg(s)
            seg.truePhiB = math.atan((xout-xin)/(3.0*(maxLayer-minLayer)))
            assocSeg.append(seg)
    return assocSeg



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
ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,2):[0,36,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80],(2,1,3):[0,36,64],(2,1,2):[0,36,80],(2,2,2):[0,36,80],(2,3,2):[0,36,80],(2,4,2):[0,36,80]}

detector_z_dict = {(1,1,3):-7.0,(1,1,2):-7.0,(1,2,2):-8.3,(1,3,2):-9.3,(1,4,2):-10.3,(2,1,3):7.0,(2,1,2):7.0,(2,2,2):8.3,(2,3,2):9.3,(2,4,2):10.3}

#tan_dict: (ec,station,ring):(yint,slope)--yint,slope of tan(theta) (y) vs. wire number (x)
tan_dict = {(1,1,3):(-0.740587,-0.00737487),(1,1,2):(-17.6517,0.301504),(1,2,2):(-0.468672,-0.00565595),(1,3,2):(-0.446212,-0.00449382),(1,4,2):(-0.451912,-0.00352335),(2,1,3):(0.737727,0.0075506),(2,1,2):(0.602858,0.00102775),(2,2,2):(0.460577,0.00586943),(2,3,2):(0.437834,0.00469009),(2,4,2):(0.538971,0.00215101)}

#ab_dict: ((ec2,station2,ring2),(ec1,station1,ring1)):(a,b) where a,b are propogation coefficients found by estimating in aFit & bFit (waiting on a better bending angle from M
ab_dict = {((2,4,2),(2,3,2)):(0,0),((2,4,2),(2,2,2)):(0,0),((2,4,2),(2,1,2)):(0,0),((2,3,2),(2,2,2)):(0,0),((2,3,2),(2,1,3)):(0,0),((2,3,2),(2,1,2)):(0,0),((2,2,2),(2,1,3)):(0,0),((2,2,2),(2,1,2)):(0,0),((1,4,2),(1,3,2)):(0,0),((1,4,2),(1,2,2)):(0,0),((1,4,2),(1,1,2)):(0,0),((1,3,2),(1,2,2)):(0,0),((1,3,2),(1,1,3)):(0,0),((1,3,2),(1,1,2)):(0,0),((1,2,2),(1,1,3)):(0,0),((1,2,2),(1,1,2)):(0,0)}

def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    if (stub.first.endcap(),stub.first.station(),stub.first.ring()) in ec_s_r_dict:
	offset, chambers, strips_per_chamber = ec_s_r_dict[stub.first.endcap(), stub.first.station(),stub.first.ring()]
	num_chamber = float(stub.first.chamber())
        num_half_strip = float(stub.second.getStrip())
        delta_phi_chamber = float(360/chambers)
        half_strips_per_chamber = strips_per_chamber*2
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
        phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
	return math.radians(phi) - math.pi

def getBend(s):
	#dict{id number : number halfstrips from middle}
	dict = {2:-4,3:4,4:-3,5:3,6:-2,7:2,8:-1,9:1,10:0}
	offset, chambers, strips_per_chamber = ec_s_r_dict[s.first.endcap(), s.first.station(),s.first.ring()]
	delta_phi_chamber = float(360/chambers)
	half_strips_per_chamber = float(strips_per_chamber*2)
	delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
	pattern = s.second.getPattern()
	#note this variable deltb_phi_bend is defined wrt the center of the detector
	delta_phi_bend = float(dict[pattern]*math.radians(delta_phi_half_strip))
	z = detector_z_dict[(s.first.endcap(),s.first.station(),s.first.ring())]
	return math.atan((2*z/0.1) * math.tan(delta_phi_bend))

# def getThetaFromEta(eta):
#    ''' takes in value of pseudorapidity and converts to theta'''
#    theta = 2*math.atan(1/math.exp(eta))
#    return theta

#def tanThetaFromEta(theta):
#    ''' take in value of theta and returns tangent of that value in radians '''
#    return math.tan(theta)

def getTanThetaCalib(a,b,wire_num):
    ''' takes in a, b (parameters found in calibration), and wire number; returns value of tan(theta) '''
    return a+b*wire_num 

##def calcPhiBendPrime(z1,z2,bend):
#    '''input (z1,z2,phi) & calculates the bending angle at location of second detector (after first/outwardmost)'''
#    c = abs((z1-z2)/z1)
#    return c*((((z1-z2)**2)/abs(z1)) + 2*(abs(z1-z2))) + bend*(1+(abs(z1-z2)/abs(z1))) 

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0


for event in events:
    N=N+1
    #real muons
    genMuons=fetchGEN(event)
    segmentsPhiDT=fetchDTSegmentsPhi(event)
    geantDT=fetchDTGEANT(event)
    geantCSC=fetchCSCGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)

    for g in genMuons:
        # drop real muons that do not land in our detectors
        if abs(g.eta())>0.9 and abs(g.eta())<1.2:
	#find the matched CSC Segments on this muon
	    cscMatchedSegments = matchedCSCStubs(g,segmentsCSC,geantCSC)
            #loop on matched segments
	    for s in cscMatchedSegments:
         	for q in cscMatchedSegments:
	            if s.first.station() != q.first.station() and s.first.endcap() == q.first.endcap():
			ec2, station2, ring2 = (s.first.endcap(),s.first.station(),s.first.ring())                
                	ec1, station1, ring1 = (q.first.endcap(),q.first.station(),q.first.ring())
			if ((ec2,station2,ring2),(ec1,station1,ring1)) in ab_dict:	
			    k = float(g.charge()/g.pt())
		    	    #get a&b values for pair of hits
			    a,b = ab_dict[(ec2,station2,ring2),(ec1,station1,ring1)]
			    
				 
			   

#f=ROOT.TFile("plots.root","RECREATE")
#f.cd()
#f.Close()




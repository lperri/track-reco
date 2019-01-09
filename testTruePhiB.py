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
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
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
    #note this variable delta_phi_bend is defined wrt the center of the detector
    delta_phi_bend = float(dict[pattern]*math.radians(delta_phi_half_strip))
    z = detector_z_dict[(s.first.endcap(),s.first.station(),s.first.ring())]
    return math.atan((2*z/0.1) * math.tan(delta_phi_bend))

    

def calcPhiBendPrime(z1,z2,bend):
    '''input (z1,z2,phi) & calculates the bending angle at location of second detector (after first/outwardmost)'''
    c = abs((z1-z2)/z1)
    return c*((((z1-z2)**2)/abs(z1)) + 2*(abs(z1-z2))) + bend*(1+(abs(z1-z2)/abs(z1))) 

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0
#histograms of true bend (from Michalis) vs calculated bend (by Leah)
# 222 (segment s not q) because thats the outermost corresponding to phi_b not phi_b_prime
#bendvsk_222_calc = ROOT.TH2D("bend_222_calc","$\Phi_{b}$ (calc) vs. k; k; $\Phi_{b}$ (calc)",100,-0.35,0.35,100,-1,1)
#bendvsk_222_true = ROOT.TH2D("bend_222_true","$\Phi_{b}$ (true) vs. k; k; $\Phi_{b}$ (true)",100,-0.35,0.35,100,-1,1)

muon_number = []
patterns = []
pt_vals = []
curvature = []
bend_calc_vals = []
bend_true_vals = []

muon_number.append('muon number:')
for x in range(1,51):
    muon_number.append(x)
patterns.append('pattern:')
pt_vals.append('pt:')
curvature.append('k:')
bend_calc_vals.append('bend Leah:')
bend_true_vals.append('bend M:')

 
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
			if ((station2,ring2),(station1,ring1)) in [((4,2),(3,2)),((4,2),(2,2)),((4,2),(1,3)),((4,2),(1,2)),((3,2),(2,2)),((3,2),(1,3)),((3,2),(1,2)),((2,2),(1,3)),((2,2),(1,2))]:
			    #z2 = detector_z_dict[(ec2, station2, ring2)]
		    	    #z1 = detector_z_dict[(ec1, station1, ring1)] 
			    k = float(g.charge()/g.pt())
			    #bend is the bending angle at the outer detector that I calculate with my function getBend(s) 
			    #compare bend to the 'true bending angle'
#			    bend_true = q.truePhiB
	            	    #deltaphi = getPhi(q) - getPhi(s)
		    	    #c = abs(z1-z2)/abs(z1)
		    	    if ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,3)):
				bend = getPhi(s) - (math.pi/2) + getBend(s)
			        #bendvsk_222_calc.Fill(k, bend)
			        bend_true = getPhi(s) - (math.pi/2) + s.truePhiB
			        #bendvsk_222_true.Fill(k, bend_true)
			        pattern = s.second.getPattern()
				pt = g.pt()
				if len(pt_vals)<51:
				    bend_calc_vals.append(bend)
				    bend_true_vals.append(bend_true)
				    patterns.append(pattern)
				    pt_vals.append(g.pt())
				    curvature.append(k)
				    				    
with open('testing_fix.txt','w') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerows(zip(*[muon_number,pt_vals,patterns,curvature,bend_calc_vals,bend_true_vals]))



#file = open('testing_fix.txt','w')
#for i in range(len(pt_vals)):
#    file.write(('muon {}: pt={},k={},pattern={},bend(Leah)={},bend(M)={}').format(i,pt_vals[i],curvature[i],patterns[i],bend_calc_vals[i],bend_true_vals[i]))
#file.close()

#f=ROOT.TFile("plots.root","RECREATE")
#f.cd()
#bendvsk_222_calc.Write()
#bendvsk_222_true.Write()
#f.Close()



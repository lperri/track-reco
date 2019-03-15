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
ec_s_r_dict = {(1,1,3):[0.086,36,64],(1,1,2):[0.082,36,80],(1,2,2):[0.086,36,80],(1,3,2):[0.069,36,80],(1,4,2):[0.073,36,80],(2,1,3):[0.095,36,64],(2,1,2):[0.079,36,80],(2,2,2):[0.074,36,80],(2,3,2):[0.087,36,80],(2,4,2):[0.086,36,80]}

detector_z_dict = {(1,1,3):-7.0,(1,1,2):-7.0,(1,2,2):-8.3,(1,3,2):-9.3,(1,4,2):-10.3,(2,1,3):7.0,(2,1,2):7.0,(2,2,2):8.3,(2,3,2):9.3,(2,4,2):10.3}

def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    if (stub.first.endcap(),stub.first.station(),stub.first.ring()) in ec_s_r_dict:
	offset, chambers, strips_per_chamber = ec_s_r_dict[stub.first.endcap(), stub.first.station(),stub.first.ring()]
	num_chamber = float(stub.first.chamber())
        num_half_strip = float(stub.second.getStrip())
        delta_phi_chamber = float((2*math.pi)/chambers)
        half_strips_per_chamber = strips_per_chamber*2
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
        phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
	while phi<(-1*math.pi):
	    phi += 2*math.pi
	while phi>math.pi:
	    phi -= 2*math.pi
	return phi


def getPsi(s):
    #dict{id number : number halfstrips from middle}
    dict = {2:4,3:-4,4:3,5:-3,6:2,7:-2,8:1,9:-1,10:0}
    offset, chambers, strips_per_chamber = ec_s_r_dict[s.first.endcap(), s.first.station(),s.first.ring()]
    delta_phi_chamber = float(2*math.pi/chambers)
    half_strips_per_chamber = float(strips_per_chamber*2)
    delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
    pattern = s.second.getPattern()
    #note this variable delta_phi_psi is defined wrt the center of the detector
    delta_phi_psi = float(dict[pattern]*(delta_phi_half_strip))
    z = detector_z_dict[(s.first.endcap(),s.first.station(),s.first.ring())]
    psi = math.atan((2*z/0.1) * math.tan(delta_phi_psi))
    return psi
    
def phiMiddleChamber(endcap,station,ring,chamber):
#   ''' takes in (endcap,station,ring) as well as the chamber number and gets global phi value to the middle of that chamber '''
    offset,chambers,strips_per_chamber = ec_s_r_dict[endcap,station,ring]
    delta_phi_chamber = float(2*math.pi/chambers)
    half_strips_per_chamber = strips_per_chamber*2
    delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
    middle_strip_num = half_strips_per_chamber/2
    phi_middle = float(offset + delta_phi_chamber*(chamber-1) + (middle_strip_num*delta_phi_half_strip))    
    while phi_middle < (-1*math.pi):
	phi_middle += 2*math.pi
    while phi_middle > math.pi:
	phi_middle -= 2*math.pi
    return phi_middle

def getBend(s,phi_middle_chamber):
    phiprime = getPhi(s) - phi_middle_chamber
    bend = getPsi(s) - phiprime
    while bend < (-1*math.pi):
	bend += 2*math.pi
    while bend > math.pi:
	bend -= 2*math.pi
    return bend
    
def calcPhiBendPrime(z1,z2,bend):
    '''input (z1,z2,phi) & calculates the bending angle at location of second detector (after first/outwardmost)'''
    c = abs((z1-z2)/z1)
    phibp = c*((((z1-z2)**2)/abs(z1)) + 2*(abs(z1-z2))) + bend*(1+(abs(z1-z2)/abs(z1))) 
    while phibp < (-1*math.pi):
	phibp += 2*math.pi
    while phibp > math.pi:
	phibp -= 2*math.pi
    return phibp

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
strip_num = []
wire_num = []
chamber_num = []
phi_vals = []
fixed_bend_l = []
fixed_bend_m = []
muon_number.append('muon number:')
for x in range(1,51):
    muon_number.append(x)
patterns.append('pattern:')
pt_vals.append('pt:')
curvature.append('k:')
strip_num.append('strip number:')
wire_num.append('wire number:')
chamber_num.append('chamber number:')
phi_vals.append('phiL:')
fixed_bend_l.append('fixed bend L:')
fixed_bend_m.append('fixed bend M:')

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
			    #deltaphi = getPhi(q) - getPhi(s)
		    	    #c = abs(z1-z2)/abs(z1)
		    	    if ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,3)):
				phi_middle_chamber = phiMiddleChamber(ec2,station2,ring2,s.first.chamber())
	#			bend = getPhi(s) - (math.pi/2) + getPsi(s)
			        #bendvsk_222_calc.Fill(k, bend)
	#		        bend_true = getPhi(s) - (math.pi/2) + s.truePhiB
			        #bendvsk_222_true.Fill(k, bend_true)
			        pattern = s.second.getPattern()
				pt = g.pt()
				fixed_bend = getBend(s, phi_middle_chamber)
				fixed_bend_true = s.truePhiB	
				while fixed_bend_true < (-1*math.pi):
        			    fixed_bend_true += 2*math.pi
    				while fixed_bend_true > math.pi:
        			    fixed_bend_true -= 2*math.pi
			        if len(pt_vals)<51:
				    patterns.append(pattern)
				    pt_vals.append(g.pt())
				    curvature.append(k)
				    strip_num.append(s.second.getStrip())  				    				    
				    wire_num.append(s.second.getKeyWG())
				    chamber_num.append(s.first.chamber())			
				    phi_vals.append(getPhi(s))
				    fixed_bend_l.append(fixed_bend)
				    fixed_bend_m.append(fixed_bend_true)
	
with open('testing_fix.txt','w') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerows(zip(*[muon_number,pt_vals,patterns,curvature,strip_num,wire_num,chamber_num,phi_vals,fixed_bend_l,fixed_bend_m]))




#f=ROOT.TFile("plots.root","RECREATE")
#f.cd()
#bendvsk_222_calc.Write()
#bendvsk_222_true.Write()
#f.Close()



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


#dictionary that links the (encap,station,ring) to the corresponding [offset, chambers, strips_per_chamber] 
#ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,1):[0,18,80],(1,2,2):[0,36,80],(1,3,1):[0,18,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80],(2,1,3):[0,36,64],(2,1,2):[0,36,80],(2,2,1):[0,18,80],(2,2,2):[0,36,80],(2,3,1):[0,18,80],(2,3,2):[0,36,80],(2,4,2):[0,36,80],(1,4,1):[0,18,80],(2,4,1):[0,18,80]}

ec_s_r_dict = {(1,1,3):[0.086,36,64],(1,1,2):[0.082,36,80],(1,2,2):[0.086,36,80],(1,3,2):[0.069,36,80],(1,4,2):[0.073,36,80],(2,1,3):[0.095,36,64],(2,1,2):[0.079,36,80],(2,2,2):[0.074,36,80],(2,3,2):[0.087,36,80],(2,4,2):[0.086,36,80]}


#need to update detector_z_dict to include new stations
detector_z_dict = {(1,1,3):-7.0,(1,1,2):-7.0,(1,2,1):-8.3,(1,2,2):-8.3,(1,3,1):-9.3,(1,3,2):-9.3,(1,4,1):-10.3,(1,4,2):-10.3,(2,1,3):7.0,(2,1,2):7.0,(2,2,1):8.3,(2,2,2):8.3,(2,3,1):9.3,(2,3,2):9.3,(2,4,1):10.3,(2,4,2):10.3}

def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    if (stub.first.endcap(),stub.first.station(),stub.first.ring()) in ec_s_r_dict:
	#print 'now were in the getPhi(s) function...'
	offset, chambers, strips_per_chamber = ec_s_r_dict[stub.first.endcap(), stub.first.station(),stub.first.ring()]
	#print 'offset, chambers, strips_per_chamber = ec_s_r_dict[stub.first.endcap(), stub.first.station(),stub.first.ring()] ',offset,chambers,strips_per_chamber
	num_chamber = float(stub.first.chamber())
	#print 'num_chamber = float(stub.first.chamber()) ',num_chamber
        num_half_strip = float(stub.second.getStrip())
	#print 'num_half_strip = float(stub.second.getStrip()) ',num_half_strip
        delta_phi_chamber = float((2*math.pi)/chambers)
	#print 'delta_phi_chamber = float((2*math.pi)/chambers)',delta_phi_chamber
        half_strips_per_chamber = strips_per_chamber*2
	#print 'half_strips_per_chamber ',half_strips_per_chamber
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
	#print 'delta_phi_half_strip ',delta_phi_half_strip
        phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
	#print 'phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip)) ',phi
	while phi<(-1*math.pi):
	    phi += 2*math.pi
	    #print 'adding 2pi becuase phi is less than -pi'
	while phi>math.pi:
	    phi -= 2*math.pi
	    #print 'subtracting 2pi because phi is greater than pi'
	#print 'the final phi is ',phi
	return phi

def getTheta(eta):
    ''' takes in value of pseudorapidity and converts to theta'''
    theta = 2*math.atan(1/math.exp(eta))
    return theta

def tanTheta(theta):
    return math.tan(theta)

def cosTheta(theta):
    return math.cos(theta)

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
phi_k_112 = ROOT.TH2D("phi_k_112","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_113 = ROOT.TH2D("phi_k_113","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_121 = ROOT.TH2D("phi_k_121","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_122 = ROOT.TH2D("phi_k_122","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_131 = ROOT.TH2D("phi_k_131","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_132 = ROOT.TH2D("phi_k_132","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_141 = ROOT.TH2D("phi_k_141","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_142 = ROOT.TH2D("phi_k_142","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)


phi_k_212 = ROOT.TH2D("phi_k_212","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_213 = ROOT.TH2D("phi_k_213","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_221 = ROOT.TH2D("phi_k_221","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_222 = ROOT.TH2D("phi_k_222","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_231 = ROOT.TH2D("phi_k_231","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_232 = ROOT.TH2D("phi_k_232","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_241 = ROOT.TH2D("phi_k_241","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)
phi_k_242 = ROOT.TH2D("phi_k_242","difference in $\phi$ vs ktan; ktan; $\Delta\Phi$",100,-1,1,100,-2,2)

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
        if abs(g.eta())>0.9 and abs(g.eta())<2.4:
	#find the matched CSC Segments on this muon
	    cscMatchedSegments = matchedCSCStubs(g,segmentsCSC,geantCSC)
            #loop on matched segments
	    for s in cscMatchedSegments:
               	ec, station, ring = (s.first.endcap(),s.first.station(),s.first.ring())
       	#	if (ec,station,ring) in ec_s_r_dict:
      		if (ec,station,ring) in ec_s_r_dict:
		    #z2 = detector_z_dict[(ec2, station2, ring2)]
       	    	    #z1 = detector_z_dict[(ec1, station1, ring1)] 
       		    k = float(g.charge()/g.pt())
		    cos_theta = cosTheta(getTheta(g.eta()))
		    diff_phi = (getPhi(s)-math.pi) - g.phi()
		    kcos = k*(cos_theta)
		    #c = abs(z1-z2)/abs(z1)
#		    while diff_phi < (-1*math.pi):
#		       diff_phi += 2*math.pi
#		    while diff_phi > math.pi:
#		       diff_phi -= 2*math.pi
    	            while diff_phi < (-1*math.pi):
   	    	        diff_phi += 2*math.pi
    	    	    while diff_phi > math.pi:
    	    	        diff_phi -= 2*math.pi	
		    #print 'ec,station,ring ',ec,station,ring
		    #print 'g.phi() ',g.phi()
		    #print 'getPhi(s) ',getPhi(s)
		    #print 'diff_phi ',diff_phi
		    #print 'k ',k

       	    	    if (ec,station,ring) == (1,1,2):
			phi_k_112.Fill(kcos,diff_phi)       			
       	    	    if (ec,station,ring) == (1,1,3):
			phi_k_113.Fill(kcos,diff_phi)       			
       	    	    if (ec,station,ring) == (1,2,1):
			phi_k_121.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (1,2,2):
			phi_k_122.Fill(kcos,diff_phi)
			print g.phi(),g.eta(),cos_theta,k,'coord=',getPhi(s),'diff=',diff_phi
                        if abs(k)<0.02 and abs(diff_phi)>0.3:
                            import pdb;pdb.set_trace()
       	    	    if (ec,station,ring) == (1,3,1):
			phi_k_131.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (1,3,2):
			phi_k_132.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (1,4,1):
			phi_k_141.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (1,4,2):
			phi_k_142.Fill(kcos,diff_phi)

       	    	    if (ec,station,ring) == (2,1,2):
			phi_k_212.Fill(kcos,diff_phi)       			
       	    	    if (ec,station,ring) == (2,1,3):
			phi_k_213.Fill(kcos,diff_phi)       			
       	    	    if (ec,station,ring) == (2,2,1):
			phi_k_221.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (2,2,2):
			phi_k_222.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (2,3,1):
			phi_k_231.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (2,3,2):
			phi_k_232.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (2,4,1):
			phi_k_241.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (2,4,2):
			phi_k_242.Fill(kcos,diff_phi)
			#phi_middle_chamber = phiMiddleChamber(ec2,station2,ring2,s.first.chamber())
       #			bend = getPhi(s) - (math.pi/2) + getPsi(s)
       		        #bendvsk_222_calc.Fill(k, bend)
       #		        bend_true = getPhi(s) - (math.pi/2) + s.truePhiB
       		        #bendvsk_222_true.Fill(k, bend_true)
       		        #pattern = s.second.getPattern()
       			#pt = g.pt()
       			#fixed_bend = getBend(s, phi_middle_chamber)
       			#fixed_bend_true = s.truePhiB	
       		#	while fixed_bend_true < (-1*math.pi):
       		#	    fixed_bend_true += 2*math.pi
   		#		while fixed_bend_true > math.pi:
       		#	    fixed_bend_true -= 2*math.pi


f=ROOT.TFile("plots.root","RECREATE")
f.cd()
phi_k_112.Write()
phi_k_113.Write()
phi_k_121.Write()
phi_k_122.Write()
phi_k_131.Write()
phi_k_132.Write()
phi_k_141.Write()
phi_k_142.Write()
phi_k_212.Write()
phi_k_213.Write()
phi_k_221.Write()
phi_k_222.Write()
phi_k_231.Write()
phi_k_232.Write()
phi_k_241.Write()
phi_k_242.Write()
f.Close()
#
#

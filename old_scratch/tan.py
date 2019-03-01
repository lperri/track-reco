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

def matchedCSCStubs(muon,segments,geant):
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
    chambers=[]
    for p in thisMuonGEANT:        
        detid=ROOT.CSCDetId(p.detUnitId())
        chambers.append(p.detUnitId())
        
    chambers=list(set(chambers))
    assocSeg=[]   
    for s in segments:
        for c in chambers:
            detid=ROOT.CSCDetId(c)
            if detid.endcap()==s.first.endcap() and detid.station()==s.first.station() and detid.ring()==s.first.ring() and detid.chamber()==s.first.chamber():
                if not (s in assocSeg):
                    assocSeg.append(s)

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

#tan_yint_slope_dict: (ec,station,ring):(yint,slope)-yint, slope of tan(theta) (y) vs. wire number (x)
tan_yint_slope_dict = {(1,1,3):(-0.740587,-0.00737487),(1,1,2):(-17.6517,0.301504),(1,2,2):(-0.468672,-0.00565595),(1,3,2):(-0.446212,-0.00449382),(1,4,2):(-0.451912,-0.00352335),(2,1,3):(0.737727,0.0075506),(2,1,2):(0.602858,0.00102775),(2,2,2):(0.460577,0.00586943),(2,3,2):(0.437834,0.00469009),(2,4,2):(0.538971,0.00215101)}


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

def getTanThetaCalib(yint,slope,wire_num):
    ''' takes in a, b (parameters found in calibration), and wire number; returns value of tan(theta) '''
    return yint+slope*wire_num 

##def calcPhiBendPrime(z1,z2,bend):
#    '''input (z1,z2,phi) & calculates the bending angle at location of second detector (after first/outwardmost)'''
#    c = abs((z1-z2)/z1)
#    return c*((((z1-z2)**2)/abs(z1)) + 2*(abs(z1-z2))) + bend*(1+(abs(z1-z2)/abs(z1))) 

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0

#tan_222_213 = ROOT.TH1D("a_222_213","($\Delta Tan(\Theta)$) vs. Curvature; $\Kappa$; ($\Delta Tan(\Theta)$)",100,-0.3,0.3,100,-0.3,0.3)
tan_222 = ROOT.TH1D("tan_222","tan 222",100,0.6,0.9)
tan_213 = ROOT.TH1D("tan_213","tan 213",100,0.6,0.9)
#tank_242_232 = ROOT.TH1D("a_242_232","($\Delta Tan(\Theta)$) vs. Curvature; $\Kappa$; ($\Delta Tan(\Theta)$)",100,-0.3,0.3,100,-0.3,0.3)
#tank_232_222 = ROOT.TH1D("a_232_222","($\Delta Tan(\Theta)$) vs. Curvature; $\Kappa$; ($\Delta Tan(\Theta)$)",100,-0.3,0.3,100,-0.3,0.3)
#tank_142_132 = ROOT.TH2D("a_142_132","($\Delta Tan(\Theta)$) vs. Curvature; $\Kappa$; ($\Delta Tan(\Theta)$)",100,-0.3,0.3,100,-0.3,0.3)
#tank_132_122 = ROOT.TH2D("a_132_122","($\Delta Tan(\Theta)$) vs. Curvature; $\Kappa$; ($\Delta Tan(\Theta)$)",100,-0.3,0.3,100,-0.3,0.3)
#tank_122_113 = ROOT.TH2D("a_122_113","($\Delta Tan(\Theta)$) vs. Curvature; $\Kappa$; ($\Delta Tan(\Theta)$)",100,-0.3,0.3,100,-0.3,0.3)

ltan_222 =[]
ltan_213 =[]

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
			#(ec,2,2),(ec,3,2),(ec,1,3) are used b/c we are propogating inward
			ec2, station2, ring2 = (s.first.endcap(),s.first.station(),s.first.ring())                
                	ec1, station1, ring1 = (q.first.endcap(),q.first.station(),q.first.ring())
			if ((station2,ring2),(station1,ring1)) in [((4,2),(3,2)),((3,2),(2,2)),((2,2),(1,3))]:
			    k = float(g.charge()/g.pt())
		    	    a2, b2 = tan_yint_slope_dict[(ec2,station2,ring2)]
			    a1, b1 = tan_yint_slope_dict[(ec1,station1,ring1)]
			    wire2 = s.second.getKeyWG()
			    wire1 = q.second.getKeyWG()    	     
			    tan_theta2 = getTanThetaCalib(a2,b2,wire2)
			    tan_theta1 = getTanThetaCalib(a1,b1,wire1)
			    #delta_tan_theta = tan_theta1 - tan_theta2
			    if ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,3)):
				tan_222.Fill(tan_theta2)
			   	tan_213.Fill(tan_theta1)
				ltan_222.append(tan_theta2)
				ltan_213.append(tan_theta1)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,2,2)):
#				tank_232_222.Fill(k,delta_tan_theta)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,3,2)):
#				tank_242_232.Fill(k,delta_tan_theta)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,2,2),(1,1,3)):
#				tank_122_113.Fill(k,delta_tan_theta)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,2,2)):
#				tank_132_122.Fill(k,delta_tan_theta)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,3,2)):
#				tank_142_132.Fill(k,delta_tan_theta)
 
			   
print "ltan_222",ltan_222[-15:]
print "ltan_213",ltan_213[-15:]
print "max ltan222",max(ltan_222)
print "max ltan213",max(ltan_213)

f=ROOT.TFile("plots.root","RECREATE")
f.cd()
tan_222.Write()
tan_213.Write()
f.Close()




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

def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    if (stub.first.endcap(),stub.first.station(),stub.first.ring()) in ec_s_r_dict:
	offset, chambers, strips_per_chamber = ec_s_r_dict[stub.first.endcap(), stub.first.station(),stub.first.ring()]
	num_chamber = float(stub.first.chamber())
        num_half_strip = float(stub.second.getStrip())
        deltb_phi_chamber = float(360/chambers)
        half_strips_per_chamber = strips_per_chamber*2
        deltb_phi_half_strip = float(deltb_phi_chamber/half_strips_per_chamber)
        phi = float(offset + (num_chamber -1)*(deltb_phi_chamber) + (num_half_strip*deltb_phi_half_strip))
	return math.radians(phi) - math.pi

def getBend(s):
    #dict{id number : number halfstrips from middle}
    dict = {2:-4,3:4,4:-3,5:3,6:-2,7:2,8:-1,9:1,10:0}
    offset, chambers, strips_per_chamber = ec_s_r_dict[s.first.endcap(), s.first.station(),s.first.ring()]
    deltb_phi_chamber = float(360/chambers)
    half_strips_per_chamber = float(strips_per_chamber*2)
    deltb_phi_half_strip = float(deltb_phi_chamber/half_strips_per_chamber)
    pattern = s.second.getPattern()
    #note this variable deltb_phi_bend is defined wrt the center of the detector
    deltb_phi_bend = float(dict[pattern]*math.radians(deltb_phi_half_strip))
    z = detector_z_dict[(s.first.endcap(),s.first.station(),s.first.ring())]
    return math.atan((2*z/0.1) * math.tan(deltb_phi_bend))


def calcPhiBendPrime(z1,z2,bend):
    '''input (z1,z2,phi) & calculates the bending angle at location of second detector (after first/outwardmost)'''
    c = abs((z1-z2)/z1)
    return (c*(((z1-z2)**2)/z1 + 2*(abs(z1-z2))) + bend*(1+abs((z1-z2)/z1)))

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0
#make histograms: hist_242_232 represents a track that goes from ME_2_3_2 to ME_2_4_2 (inverse propogation) -- to fit a
#b_242_232 = ROOT.TH2D("b_242_232","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,0.1,0.3)
#b_242_222 = ROOT.TH2D("b_242_222","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_242_212 = ROOT.TH2D("b_242_212","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_232_222 = ROOT.TH2D("b_232_222","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_232_213 = ROOT.TH2D("b_232_213","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_232_212 = ROOT.TH2D("b_232_212","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_222_213 = ROOT.TH2D("b_222_213","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_222_212 = ROOT.TH2D("b_222_212","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#
#b_142_132 = ROOT.TH2D("b_142_132","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_142_122 = ROOT.TH2D("b_142_122","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_142_112 = ROOT.TH2D("b_142_112","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_132_122 = ROOT.TH2D("b_132_122","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_132_113 = ROOT.TH2D("b_132_113","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_132_112 = ROOT.TH2D("b_132_112","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_122_113 = ROOT.TH2D("b_122_113","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)
#b_122_112 = ROOT.TH2D("b_122_112","($\Phi_{b}'$-d$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Phi_{b}'$-d$\Phi_{b}$)",100,-0.3,0.3,100,-6.5,6.5)

#make histograms to fit b

#b_232_222_y = []
#b_232_222_x = []
#b_222_213_y = []
#b_222_213_x = []
#b_242_232_y = []
#b_242_232_x = []
#b_222_212_x = []
#b_222_212_y = []
#b_232_212_x = []
#b_232_212_y = []
#b_232_213_x = []
#b_232_213_y = []
#b_242_212_x = []
#b_242_212_y = []
#b_242_213_x = []
#b_242_213_y = []
#b_242_222_x = []
#b_242_222_y = []
#
#b_122_113_y = []
#b_122_113_x = []
#b_132_122_y = []
#b_132_122_x = []
#b_142_132_y = []
#b_142_132_x = []
#b_122_112_x = []
#b_122_112_y = []
#b_132_112_x = []
#b_132_112_y = []
#b_132_113_x = []
#b_132_113_y = []
#b_142_112_x = []
#b_142_112_y = []
#b_142_122_x = []
#b_142_122_y = []

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
			if ((station2,ring2),(station1,ring1)) in [((4,2),(3,2)),((4,2),(2,2)),((4,2),(1,3)),((4,2),(1,2)),((3,2),(2,2)),((3,2),(1,3)),((3,2),(1,2)),((2,2),(1,3)),((2,2),(1,2))]:
			    z2 = detector_z_dict[(ec2, station2, ring2)]
		    	    z1 = detector_z_dict[(ec1, station1, ring1)] 
			    k = float(g.charge()/g.pt())
			    bend2 = getBend(s)
		    	    c = abs(z1-z2)/abs(z1)
			    phibendprime = calcPhiBendPrime(z1,z2,bend2) 
			    d = float(1+(abs(z1-z2)/abs(z1)))
			    if ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,3,2)):
				import pdb;pdb.set_trace()
				print 'we have two hits 242->232.'
				print 'k ',k
				print 'phibendprime ',phibendprime
				print 'd ',d
				print 'bend2 ',bend2
				print 'phibendprime - d*bend2 = ',phibendprime-d*bend2
#    			        b_242_232_x.append(k)
#			        b_242_232_y.append(phibendprime - d*bend2)
			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,2,2)):
				import pdb;pdb.set_trace()
                                print 'we have two hits 242->222.'
                                print 'k ',k
                                print 'phibendprime ',phibendprime
                                print 'd ',d
                                print 'bend2 ',bend2
                                print 'phibendprime - d*bend2 = ',phibendprime-d*bend2    			        
#				b_242_222_x.append(k)
#    				b_242_222_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,1,2)):
#    				b_242_212_x.append(k)
#    				b_242_212_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,2,2)):
#    				b_232_222_x.append(k)
#    				b_232_222_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,1,3)):
#    				b_232_213_x.append(k)
#    				b_232_213_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,1,2)):
#    				b_232_212_x.append(k)
#    				b_232_212_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,3)):
#    				b_222_213_x.append(k)
#    				b_222_213_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,2)):
#    				b_222_212_x.append(k)
#    				b_222_212_y.append(phibendprime-d*bend2)
#
#
#
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,3,2)):
#    				b_142_132_x.append(k)
#    				b_142_132_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,2,2)):
#    				b_142_122_x.append(k)
#    				b_142_122_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,1,2)):
#    				b_142_112_x.append(k)
#    				b_142_112_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,2,2)):
#    				b_132_122_x.append(k)
#    				b_132_122_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,1,3)):
#    				b_132_113_x.append(k)
#    				b_132_113_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,1,2)):
#    				b_132_112_x.append(k)
#    				b_132_112_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,2,2),(1,1,3)):
#    				b_122_113_x.append(k)
#    				b_122_113_y.append(phibendprime-d*bend2)
#			    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,2,2),(1,1,2)):
#    				b_122_112_x.append(k)
#    				b_122_112_y.append(phibendprime-d*bend2)
 


#f=ROOT.TFile("plots.root","RECREATE")
#f.cd()
#b_242_232.Write()
#b_232_222.Write()
#b_222_213.Write()
#b_142_132.Write()
#b_132_122.Write()
#b_122_113.Write()
#f.Close()


#get efficient binning
#print 'range b_242_232_y:',[min(b_242_232_y),max(b_242_232_y)]
#print 'range b_242_232_x:',[min(b_242_232_x),max(b_242_232_x)]
#print 'range b_242_222_y:',[min(b_242_222_y),max(b_242_222_y)]
#print 'range b_242_222_x:',[min(b_242_222_x),max(b_242_222_x)]
#print 'range b_242_212_y:',[min(b_242_212_y),max(b_242_212_y)]
#print 'range b_242_212_x:',[min(b_242_212_x),max(b_242_212_x)]
#print 'range b_232_222_y:',[min(b_232_222_y),max(b_232_222_y)]
#print 'range b_232_222_x:',[min(b_232_222_x),max(b_232_222_x)]
#print 'range b_232_213_y:',[min(b_232_213_y),max(b_232_213_y)]
#print 'range b_232_213_x:',[min(b_232_213_x),max(b_232_213_x)]
#print 'range b_232_212_y:',[min(b_232_212_y),max(b_232_212_y)]
#print 'range b_232_212_x:',[min(b_232_212_x),max(b_232_212_x)]
#print 'range b_222_213_y:',[min(b_222_213_y),max(b_222_213_y)]
#print 'range b_222_213_x:',[min(b_222_213_x),max(b_222_213_x)]
#print 'range b_222_212_y:',[min(b_222_212_y),max(b_222_212_y)]
#print 'range b_222_212_x:',[min(b_222_212_x),max(b_222_212_x)]
#
#
#print 'range b_142_132_y:',[min(b_142_132_y),max(b_142_132_y)]
#print 'range b_142_132_x:',[min(b_142_132_x),max(b_142_132_x)]
#print 'range b_142_122_y:',[min(b_142_122_y),max(b_142_122_y)]
#print 'range b_142_122_x:',[min(b_142_122_x),max(b_142_122_x)]
#print 'range b_142_112_y:',[min(b_142_112_y),max(b_142_112_y)]
#print 'range b_142_112_x:',[min(b_142_112_x),max(b_142_112_x)]
#print 'range b_132_122_y:',[min(b_132_122_y),max(b_132_122_y)]
#print 'range b_132_122_x:',[min(b_132_122_x),max(b_132_122_x)]
#print 'range b_132_113_y:',[min(b_132_113_y),max(b_132_113_y)]
#print 'range b_132_113_x:',[min(b_132_113_x),max(b_132_113_x)]
#print 'range b_132_112_y:',[min(b_132_112_y),max(b_132_112_y)]
#print 'range b_132_112_x:',[min(b_132_112_x),max(b_132_112_x)]
#print 'range b_122_113_y:',[min(b_122_113_y),max(b_122_113_y)]
#print 'range b_122_113_x:',[min(b_122_113_x),max(b_122_113_x)]
#print 'range b_122_112_y:',[min(b_122_112_y),max(b_122_112_y)]
#print 'range b_122_112_x:',[min(b_122_112_x),max(b_122_112_x)]

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

    

def calcPhiBendPrime(z1,z2,bend):
    '''input (z1,z2,phi) & calculates the bending angle at location of second detector (after first/outwardmost)'''
    c = abs((z1-z2)/z1)
    return c*((((z1-z2)**2)/abs(z1)) + 2*(abs(z1-z2))) + bend*(1+(abs(z1-z2)/abs(z1))) 

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0
#make histograms: hist_242_232 represents a track that goes from ME_2_3_2 to ME_2_4_2 (inverse propogation) -- to fit a
a_242_232 = ROOT.TH2D("a_242_232","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-6.3,6.3)
a_242_222 = ROOT.TH2D("a_242_222","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-6.3,6.3)
a_242_212 = ROOT.TH2D("a_242_212","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-0.3,0.3)
a_232_222 = ROOT.TH2D("a_232_222","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-6.3,6.3)
a_232_213 = ROOT.TH2D("a_232_213","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.2,0.2,100,-0.3,0.3)
a_232_212 = ROOT.TH2D("a_232_212","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.3,0.3,100,-0.3,0.3)
a_222_213 = ROOT.TH2D("a_222_213","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-0.2,0.2)
a_222_212 = ROOT.TH2D("a_222_212","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-0.4,0.4)

a_142_132 = ROOT.TH2D("a_142_132","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.3,0.3,100,-6.0,6.0)
a_142_122 = ROOT.TH2D("a_142_122","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.3,0.3,100,-6.3,6.3)
a_142_112 = ROOT.TH2D("a_142_112","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.3,0.3,100,-0.3,6.3)
a_132_122 = ROOT.TH2D("a_132_122","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-6.3,6.3)
a_132_113 = ROOT.TH2D("a_132_113","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.2,0.2,100,-0.3,0.3)
a_132_112 = ROOT.TH2D("a_132_112","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-0.3,6.3)
a_122_113 = ROOT.TH2D("a_122_113","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.3,0.3,100,-0.2,0.1)
a_122_112 = ROOT.TH2D("a_122_112","($\Delta\Phi$-c$\Phi_{b}$) vs. Curvature; $\Kappa$; ($\Delta\Phi$-c$\Phi_{b}$)",100,-0.4,0.4,100,-0.2,6.3)







#make histograms to fit a

#a_232_222_y = []
#a_232_222_x = []
#a_222_213_y = []
#a_222_213_x = []
#a_242_232_y = []
#a_242_232_x = []
#a_222_212_x = []
#a_222_212_y = []
#a_232_212_x = []
#a_232_212_y = []
#a_232_213_x = []
#a_232_213_y = []
#a_242_212_x = []
#a_242_212_y = []
#a_242_222_x = []
#a_242_222_y = []
#
#a_122_113_y = []
#a_122_113_x = []
#a_132_122_y = []
#a_132_122_x = []
#a_142_132_y = []
#a_142_132_x = []
#a_122_112_x = []
#a_122_112_y = []
#a_132_112_x = []
#a_132_112_y = []
#a_132_113_x = []
#a_132_113_y = []
#a_142_112_x = []
#a_142_112_y = []
#a_142_122_x = []
#a_142_122_y = []

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
			    z2 = detector_z_dict[(ec2, station2, ring2)]
		    	    z1 = detector_z_dict[(ec1, station1, ring1)] 
			    k = float(g.charge()/g.pt())
			    bend2 = getBend(s)
	            	    deltaphi = getPhi(q) - getPhi(s)
		    	    c = abs(z1-z2)/abs(z1)
		    	    if ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,3,2)):
				a_242_232.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,2,2)):
				a_242_222.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,4,2),(2,1,2)):
				a_242_212.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,2,2)):
				a_232_222.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,1,3)):
				a_232_213.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,3,2),(2,1,2)):
				a_232_212.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,3)):
				a_222_213.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((2,2,2),(2,1,2)):
				a_222_212.Fill(k,deltaphi-c*bend2)
				#import pdb;pdb.set_trace()
				#print "curvature,k",k
				#print "phiBouter",bend2
				#print "delta phi",deltaphi
				#print "phiOUTER",getPhi(s)
				#print "phiInner",getPhi(q)
				#print "c",c

		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,3,2)):
				a_142_132.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,2,2)):
				a_142_122.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,4,2),(1,1,2)):
				a_142_112.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,2,2)):
				a_132_122.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,1,3)):
				a_132_113.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,3,2),(1,1,2)):
				a_132_112.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,2,2),(1,1,3)):
				a_122_113.Fill(k,deltaphi-c*bend2)
		    	    elif ((ec2,station2,ring2),(ec1,station1,ring1)) == ((1,2,2),(1,1,2)):
				a_122_112.Fill(k,deltaphi-c*bend2)





f=ROOT.TFile("plots.root","RECREATE")
f.cd()
a_242_232.Write()
a_242_222.Write()
a_242_212.Write()
a_232_222.Write()
a_232_213.Write()
a_232_212.Write()
a_222_213.Write()
a_222_212.Write()
a_142_132.Write()
a_142_122.Write()
a_142_112.Write()
a_132_122.Write()
a_132_113.Write()
a_132_112.Write()
a_122_113.Write()
a_122_112.Write()
f.Close()


#get efficient binning
#print 'range a_242_232_y:',[min(a_242_232_y),max(a_242_232_y)]
#print 'range a_242_232_x:',[min(a_242_232_x),max(a_242_232_x)]
#print 'range a_242_222_y:',[min(a_242_222_y),max(a_242_222_y)]
#print 'range a_242_222_x:',[min(a_242_222_x),max(a_242_222_x)]
#print 'range a_242_212_y:',[min(a_242_212_y),max(a_242_212_y)]
#print 'range a_242_212_x:',[min(a_242_212_x),max(a_242_212_x)]
#print 'range a_232_222_y:',[min(a_232_222_y),max(a_232_222_y)]
#print 'range a_232_222_x:',[min(a_232_222_x),max(a_232_222_x)]
#print 'range a_232_213_y:',[min(a_232_213_y),max(a_232_213_y)]
#print 'range a_232_213_x:',[min(a_232_213_x),max(a_232_213_x)]
#print 'range a_232_212_y:',[min(a_232_212_y),max(a_232_212_y)]
#print 'range a_232_212_x:',[min(a_232_212_x),max(a_232_212_x)]
#print 'range a_222_213_y:',[min(a_222_213_y),max(a_222_213_y)]
#print 'range a_222_213_x:',[min(a_222_213_x),max(a_222_213_x)]
#print 'range a_222_212_y:',[min(a_222_212_y),max(a_222_212_y)]
#print 'range a_222_212_x:',[min(a_222_212_x),max(a_222_212_x)]
#
#
#
#print 'range a_142_132_y:',[min(a_142_132_y),max(a_142_132_y)]
#print 'range a_142_132_x:',[min(a_142_132_x),max(a_142_132_x)]
#print 'range a_142_122_y:',[min(a_142_122_y),max(a_142_122_y)]
#print 'range a_142_122_x:',[min(a_142_122_x),max(a_142_122_x)]
#print 'range a_142_112_y:',[min(a_142_112_y),max(a_142_112_y)]
#print 'range a_142_112_x:',[min(a_142_112_x),max(a_142_112_x)]
#print 'range a_132_122_y:',[min(a_132_122_y),max(a_132_122_y)]
#print 'range a_132_122_x:',[min(a_132_122_x),max(a_132_122_x)]
#print 'range a_132_113_y:',[min(a_132_113_y),max(a_132_113_y)]
#print 'range a_132_113_x:',[min(a_132_113_x),max(a_132_113_x)]
#print 'range a_132_112_y:',[min(a_132_112_y),max(a_132_112_y)]
#print 'range a_132_112_x:',[min(a_132_112_x),max(a_132_112_x)]
#print 'range a_122_113_y:',[min(a_122_113_y),max(a_122_113_y)]
#print 'range a_122_113_x:',[min(a_122_113_x),max(a_122_113_x)]
#print 'range a_122_112_y:',[min(a_122_112_y),max(a_122_112_y)]
#print 'range a_122_112_x:',[min(a_122_112_x),max(a_122_112_x)]

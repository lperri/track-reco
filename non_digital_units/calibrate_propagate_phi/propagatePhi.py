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


#dictionary that links the (encap,station,ring) to the corresponding [offset, chambers, strips_per_chamber] 
ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,1):[0,18,80],(1,2,2):[0,36,80],(1,3,1):[0,18,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80],(2,1,3):[0,36,64],(2,1,2):[0,36,80],(2,2,1):[0,18,80],(2,2,2):[0,36,80],(2,3,1):[0,18,80],(2,3,2):[0,36,80],(2,4,2):[0,36,80],(1,4,1):[0,18,80],(2,4,1):[0,18,80]}
#ec_s_r_dict = {(1,1,3):[0.086,36,64],(1,1,2):[0.082,36,80],(1,2,2):[0.086,36,80],(1,3,2):[0.069,36,80],(1,4,2):[0.073,36,80],(2,1,3):[0.095,36,64],(2,1,2):[0.079,36,80],(2,2,2):[0.074,36,80],(2,3,2):[0.087,36,80],(2,4,2):[0.086,36,80]}
#from fitPhiK.py -- parameters from the histogram slices of phi vs k/cosTheta (ec,station,ring):[p0,p1,p2] where they are from the function [p0]*x/(1+[p1]*abs(x)) + [p2]
# 
phi_k_dict = {(1,1,2):[-1.58,-0.33,-0.0855],(1,1,3):[-1.84,-0.48,0.085],(1,2,1):[-1.06,-0.6,0.084],(1,2,2):[-1.8,-0.28,0.0857],(1,3,1):[-0.9,-0.3,0.084],(1,3,2):[-1.8,0.02,0.086],(1,4,1):[-1,0.04,0.084],(1,4,2):[-1.7,0.036,0.086],(2,1,2):[1.6,-0.18,0.085],(2,1,3):[1.8,-0.5,0.084],(2,2,1):[1,-0.7,0.083],(2,2,2):[1.8,-0.24,0.086],(2,3,1):[1,-0.3,0.085],(2,3,2):[1.76,-0.11,0.0855],(2,4,1):[0.99,-0.09,0.085],(2,4,2):[1.69,-0.038,0.0852]}

rms_dict = {(1,1,2):[0.0134,0.2],(1,1,3):[0.016,0.11],(1,2,1):[0.055,0.22],(1,2,2):[0.0122,0.19],(1,3,1):[0.204,0],(1,3,2):[0.103,0],(1,4,1):[0.204,0],(1,4,2):[0.103,0],(2,1,2):[0.102,0],(2,1,3):[0.085,0],(2,2,1):[0.199,0],(2,2,2):[0.102,0.15],(2,3,1):[0.0288,0.23],(2,3,2):[0.011,0.22],(2,4,1):[0.0275,0.2],(2,4,2):[0.0105,0.24]}

#rms_dict ={(1,1,2):[0.0134,0.2],(1,1,3):[0.016,0.11],(1,2,1):[0.055,0.22],(1,2,2):[0.0122,0.19],(1,3,1):[0.0288,0.23],(1,3,2):[0.103,0],(1,4,1):[0.204,0],(1,4,2):[0.103,0],(2,1,2):[0.102,0],(2,1,3):[0.085,0],(2,2,1):[0.199,0],(2,2,2):[0.102,0.15],(2,3,1):[0.0288,0.23],(2,3,2):[0.011,0.22],(2,4,1):[0.0275,0.2],(2,4,2):[0.0105,0.24]} 

#need to update detector_z_dict to include new stations
#detector_z_dict = {(1,1,3):-7.0,(1,1,2):-7.0,(1,2,2):-8.3,(1,3,2):-9.3,(1,4,2):-10.3,(2,1,3):7.0,(2,1,2):7.0,(2,2,2):8.3,(2,3,2):9.3,(2,4,2):10.3}

def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
    if (ec,station,ring) in ec_s_r_dict:
	offset, chambers, strips_per_chamber = ec_s_r_dict[ec,station,ring]
	num_chamber = float(stub.first.chamber())
        num_half_strip = float(stub.second.getStrip())
        delta_phi_chamber = float((2*math.pi)/chambers)
        half_strips_per_chamber = float(strips_per_chamber)*2
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
	#print 'offset, chambers, strips_per_chamber = ec_s_r_dict[ec,station,ring] ',offset,chambers,strips_per_chamber
	#print 'num_chamber = float(stub.first.chamber()) ',num_chamber
        #print 'num_half_strip = float(stub.second.getStrip()) ',num_half_strip
        #print 'delta_phi_chamber = float((2*math.pi)/chambers) ',delta_phi_chamber
        #print 'half_strips_per_chamber = strips_per_chamber*2 ',half_strips_per_chamber
        #print 'delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber) ',delta_phi_half_strip
	if (ec,station) == (2,1) or (ec,station) == (2,2) or (ec,station) == (1,3) or (ec,station) == (1,4): 
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
	elif (ec,station) == (2,3) or (ec,station) == (2,4) or (ec,station) == (1,1) or (ec,station) == (1,2):
	#    print 'ok in this part.'
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip)) - 2*math.pi
         #   print 'phi = float(offset + (chambers - num_chamber -1)*deltaphichamber + (hs_per_cham - num_hs)*deltaphihs) ',phi 
        while phi<(-1*math.pi):
	    phi += 2*math.pi
	while phi>math.pi:
	    phi -= 2*math.pi
        #import pdb; pdb.set_trace()
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

def graphRMSofProfileY(hist_2d):
    ''' must be 2d hist '''
    graph_error = ROOT.TGraphErrors()
    for i in range(1,hist_2d.GetNbinsX()+1):
        proj_y = hist_2d.ProjectionY("",i,i)
	k = (hist_2d.GetXaxis()).GetBinCenter(i)
	rms = proj_y.GetRMS()
	rms_error = proj_y.GetRMSError()
	graph_error.SetPoint(i,k,rms)
	graph_error.SetPointError(i,0,rms_error)	
    return graph_error


def resolution(a,b,k):
    return math.sqrt(a**2+(b**2)*k**2)

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0
#bendvsk_222_calc = ROOT.TH2D("bend_222_calc","$\Phi_{b}$ (calc) vs. k; k; $\Phi_{b}$ (calc)",100,-0.35,0.35,100,-1,1)
#bendvsk_222_true = ROOT.TH2D("bend_222_true","$\Phi_{b}$ (true) vs. k; k; $\Phi_{b}$ (true)",100,-0.35,0.35,100,-1,1)


#possible_paths = [((2,4,2),(2,3,2)),((2,4,2),(2,2,2)),((2,4,2),(2,1,2)),((2,4,2),(2,2,1)),((2,3,2),(2,2,2)),((2,3,2),(2,1,3)),((2,3,2),(2,1,2)),((2,3,2),(2,2,1)),((2,2,2),(2,1,3)),((2,2,2),(2,1,2)),((1,4,2),(1,3,2)),((1,4,2),(1,2,2)),((1,4,2),(1,1,2)),((1,4,2),(1,2,1)),((1,3,2),(1,2,2)),((1,3,2),(1,1,3)),((1,3,2),(1,1,2)),((1,3,2),(1,2,1)),((1,2,2),(1,1,3)),((1,2,2),(1,1,2)),((1,4,1),(1,3,1)),((1,4,1),(1,2,1)),((2,4,1),(2,3,1)),((2,4,1),(2,2,1))]

prop_phi_112 = ROOT.TH1D("prop_phi_112","prop_phi_112",100,-math.pi,math.pi)
prop_phi_113 = ROOT.TH1D("prop_phi_113","prop_phi_113",100,-math.pi,math.pi)
prop_phi_122 = ROOT.TH1D("prop_phi_122","prop_phi_122",100,-math.pi,math.pi)
prop_phi_132 = ROOT.TH1D("prop_phi_132","prop_phi_132",100,-math.pi,math.pi)
prop_phi_142 = ROOT.TH1D("prop_phi_142","prop_phi_142",100,-math.pi,math.pi)
prop_phi_212 = ROOT.TH1D("prop_phi_212","prop_phi_212",100,-math.pi,math.pi)
prop_phi_213 = ROOT.TH1D("prop_phi_213","prop_phi_213",100,-math.pi,math.pi)
prop_phi_222 = ROOT.TH1D("prop_phi_222","prop_phi_222",100,-math.pi,math.pi)
prop_phi_232 = ROOT.TH1D("prop_phi_232","prop_phi_232",100,-math.pi,math.pi)
prop_phi_242 = ROOT.TH1D("prop_phi_242","prop_phi_242",100,-math.pi,math.pi)
prop_phi_221 = ROOT.TH1D("prop_phi_221","prop_phi_221",100,-math.pi,math.pi)
prop_phi_121 = ROOT.TH1D("prop_phi_121","prop_phi_121",100,-math.pi,math.pi)
prop_phi_131 = ROOT.TH1D("prop_phi_131","prop_phi_131",100,-math.pi,math.pi)
prop_phi_141 = ROOT.TH1D("prop_phi_141","prop_phi_141",100,-math.pi,math.pi)
prop_phi_231 = ROOT.TH1D("prop_phi_231","prop_phi_231",100,-math.pi,math.pi)
prop_phi_241 = ROOT.TH1D("prop_phi_241","prop_phi_241",100,-math.pi,math.pi)


for event in events:
    N=N+1
    if N > 10000:
	break
#    if N in [50000,100000,150000,200000]:
#	print ('analyzing event {}').format(N)
#    if N > 10000:
#	break
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
  		#if (ec,station,ring) == (1,3,1) or (ec,station,ring) == (2,3,1):
  	 	    #print 'ec,station,ring and ec1,station1,ring1: ',(ec,station,ring) 
  	            k = float(g.charge()/g.pt())
  	            cos_theta = cosTheta(getTheta(g.eta()))
  	            kcos = k/(cos_theta)
  	            p0_s,p1_s,p2_s = phi_k_dict[ec,station,ring]
 		    #print 'p0_s, p1_s, p2_s: ',(p0_s,p1_s,p2_s)
  	           
		    propped_phi_s = g.phi() + (kcos*p0_s)/(1+p1_s*abs(kcos)) + p2_s
      		    
  	            a_s,b_s = rms_dict[ec,station,ring]
  		    #print 'a_s, b_s:  ',(a_s,b_s)
  	            
  	            res_s = resolution(a_s,b_s,kcos)
  	
      	            while propped_phi_s < (-1*math.pi):
     	                propped_phi_s += 2*math.pi
      	            while propped_phi_s > math.pi:
      	                propped_phi_s -= 2*math.pi	
  
  
  	            #print 'propped_phi_s = g.phi() + math.pi + (kcos*p0_s)/(1+p1_s*kcos) + p2_s:  ',propped_phi_s
      	            phi_diff_s = abs(propped_phi_s) - abs(getPhi(s))

                    while phi_diff_s < (-1*math.pi):
     	                phi_diff_s += 2*math.pi
      	            while phi_diff_s > math.pi:
      	                phi_diff_s -= 2*math.pi	
  	            
 		     
  		    phi_diff_res = (phi_diff_s)/res_s
		    if propped_phi_s*getPhi(s) < 0:
                    
                        print 'ec,station,ring',ec,station,ring
                        print 'chamber',s.first.chamber()
		        print 'strip',s.second.getStrip()
  	                print 'g.phi(), propped_phi, getPhi(s):  ',(g.phi(),propped_phi_s,getPhi(s))
                        print 'phi_diff_s ',phi_diff_s
		        print 'res_s ',res_s 
                        print 'phi_diff_res ',phi_diff_res
                        import pdb;pdb.set_trace()
                    if (ec,station,ring) == (2,4,2):
  	                prop_phi_242.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,2,2):
  	                prop_phi_222.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,1,2):
  	                prop_phi_212.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,2,1):
  	                prop_phi_221.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,3,2):
  	                prop_phi_232.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,1,3):
  	                prop_phi_213.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,4,1):
  	                prop_phi_241.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (2,3,1):
  	                prop_phi_231.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,4,2):
  	                prop_phi_142.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,2,2):
  	                prop_phi_122.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,1,2):
  	                prop_phi_112.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,2,1):
  	                prop_phi_121.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,3,2):
  	                prop_phi_132.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,1,3):
  	                prop_phi_113.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,4,1):
  	                prop_phi_141.Fill(phi_diff_res)
  	            elif (ec,station,ring) == (1,3,1):
  	                prop_phi_131.Fill(phi_diff_res)


g=ROOT.TFile("propagatePhi.root","RECREATE")
g.cd()
prop_phi_221.Write()
prop_phi_231.Write()
prop_phi_241.Write()
prop_phi_212.Write()
prop_phi_213.Write()
prop_phi_222.Write()
prop_phi_232.Write()
prop_phi_242.Write()
prop_phi_121.Write()
prop_phi_131.Write()
prop_phi_141.Write()
prop_phi_112.Write()
prop_phi_113.Write()
prop_phi_122.Write()
prop_phi_132.Write()
prop_phi_142.Write()
g.Close()


#
#

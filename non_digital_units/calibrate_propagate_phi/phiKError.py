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


#dictionary that links the (encap,station,ring) to the corresponding [offset, chambers, strips_per_chamber] 
#ec_s_r_dict = {(1,1,3):[0,36,64],(1,1,2):[0,36,80],(1,2,1):[0,18,80],(1,2,2):[0,36,80],(1,3,1):[0,18,80],(1,3,2):[0,36,80],(1,4,2):[0,36,80],(2,1,3):[0,36,64],(2,1,2):[0,36,80],(2,2,1):[0,18,80],(2,2,2):[0,36,80],(2,3,1):[0,18,80],(2,3,2):[0,36,80],(2,4,2):[0,36,80],(1,4,1):[0,18,80],(2,4,1):[0,18,80]}
#ec_s_r_dict = {(1,1,3):[0.086,36,64],(1,1,2):[0.082,36,80],(1,2,2):[0.086,36,80],(1,3,2):[0.069,36,80],(1,4,2):[0.073,36,80],(2,1,3):[0.095,36,64],(2,1,2):[0.079,36,80],(2,2,2):[0.074,36,80],(2,3,2):[0.087,36,80],(2,4,2):[0.086,36,80]}


#ec_s_r_dict [CORRESPONDING TO 'NEWEST_PLOTS',PRE-ECSWITCH] = {(1,1,3):[0.096,36,64],(1,1,2):[0.094,36,80],(1,2,1):[0.12,18,80],(1,2,2):[0.098,36,80],(1,3,1):[0.091,18,80],(1,3,2):[0.084,36,80],(1,4,2):[0.084,36,80],(2,1,3):[0.089,36,64],(2,1,2):[0.086,36,80],(2,2,1):[0.093,18,80],(2,2,2):[0.074,36,80],(2,3,1):[0.11,18,80],(2,3,2):[0.081,36,80],(2,4,2):[0.081,36,80],(1,4,1):[0.080,18,80],(2,4,1):[0.12,18,80]}

#post-endcap switch!!! looks correct:
ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}



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
        if (ec,station) == (2,1) or (ec,station) == (2,2) or (ec,station) == (1,3) or (ec,station) == (1,4):
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
        elif (ec,station) == (2,3) or (ec,station) == (2,4) or (ec,station) == (1,1) or (ec,station) == (1,2):
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (half_strips_per_chamber+1-num_half_strip)*delta_phi_half_strip)
        while phi<(-1*math.pi):
            phi += 2*math.pi
        while phi>math.pi:
            phi -= 2*math.pi
        
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
    graph_error.SetTitle(("RMS for {}").format(hist_2d))
    for i in range(1,hist_2d.GetNbinsX()+1):
        proj_y = hist_2d.ProjectionY("",i,i)
	k = (hist_2d.GetXaxis()).GetBinCenter(i)
	rms = proj_y.GetRMS()
	rms_error = proj_y.GetRMSError()
	graph_error.SetPoint(i,k,rms)
	graph_error.SetPointError(i,0,rms_error)	
    return graph_error


events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0
#bendvsk_222_calc = ROOT.TH2D("bend_222_calc","$\Phi_{b}$ (calc) vs. k; k; $\Phi_{b}$ (calc)",100,-0.35,0.35,100,-1,1)
#bendvsk_222_true = ROOT.TH2D("bend_222_true","$\Phi_{b}$ (true) vs. k; k; $\Phi_{b}$ (true)",100,-0.35,0.35,100,-1,1)


phi_k_112 = ROOT.TH2D("phi_k_112","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_113 = ROOT.TH2D("phi_k_113","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_121 = ROOT.TH2D("phi_k_121","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_122 = ROOT.TH2D("phi_k_122","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_131 = ROOT.TH2D("phi_k_131","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_132 = ROOT.TH2D("phi_k_132","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_141 = ROOT.TH2D("phi_k_141","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_142 = ROOT.TH2D("phi_k_142","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)


phi_k_212 = ROOT.TH2D("phi_k_212","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_213 = ROOT.TH2D("phi_k_213","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_221 = ROOT.TH2D("phi_k_221","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_222 = ROOT.TH2D("phi_k_222","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_231 = ROOT.TH2D("phi_k_231","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_232 = ROOT.TH2D("phi_k_232","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_241 = ROOT.TH2D("phi_k_241","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)
phi_k_242 = ROOT.TH2D("phi_k_242","$\Delta$$\phi$ vs  k/cos; k/cos; $\Delta\Phi$",100,-0.7,0.7,100,-1.4,1.4)


for event in events:
    N=N+1
    if N in [100,1000,50000,100000,150000,200000]:
	print ('analyzing event {}').format(N)
    #if N>10000:
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
#		if (ec,station,ring) == (2,2,2):
       		    k = float(g.charge()/g.pt())
		    cos_theta = cosTheta(getTheta(g.eta()))
		    diff_phi = (getPhi(s)) - g.phi()
		    kcos = k/(cos_theta)
    	            while diff_phi < (-1*math.pi):
   	    	        diff_phi += 2*math.pi
    	    	    while diff_phi > math.pi:
    	    	        diff_phi -= 2*math.pi	
     #               print 'getPhi(s) ',getPhi(s)
     #               print 'g.phi() ',g.phi()
     #               print 'diff_phi ',diff_phi
     #  	    	    import pdb;pdb.set_trace()
                    if (ec,station,ring) == (1,1,2):
			phi_k_112.Fill(kcos,diff_phi)       			
       	    	    if (ec,station,ring) == (1,1,3):
			phi_k_113.Fill(kcos,diff_phi)       			
       	    	    if (ec,station,ring) == (1,2,1):
			phi_k_121.Fill(kcos,diff_phi)
       	    	    if (ec,station,ring) == (1,2,2):
			phi_k_122.Fill(kcos,diff_phi)
                       #if abs(k)<0.02 and abs(diff_phi)>0.1:
                       #    import pdb;pdb.set_trace()
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


print 'computing RMS...'
phi_k_error_112 = graphRMSofProfileY(phi_k_112)
phi_k_error_112.SetNameTitle("rms_112","rms_112")
phi_k_error_113 = graphRMSofProfileY(phi_k_113)
phi_k_error_113.SetNameTitle("rms_113","rms_113")
phi_k_error_121 = graphRMSofProfileY(phi_k_121)
phi_k_error_121.SetNameTitle("rms_121","rms_121")
phi_k_error_122 = graphRMSofProfileY(phi_k_122)
phi_k_error_122.SetNameTitle("rms_122","rms_122")
phi_k_error_131 = graphRMSofProfileY(phi_k_131)
phi_k_error_131.SetNameTitle("rms_131","rms_131")
phi_k_error_132 = graphRMSofProfileY(phi_k_132)
phi_k_error_132.SetNameTitle("rms_132","rms_132")
phi_k_error_141 = graphRMSofProfileY(phi_k_141)
phi_k_error_141.SetNameTitle("rms_141","rms_141")
phi_k_error_142 = graphRMSofProfileY(phi_k_142)
phi_k_error_142.SetNameTitle("rms_142","rms_142")
phi_k_error_212 = graphRMSofProfileY(phi_k_212)
phi_k_error_212.SetNameTitle("rms_212","rms_212")
phi_k_error_213 = graphRMSofProfileY(phi_k_213)
phi_k_error_213.SetNameTitle("rms_213","rms_213")
phi_k_error_221 = graphRMSofProfileY(phi_k_221)
phi_k_error_221.SetNameTitle("rms_221","rms_221")
phi_k_error_222 = graphRMSofProfileY(phi_k_222)
phi_k_error_222.SetNameTitle("rms_222","rms_222")
phi_k_error_231 = graphRMSofProfileY(phi_k_231)
phi_k_error_231.SetNameTitle("rms_231","rms_231")
phi_k_error_232 = graphRMSofProfileY(phi_k_232)
phi_k_error_232.SetNameTitle("rms_232","rms_232")
phi_k_error_241 = graphRMSofProfileY(phi_k_241)
phi_k_error_241.SetNameTitle("rms_241","rms_241")
phi_k_error_242 = graphRMSofProfileY(phi_k_242)
phi_k_error_242.SetNameTitle("rms_242","rms_242")

print 'writing files...'
f=ROOT.TFile("phi_k_plots_ecSwitch.root","RECREATE")
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

g=ROOT.TFile("rms_plots_ecSwitch.root","RECREATE")
g.cd()
phi_k_error_112.Write()
phi_k_error_113.Write()
phi_k_error_121.Write()
phi_k_error_122.Write()
phi_k_error_131.Write()
phi_k_error_132.Write()
phi_k_error_141.Write()
phi_k_error_142.Write()
phi_k_error_212.Write()
phi_k_error_213.Write()
phi_k_error_221.Write()
phi_k_error_222.Write()
phi_k_error_231.Write()
phi_k_error_232.Write()
phi_k_error_241.Write()
phi_k_error_242.Write()
g.Close()


#
#

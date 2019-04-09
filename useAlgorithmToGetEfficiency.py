import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle
from array import array
import csv
import numpy

def median(lst):
    """ returns median out of a list """
    return numpy.median(numpy.array(lst))

def cscSegmentsPhi(event,ontime=True,twinMux=True):
    """ for a given event, cscSegmentsPhi returns """ 
    # CSC local trigger measures phi coordinate; Local Charged Tracks (LCTs) are synonymous with muon segments
    phiSeg = Handle('std::vector<std::pair<CSCDetId,CSCCorrelatedLCTDigi> >')
    # obtain trigger hits for a given event
    event.getByLabel('cscTriggerHits',phiSeg)
    #
    digis = phiSeg.product()
    return digis

def fetchCSCGEANT(event):
    """ obtain simulated particles for a given event """
    geantH  = Handle  ('vector<PSimHit>')
    event.getByLabel('g4SimHits:MuonCSCHits',geantH)
    geant=filter(lambda x: x.pabs()>0.5 and abs(x.particleType())==13,geantH.product())
    return geant

class CSCSeg(object):
    """ CSC segment object """
    def __init__(self,s):
        self.seg=s
        self.truePhiB=-1
    def __getattr__(self,name):
        return getattr(self.seg,name)

def matchedCSCStubs(muon,segments,geant):
    """ return matched stubs in CSC chambers """
    
    # muons will appear in oppositely charged pairs, so thisMuonGEANT reflects this -- positively charged muons correspond to particle ID of -13 (CMS convention) and negatively charged muons correspond to particle ID of 13.
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==-13) or ((muon.charge()<0) and x.particleType()==13),geant)
    # thisMuonGEANT is a PSimHit Object
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
    # collect all simulated partilces for an event
    event.getByLabel('genParticles',genH)
    # only keep the muons that are in the correct eta range
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





gen_muons_eta = ROOT.TH1D('gen_muons','gen_muons;phi;muons (>= 2 assoc. stubs)',100,1.2,2.4)
gen_muons_pt = ROOT.TH1D('gen_muons','gen_muons;phi;muons (>= 2 assoc. stubs)',100,0,100)
gen_muons_phi = ROOT.TH1D('gen_muons','gen_muons;phi;muons (>= 2 assoc. stubs)',100,-3.14,3.14)
my_muons_eta = ROOT.TH1D('my_muons','my_muons;phi;muons (>= 2 assoc. stubs)',100,1.2,2.4)
my_muons_pt = ROOT.TH1D('my_muons','my_muons;phi;muons (>= 2 assoc. stubs)',100,0,100)
my_muons_phi = ROOT.TH1D('my_muons','my_muons;phi;muons (>= 2 assoc. stubs)',100,-3.14,3.14)



    
#eta_ranges has key (min_eta, max_eta) and value (ec,station,ring)
eta_ranges = {(1,1,3):(0.90,1.20), (1,1,2):(1.18,1.70), (1,2,2):(1.00,1.60), (1,2,1):(1.60,2.40), (1,3,2):(1.10,1.70), (1,3,1):(1.70,2.40), (1,4,1):(1.80,2.40), (1,4,2):(1.16,1.80), (2,1,3):(-1.20,-0.90), (2,1,2):(-1.70,-1.18), (2,2,2):(-1.60,-1.00), (2,2,1):(-2.40,-1.60), (2,3,2):(-1.70,-1.10), (2,3,1):(-2.40,-1.70),(2,4,1):(-2.40,-1.80),(2,4,2):(-1.80,-1.16)}

#phi_prop_coefs gives you the a,b,c i.e. propagation coefficients -- used in the muon.phiProp method
phi_prop_coefs = {(1,1,2):[-1.58,-0.33,-0.0855],(1,1,3):[-1.84,-0.48,0.085],(1,2,1):[-1.06,-0.6,0.084],(1,2,2):[-1.8,-0.28,0.0857],(1,3,1):[-0.9,-0.3,0.084],(1,3,2):[-1.8,0.02,0.086],(1,4,1):[-1,0.04,0.084],(1,4,2):[-1.7,0.036,0.086],(2,1,2):[1.6,-0.18,0.085],(2,1,3):[1.8,-0.5,0.084],(2,2,1):[1,-0.7,0.083],(2,2,2):[1.8,-0.24,0.086],(2,3,1):[1,-0.3,0.085],(2,3,2):[1.76,-0.11,0.0855],(2,4,1):[0.99,-0.09,0.085],(2,4,2):[1.69,-0.038,0.0852]}

#wire_prop_coefs gives you x,y,w i.e. propagation coefficients -- used in the muon.wireProp method 
wire_prop_coefs = {(1, 2, 1): [499.9, -303.3, 39.59], (2, 1, 3): [273.5, 364.9, 106.9], (2, 3, 2): [336.7, 334.3, 79.78], (1, 3, 2): [332.3, -327.8, 77.39], (2, 2, 2): [312, 330.2, 83.64], (1, 1, 3): [285, -388, 118], (1, 4, 2): [339.5, -315.5, 69.73], (1, 4, 1): [812.3, -564.8, 94.43], (1, 3, 1): [524.1, -324.8, 44.81], (1, 1, 2): [197.2, -86.25, -21.13], (1, 2, 2): [310, -327, 82.4], (2, 1, 2): [360.4, 316.2, 59.59], (2, 4, 1): [803, 555.9, 92.31], (2, 4, 2): [349, 328.5, 74.13], (2, 3, 1): [473.1, 276.8, 33.59], (2, 2, 1): [470, 270.9, 30.84]}

#wire_rms_coefs gives you the coefficients to calculate the rms for the variable wire -- used in the muon.stubThresholdTest method
wire_rms_coefs = {(2, 1, 2): [3.78391, 14.105], (2, 1, 3): [1.3635, 16.8427], (2, 2, 2): [1.85818, 16.8152], (2, 3, 2): [1.7054, 16.2489], (1, 3, 2): [1.69797, 17.1057], (1, 2, 1): [4.2155, 18.7755], (1, 2, 2): [1.85454, 16.7872], (1, 4, 1): [3.21621, 21.5728], (1, 3, 1): [3.32301, 17.47], (2, 4, 2): [1.65894, 16.0229], (2, 2, 1): [4.41787, 14.9874], (1, 4, 2): [1.61513, 17.1709], (1, 1, 3): [1.40545, 16.968], (2, 4, 1): [3.25701, 22.0568], (1, 1, 2): [3.7663, 12.4721], (2, 3, 1): [3.48749, 15.1532]}


#phi_rms_coefs gives you the coefficients to calculate the rms for the variable phi -- used in the muon.stubThresholdTest method 
phi_rms_coefs = {(1,1,2):[0.0134,0.2],(1,1,3):[0.016,0.11],(1,2,1):[0.055,0.22],(1,2,2):[0.0122,0.19],(1,3,1):[0.204,0],(1,3,2):[0.103,0],(1,4,1):[0.204,0],(1,4,2):[0.103,0],(2,1,2):[0.102,0],(2,1,3):[0.085,0],(2,2,1):[0.199,0],(2,2,2):[0.102,0.15],(2,3,1):[0.0288,0.23],(2,3,2):[0.011,0.22],(2,4,1):[0.0275,0.2],(2,4,2):[0.0105,0.24]}

#phi_pull_3_sig gives you the threshold at 3 sigma for the variable phi -- used in the muon.stubThresholdTest method
phi_pull_3_sig = {(2,2,1):1.833, (2,3,1):1.440, (2,4,1):1.424, (2,1,2):1.780, (2,1,3):2.191, (2,2,2):1.835, (2,3,2):1.882, (2,4,2):1.938, (1,2,1):1.964, (1,3,1):1.424, (1,4,1):1.430, (1,1,2):1.647, (1,1,3):2.104, (1,2,2):1.734, (1,3,2):2.157, (1,4,2):2.021}

wire_pull_3_sig = {(2,1,2):2.661,(2,3,2):2.848,(1,3,2):2.835,(2,2,2):2.949,(1,1,3):2.524,(1,4,1):2.469,(1,3,1):2.282,(1,2,1):2.267,(1,1,2):2.672,(1,2,2):2.931,(1,4,2):2.791,(2,4,1):2.451,(2,4,2):2.766,(2,3,1):2.233,(2,1,3):2.468,(2,2,1):2.236}



# create dictionaries for 4 sigma as well
phi_pull_4_sig = {}
wire_pull_4_sig = {}
for key in phi_pull_3_sig:
    pull_3 = phi_pull_3_sig[key] 
    phi_pull_4_sig[key] = (pull_3/(3.0))*4
for key in wire_pull_3_sig:
    pull_3 = wire_pull_3_sig[key]
    wire_pull_4_sig[key] = (pull_3/(3.0))*4



class muon(object):

    """ reconstructed muon object """

    def __init__(self,gen_eta,gen_phi,k):
        """ muon has gen variables measured in tracker; true muon defined by number of associated stubs """
        self.gen_eta = gen_eta
        self.gen_phi = gen_phi
        self.k = k
    
    def possibleEcStationRing(self):
        """ obtain a list of possible (EC,Station,Ring) given gen eta """
        possible_esr = []
        for key in eta_ranges:
            if self.gen_eta >= (eta_ranges[key][0]) and self.gen_eta <= (eta_ranges[key][1]):
                if key in ec_s_r_dict:
                    possible_esr.append(key)
        return possible_esr
   
    def angleInRange(self,angle):
        """ makes sure angle is between -pi,pi """
        while angle < (-1*math.pi):
            angle += 2*math.pi
        while angle > math.pi:
            angle -= 2*math.pi
        return angle

    def getTheta(self):
        """ takes in value of pseudorapidity and converts to theta"""
        return 2*math.atan(1/math.exp(self.gen_eta))

    def phiPropagated(self,(ec,station,ring)):
        """ returns propogated phi """
        a,b,c = phi_prop_coefs[(ec,station,ring)]
        theta = self.getTheta()
        kcos = k/math.cos(theta)
        phi = self.gen_phi + (kcos*a)/(1+b*abs(kcos)) + c
        return self.angleInRange(phi)

    def wirePropagated(self,(ec,station,ring)):
        """ returns propogated wire number  """
        x,y,w = wire_prop_coefs[(ec,station,ring)] 
        return x+(y*self.gen_eta)+(w*(self.gen_eta**2))

    def getPhi(self,stub):
        """ returns phi value based on half strip number """
        # which endcap, station, ring got hit
        ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
        if (ec,station,ring) in ec_s_r_dict:
            offset, chambers, strips_per_chamber = ec_s_r_dict[ec,station,ring]
            # num_chamber gives the chamber number that detected a hit
            num_chamber = float(stub.first.chamber())
            # num_half_strip gives the half strip that was hit
            num_half_strip = float(stub.second.getStrip())
            # delta_phi_chamber gives the angle subtended by one chamber in radians -- depends on number of chambers in a ring (which is just based on geometry of the detector)
            delta_phi_chamber = float((2*math.pi)/chambers)
            # the CSCs have overlayed individual strips in layers, such that the resolution ends up being .5 strip
            half_strips_per_chamber = float(strips_per_chamber)*2
            # delta_phi_half_strip is the angle subtended by a half strip in radians
            delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
            # the way phi is calculated varies depending on the endcap & station because of the geometry of the detector and which direction the half strip numbers increase
            if (ec,station) == (2,1) or (ec,station) == (2,2) or (ec,station) == (1,3) or (ec,station) == (1,4):
                # for these (ec,station) combinations, phi increases in the direction of increasing strip number
                phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
            elif (ec,station) == (2,3) or (ec,station) == (2,4) or (ec,station) == (1,1) or (ec,station) == (1,2):
                # for these (ec,station) combinations, phi decreases in the direction of increasing strip number
                phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip)) - 2*math.pi
            return self.angleInRange(phi)


    def resolution(self,coefs):
        """ returns resolution -- constant term and k-dependent term """
        a,b = coefs
        return math.sqrt(a**2+(b**2)*(self.k)**2)


    def stubThresholdTest(self,stub):
        
        """ returns stub if it passes both thresholds -- i.e. if it's possible that it started at tracker and propagated there """ 
        ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
        if (ec,station,ring) in ec_s_r_dict:
            #determine phi threshold
            phi_threshold = phi_pull_4_sig[(ec,station,ring)]
            phi_res = self.resolution(phi_rms_coefs[(ec,station,ring)])
            phi_diff = self.angleInRange(self.getPhi(stub) - self.phiPropagated((ec,station,ring))) 
            #determine wire threshold
            wire_threshold = wire_pull_4_sig[(ec,station,ring)]
            wire_res = self.resolution(wire_rms_coefs[(ec,station,ring)])
            wire_diff = stub.second.getKeyWG() - int(self.wirePropagated((ec,station,ring)))
            #test the stub against both
            if abs(phi_diff/phi_res) < phi_threshold and abs(wire_diff/wire_res) < wire_threshold:
                return stub 


    def collectAcceptedStubsPerECSR(self,segments,(ec,station,ring)):
        """ return all accepted stubs (those that pass thresholds) within a given ECSR -- not yet associated """
        accepted_stubs = []
        for stub in segments:
            if stub.first.endcap() == ec and stub.first.station() == station and stub.first.ring() == ring:
                if self.stubThresholdTest(stub):
                    accepted_stubs.append((self.stubThresholdTest(stub)))
        return accepted_stubs


    def closestStubPerECSR(self,accepted_stubs,(ec,station,ring)):
        """ compare each accepted stub's phi to the tracker values -- take only closest one as associated PER EC,STATION,RING"""
        zipped_stub_distance = [(stub,self.angleInRange(self.getPhi(stub) - self.phiPropagated((ec,station,ring)))) for stub in accepted_stubs]
        #minimize zipped list of stubs and their respective distance to the propagated track
        return min(zipped_stub_distance, key=lambda x: x[1])[0]

    def main(self,segments):
        """ collect all associated stubs in a given ECSR and return the number of associated stubs """
        associated_stubs = []
        for combo in self.possibleEcStationRing():
            ec,station,ring = combo
            accepted_stubs = self.collectAcceptedStubsPerECSR(segments,combo)
            if accepted_stubs:
                associated_stubs.append(self.closestStubPerECSR(accepted_stubs,combo))
        #print 'my algorithm gives this many stubs: ',len(associated_stubs)
        alg_ecsr = []
        for stub in associated_stubs:
            alg_ecsr.append((stub.first.endcap(),stub.first.station(),stub.first.ring()))
        #print 'for my algorithm, these assoc. stubs have ec,station,ring of ',alg_ecsr
        return len(associated_stubs) 




ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0

for event in events:
    N=N+1
    if N%10000 == 0:
        print ('--- analyzing event {} ---').format(N)
    genMuons=fetchGEN(event)
    geantCSC=fetchCSCGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)
    i=0
    for g in genMuons:
        i+=1
        #print 'muon {}'.format(i)
        # drop real muons that do not land in our detectors
        if abs(g.eta())>1.2 and abs(g.eta())<2.4:
            
            cscMatchedSegments = matchedCSCStubs(g,segmentsCSC,geantCSC)
            segments = segmentsCSC
            pt = g.pt()
            phi = g.phi()
            eta = g.eta()
            k = g.charge()/g.pt() 
            mu = muon(g.eta(),g.phi(),k)
            
            ecsr = []
            #here we get the # stubs found by THEIR Alg
            lenCSCMatchedSeg = 0
            for stub in cscMatchedSegments:
                if (stub.first.station(),stub.first.ring()) != (1,1): 
                    ecsr.append((stub.first.endcap(),stub.first.station(),stub.first.ring()))
                    lenCSCMatchedSeg += 1
                
            num_stubs_my_alg = mu.main(segments)
            if lenCSCMatchedSeg >= 2:
                gen_muons_eta.Fill(eta)
                gen_muons_phi.Fill(phi)
                gen_muons_pt.Fill(pt)
            if num_stubs_my_alg >= 2:
                my_muons_eta.Fill(eta)
                my_muons_phi.Fill(phi)
                my_muons_pt.Fill(pt)





f = ROOT.TFile('efficiency.root','RECREATE')
f.cd()
gen_muons_eta.Write('genMuonsEta')
gen_muons_pt.Write('genMuonsPT')
gen_muons_phi.Write('genMuonsPhi')
my_muons_eta.Write('myMuonsEta')
my_muons_pt.Write('myMuonsPT')
my_muons_phi.Write('myMuonsPhi')
eff_eta = ROOT.TGraphAsymmErrors(my_muons_eta,gen_muons_eta)
eff_pt = ROOT.TGraphAsymmErrors(my_muons_pt,gen_muons_pt)
eff_phi = ROOT.TGraphAsymmErrors(my_muons_phi,gen_muons_phi)
eff_eta.Write('efficiency_eta')
eff_pt.Write('efficiency_pt')
eff_phi.Write('efficiency_phi')

f.Close()

#
#

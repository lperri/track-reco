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





    
#eta_ranges is a dictionary with key (min_eta, max_eta) and value (ec,station,ring)
eta_ranges = {(1,1,3):(0.9,1.4), (1,1,2):(1.18,1.7), (1,2,2):(1,1.6), (1,2,1):(1.6,2.6), (1,3,2):(1.1,1.7), (1,3,1):(1.7,2.6), (1,4,1):(1.8,2.6), (1,4,2):(1.16,1.8), (2,1,3):(-1.4,-0.9), (2,1,2):(-1.7,-1.18), (2,2,2):(-1.6,-1), (2,2,1):(-2.6,-1.6), (2,3,2):(-1.7,-1.1), (2,3,1):(-2.6,-1.7),(2,4,1):(-2.6,-1.8),(2,4,2):(-1.8,-1.16)}

#phi_prop_coefs gives you the a,b,c in the phiProp function
phi_prop_coefs = {(1,1,2):[-1.58,-0.33,-0.0855],(1,1,3):[-1.84,-0.48,0.085],(1,2,1):[-1.06,-0.6,0.084],(1,2,2):[-1.8,-0.28,0.0857],(1,3,1):[-0.9,-0.3,0.084],(1,3,2):[-1.8,0.02,0.086],(1,4,1):[-1,0.04,0.084],(1,4,2):[-1.7,0.036,0.086],(2,1,2):[1.6,-0.18,0.085],(2,1,3):[1.8,-0.5,0.084],(2,2,1):[1,-0.7,0.083],(2,2,2):[1.8,-0.24,0.086],(2,3,1):[1,-0.3,0.085],(2,3,2):[1.76,-0.11,0.0855],(2,4,1):[0.99,-0.09,0.085],(2,4,2):[1.69,-0.038,0.0852]}

#wire_prop_coefs gives you x,y,w in the wireProp function 
wire_prop_coefs = {(1, 2, 1): [499.9, -303.3, 39.59], (2, 1, 3): [273.5, 364.9, 106.9], (2, 3, 2): [336.7, 334.3, 79.78], (1, 3, 2): [332.3, -327.8, 77.39], (2, 2, 2): [312, 330.2, 83.64], (1, 1, 3): [285, -388, 118], (1, 4, 2): [339.5, -315.5, 69.73], (1, 4, 1): [812.3, -564.8, 94.43], (1, 3, 1): [524.1, -324.8, 44.81], (1, 1, 2): [197.2, -86.25, -21.13], (1, 2, 2): [310, -327, 82.4], (2, 1, 2): [360.4, 316.2, 59.59], (2, 4, 1): [803, 555.9, 92.31], (2, 4, 2): [349, 328.5, 74.13], (2, 3, 1): [473.1, 276.8, 33.59], (2, 2, 1): [460, 270.9, 30.84]}


wire_rms_coefs = {(2, 1, 2): [3.78391, 14.105], (2, 1, 3): [1.3635, 16.8427], (2, 2, 2): [1.85818, 16.8152], (2, 3, 2): [1.7054, 16.2489], (1, 3, 2): [1.69797, 17.1057], (1, 2, 1): [4.2155, 18.7755], (1, 2, 2): [1.85454, 16.7872], (1, 4, 1): [3.21621, 21.5728], (1, 3, 1): [3.32301, 17.47], (2, 4, 2): [1.65894, 16.0229], (2, 2, 1): [4.41787, 14.9874], (1, 4, 2): [1.61513, 17.1709], (1, 1, 3): [1.40545, 16.968], (2, 4, 1): [3.25701, 22.0568], (1, 1, 2): [3.7663, 12.4721], (2, 3, 1): [3.48749, 15.1532]}


#phi_rms_coefs gives you the coefficients to calculates the rms within the stubThresholdTest function 
phi_rms_coefs = {(1,1,2):[0.0134,0.2],(1,1,3):[0.016,0.11],(1,2,1):[0.055,0.22],(1,2,2):[0.0122,0.19],(1,3,1):[0.204,0],(1,3,2):[0.103,0],(1,4,1):[0.204,0],(1,4,2):[0.103,0],(2,1,2):[0.102,0],(2,1,3):[0.085,0],(2,2,1):[0.199,0],(2,2,2):[0.102,0.15],(2,3,1):[0.0288,0.23],(2,3,2):[0.011,0.22],(2,4,1):[0.0275,0.2],(2,4,2):[0.0105,0.24]}


phi_pull_3_sig = {(2,2,1):1.833, (2,3,1):1.440, (2,4,1):1.424, (2,1,2):1.780, (2,1,3):2.191, (2,2,2):1.835, (2,3,2):1.882, (2,4,2):1.938, (1,2,1):1.964, (1,3,1):1.424, (1,4,1):1.430, (1,1,2):1.647, (1,1,3):2.104, (1,2,2):1.734, (1,3,2):2.157, (1,4,2):2.021}

wire_pull_3_sig = {(2,1,2):2.661,(2,3,2):2.848,(1,3,2):2.835,(2,2,2):2.949,(1,1,3):2.524,(1,4,1):2.469,(1,3,1):2.282,(1,2,1):2.267,(1,1,2):2.672,(1,2,2):2.931,(1,4,2):2.791,(2,4,1):2.451,(2,4,2):2.766,(2,3,1):2.233,(2,1,3):2.468,(2,2,1):0}

# FIX 2,2,1 SOMETHING WENT WRONG DURING CALIBRATION



class muon(object):

    ''' reconstructed muon object '''

    def __init__(self,gen_eta,gen_phi,k):
        ''' muon has gen variables measured in tracker; true muon defined by number of associated stubs '''
        self.gen_eta = gen_eta
        self.gen_phi = gen_phi
        self.k = k
    
    def possibleEcStationRing(self):
        ''' obtain a list of possible (EC,Station,Ring) given gen eta '''
        possible_esr = []
        for key in eta_ranges:
            if self.gen_eta >= key[0] and self.gen_eta <= key[1]:
                if key in ec_s_r_dict:
                    possible_esr.append(key)
        return possible_esr
    
    def angleInRange(self,angle):
        ''' makes sure angle is between -pi,pi '''
        while angle < (-1*math.pi):
            angle += 2*math.pi
        while angle > math.pi:
            angle -= 2*math.pi
        return angle

    def getTheta(self):
        ''' takes in value of pseudorapidity and converts to theta'''
        return 2*math.atan(1/math.exp(self.gen_eta))

    def phiPropagated(self,(ec,station,ring)):
        ''' returns propogated phi '''
        a,b,c = phi_prop_coefs[(ec,station,ring)]
        theta = self.getTheta()
        kcos = k/math.cos(theta)
        phi = self.gen_phi + (kcos*a)/(1+b*abs(kcos)) + c
        return self.angleInRange(phi)

    def wirePropagated(self,(ec,station,ring)):
        ''' returns propogated wire number  '''
        x,y,w = wire_prop_coefs[(ec,station,ring)] 
        return x+(y*self.gen_eta)+(w*(self.gen_eta**2))

    def getPhi(self,stub):
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
                phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip)) - 2*math.pi
            return self.angleInRange(phi)


    def resolution(self,coefs):
        a,b = coefs
        return math.sqrt(a**2+(b**2)*(self.k)**2)


    def stubThresholdTest(self,stub):
        
        ''' returns stub if it passes both thresholds -- i.e. if it's possible that it started at tracker and propagated there '''
        ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
        if (ec,station,ring) in ec_s_r_dict:
            #determine phi threshold
            phi_threshold = phi_pull_3_sig[(ec,station,ring)]
            phi_res = self.resolution(phi_rms_coefs[(ec,station,ring)])
            
            phi_diff = self.angleInRange(self.getPhi(stub) - self.gen_phi) 
            #determine wire threshold
            wire_threshold = wire_pull_3_sig[(ec,station,ring)]
            wire_res = self.resolution(wire_rms_coefs[(ec,station,ring)])
            wire_diff = stub.second.getKeyWG() - int(self.wirePropagated((ec,station,ring)))
            #test the stub against both
            if abs(phi_diff/phi_res) < phi_threshold and abs(wire_diff/wire_res) < wire_threshold:
                return stub 

    def collectAcceptedStubs(self,segments):
        ''' return all accepted stubs -- not yet associated '''
        accepted_stubs = []
        for s in segments:
            #run each function
            if self.stubThresholdTest(s):
                accepted_stubs.append(self.stubThresholdTest(s))
        return accepted_stubs
    
    
    def associatedStub(self,accepted_stubs,(ec,station,ring)):
        ''' compare each accepted stub's phi to the tracker values -- take only closest one as associated PER EC,STATION,RING'''
        zipped_stub_distance = [(s,self.angleInRange(self.getPhi(s) - self.phiPropagated((ec,station,ring)))) for s in accepted_stubs]
        #minimize zipped list of stubs and their respective distance to the propagated track
        return min(zipped_stub_distance, key=lambda x: x[1])

    def main(self,segments):
        associated_stubs = []
        for combo in self.possibleEcStationRing():
            ec,station,ring = combo
        # return the number of associated
            accepted_stubs = self.collectAcceptedStubs(segments)
            if accepted_stubs:
                associated_stubs.append(self.associatedStub(accepted_stubs,(ec,station,ring)))
        return len(associated_stubs) 




ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}

events = Events(['file:/scratch2/Leah/CMSSW_10_1_7/src/cscHits.root'])

N=0


for event in events:
    N=N+1
    if N > 50:
	break
    print ('analyzing event {}').format(N)
    #real muons
    genMuons=fetchGEN(event)
    segmentsPhiDT=fetchDTSegmentsPhi(event)
    geantDT=fetchDTGEANT(event)
    geantCSC=fetchCSCGEANT(event)
    segmentsCSC=cscSegmentsPhi(event)
    i=0
    for g in genMuons:
        i+=1
        print
        print
        print "muon {}".format(i)
        # drop real muons that do not land in our detectors
        if abs(g.eta())>0.9 and abs(g.eta())<2.4:
            cscMatchedSegments = matchedCSCStubs(g,segmentsCSC,geantCSC)
            segments = segmentsCSC

            k = g.charge()/g.pt() 
            mu = muon(g.eta(),g.phi(),k)
            print 'all segs: ',mu.main(segments)
            print 'matched segs: ',mu.main(cscMatchedSegments)
#
#

import ROOT
from DataFormats.FWLite import Events, Handle
ROOT.FWLiteEnabler.enable()
from makeHistograms import ec_s_r
import math

tag='/scratch3/MTF/data/190401/singleMu0'
isData=False

def fetchNewHits(event):
    """ obtain events """
    emtfHitH    = Handle  ('std::vector<l1t::EMTFHit>')
    event.getByLabel('simEmtfDigis',emtfHitH)
    #BX = bunch crossing
    return filter(lambda x: x.BX()==0,emtfHitH.product())

def fetchGEN(event,etaMin=0.8,etaMax=2.5):
    """ obtain gen muons for a given eta range """
    genH  = Handle  ('vector<reco::GenParticle>')
    event.getByLabel('genParticles',genH)
    genMuons=filter(lambda x: abs(x.pdgId())==13 and x.status()==1 and abs(x.eta())<etaMax and abs(x.eta())>etaMin,genH.product())
    return genMuons

def matchedStubs(g,hits):
    """ returns matched stubs based on proximity (in eta) to tracker value """
    matched=[]
    for hit in hits:
        if abs(g.eta()-hit.Eta())<0.3:
            matched.append(hit)
    return matched


#def matchedStubs(g,all_hits):
#    """ returns matched stubs for a muon based on hit proximity to gen phi """
#    matched_stubs = []
#    for sector in ec_s_r:
#        stubs_in_sector = []
#        for hit in all_hits:
#            if (hit.Endcap(),hit.Station(),hit.Ring()) == sector:
#                stubs_in_sector.append(hit)
#        if stubs_in_sector:
#            phiDiff_etaDiff_Stub = [(abs(g.phi() - math.radians(stub.Phi_glob())),(abs(g.eta()-stub.Eta())),stub) for stub in stubs_in_sector]
#            #print(sector, phiDiff_etaDiff_Stub)
#            min_delta_phi,delta_eta,stub = min(phiDiff_etaDiff_Stub, key = lambda t: t[0])
#            #print('min of phi values of stubs in sector ',min_delta_phi)
#            if delta_eta < 0.3 and g.phi()*(math.radians(stub.Phi_glob())) > 0:
#                 matched_stubs.append(stub)
#            #print('matched stubs for this muon',matched_stubs)
#            #print('we have matched stubs in ',[(s.Endcap(),s.Station(),s.Ring()) for s in matched_stubs])
#    return matched_stubs




events = Events([tag+'.root'])


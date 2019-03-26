import ROOT
from DataFormats.FWLite import Events, Handle
ROOT.FWLiteEnabler.enable()

tag='/scratch3/MTF/data/190220/singleMu200_Aging'
isData=False

def fetchNewHits(event):
    """ obtain events """
    emtfHitH    = Handle  ('std::vector<l1t::EMTFHit>')
    event.getByLabel('simEmtfDigis',emtfHitH)
    #BX = bunch crossing
    return filter(lambda x: x.BX()==0,emtfHitH.product())

def fetchGEN(event,etaMin=1.2,etaMax=2.5):
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

events = Events([tag+'.root'])


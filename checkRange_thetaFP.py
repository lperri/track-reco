from angleManipulations import *
from retrieveDataAsEvents import *

histo = ROOT.TH1D('fp','fp',128,0,127)

#list ec_s_r gives all possible (endcap,station,ring) combinations that we care about for this project
ec_s_r = [(1,1,2),(1,1,3),(1,2,1),(1,2,2),(1,3,1),(1,3,2),(1,4,2),(-1,1,3),(-1,1,2),(-1,2,1),(-1,2,2),(-1,3,1),(-1,3,2),(-1,4,2),(1,4,1),(-1,4,1)]

for event in events:
    gen_muons = fetchGEN(event,0.8,2.5)
    emtf_hits = fetchNewHits(events)
    for gen_muon in gen_muons:
        matched_stubs = matchedStubs(gen_muon, emtf_hits)
        for hit in matched_stubs:
            ec,station,ring = hit.Endcap(),hit.Station(),hit.Ring()
            if (ec,station,ring) not in ec_s_r:
                break;
            histo.Fill(hit.Theta_fp())

f = ROOT.TFile('checkRangeThetaFP_newHits.root','RECREATE')
f.cd()
histo.Write()
f.Close()


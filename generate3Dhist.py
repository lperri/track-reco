from angleManipulations import *
from retrieveDataAsEvents import *
import math
import ROOT
ROOT.FWLiteEnabler.enable()

#list ec_s_r gives all possible (endcap,station,ring) combinations that we care about for this project
ec_s_r = [(1,1,2),(1,1,3),(1,2,1),(1,2,2),(1,3,1),(1,3,2),(1,4,2),(-1,1,3),(-1,1,2),(-1,2,1),(-1,2,2),(-1,3,1),(-1,3,2),(-1,4,2),(1,4,1),(-1,4,1)]

#list ec_s gives all possible (endcap,station) combinations that we care about
ec_s = [(1,1),(1,2),(1,3),(1,4),(-1,1),(-1,2),(-1,3),(-1,4)]

hists = {}
for combo in ec_s:
    if combo[0] == 1:
        ec_s__string = 'hist_ec'+str(combo[0])+'_st'+str(combo[1])
    else:
        ec_s__string = 'hist_ec'+'n1'+'_st'+str(combo[1])
    hists[(combo)] = ROOT.TH3D(ec_s__string,ec_s__string,100,-4000,4000,100,-5000,5000,128,0,127)

hists_separate_chambers = {}
for combo in ec_s_r:
    if combo[0] == 1:
        ec_s_r__string = 'hist_ec'+str(combo[0])+'_st'+str(combo[1])+'_r'+str(combo[2])
    else:
        ec_s_r__string = 'hist_ec'+str(combo[0])+'_st'+str(combo[1])+'_r'+str(combo[2])
    hists_separate_chambers[(combo)] = ROOT.TH3D(ec_s_r__string,ec_s_r__string,100,-4000,4000,100,-5000,5000,128,0,127)

counter = -1
for event in events:
    gen_muons = fetchGEN(events,0.8,2.5)
    if counter%10000==0:
        print 'analyzing event ',counter
    counter += 1
    emtf_hits = fetchNewHits(events)
    for gen_muon in gen_muons:
        matched_stubs = matchedStubs(gen_muon,emtf_hits)
        k = gen_muon.charge()/gen_muon.pt()
        gen_theta = (etaToTheta(gen_muon.eta()))

        #gen_theta = angleInRadianRange(etaToTheta(gen_muon.eta()),(0,math.pi))
        k_digi = int(8192*(k/math.cos(gen_theta)))
        for hit in matched_stubs:
            ec,station,ring = hit.Endcap(),hit.Station(),hit.Ring()
            if (ec,station,ring) not in ec_s_r:
                break;
            phi_diff = hit.Phi_fp() - genPhiToDigi(gen_muon,hit)
            theta_fp = hit.Theta_fp()
            if counter == 17:
                print 'at event 17'
                print 'theta fp',hit.Theta_fp()
                try:
                    hists[(ec,station)].Fill(k_digi,phi_diff,theta_fp)
                except:
                    print 'did not work'
                hists_separate_chambers[(ec,station,ring)].Fill(k_digi,phi_diff,theta_fp)



f = ROOT.TFile('hists_3d_402.root','RECREATE')
f.cd()
for combo in ec_s:
    hists[combo].Write()
    #proj_yx = hists[combo].Project3D('yx')
    #proj_yx.Write('proj_k_phi_ME'+str(combo))
    proj_z = hists[combo].ProjectionZ(str(combo)+'_ProjectionZ')
    proj_z.Write('proj_z_ME'+str(combo))
f.Close()

g = ROOT.TFile('hists_3d_separate_chambers_402.root','RECREATE')
g.cd()
for combo in ec_s_r:
    hists_separate_chambers[combo].Write()
    #proj_yx = hists[combo].Project3D('yx')
    #proj_yx.Write('proj_k_phi_ME'+str(combo))
    proj_z = hists_separate_chambers[combo].ProjectionZ(str(combo)+'_ProjectionZ')
    proj_z.Write('proj_z_ME'+str(combo))
g.Close()

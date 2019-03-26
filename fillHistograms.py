from retrieveDataAsEvents import *
from angleManipulations import *
from makeHistograms import *

counter = -1
for event in events:
    
    if counter%10000 == 0:
        print 'analyzing event {}'.format(counter)
    if counter == 50000:
        break;
    counter += 1
    gen_muons = fetchGEN(event,1.2,2.5)
    emtf_hits = fetchNewHits(events)

    for gen_muon in gen_muons:
        
        matched_stubs = matchedStubs(gen_muon,emtf_hits)
        gen_theta = angleInRadianRange(etaToTheta(gen_muon.eta()),(0,math.pi))
        gen_theta_digi = genThetaToDigi(gen_theta)
        k = gen_muon.charge()/gen_muon.pt()
        k_digi = int(8192*(k/math.cos(gen_theta)))

        for hit in matched_stubs:
            
            ec, station, ring = hit.Endcap(), hit.Station(), hit.Ring()
            if (ec, station, ring) not in ec_s_r:
                break;
            
            # note that gen_phi as defined below first gets put into global coordinates WRT the center of CMS then gets put in range (0,2PI) in preparation for conversion to digital coordinates
            gen_phi_digi = genPhiToDigi(gen_muon,hit)
            phi_difference = hit.Phi_fp() - gen_phi_digi
            hist_phi_calibration[(ec,station,ring)].Fill(k_digi,phi_difference)
            hist_theta_calibration[(ec,station,ring)].Fill(gen_theta,hit.Theta_fp())


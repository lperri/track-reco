from retrieveDataAsEvents import *
from angleManipulations import *
from makeHistograms import *

counter = -1
for event in events:
    
    if counter%10000 == 0:
        print 'analyzing event {}'.format(counter)
    #if counter == 100000:
    #    break;
    counter += 1
    gen_muons = fetchGEN(events,0.8,2.5)
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
            if phi_difference > -3800 and phi_difference < 4200:
            print 'k',k
            print 'gen_phi',gen_muon.phi()
            print 'gen_phi_digi',gen_phi_digi
            print 
            hist_phi_calibration[(ec,station,ring)].Fill(k_digi,phi_difference)

#            hist_theta_calibration[(ec,station,ring)].Fill(gen_theta,hit.Theta_fp())
            #hist_theta_calibration[(ec,station)].Fill(gen_theta,hit.Theta_fp())
            #hist_3d[(ec,station)].Fill(k_digi,phi_difference,hit.Theta_fp()) 
            #hist_1d[(ec,station)].Fill(hit.Theta_fp()) 
            

from __future__ import print_function
import ROOT,itertools,math      #
from array import array         # 
from DataFormats.FWLite import Events, Handle
ROOT.FWLiteEnabler.enable()

#tag='singleMuonOfficial'
#isData=False
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

def phiRangePi(phi):
    """ puts phi in CMS range, -pi to pi """
    while phi < (-1*math.pi):
        phi += 2*math.pi
    while phi > math.pi:
        phi -= 2*math.pi
    return phi

def phiRange2Pi(phi): 
    """ puts phi in normal range, 0 to 2pi """
    while phi < 0:
        phi += 2*math.pi
    while phi > 2*math.pi:
        phi -= 2*math.pi
    return phi

def etaToTheta(eta):
    """ converts eta to theta """
    return 2*math.atan(1/math.exp(eta))

def thetaRangePi(theta):
    """ puts theta in range, 0 to pi """
    while theta < 0:
        theta += 2*math.pi
    while theta > math.pi:
        theta -= 2*math.pi
    return theta

def convertToDigital(which_angle,angle_value):
    """ which_angle = 'phi' or 'theta'; angle_value = angle in radians, must be in correct range (for phi must be in range 0 to 2pi, theta must be 0 to pi); digital range -- (from 0 to) 4920 for phi and (from 0 to) 127 for theta -- will return the angle digitized """
    if which_angle == 'phi':
        max_digital_range = 4920
        max_radian_range = 2*math.pi
    elif which_angle == 'theta':
        max_digital_range = 127
        max_radian_range = math.pi
    return int(angle_value)*(max_digital_range/max_radian_range)

def noNegativesInRootFileName(name):
    """ To avoid syntax error -- takes in the combiantion, (ec, station, ring) and changes negative endcap (-1) to (n1) -- returns 
    n1' string """
    if name[0] < 0:
        return 'n1'


########### set up histograms ###########
events=Events([tag+'.root'])

#list ec_s_r gives all possible (endcap,station,ring) combinations that we care about for this project
ec_s_r = [(1,1,2),(1,1,3),(1,2,1),(1,2,2),(1,3,1),(1,3,2),(1,4,2),(-1,1,3),(-1,1,2),(-1,2,1),(-1,2,2),(-1,3,1),(-1,3,2),(-1,4,2),(1,4,1),(-1,4,1)]

# hist_phi_offset is dictionary of histograms with (endcap,station,ring) as the key; value is phi_offset_ECSR histgoram and 
# hist_theta_offset is dictionary of histograms with (endcap,station,ring) as the key; value is theta_offset_ECSR histgoram
hist_phi_offset = {}
hist_theta_offset = {}

#fill both dictionaries with the correct histogram objects
for combo in ec_s_r:
    ec_s_r__string = str(combo[0]) + str(combo[1]) + str(combo[2])
    hist_phi_offset[(combo)] = ROOT.TH2D('phi_offset_'+str(ec_s_r__string),'phi_offset_'+str(ec_s_r__string)+';curvature;'+'phi_fp - gen_phi',100,-4000,4000,100,0,4920)
    hist_theta_offset[(combo)] = ROOT.TH2D('theta_offset_'+str(ec_s_r__string),'theta_offset_'+str(ec_s_r__string)+';gen_theta;'+'theta_fp',100,0,127,100,0,127)

# phi_offset_slope_loss is (endcap,station,ring):[offset,slope,loss]
phi_offset_slope_loss = {}
# theta_offset_slope_loss is (endcap,station,ring):[offset,slope,loss]
theta_offset_slope_loss = {}


########### main loop ###########
counter=-1
for event in events:
    if counter%10000 == 0:
        print ('analyzing event {}'.format(counter))
    #if counter > 10000:
    #    break;
    counter+=1
    gen=fetchGEN(event,1.2,2.5)
    emtf=fetchNewHits(event)

    for g in gen:
        matched=matchedStubs(g,emtf)
        gen_theta = thetaRangePi(etaToTheta(g.eta()))  
        k_digi_cos = int(((g.charge()/g.pt())/(math.cos(gen_theta)))*8192)
        
        for hit in matched:
            ec,station,ring = hit.Endcap(),hit.Station(),hit.Ring()
            if (ec,station,ring) not in ec_s_r:
                break;
            gen_phi_digi = convertToDigital('phi',phiRange2Pi(g.phi()))
            phi_difference = hit.Phi_fp() - g.phi()
            hist_phi_offset[(ec,station,ring)].Fill(k_digi_cos,phi_difference)
            hist_theta_offset[(ec,station,ring)].Fill(gen_theta,hit.Theta_fp())





#now fit all histograms and store parameters in dictionary phi_offset_slope, theta_offset_slope
#also save root file of all histograms

f = ROOT.TFile('phi_theta_offset.root', 'RECREATE')
f.cd()

# the 'combo' is jsut the (endcap, station ring) and this is the same for hist_phi_offset & hist_theta_offset, therefore only need to loop once
for combo in hist_phi_offset:
    
    # protect against root throwing syntax error
    ec = noNegativesInRootFileName(combo)
    if combo[0] == -1:
        name = ec+str(combo[1])+str(combo[2])
    else:
        name = str(combo[0])+str(combo[1])+str(combo[2])
    
    # write 2d histograms
    obj_phi = hist_phi_offset[combo]
    c1 = ROOT.TCanvas('phi_offset_'+str(name), 'phi_offset_'+str(name))
    c1.SetGrid()
    c1.cd()
    obj_phi.Write('phi_offset_'+str(name))
    obj_phi.Draw('AP')
    c1.Update()
    
    obj_theta = hist_theta_offset[combo]
    c2 = ROOT.TCanvas('theta_offset_'+str(name), 'theta_offset_'+str(name))
    c2.SetGrid()
    c2.cd()
    obj_theta.Write('theta_offset_'+str(name))
    obj_theta.Draw('AP')
    c2.Update()

    # fitting sliced 2d histograms
    c3 = ROOT.TCanvas('fitted'+str(name),'fitted'+str(name))
    c2.cd()
    sliced = obj_phi.ProfileX()
    fit_function = ROOT.TF1('fit_function','x*[1]/(1+[2]*abs(x)) + [0]',sliced.GetXaxis().GetXmin(),sliced.GetXaxis().GetXmax())
    sliced.Fit('fit_function','S')
    sliced.Draw()
    sliced.Write('phi_offset_'+str(name)+'fitted')
    pointer = sliced.GetFunction('fit_function')
    offset  = pointer.GetParameter(0)
    slope = pointer.GetParameter(1)
    energy_loss = pointer.GetParameter(2) 
    phi_offset_slope_loss[combo] = [offset,slope,energy_loss]
    c2.Update()

#for combo in hist_theta_offset:
#    c3 = ROOT.TCanvas('theta_offset_'+str(combo), 'theta_offset_'+str(combo))
#    hist_theta_offset[combo].Write('theta_offset_'+str(combo))
#    hist_theta_offset[combo].Draw('theta_offset_'+str(combo))

f.Close()
print(phi_offset_slope_loss)



###### note to self: these attributes are defined as follows ######
#print('Gen Muon pt={pt},eta={eta},phi={phi},charge={charge}'.format(pt=g.pt(),eta=g.eta(),phi=g.phi(),charge=g.charge()))
#print('----New hit BX={bx} RPC={rpc} GEM={gem} CSC={csc} endcap={endcap} station={station} ring={ring} chamber={chamber} phi (0-4920) ={phi} theta (0-127)={theta}'.format(bx=hit.BX(),rpc=hit.Is_RPC(), gem=hit.Is_GEM(),csc=hit.Is_CSC(),endcap=hit.Endcap(), station=hit.Station(), ring=hit.Ring(),chamber= hit.Chamber(),phi=hit.Phi_fp(), theta=hit.Theta_fp()))

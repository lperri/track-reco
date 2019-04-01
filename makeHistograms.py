import ROOT
ROOT.FWLiteEnabler.enable()
import math

#list ec_s_r gives all possible (endcap,station,ring) combinations that we care about for this project
ec_s_r = [(1,1,2),(1,1,3),(1,2,1),(1,2,2),(1,3,1),(1,3,2),(1,4,2),(-1,1,3),(-1,1,2),(-1,2,1),(-1,2,2),(-1,3,1),(-1,3,2),(-1,4,2),(1,4,1),(-1,4,1)]
ec_s = [(1,1),(1,2),(1,3),(1,4),(-1,1),(-1,2),(-1,3),(-1,4)]
# hist_phi_offset is dictionary of histograms with (endcap,station,ring) as the key; value is phi_offset_ECSR histgoram and 
# hist_theta_offset is dictionary of histograms with (endcap,station,ring) as the key; value is theta_offset_ECSR histgoram

hist_phi_calibration = {}
hist_theta_calibration = {}

#fill both dictionaries with the correct histogram objects

#for combo in ec_s_r:
for combo in ec_s:
    ec_s_r__string = str(combo[0]) + str(combo[1]) + str(combo[2])
    ec_s__string = str(combo[0]) + str(combo[1])
    #hist_phi_calibration[(combo)] = ROOT.TH2D('phi_offset_'+str(ec_s_r__string),'phi_offset_'+str(ec_s_r__string)+';curvature;'+'phi_fp - gen_phi',100,-4000,4000,100,-4920,4920)
    if combo[0] == 1:
    # positive endcap should only go up to pi/2
        hist_theta_calibration[(combo)] = ROOT.TH2D('theta_offset_'+str(ec_s_r__string),'theta_offset_'+str(ec_s_r__string)+';gen_theta;'+'theta_fp',100,0,math.pi/2,128,0,127)
    elif combo[0] == -1:
    # negative endcap should go from pi/2 --> pi
        hist_theta_calibration[(combo)] = ROOT.TH2D('theta_offset_'+str(ec_s_r__string),'theta_offset_'+str(ec_s_r__string)+';gen_theta;'+'theta_fp',100,math.pi/2,math.pi,128,0,127)



# phi_offset_slope_loss is (endcap,station,ring):[offset,slope,loss]
phi_offset_slope_loss = {}
# theta_offset_slope_loss is (endcap,station,ring):[offset,slope,loss]
theta_offset_slope = {}


import ROOT
from fillHistograms import *

# fit all histograms and store parameters in dictionary phi_offset_slope_loss, theta_offset_slope
# save root file of all histograms

root_file_name_phi = 'phi_calibration'
root_file_name_theta = 'theta_calibration'
fit_function_phi = 'x*[1]/(1+[2]*abs(x)) + [0]'
fit_function_theta = '[0] + [1]*x'


def drawHistogram(root_file_name,dictionary,combo,name):
    """ inputs: dictionary (which dictionary containing the histogram objects do you want to use), combo ((ec,station,ring)) i.e. the key of the dictionary letting it know which histogram you want, and name (what you want to call histogram); will draw histogram on canvas and write object to root file """
    obj = dictionary[combo]
    c1 = ROOT.TCanvas(root_file_name+'_'+name,root_file_name+'_'+name)
    c1.SetGrid()
    c1.cd()
    obj.Write(root_file_name+name)
    obj.Draw('AP')
    c1.Update()
    return obj
    
def fitSliced2DHistogram(root_file_name,obj,params_dictionary,combo,name,function):
    """ inputs: obj (histogram object), params_dictionary (which dictionary containing the fitting parameters do you want to update), combo ((ec,station,ring)) i.e. the key of the dictionary letting it know which histogram you want, name (what you want to call histogram) and function (function to fit sliced 2d histogram to) """
    c2 = ROOT.TCanvas(root_file_name+'_'+name+'_fitted',root_file_name+'_'+name+'_fitted')
    c2.cd()
    sliced_obj = obj.ProfileX()
    fit_function = ROOT.TF1('fit_function',function,sliced_obj.GetXaxis().GetXmin(),sliced_obj.GetXaxis().GetXmax())
    sliced_obj.Fit('fit_function','S','Q')
    ROOT.gStyle.SetOptStat()
    sliced_obj.Draw('AP')
    sliced_obj.Write(root_file_name+'_'+name+'_fitted')
    pointer = sliced_obj.GetFunction('fit_function')
    offset = pointer.GetParameter(0)
    slope = pointer.GetParameter(1)
    try:
        energy_loss = pointer.GetParameter(2)
        params_dictionary[combo] = [offset,slope,energy_loss]
    except:
        params_dictionary[combo] = [offset, slope]
    c2.Update()
    
    
def generateRootName(combo):
    """ input: combo (ec,station,ring), generates name """
    if combo[0] == -1:
        name = 'n1'+str(combo[1])+str(combo[2])
    else:
        name = str(combo[0])+str(combo[1])+str(combo[2])
    return name


f = ROOT.TFile(root_file_name_phi+'.root','RECREATE')
f.cd()
for combo in hist_phi_calibration:
    # construct histogram names based on (ec,station,ring)
    name = generateRootName(combo)
    #draw abd write 2D histograms
    drawn_obj = drawHistogram(root_file_name_phi,hist_phi_calibration,combo,name)
    # fit sliced 2D histograms
    fitSliced2DHistogram(root_file_name_phi,drawn_obj,phi_offset_slope_loss,combo,name,fit_function_phi)
f.Close()


g = ROOT.TFile(root_file_name_theta+'.root','RECREATE')
g.cd()
# do same for theta
for combo in hist_theta_calibration:
    name = generateRootName(combo)
    drawn_obj = drawHistogram(root_file_name_theta,hist_theta_calibration,combo,name)
    fitSliced2DHistogram(root_file_name_theta,drawn_obj,theta_offset_slope,combo,name,fit_function_theta)
g.Close()

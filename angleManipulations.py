import math

def angleInRadianRange(angle,(minimum,maximum)):
    """ takes an angle in radians and returns it in the specified range """
    while angle < minimum:
        angle += 2*math.pi
    while angle > maximum:
        angle += 2*math.pi
    return angle

def etaToTheta(eta):
    """ converts eta to theta """
    return 2*math.atan(1/math.exp(eta))

def genThetaToDigi(angle_value):
    """ angle_value = angle in radians, must be in correct range (0 to pi); digital range is from 0 to 127 for theta -- will return the angle digitized """
    max_digital_range = 127
    max_radian_range = math.pi
    return int(angle_value)*(max_digital_range/max_radian_range)

def genPhiToDigi(gen_muon,hit):
    """ takes gen phi and puts into global digital coordinates """
    # the offset of pi/12 appears because the first sector starts at 15 degrees or pi/12 rad
    sector_size = math.pi/3
    #return int((angleInRadianRange(gen_muon.phi(),(0,2*math.pi)) - (hit.Sector()-1)*math.pi/3.0 - math.pi/12)*4920*3/math.pi)
    return int((angleInRadianRange(gen_muon.phi(),(0,2*math.pi)) - (hit.Sector()-1)*sector_size - math.pi/12)*4920/sector_size)

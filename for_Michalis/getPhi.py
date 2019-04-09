ec_s_r_dict = {(1,1,3):[0.085,36,64],(1,1,2):[0.086,36,80],(1,2,1):[0.085,18,80],(1,2,2):[0.086,36,80],(1,3,1):[0.092,18,80],(1,3,2):[0.090,36,80],(1,4,2):[0.090,36,80],(2,1,3):[0.090,36,64],(2,1,2):[0.090,36,80],(2,2,1):[0.092,18,80],(2,2,2):[0.090,36,80],(2,3,1):[0.085,18,80],(2,3,2):[0.086,36,80],(2,4,2):[0.086,36,80],(1,4,1):[0.091,18,80],(2,4,1):[0.085,18,80]}

def getPhi(stub):
    '''takes in hit and returns global phi between -pi and pi radians'''
    ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
    if (ec,station,ring) in ec_s_r_dict:
        offset, chambers, strips_per_chamber = ec_s_r_dict[ec,station,ring]
        num_chamber = float(stub.first.chamber())
        num_half_strip = float(stub.second.getStrip())
        delta_phi_chamber = float((2*math.pi)/chambers)
        half_strips_per_chamber = float(strips_per_chamber)*2
        delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
        if (ec,station) == (2,1) or (ec,station) == (2,2) or (ec,station) == (1,3) or (ec,station) == (1,4):
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + ((half_strips_per_chamber+1-num_half_strip)*delta_phi_half_strip))
        elif (ec,station) == (2,3) or (ec,station) == (2,4) or (ec,station) == (1,1) or (ec,station) == (1,2):
            phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
        while phi<(-1*math.pi):
            phi += 2*math.pi
        while phi>math.pi:
            phi -= 2*math.pi
        return phi


class muon(object):

    """ reconstructed muon object """

    def __init__(self,gen_eta,gen_phi,k):
        """ muon has gen variables measured in tracker; true muon defined by number of associated stubs """
        self.gen_eta = gen_eta
        self.gen_phi = gen_phi
        self.k = k
    
    def possibleEcStationRing(self):
        """ obtain a list of possible (EC,Station,Ring) given gen eta """
        possible_esr = []
        for key in eta_ranges:
            if self.gen_eta >= (eta_ranges[key][0]) and self.gen_eta <= (eta_ranges[key][1]):
                if key in ec_s_r_dict:
                    possible_esr.append(key)
        return possible_esr
   
    def angleInRange(self,angle):
        """ makes sure angle is between -pi,pi """
        while angle < (-1*math.pi):
            angle += 2*math.pi
        while angle > math.pi:
            angle -= 2*math.pi
        return angle

    def getTheta(self):
        """ takes in value of pseudorapidity and converts to theta"""
        return 2*math.atan(1/math.exp(self.gen_eta))

    def phiPropagated(self,(ec,station,ring)):
        """ returns propogated phi """
        a,b,c = phi_prop_coefs[(ec,station,ring)]
        theta = self.getTheta()
        kcos = k/math.cos(theta)
        phi = self.gen_phi + (kcos*a)/(1+b*abs(kcos)) + c
        return self.angleInRange(phi)

    def wirePropagated(self,(ec,station,ring)):
        """ returns propogated wire number  """
        x,y,w = wire_prop_coefs[(ec,station,ring)] 
        return x+(y*self.gen_eta)+(w*(self.gen_eta**2))

    def getPhi(self,stub):
        """ returns phi value based on half strip number """
        # which endcap, station, ring got hit
        ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
        if (ec,station,ring) in ec_s_r_dict:
            offset, chambers, strips_per_chamber = ec_s_r_dict[ec,station,ring]
            # num_chamber gives the chamber number that detected a hit
            num_chamber = float(stub.first.chamber())
            # num_half_strip gives the half strip that was hit
            num_half_strip = float(stub.second.getStrip())
            # delta_phi_chamber gives the angle subtended by one chamber in radians -- depends on number of chambers in a ring (which is just based on geometry of the detector)
            delta_phi_chamber = float((2*math.pi)/chambers)
            # the CSCs have overlayed individual strips in layers, such that the resolution ends up being .5 strip
            half_strips_per_chamber = float(strips_per_chamber)*2
            # delta_phi_half_strip is the angle subtended by a half strip in radians
            delta_phi_half_strip = float(delta_phi_chamber/half_strips_per_chamber)
            # the way phi is calculated varies depending on the endcap & station because of the geometry of the detector and which direction the half strip numbers increase
            if (ec,station) == (2,1) or (ec,station) == (2,2) or (ec,station) == (1,3) or (ec,station) == (1,4):
                # for these (ec,station) combinations, phi increases in the direction of increasing strip number
                phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip))
            elif (ec,station) == (2,3) or (ec,station) == (2,4) or (ec,station) == (1,1) or (ec,station) == (1,2):
                # for these (ec,station) combinations, phi decreases in the direction of increasing strip number
                phi = float(offset + (num_chamber -1)*(delta_phi_chamber) + (num_half_strip*delta_phi_half_strip)) - 2*math.pi
            return self.angleInRange(phi)


    def resolution(self,coefs):
        """ returns resolution -- constant term and k-dependent term """
        a,b = coefs
        return math.sqrt(a**2+(b**2)*(self.k)**2)


    def stubThresholdTest(self,stub):
        
        """ returns stub if it passes both thresholds -- i.e. if it's possible that it started at tracker and propagated there """ 
        ec,station,ring = (stub.first.endcap(),stub.first.station(),stub.first.ring())
        if (ec,station,ring) in ec_s_r_dict:
            #determine phi threshold
            phi_threshold = phi_pull_4_sig[(ec,station,ring)]
            phi_res = self.resolution(phi_rms_coefs[(ec,station,ring)])
            phi_diff = self.angleInRange(self.getPhi(stub) - self.phiPropagated((ec,station,ring))) 
            #determine wire threshold
            wire_threshold = wire_pull_4_sig[(ec,station,ring)]
            wire_res = self.resolution(wire_rms_coefs[(ec,station,ring)])
            wire_diff = stub.second.getKeyWG() - int(self.wirePropagated((ec,station,ring)))
            #test the stub against both
            if abs(phi_diff/phi_res) < phi_threshold and abs(wire_diff/wire_res) < wire_threshold:
                return stub 


    def collectAcceptedStubsPerECSR(self,segments,(ec,station,ring)):
        """ return all accepted stubs (those that pass thresholds) within a given ECSR -- not yet associated """
        accepted_stubs = []
        for stub in segments:
            if stub.first.endcap() == ec and stub.first.station() == station and stub.first.ring() == ring:
                if self.stubThresholdTest(stub):
                    accepted_stubs.append((self.stubThresholdTest(stub)))
        return accepted_stubs


    def closestStubPerECSR(self,accepted_stubs,(ec,station,ring)):
        """ compare each accepted stub's phi to the tracker values -- take only closest one as associated PER EC,STATION,RING"""
        zipped_stub_distance = [(stub,self.angleInRange(self.getPhi(stub) - self.phiPropagated((ec,station,ring)))) for stub in accepted_stubs]
        #minimize zipped list of stubs and their respective distance to the propagated track
        return min(zipped_stub_distance, key=lambda x: x[1])[0]

    def main(self,segments):
        """ collect all associated stubs in a given ECSR and return the number of associated stubs """
        associated_stubs = []
        for combo in self.possibleEcStationRing():
            ec,station,ring = combo
            accepted_stubs = self.collectAcceptedStubsPerECSR(segments,combo)
            if accepted_stubs:
                associated_stubs.append(self.closestStubPerECSR(accepted_stubs,combo))
        #print 'my algorithm gives this many stubs: ',len(associated_stubs)
        alg_ecsr = []
        for stub in associated_stubs:
            alg_ecsr.append((stub.first.endcap(),stub.first.station(),stub.first.ring()))
        #print 'for my algorithm, these assoc. stubs have ec,station,ring of ',alg_ecsr
        return len(associated_stubs) 
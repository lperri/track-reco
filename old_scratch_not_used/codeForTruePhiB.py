class CSCSeg(object):
    def __init__(self,s):
        self.seg=s
        self.truePhiB=-1
    def __getattr__(self,name):
        return getattr(self.seg,name)
    


def matchedCSCStubs(muon,segments,geant):
    thisMuonGEANT = filter(lambda x: (muon.charge()>0 and x.particleType()==13) or ((muon.charge()<0) and x.particleType()==-13),geant)
    chambers=[]
    for p in thisMuonGEANT:        
        detid=ROOT.CSCDetId(p.detUnitId())
        chambers.append(p.detUnitId())
      
    assocSeg=[]   
    for s in segments:
        minLayer=8;
        maxLayer=0;
        xin=0;
        yin=0;
        xout=0;
        yout=0;

        associated=False
        for c in thisMuonGEANT:
            detid=ROOT.CSCDetId(c.detUnitId())
#            print 'Det=',detid.station(),detid.ring(),detid.chamber(),detid.layer(),'Angle=',c.entryPoint().phi().value(),'entry=',c.entryPoint().x(),c.entryPoint().y(),c.entryPoint().z()
            if detid.layer()<minLayer:
                minLayer=detid.layer()
                xin =c.entryPoint().x()
                yin =c.entryPoint().y()
            if detid.layer()>maxLayer:
                maxLayer=detid.layer()
                xout =c.entryPoint().x()
                yout =c.entryPoint().y()
            if detid.endcap()==s.first.endcap() and detid.station()==s.first.station() and detid.ring()==s.first.ring() and detid.chamber()==s.first.chamber():
                associated=True
        if associated:
            seg = CSCSeg(s)
            seg.truePhiB = math.atan((xout-xin)/(3.0*(maxLayer-minLayer)))
            assocSeg.append(seg)
    return assocSeg


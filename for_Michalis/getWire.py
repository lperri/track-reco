def wireNumCalib(stub):
    dict_wire_calib = {(1, 2, 1): [499.9, -303.3, 39.59], (2, 1, 3): [273.5, 364.9, 106.9], (2, 3, 2): [336.7, 334.3, 79.78], (1, 3, 2): [332.3, -327.8, 77.39], (2, 2, 2): [312, 330.2, 83.64], (1, 1, 3): [285, -388, 118], (1, 4, 2): [339.5, -315.5, 69.73], (1, 4, 1): [812.3, -564.8, 94.43], (1, 3, 1): [524.1, -324.8, 44.81], (1, 1, 2): [197.2, -86.25, -21.13], (1, 2, 2): [310, -327, 82.4], (2, 1, 2): [360.4, 316.2, 59.59], (2, 4, 1): [803, 555.9, 92.31], (2, 4, 2): [349, 328.5, 74.13], (2, 3, 1): [473.1, 276.8, 33.59], (2, 2, 1): [460, 270.9, 30.84]}
    ec,station,ring = stub.first.endcap(),stub.first.station(),stub.first.ring()
    x,y,w = dict_wire_calib[(ec,station,ring)]
    eta = stub.eta()
    return x + (y*eta) + w*(eta**2)

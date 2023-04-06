# -*- coding: utf-8 -*-
"""
Created on Mon May 10 17:12:14 2021
@author: Evangelos
"""

from pyvisa import ResourceManager

class LO():

    def __init__(self,address,reset=False):
        #creates instance of the class, shorthand for commands, and optionally resets the device
        self.inst = ResourceManager().open_resource(address)
        self.w = self.inst.write
        self.r = self.inst.read
        self.q = self.inst.query
        if reset:
            self.w("*RST")
        #self.w("DISP:PSAV ON")
        #self.w("DISP:UPD OFF")

    def id(self):
        # queries instrument identification info
        return self.q("*IDN?")

    def RF_ON(self):
        #enables RF output
        self.w("OUTP ON")

    def RF_OFF(self):
        # disables RF output
        self.w("OUTP OFF")

    def set_freq(self, value):
        # sets the RF frequency in GHz
        comm = "FREQ %s GHz" % value
        self.w(comm)

    def freq_mode(self,mode):
        #sets the frequency mode (sweep or CW)
        self.w(f"FREQ:MODE {mode}")
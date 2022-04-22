 # $Id: geom.py,v 1.61 2012/10/08 16:14:09 samn Exp $ 

# from pyinit import *
# from labels import *
# h.celsius = 37
from neuron import h
		
###############################################################################
# General Cell
###############################################################################
class Cell:
  "General cell"	
  def __init__ (self,x=0,y=0,z=0,ID=0,ty=0):
    self.x=x
    self.y=y
    self.z=z
    self.ID=ID
    self.ty = ty
    self.snames = [] # list of section names
    self.all_sec = []
    self.add_comp('soma', False)
    self.set_morphology()
    self.set_conductances()
    self.set_synapses()
    self.calc_area()
    self.infod = {} # dictionary for storing indices into nsl,ncl
    self.poID = [] # list of postsynaptic IDs (indices into Network's ce)
    self.poNC = [] # list of pointers to postsynaptic NetCons
    self.poSY = [] # synapse type (code from labels.py)
    #self.movetopos()

  # saves information on a synapse to another cell
  #  poid is postsynaptic id, nc is NetCon, syty is synapse code (from labels.py)
  def savesyinfo (self,poid,nc,syty):
    self.poID.append(poid)
    self.poNC.append(nc)
    self.poSY.append(syty)

  # get number of outgoing connections
  def getdvi (self): return len(self.poID)
  def set_morphology (self): pass			
  def set_conductances (self): pass		
  def set_synapses (self): pass		
          
  def add_comp (self, name, rec):
    self.snames.append( name )
    self.__dict__[name] = h.Section(cell=self, name=name)
    self.all_sec.append(self.__dict__[name])    
    if rec: # Record voltage
      self.__dict__[name+"_volt"] = h.Vector(int(h.tstop/h.dt)+1)
      self.__dict__[name+"_volt"].record(self.__dict__[name](0.5)._ref_v)
      self.__dict__[name+"_volt"].label(name+"_volt")
	
  def plot_volt (self, name, fig=1):
    figure(fig)
    volt = self.__dict__[name+"_volt"].to_python()
    plot(arange(len(volt))*h.dt, volt)
		
  def calc_area (self):
    self.total_area = 0
    self.n = 0
    for sect in self.all_sec:
      self.total_area += h.area(0.5,sec=sect)
      self.n+=1



#-------------------------------
# Thalamic reticular cell (sRE)
class sRE (Cell):
  def __init__ (self,x=0,y=0,z=0,ID=0,ty=0):
    Cell.__init__(self,x,y,z,ID,ty)
    self.soma.insert('k_ion')
    self.soma.insert('na_ion')
    self.soma.insert('ca_ion')
    # using an ohmic current rather than GHK equation
    h.ion_style("ca_ion",0,1,0,0,0)
    self.soma.Ra = 100.0
    self.soma.nseg = 1
    self.soma.diam = 70.0 # one-compartment of 14260 um2
    self.soma.L = 64.86
    self.soma.insert('pas') # leak current 
    self.soma.insert('hh2ad') # Hodgin-Huxley INa and IK -- HH2.mod
    self.soma.insert('itre') # reticular IT current -- IT2.mod
    self.soma.insert('kl') # K leak -- kl.mod
    self.soma.insert('cadad') # calcium decay    
    self.soma.e_pas = -77
    self.soma.g_pas = 5e-5
    self.soma.ena= 50
    self.soma.ek = -95
    self.soma.gnabar_hh2ad = 90e-3 # 0.1
    self.soma.vtraub_hh2ad = -50.0 # adjusts threshold
    self.soma.gkbar_hh2ad = 0.01
    self.soma.gmax_itre = 2e-3
    h.erev_kl = self.soma.ek
    self.soma.gmax_kl=3e-6
    h.q10h_itre = 3.0
    #shape_soma(self)

#   def set_synapses (self):
#     self.somaGABAf 	= GABAAFast(sect=self.soma, loc=0.5)
#     self.somaGABAss	= GABAASlow(sect=self.soma, loc=0.5)
#     self.somaNMDA 	= SynapseNMDA(sect=self.soma, loc=0.5)
#     if STDP: self.somaAMPAf = SynapseSTDP(sect=self.soma,loc=0.5,tau=5.35,e=0,dtau=34,ptau=17,d=0.5,p=0.5)
#     else: self.somaAMPAf = AMPAFast(sect=self.soma, loc=0.5)
  
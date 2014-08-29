#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('TkAgg')

from numpy import arange, ndarray
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import Tkinter as Tk
import sys

#from pylab import *
#from Tkinter import *

class Simulations(): # class class class class class class class class
  def __init__(self):
    c   = 3.0e+11        # speed of light, um/ms
    cs  = 2.8e-11        # cross section, um^2 
    n   = 1.5
    self.WI  = (c/n)*cs       # 18.0,  um^3/ms   
    self.Nv  = 5.0e+7         # 5.0e+7 1/um^3 - density of active particles
    self.W2  = 1./0.23        # 1/ms - 1/spontaneous emission lifetime (2->1 transitions)
    self.W3  = 1.e-3*self.W2  # 1/ms - 1/spontaneous emission lifetime (3->2 transitions)
    self.eta = 1.e-16
    tb, te, self.ts = 0.0, 1.0, 1.0e-4

    self.x  = arange(tb,te,self.ts)
    self.Np = len(self.x)
    self.y1 = ndarray(shape=(self.Np), dtype='float')
    self.y2 = ndarray(shape=(self.Np), dtype='float')
#    self.y3 = ndarray(shape=(self.Np), dtype='float')
    self.y4 = ndarray(shape=(self.Np), dtype='float')

    self.PlotWindow = Tk.Toplevel()
    self.PlotWindow.title('numerical simulation of kinetic equations')
    fig = Figure(figsize=(10,6), dpi=100)
    self.g1 = fig.add_subplot(211)
    self.g2 = fig.add_subplot(212)
    self.canvas = FigureCanvasTkAgg(fig, self.PlotWindow)
    self.canvas.show()
    self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.PlotWindow)
    self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

  def Plot(self,event):
    
    def Eq1(N1,N):     
      return -N1*(Wp+0.5*self.W2) + 0.5*self.W2*(self.Nv-N)
    def Eq2(N1,N,U):   
      return  N1*(Wp+self.W3-0.5*self.W2) + 0.5*self.W2*(self.Nv-N) - self.W3*(self.Nv+N) - 2.0*N*U*self.WI
    def Eq3(N,U):      
      if Tc>0.0: return  self.WI*N*U - U/Tc 
      else:      return  0.0 
    
    Wp = 1.e-3*self.WPv.get()  # pump power 1/ms
    Tc = 1.e-6*self.TLv.get()  # photon life time, ms
    Zn = Wp*self.W2 + Wp*self.W3 + self.W2*self.W3
    if Wp==0.0:  Ne = self.Nv
    else:        Ne = self.Nv*Wp*(self.W2-self.W3)/Zn
    To  = (2.*Wp+self.W2)/Zn
    N1, N, U  = self.Nv, 0.0, 0.0
    if Tc==0.0:  No = 0.0;               Uo = 1.0;            Wt = 1.0
    else:        No = 1.0/(self.WI*Tc);  Uo = Tc/To*(Ne-No);  Wt = self.W3*No/(self.Nv-No)
#    print 'Ne/Nv=%.2e  No/Nv=%.2e  Uo=%.2e  Wp/Wt=%2e' % (Ne/self.Nv, No/self.Nv, Uo, Wp/Wt)
    dt = self.ts
    p  = 0
    for t in self.x:
      self.y1[p] =  N/Ne
      self.y2[p] = No/Ne
      if U/Uo<50.0: self.y4[p] = U/Uo
      else: self.y4[p] = 0.0
      p+=1
      n11 = dt*Eq1(N1,           N         );   n1 = dt*Eq2(N1,           N,          U         );   u1 = dt*Eq3(N,          U         )
      n12 = dt*Eq1(N1 + 0.5*n11, N + 0.5*n1);   n2 = dt*Eq2(N1 + 0.5*n11, N + 0.5*n1, U + 0.5*u1);   u2 = dt*Eq3(N + 0.5*n1, U + 0.5*u1)
      n13 = dt*Eq1(N1 + 0.5*n12, N + 0.5*n2);   n3 = dt*Eq2(N1 + 0.5*n12, N + 0.5*n2, U + 0.5*u2);   u3 = dt*Eq3(N + 0.5*n2, U + 0.5*u2)
      n14 = dt*Eq1(N1 +     n13, N +     n3);   n4 = dt*Eq2(N1 +     n13, N +     n3, U +     u3);   u4 = dt*Eq3(N +     n3, U +     u3)
      dN1 = (n11 + 2.*n12 + 2.*n13 + n14)/6.;   dN = (n1 + 2.*n2 + 2.*n3 + n4)/6.;                   dU = (u1 + 2.*u2 + 2.*u3 + u4)/6.
      N1 += dN1
      N  += dN
      if N<No and U<=self.eta*Uo:  U  = self.eta*Uo
      else:                        U += dU

    self.g1.clear()
    self.g1.grid(True)
    self.g1.plot(self.x,100.*self.y1,color='b')
    self.g1.plot(self.x,100.*self.y2,color='r')
    self.g1.set_ylabel('inverse population') 
    self.g1.legend(["N/Ne, %", "No/Ne, %"],loc=2)

    self.g2.clear()
    self.g2.grid(True)
    self.g2.set_xlabel('time, ms') 
    self.g2.set_ylabel('lasing power') 
    self.g2.plot(self.x,self.y4,color='b')
    self.g2.legend(["U/Uo"],loc=2)
    self.canvas.show()
    self.toolbar.update()
          
class Application(Tk.Frame,Simulations): # class class class class class class class class
  def __init__(self, master=None):
    Tk.Frame.__init__(self, master)
    self.grid()
    Simulations.__init__(self)

    self.WPv, self.TLv = Tk.DoubleVar(), Tk.DoubleVar()
    self.S1 = Tk.Scale(self, orient=Tk.HORIZONTAL, length=400, label='Pumping speed, 1/us', 
                    variable=self.WPv, from_=0.0,   to_=5.0,  resolution=0.05, tickinterval=2.5,   command=self.Plot)
    self.S2 = Tk.Scale(self, orient=Tk.HORIZONTAL, length=400, label='Photon lifetime in the laser cavity, ns', 
                    variable=self.TLv, from_=0.0,   to_=50.0,  resolution=1.0, tickinterval=25.0,    command=self.Plot)
    self.S1.set(1.0);                self.S2.set(0.0)
    self.S1.grid( row=0, column=0);  self.S2.grid( row=1, column=0)

app = Application()
app.master.title("Nd:YAG laser")
app.mainloop()


#!/usr/bin/env python3
import sys
from matplotlib.backends.qt_compat import QtCore, QtWidgets
if QtCore.qVersion() >= "5.": from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
else:                         from matplotlib.backends.backend_qt4agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class KineticEquations:
  def __init__(self):
    c   = 3.0e+11             # speed of light  [μm/ms]
    cs  = 2.8e-11             # induced radiation cross section [μm^2] 
    n   = 1.82                # refraction index
    self.Wi  = (c/n)*cs       # induced radiation constant [μm^3/ms]   
    self.Nv  = 1.38e+8        # density of active particles at 1% doping [1/μm^3]
    self.Nv *= 1.0            # doping concentration acount
    self.W32 = 1./0.23        # 1/ms - 1/spontaneous emission lifetime (3=>2 transitions)
    self.W21 = 1e+4*self.W32  # 1/ms - 1/spontaneous emission lifetime (2=>1 transitions)
    self.eta = 1.e-16         # U/Ucw ratio that is treated as zero
    self.dt  = 5e-5           # time step for calculations [ms]
    self.Np  = 21000          # number of points in calculations
    self.x  = [self.dt*i for i in range(self.Np)]
    self.y1 = [0.0       for i in range(self.Np)]
    self.y2 = [0.0       for i in range(self.Np)]
    self.y3 = [0.0       for i in range(self.Np)]
  
  def Solve(self, Wp, Tc):
    def Eq1(N1,N):     
      return -N1*(Wp+0.5*self.W21) + 0.5*self.W21*(self.Nv-N)
    def Eq2(N1,N,U):   
      return  N1*(Wp+self.W32-0.5*self.W21) + 0.5*self.W21*(self.Nv-N) - self.W32*(self.Nv+N) - 2.0*N*U*self.Wi
    def Eq3(N,U):      
      if Tc>0.0: return  self.Wi*N*U - U/Tc 
      else:      return  0.0 
    
    Zn = Wp*(self.W21 + self.W32) + self.W21*self.W32
    No = Wp*(self.W21 - self.W32)/Zn*self.Nv
    To = (2.*Wp + self.W21)/Zn
    if Tc==0.0:  Ncw = 0.0;              Ucw = 0.0            
    else:        Ncw = 1.0/(self.Wi*Tc); Ucw = Tc/To*(No-Ncw)
    N1, N, U  = self.Nv, 0.0, 0.0
    dt = self.dt
    for p in range(self.Np):
      if p>20000: Wp = 0.0
      if Ncw > 0: self.y1[p] = N/Ncw; self.y2[p] = Ncw/No
      else:       self.y1[p] = N/No;  self.y2[p] = None
      if Ucw > 0: self.y3[p] = U/Ucw
      else:       self.y3[p] = 0.0
      n11 = dt*Eq1(N1,           N         );   n1 = dt*Eq2(N1,           N,          U         );   u1 = dt*Eq3(N,          U         )
      n12 = dt*Eq1(N1 + 0.5*n11, N + 0.5*n1);   n2 = dt*Eq2(N1 + 0.5*n11, N + 0.5*n1, U + 0.5*u1);   u2 = dt*Eq3(N + 0.5*n1, U + 0.5*u1)
      n13 = dt*Eq1(N1 + 0.5*n12, N + 0.5*n2);   n3 = dt*Eq2(N1 + 0.5*n12, N + 0.5*n2, U + 0.5*u2);   u3 = dt*Eq3(N + 0.5*n2, U + 0.5*u2)
      n14 = dt*Eq1(N1 +     n13, N +     n3);   n4 = dt*Eq2(N1 +     n13, N +     n3, U +     u3);   u4 = dt*Eq3(N +     n3, U +     u3)
      dN1 = (n11 + 2.*n12 + 2.*n13 + n14)/6.;   dN = (n1 + 2.*n2 + 2.*n3 + n4)/6.;                   dU = (u1 + 2.*u2 + 2.*u3 + u4)/6.
      N1 += dN1
      N  += dN
      if N<Ncw and U<=self.eta*Ucw: U  = self.eta*Ucw
      else:                         U += dU

class ApplicationWindow(QtWidgets.QMainWindow):
  lbls = {'ru': ['скорость накачки: {:5.2f} [1/мкс]', 'время жизни фотона в резонаторе: {:5.2f} [нс]',
                 'инверсная населенность', 'время, мс', 'мощность излучения'],
          'en': ['pump speed: {:5.2f} [1/μs]',        'photon lifetime: {:5.2f} [ns]',
                 'inverse population',     'time, ms',  'radiation power']}        
  def __init__(self, lg='en'):
    super().__init__()
    self._main = QtWidgets.QWidget()
    self.lbls = self.lbls[lg]
    self.setCentralWidget(self._main)
    layout   = QtWidgets.QVBoxLayout(self._main)
    self.canvas   = FigureCanvas(Figure(figsize=(8, 6)))
    self.controls = [QtWidgets.QLabel(), QtWidgets.QSlider(QtCore.Qt.Horizontal),
                     QtWidgets.QLabel(), QtWidgets.QSlider(QtCore.Qt.Horizontal)]
    layout.addWidget(self.canvas)
    self.addToolBar(NavigationToolbar(self.canvas, self))
    self.SaveAs = QtWidgets.QPushButton('Save raw data')
    self.canvas.toolbar.addWidget(self.SaveAs)
    for el in self.controls: 
      layout.addWidget(el)
      el.setMaximumHeight(25)
    for w in [1,3]:
      self.controls[w].setMinimum(0)
      self.controls[w].setMaximum(199)
      self.controls[w].setSingleStep(1)
    self.controls[1].valueChanged.connect(self._Pump)
    self.controls[3].valueChanged.connect(self._Time)
    self.SaveAs.pressed.connect(self._SaveAs)
    self.SdM = KineticEquations()
    self.population = self.canvas.figure.add_subplot(211)
    self._radiation = self.canvas.figure.add_subplot(212)
    self.controls[1].setValue(10)
    self.controls[3].setValue(10)

  def _SaveAs(self):
    fileName = QtWidgets.QFileDialog.getSaveFileName(self, "Save raw data", "./", "Data Files (*.txt *.dat)")
    Np = len(self.SdM.x)
    try:
      with open(fileName[0],'w') as fp:
        if self.SdM.y2[0] is None:
          fp.write('# t[ms]   N/No\n')
          for i in range(Np):
            fp.write("{:08.5f}  {:08.5f}\n".format(self.SdM.x[i], self.SdM.y1[i]))
        else:
          fp.write('# t[ms]   N/Ncw     Ncw/No    U/Ucw\n')
          for i in range(Np):
            fp.write("{:08.5f}  {:08.5f}  {:08.5f}  {:08.5f}\n".format(self.SdM.x[i], self.SdM.y1[i], self.SdM.y2[i], self.SdM.y3[i]))
    except FileNotFounError:
      pass
     
  def _Pump(self, value):
    self.controls[0].setText(self.lbls[0].format(0.025*value))
    self.Solve()

  def _Time(self, value):
    self.controls[2].setText(self.lbls[1].format(0.05*value))
    self.Solve()

  def Solve(self):
    Wp = self.controls[1].value()*5e-5  # pump power 1/ms
    Tc = self.controls[3].value()*1e-7  # photon life time, ms
    self.SdM.Solve(Wp, Tc)
    self.population.clear()
    self.population.grid(True)
    self.population.set_ylabel(self.lbls[2]) 
    if self.SdM.y2[0]==None:
      self.population.plot(self.SdM.x, self.SdM.y1, color='b')
      self.population.legend(["$N/N_0$"], loc=2)
    else:
      self.population.plot(self.SdM.x, self.SdM.y1, color='b')
      self.population.plot(self.SdM.x, self.SdM.y2, color='r')
      self.population.legend(["$N/N_{cw}$", "$N_{cw}/N_0$"], loc=2)

    self._radiation.clear()
    self._radiation.grid(True)
    self._radiation.set_xlabel(self.lbls[3]) 
    self._radiation.set_ylabel(self.lbls[4]) 
    self._radiation.plot(self.SdM.x, self.SdM.y3, color='b')
    self._radiation.legend(["$U/U_{cw}$"],loc=2)
    self.canvas.draw()
  
  
if __name__ == "__main__":
  # Check whether there is already a running QApplication (e.g., if running from an IDE).
  qapp = QtWidgets.QApplication.instance()
  lang = 'en'
  if not qapp:
    qapp = QtWidgets.QApplication(sys.argv)
    if len(sys.argv)>1 and 'ru' in sys.argv[1]:
      lang = 'ru' 
  app = ApplicationWindow(lg=lang)
  app.show()
  app.activateWindow()
  app.raise_()
  qapp.exec_()
 

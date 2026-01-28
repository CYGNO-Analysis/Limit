import uproot
import array as array
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from ROOT import *

file = uproot.open("results/test/test_mcmc.root")

phase = file['test_mcmc']['Phase'].array(library = 'np')
mub = file['test_mcmc']['mub'].array(library = 'np')
mus = file['test_mcmc']['mus'].array(library = 'np')

mub = mub[np.where(phase>0)]
mus = mus[np.where(phase>0)]

#ci = stats.t.interval(0.90, len(mus)-1)
ci = np.quantile(mus, 0.90)
cib = np.quantile(mub, 0.90)
print('quantile 90% :',ci)
print('quantile 90% :',cib)

plt.hist(mub, bins = 50,  density=True)
plt.xlabel('mub')
plt.savefig('mub.pdf')
plt.show()

plt.hist(mus, bins = 50 ,  density=True)
plt.xlabel('mus')
plt.savefig('mus.pdf')
plt.show()
"""
c = TCanvas()
h = TH1F("h","h",50,0.,2000.)
h1 = TH1F("h1","h1",50,105000.,106500.)
for i in mus:
    h.Fill(i)
for i in mub:
    h1.Fill(i)

c.Draw()
h.DrawNormalized()

c1=TCanvas()
c1.Draw()
h1.DrawNormalized()
"""

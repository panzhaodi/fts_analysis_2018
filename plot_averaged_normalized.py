import pylab as pl
from scipy import interpolate
from glob import glob
from pylab import *
import pylab as pl
import pickle
import glob
from scipy import interpolate
wafer=['w174', 'w176', 'w177', 'w187', 'w188', 'w201', 'w180', 'w203', 'w172', 'w181']  
#colorlist= [ 'red', 'hotpink','blue', 'teal','cyan', 'green', 'lime','orange', 'lightsalmon','purple', 'mediumorchid','brown',]
colorlist= [ 'chocolate', 'hotpink','gold', 'red', 'blue', 'green','orchid', 'burlywood','darksalmon', 'gray','darkviolet', 'yellowgreen', 'whitesmoke']

wafer=wafer[9]

'''
This code plot the averaged normalized band for every wafer. 
'''

#This will remove the dark bolometers from the plotting list
def filter_spectra(datafile):
      data=pickle.load(open(datafile, 'rb'))
      boloname= datafile.split('/')[-1].split('center_pixel')[0]
      wafer=boloname.split('_')[0].split('w')[1]
      pixel=boloname.split('.')[0].split('_')[1]
      band= boloname.split('.')[-2]
      pol=  boloname.split('.')[-1]
      fts_data= data['avg']
      #print wafer, pixel, band ,polavg_
      ind=pl.where((fts_data['freq']>70) & (fts_data['freq']<275))[0]
      low=pl.where((fts_data['freq']>10) & (fts_data['freq']<50))[0]
      harmonic=pl.where((fts_data['freq']>(2.*(int(band)-30))) & (fts_data['freq']<(2.*(int(band)+30))))[0]
      out_ind=pl.concatenate((low,harmonic))
      passed=True
      if pl.mean(abs(fts_data['real'][out_ind])/max(fts_data['real'][ind]))>0.1 or  pl.mean(abs(fts_data['real'][ind]))/pl.mean(abs(fts_data['stdev'][ind]))<5:
         passed=False

      if int(wafer)>=172 and int(wafer)<=179 and pixel in ['266', '225', '75', '6', '47', '197']:
         passed=False
      if int(wafer)>=180 and int(wafer)<=184 :
         if pixel in ['266', '225', '75', '6', '47', '197']:
            passed=False
         if pixel in ['261', '238', '211', '180', '145'] and band+'.'+pol in ['220.x', '220.y', '150.y', '90.y']:
            passed=False
      if int(wafer)>=185 : 
         if pixel in ['271', '250', '225', '196', '163']:
            passed=False
         if pixel in ['261', '238', '211', '180', '145'] and band+'.'+pol in ['220.x', '220.y', '150.y', '90.y']:
            passed=False
      if passed== False:
        print 'filtered'+datafile
      return passed    




      





# simulation files
sim90=np.loadtxt('./spectra/500nm_90')
sim150=np.loadtxt('./spectra/500nm_150')
sim220=np.loadtxt('./spectra/500nm_220')
pl.figure()
pl.clf()
pl.plot(sim150[:,0],(sim150[:,1]**2)/max(sim150[:,1]**2),'k--',label='Simulation(500nm ${SiO_2}$)')
pl.plot(sim220[:,0],(sim220[:,1]**2)/max(sim220[:,1]**2),'k--')
pl.plot(sim90[:,0],(sim90[:,1]**2)/max(sim90[:,1]**2),'k--')

#All the spectra files should have been generated. 
#Firstly plot the centered pixels. 
data_90ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/zero/'+wafer+'*.90.*.pkl')
data_150ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/zero/'+wafer+'*.150.*.pkl')
data_220ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/zero/'+wafer+'*.220.*.pkl')


#print data_90ghz

#plot averaged spectra


data_file=data_90ghz[0]
data=pickle.load(open(data_file,'rb'))['avg']
avg_90=np.zeros(len(data['freq'][np.where(data['freq']<300)]))
avg_150=np.zeros(len(data['freq'][np.where(data['freq']<300)]))
avg_220=np.zeros(len(data['freq'][np.where(data['freq']<300)]))


for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      band=  pl.where((data['freq']>70) & (data['freq']<150))[0]
      avg_90=avg_90+np.array(data['real'])/max(np.array(data['real'])[band])
      #pl.plot(data['freq'],np.array(data['real'])/max(np.array(data['real'])[band]) ,color='b',lw=0.5) 
  except EOFError:
      print data_file

pl.plot(data['freq'], avg_90/max(avg_90[band]), color='b')

for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      band=  pl.where((data['freq']>100) & (data['freq']<200))[0]
      avg_150=avg_150+np.array(data['real'])/max(np.array(data['real'])[band])
  except EOFError:
      print data_file
pl.plot(data['freq'], avg_150/max(avg_150[band]), color='b')

for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      #pl.figure()
      band=  pl.where((data['freq']>150) & (data['freq']<280))[0]
      avg_220=avg_220+np.array(data['real'])/max(np.array(data['real'])[band])
  except EOFError:
      print data_file

pl.plot(data['freq'], avg_220/max(avg_220[band]), color='b')


data_90ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/one/'+wafer+'*.90.*.pkl')
data_150ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/one/'+wafer+'*.150.*.pkl')
data_220ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/one/'+wafer+'*.220.*.pkl')

#plot averaged spectra for the neighbor. 

for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      band=  pl.where((data['freq']>70) & (data['freq']<150))[0]
      avg_90=avg_90+np.array(data['real'])/max(np.array(data['real'])[band])
      #pl.plot(data['freq'],np.array(data['real'])/max(np.array(data['real'])[band]) ,color='b',lw=0.5) 
  except EOFError:
      print data_file

pl.plot(data['freq'], avg_90/max(avg_90[band]), color='g')

for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      band=  pl.where((data['freq']>100) & (data['freq']<200))[0]
      avg_150=avg_150+np.array(data['real'])/max(np.array(data['real'])[band])
  except EOFError:
      print data_file
pl.plot(data['freq'], avg_150/max(avg_150[band]), color='g')

for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      #pl.figure()
      band=  pl.where((data['freq']>150) & (data['freq']<280))[0]
      avg_220=avg_220+np.array(data['real'])/max(np.array(data['real'])[band])
  except EOFError:
      print data_file

pl.plot(data['freq'], avg_220/max(avg_220[band]), color='g')

data_90ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/two/'+wafer+'*.90.*.pkl')
data_150ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/two/'+wafer+'*.150.*.pkl')
data_220ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/two/'+wafer+'*.220.*.pkl')

#plot averaged spectra for the second neighbor. 

for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      band=  pl.where((data['freq']>70) & (data['freq']<150))[0]
      avg_90=avg_90+np.array(data['real'])/max(np.array(data['real'])[band])
      #pl.plot(data['freq'],np.array(data['real'])/max(np.array(data['real'])[band]) ,color='b',lw=0.5) 
  except EOFError:
      print data_file

pl.plot(data['freq'], avg_90/max(avg_90[band]), color='r')

for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      band=  pl.where((data['freq']>100) & (data['freq']<200))[0]
      avg_150=avg_150+np.array(data['real'])/max(np.array(data['real'])[band])
  except EOFError:
      print data_file
pl.plot(data['freq'], avg_150/max(avg_150[band]), color='r')

for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<300)]
      data['real']=data['real'][np.where(data['freq']<300)]
      #pl.figure()
      band=  pl.where((data['freq']>150) & (data['freq']<280))[0]
      avg_220=avg_220+np.array(data['real'])/max(np.array(data['real'])[band])
  except EOFError:
      print data_file

pl.plot(data['freq'], avg_220/max(avg_220[band]), color='r')


pl.plot([],color='b', label='Center')
pl.plot([],color='g', label='Neighbor')
pl.plot([],color='r', label='Second Neighbor')
pl.axis([50,300, -.2, 1])
pl.legend(loc=3, ncol=2)
pl.grid()
pl.title(wafer+ ' normalized averaged spectra')
#
pl.xlabel('Frequency(GHz)')
pl.ylabel('Normalized unit')
pl.savefig(wafer+'averagednormalized')

pl.show()



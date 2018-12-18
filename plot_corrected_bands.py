import pylab as pl
from scipy import interpolate
from glob import glob
from pylab import *
import pylab as pl
import pickle
import glob
from scipy import interpolate
wafer_lst=['w174', 'w176', 'w177', 'w187', 'w188', 'w201', 'w180', 'w203', 'w172', 'w181']  
wafer=wafer_lst[8]
#colorlist= [ 'red', 'hotpink','blue', 'teal','cyan', 'green', 'lime','orange', 'lightsalmon','purple', 'mediumorchid','brown',]
colorlist= [ 'chocolate', 'hotpink','gold', 'red', 'blue', 'green','orchid', 'burlywood','darksalmon', 'gray','darkviolet', 'yellowgreen', 'whitesmoke']
'''
This code plots the bands corrected for the Lyot stop, Rayleigh-Jeans source, AOmega of the detector, and the nu^2 from the Mylar filter. All corrections are encapculated in fts_correction.txt. Brad generated this file (ask him if you have questions). So the corrected data should only be the spectrum of the detector (plus the lenses, lenslets, and metal mesh filters. )

This is the spectra that you should compare to simulations. 
'''

# Remove dark bolos
def filter_spectra(datafile):
    try:
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
    except:
      return False 



correction =  zip(*np.loadtxt('fts_correction.txt', usecols=[0,1]))
correction_function=  interpolate.UnivariateSpline(correction[0],correction[1], s=0) 
'''
frequency=np.linspace(8,300,2000)
pl.figure()
pl.plot(frequency, correction_function(frequency)/max(correction_function(frequency)))
pl.xlabel('Frequency')
pl.ylabel('Correction factor, normalized')
pl.title('Correction factor (blackbody, mylar, lyot stop)')
'''

#
pl.figure(figsize=[9.5,7])
pl.clf()

sim90=np.loadtxt('./spectra/500nm_90')
sim150=np.loadtxt('./spectra/500nm_150')
sim220=np.loadtxt('./spectra/500nm_220')
 
#pl.plot(sim150[:,0],(sim150[:,1]**2)/max(sim150[:,1]**2),'k--',linewidth=2,label='Simulation for 500nm $  {SiO_2}$')
#pl.plot(sim220[:,0],(sim220[:,1]**2)/max(sim220[:,1]**2),'k--',linewidth=2)
#pl.plot(sim90[:,0],(sim90[:,1]**2)/max(sim90[:,1]**2),'k--',linewidth=2)

data_90ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/zero/'+wafer+'*.90.*.pkl')
data_150ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/zero/'+wafer+'*.150.*.pkl')
data_220ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/zero/'+wafer+'*.220.*.pkl')

frequency=np.linspace(8,300,2000)
band_avg_90=np.zeros(len(frequency))
band_avg_150=np.zeros(len(frequency))
band_avg_220=np.zeros(len(frequency))
#n220=0
for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<200)]
      data['real']=data['real'][np.where(data['freq']<200)]
      #data['freq']=data['freq'][np.where(data['freq']>65)]
      #data['real']=data['real'][np.where(data['freq']>65)]
      func90=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_90=band_avg_90+func90(frequency)
      #n90=n90+1
  except EOFError:
      print data_file
for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<200)]
      data['real']=data['real'][np.where(data['freq']<200)]
      #data['freq']=data['freq'][np.where(data['freq']>75)]
      #data['real']=data['real'][np.where(data['freq']>75)]
      func150=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_150=band_avg_150+func150(frequency)
      #n150=n150+1
  except EOFError:
      print data_file
for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      #pl.figure()
      data['freq']=data['freq'][np.where(data['freq']<270)]
      data['real']=data['real'][np.where(data['freq']<270)]
      #data['freq']=data['freq'][np.where(data['freq']>150)]
      #data['real']=data['real'][np.where(data['freq']>150)]
      func220=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_220=band_avg_220+func220(frequency)
      #n220=n220+1
      #pl.title(data_file)
  except EOFError:
      print data_file


#plot averaged spectra

frequency=np.linspace(65,285,2000)
band_avg_90=np.zeros(len(frequency))
band_avg_150=np.zeros(len(frequency))
band_avg_220=np.zeros(len(frequency))
#n90=0
#n150=0
#n220=0
for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      func90=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_90=band_avg_90+func90(frequency)
      #n90=n90+1
  except EOFError:
      print data_file
for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      func150=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_150=band_avg_150+func150(frequency)
      #n150=n150+1
  except EOFError:
      print data_file
for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      #pl.figure()
      func220=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_220=band_avg_220+func220(frequency)
      #n220=n220+1
      #pl.title(data_file)
  except EOFError:
      print data_file
ndf=np.exp(-0.3*(frequency/29.9792)**1.2*1.27)
band_avg_90=band_avg_90 / correction_function(frequency)
band_avg_150=band_avg_150 / correction_function(frequency)
band_avg_220=band_avg_220 / correction_function(frequency)
band_center_90= np.sum(np.array(band_avg_90)*np.array(frequency))/np.sum(band_avg_90)
band_width_90=np.sum(band_avg_90)*(frequency[1]-frequency[0])/max(band_avg_90)
band_center_150= np.sum(np.array(band_avg_150)*np.array(frequency))/np.sum(band_avg_150)
band_width_150=np.sum(band_avg_150)*(frequency[1]-frequency[0])/max(band_avg_150)
band_center_220= np.sum(np.array(band_avg_220)*np.array(frequency))/np.sum(band_avg_220)
band_width_220=np.sum(band_avg_220)*(frequency[1]-frequency[0])/max(band_avg_220)
ind90=np.where(np.logical_and(frequency<150, frequency>50))
ind150=np.where(np.logical_and(frequency<200, frequency>100))
ind220=np.where(np.logical_and(frequency<300, frequency>150))
data_to_store={}
data_to_store['center']={}
data_to_store['neighbor']={}
data_to_store['second']={}
for key in data_to_store.keys():
 data_to_store[key]['90']={}
 data_to_store[key]['150']={}
 data_to_store[key]['220']={}
 


pl.plot(frequency[ind90], band_avg_90[ind90]/max(band_avg_90[ind90]),label='Centered pixels',color='b')
pl.plot(frequency[ind150], band_avg_150[ind150]/max(band_avg_150[ind150]),color='b')
pl.plot(frequency[ind220], band_avg_220[ind220]/max(band_avg_220[ind220]),color='b')

data_to_store['center']['90']['freq']= frequency[ind90]
data_to_store['center']['150']['freq']= frequency[ind150]
data_to_store['center']['220']['freq']= frequency[ind220]
data_to_store['center']['90']['spectrum']= band_avg_90[ind90]/max(band_avg_90[ind90])
data_to_store['center']['150']['spectrum']= band_avg_150[ind150]/max(band_avg_150[ind150])
data_to_store['center']['220']['spectrum']= band_avg_220[ind220]/max(band_avg_220[ind220])


pl.plot(sim150[:,0],(sim150[:,1]**2)/max(sim150[:,1]**2),'k--',linewidth=2 )
pl.plot(sim220[:,0],(sim220[:,1]**2)/max(sim220[:,1]**2),'k--',linewidth=2)
pl.plot(sim90[:,0],(sim90[:,1]**2)/max(sim90[:,1]**2),'k--',linewidth=2)



############Plot neighbor
data_90ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/one/'+wafer+'*.90.*.pkl')
data_150ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/one/'+wafer+'*.150.*.pkl')
data_220ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/one/'+wafer+'*.220.*.pkl')

frequency=np.linspace(8,300,2000)
band_avg_90=np.zeros(len(frequency))
band_avg_150=np.zeros(len(frequency))
band_avg_220=np.zeros(len(frequency))
#n220=0
for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<200)]
      data['real']=data['real'][np.where(data['freq']<200)]
      #data['freq']=data['freq'][np.where(data['freq']>65)]
      #data['real']=data['real'][np.where(data['freq']>65)]
      func90=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_90=band_avg_90+func90(frequency)
      #n90=n90+1
  except EOFError:
      print data_file
for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<200)]
      data['real']=data['real'][np.where(data['freq']<200)]
      #data['freq']=data['freq'][np.where(data['freq']>75)]
      #data['real']=data['real'][np.where(data['freq']>75)]
      func150=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_150=band_avg_150+func150(frequency)
      #n150=n150+1
  except EOFError:
      print data_file
for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      #pl.figure()
      data['freq']=data['freq'][np.where(data['freq']<270)]
      data['real']=data['real'][np.where(data['freq']<270)]
      #data['freq']=data['freq'][np.where(data['freq']>150)]
      #data['real']=data['real'][np.where(data['freq']>150)]
      func220=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_220=band_avg_220+func220(frequency)
      #n220=n220+1
      #pl.title(data_file)
  except EOFError:
      print data_file


#plot averaged spectra

frequency=np.linspace(65,285,2000)
band_avg_90=np.zeros(len(frequency))
band_avg_150=np.zeros(len(frequency))
band_avg_220=np.zeros(len(frequency))
#n90=0
#n150=0
#n220=0
for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      func90=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_90=band_avg_90+func90(frequency)
      #n90=n90+1
  except EOFError:
      print data_file
for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      func150=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_150=band_avg_150+func150(frequency)
      #n150=n150+1
  except EOFError:
      print data_file
for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      #pl.figure()
      func220=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_220=band_avg_220+func220(frequency)
      #n220=n220+1
      #pl.title(data_file)
  except EOFError:
      print data_file
ndf=np.exp(-0.3*(frequency/29.9792)**1.2*1.27)
band_avg_90=band_avg_90/ correction_function(frequency)
band_avg_150=band_avg_150 / correction_function(frequency)
band_avg_220=band_avg_220 / correction_function(frequency)
band_center_90= np.sum(np.array(band_avg_90)*np.array(frequency))/np.sum(band_avg_90)
band_width_90=np.sum(band_avg_90)*(frequency[1]-frequency[0])/max(band_avg_90)
band_center_150= np.sum(np.array(band_avg_150)*np.array(frequency))/np.sum(band_avg_150)
band_width_150=np.sum(band_avg_150)*(frequency[1]-frequency[0])/max(band_avg_150)
band_center_220= np.sum(np.array(band_avg_220)*np.array(frequency))/np.sum(band_avg_220)
band_width_220=np.sum(band_avg_220)*(frequency[1]-frequency[0])/max(band_avg_220)
ind90=np.where(np.logical_and(frequency<150, frequency>50))
ind150=np.where(np.logical_and(frequency<200, frequency>100))
ind220=np.where(np.logical_and(frequency<300, frequency>150))

pl.plot(frequency[ind90], band_avg_90[ind90]/max(band_avg_90[ind90]),label='Neighbor pixels',color='g')
pl.plot(frequency[ind150], band_avg_150[ind150]/max(band_avg_150[ind150]),color='g')
pl.plot(frequency[ind220], band_avg_220[ind220]/max(band_avg_220[ind220]),color='g')
pl.plot(sim150[:,0],(sim150[:,1]**2)/max(sim150[:,1]**2),'k--',linewidth=2)
pl.plot(sim220[:,0],(sim220[:,1]**2)/max(sim220[:,1]**2),'k--',linewidth=2)
pl.plot(sim90[:,0],(sim90[:,1]**2)/max(sim90[:,1]**2),'k--',linewidth=2)


data_to_store['neighbor']['90']['freq']= frequency[ind90]
data_to_store['neighbor']['150']['freq']= frequency[ind150]
data_to_store['neighbor']['220']['freq']= frequency[ind220]
data_to_store['neighbor']['90']['spectrum']= band_avg_90[ind90]/max(band_avg_90[ind90])
data_to_store['neighbor']['150']['spectrum']= band_avg_150[ind150]/max(band_avg_150[ind150])
data_to_store['neighbor']['220']['spectrum']= band_avg_220[ind220]/max(band_avg_220[ind220])

############plot second
data_90ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/two/'+wafer+'*.90.*.pkl')
data_150ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/two/'+wafer+'*.150.*.pkl')
data_220ghz=glob.glob('/home/zpan/Desktop/pole_fts_2018/zpan_tmp/two/'+wafer+'*.220.*.pkl')

frequency=np.linspace(8,300,2000)
band_avg_90=np.zeros(len(frequency))
band_avg_150=np.zeros(len(frequency))
band_avg_220=np.zeros(len(frequency))
#n220=0
for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<200)]
      data['real']=data['real'][np.where(data['freq']<200)]
      #data['freq']=data['freq'][np.where(data['freq']>65)]
      #data['real']=data['real'][np.where(data['freq']>65)]
      func90=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_90=band_avg_90+func90(frequency)
      #n90=n90+1
  except EOFError:
      print data_file
for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      data['freq']=data['freq'][np.where(data['freq']<200)]
      data['real']=data['real'][np.where(data['freq']<200)]
      #data['freq']=data['freq'][np.where(data['freq']>75)]
      #data['real']=data['real'][np.where(data['freq']>75)]
      func150=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_150=band_avg_150+func150(frequency)
      #n150=n150+1
  except EOFError:
      print data_file
for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      #pl.figure()
      data['freq']=data['freq'][np.where(data['freq']<270)]
      data['real']=data['real'][np.where(data['freq']<270)]
      #data['freq']=data['freq'][np.where(data['freq']>150)]
      #data['real']=data['real'][np.where(data['freq']>150)]
      func220=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']>50)]), s=0) 
      band_avg_220=band_avg_220+func220(frequency)
      #n220=n220+1
      #pl.title(data_file)
  except EOFError:
      print data_file


#plot averaged spectra

frequency=np.linspace(65,285,2000)
band_avg_90=np.zeros(len(frequency))
band_avg_150=np.zeros(len(frequency))
band_avg_220=np.zeros(len(frequency))
#n90=0
#n150=0
#n220=0
for data_file in data_90ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      func90=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_90=band_avg_90+func90(frequency)
      #n90=n90+1
  except EOFError:
      print data_file
for data_file in data_150ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      func150=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_150=band_avg_150+func150(frequency)
      #n150=n150+1
  except EOFError:
      print data_file
for data_file in data_220ghz:
 if filter_spectra(data_file):
  try:
      data=pickle.load(open(data_file,'rb'))['avg']
      #pl.figure()
      func220=interpolate.UnivariateSpline(data['freq'],np.array(data['real'])/max(data['real'][np.where(data['freq']<400)]), s=0) 
      band_avg_220=band_avg_220+func220(frequency)
      #n220=n220+1
      #pl.title(data_file)
  except EOFError:
      print data_file
ndf=np.exp(-0.3*(frequency/29.9792)**1.2*1.27)
band_avg_90=band_avg_90/ correction_function(frequency)
band_avg_150=band_avg_150 / correction_function(frequency)
band_avg_220=band_avg_220 / correction_function(frequency)
band_center_90= np.sum(np.array(band_avg_90)*np.array(frequency))/np.sum(band_avg_90)
band_width_90=np.sum(band_avg_90)*(frequency[1]-frequency[0])/max(band_avg_90)
band_center_150= np.sum(np.array(band_avg_150)*np.array(frequency))/np.sum(band_avg_150)
band_width_150=np.sum(band_avg_150)*(frequency[1]-frequency[0])/max(band_avg_150)
band_center_220= np.sum(np.array(band_avg_220)*np.array(frequency))/np.sum(band_avg_220)
band_width_220=np.sum(band_avg_220)*(frequency[1]-frequency[0])/max(band_avg_220)

ind90=np.where(np.logical_and(frequency<150, frequency>50))
ind150=np.where(np.logical_and(frequency<200, frequency>100))
ind220=np.where(np.logical_and(frequency<300, frequency>150))

pl.plot(frequency[ind90], band_avg_90[ind90]/max(band_avg_90[ind90]),label='Second pixels',color='r')
pl.plot(frequency[ind150], band_avg_150[ind150]/max(band_avg_150[ind150]),color='r')
pl.plot(frequency[ind220], band_avg_220[ind220]/max(band_avg_220[ind220]),color='r')

pl.plot(sim150[:,0],(sim150[:,1]**2)/max(sim150[:,1]**2),'k--',linewidth=2,label='Simulation for 500nm $  {SiO_2}$')
pl.plot(sim220[:,0],(sim220[:,1]**2)/max(sim220[:,1]**2),'k--',linewidth=2)
pl.plot(sim90[:,0],(sim90[:,1]**2)/max(sim90[:,1]**2),'k--',linewidth=2)
#end

data_to_store['second']['90']['freq']= frequency[ind90]
data_to_store['second']['150']['freq']= frequency[ind150]
data_to_store['second']['220']['freq']= frequency[ind220]
data_to_store['second']['90']['spectrum']= band_avg_90[ind90]/max(band_avg_90[ind90])
data_to_store['second']['150']['spectrum']= band_avg_150[ind150]/max(band_avg_150[ind150])
data_to_store['second']['220']['spectrum']= band_avg_220[ind220]/max(band_avg_220[ind220])

pl.axis([50,300,-.2,1.1])
pl.legend(loc=3, ncol=2)
pl.xlabel('Frequency(GHz)')
pl.ylabel('Normalized amplitude')
pl.grid()
pl.title(wafer+ ' averaged spectra vs distance, corrected for everything')
#pl.clf()
pl.savefig(wafer+'averaged_spectra_corrected')
pickle.dump(data_to_store, open(wafer+'spectra_unnormalized_averaged.pkl', 'wb'))
pl.show()



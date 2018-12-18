'''
 Zhaodi Pan  Oct 31 2015
 This code is to do  real_time analysis for the FTS spectra we generated
'''
import pylab as pl
import sptpol_software.util.files as files
import cPickle as pkl
from glob import glob
from pyfts_package import *
import numpy as np
from netCDF4 import Dataset
import pickle

# For chopped scans, you might have to shift the frequencies by 60 Hz. I'm not sure why this is yet.
global   data_dir
data_dir='./zero/'

# This part runs the FTS analysis for a pixel, given the centered pixel, band, optical speed. 
# We set the optics correction to False and correct for it later.
# The centered pixel may not be the same as the boloname's pixel, because it might be a neighbor or a second neighbor. 

def fts(netcdf='W98_114px_20161023_0442',boloname='W98_114.150.x',center_pixel='114', band=150,speed=.5,chop=False,optics_correction=False,save=True,return_data=True,plot=True, num=1):


  rate = 152.587890625   #sampling rate in Hz
  #speed = optical velocity in cm/s

  bandpass = [band-30,band+30]   #approximate band of the detector, used to find interferrorgrams in code
  
  bolo_data=glob(netcdf+'.nc')[0]
  bolo_data=Dataset(bolo_data,mode='r')
  #print bolo_data.variables.keys()
  data_i=pl.array(bolo_data.variables[boloname+'_I'])
  data_q=pl.array(bolo_data.variables[boloname+'_Q'])
    
  cts = np.sqrt(abs(abs(np.array(data_i))**2+abs(np.array(data_q))**2))
  pl.plot(cts)
  if chop:
    cts=-1*cts
  #pl.figure()
  #pl.clf()
  #pl.plot(range(len(cts)),cts)   #plot it
  #cts = cts[1.5e4:*]  # if you want to cut some data points off the beginning


  out = analyze_fts_scan(cts, speed, band=bandpass, rate=rate, chop=chop,hann=True, pband = [10, 400],absv=False,plot=plot,length=7.5)  #process the FTS data# if it's high enough signal-to-noise drop the /  abs at the end
  #print cts

  avg = add_fts_scan(out)  #average the fts scans
  if avg!=0: 
   figure=pl.figure()
   pl.clf()
   pl.plot(avg['freq'], avg['real']/max(avg['real'])) # plot the average
   pl.xlim(10, 300)
   pl.xlabel('Frequency (GHz)')
   pl.ylabel('Normalized Response')
   pl.plot(avg['freq'], avg['im']/max(avg['real']), 'b--') #overplot the imagninary component to get an estimate of the noise
   pl.grid(b=True)
   pl.title('Averaged Spectra for '+boloname)
   pl.savefig(data_dir+ boloname+'center_pixel'+center_pixel+'spectrum.png')
   figure2=pl.figure()
   pl.clf()
   pl.plot(cts )
   pl.title(boloname)
   pl.savefig(data_dir+ boloname+'center_pixel'+center_pixel+'timestream.png')
   if optics_correction==True:
    optics=pkl.load(open('/home/cryo/Data/spectra/fts_analysis/Optics_chain.pkl','r'))
    ind1=pl.where((optics[0]>60) & (optics[0]<240))
    ind2=pl.where((avg['freq']>60) & (avg['freq']<240))
    ind3=pl.where((avg['freq']>60) & (avg['freq']<300))
    
    pl.figure(4)
    pl.plot(optics[0][ind1],(avg['real'][ind2]/optics[1][ind1])/max(avg['real'][ind2]/optics[1][ind1]))
    pl.plot(optics[0][ind1],(avg['im'][ind2]/optics[1][ind1])/max(avg['real'][ind2]/optics[1][ind1]),'b--')
    pl.grid(b=True)
    pl.xlabel('Frequency (GHz)')
    pl.ylabel('Normalized Response')
    pl.title('Optics-Corrected Avg Spectra for '+sq+'Ch'+str(chan))

   if return_data:
    if save:
#      pkl.dump([avg['freq'][ind3],avg['real'][ind3]/max(avg['real'][ind3])],open(data_dir+parser_dir+'/'+parser_dir.rpartition('/')[-1]+'_'+sq+ch+'_raw_spectrum.pkl','wb'))
#      pkl.dump([optics[0][ind1],(avg['real'][ind2]/optics[1][ind1])/max(avg['real'][ind2]/optics[1][ind1])],open(data_dir+parser_dir+'/'+parser_dir.rpartition('/')[-1]+'_'+sq+ch+'_corrected_spectrum.pkl','wb'))
#      file2=open(data_dir+parser_dir+'/'+parser_dir.rpartition('/')[-1]+'_'+sq+ch+'spectrum.txt','w')
#      for i in range(len(avg['freq'])):
#             file2.write("%s  %s\n" % (avg['freq'][i], avg['real'][i]/max(avg['real'])))
#      file2.close()
#      pkl.dump(avg,open(data_dir+parser_dir+'/'+parser_dir.rpartition('/')[-1]+'_'+sq+ch+'_spectrum.pkl','wb'))
      pkl.dump({'avg':avg, 'cts':{'i':np.array(data_i), 'q':np.array(data_q)}},open(data_dir+boloname+'center_pixel'+center_pixel+'.pkl','wb'))
    if optics_correction:  
      return [optics[0][ind1],(avg['real'][ind2]/optics[1][ind1])/max(avg['real'][ind2]/optics[1][ind1])]
    else:
      return avg #[(avg['freq'], avg['real']/max(avg['real']))]


# This function plot all the netcdf files. 
def netcdf_plot(netcdf_files={'157_40_75mms.nc':'157.40', '157_88_6mms.nc':'157.88', '157_88_75mms.nc':'157.88'  },distance_range=12):
#'157_40_75mms.nc':'157.40', '157_88_6mms.nc':'157.88', '157_88_75mms.nc':'157.88' 75mms

#'139_210_6mms.nc':'139.210', '157_40_6mms.nc':'157.40', '157_88_6mms_2.nc':'157.88','157_104_6mms.nc':'157.104'

#'136_15.nc':'136.15','136_255.nc':'136.255','136_50.nc':'136.50','139_50.nc':'139.50','139_57.nc':'139.57','139_75.nc':'139.75', '142_117.nc':'142.117', '142_148.nc':'142.148', '142_216.nc':'142.216', '142_216_2.nc':'142.216', '142_242.nc':'142.242', '142_242_2.nc':'142.242', '147_186.nc':'147.186', '147_186_2.nc':'147.186', '147_241.nc':'147.241', '147_241_2.nc':'147.241', '148_41.nc':'148.41', '148_41_2.nc':'148.41', '148_58.nc':'148.58', '148_124.nc':'148.124', '148_124_2.nc':'148.124', '148_200.nc':'148.200', '148_237.nc':'148.237', '148_244.nc':'148.244', '152_209.nc':'152.209', '152_224.nc':'152.224', '153_38.nc':'153.38', '153_117.nc':'153.117', '153_117_2.nc':'153.117', '153_148.nc':'153.148', '153_148_2.nc':'153.148', '153_148_3.nc':'153.148', '157_79.nc':'157.79', '157_232.nc':'157.232', '157_232_2.nc':'157.232', '157_232_3.nc':'157.232', '158_78.nc':'158.78', '158_78_2.nc':'158.78', '158_214.nc':'158.214', '162_40.nc':'162.40', '162_112.nc':'162.112', '162_112_2.nc':'162.112', '162_131.nc':'162.131'
##########'147_186.nc':'147.186','147_186_2.nc':'147.186','147_241.nc':'147.241', '147_241_2.nc':'147.241', '148_58.nc':'148.58', '148_58_2.nc':'148.58', '148_124.nc':'148.124', '148_124_2.nc':'148.124', '148_200.nc':'148.200', '148_237.nc':'148.237', '148_244.nc':'148.244', '152_224.nc':'152.224', '153_38.nc':'153.38', '153_117_1.nc':'153.117', '153_117_2.nc':'153.117', '153_148_1.nc':'153.148', '153_148_2.nc':'153.148', '153_148_3.nc':'153.148', '158_78_1.nc':'158.78', '158_78_2.nc':'158.78', '158_214.nc':'158.214', '162_112_1.nc':'162.112', '162_112_2.nc':'162.112' 
#'fts_first_cooldown_pixel_86_87_fts_scan_7.5mm_sec.nc':'92.87', 'fts_first_cooldown_pixel_160_161_fts_scan_7.5mm_sec.nc':'92.160', 'fts_first_cooldown_pixel_160_161_fts_scan_7.5mm_sec.nc':'92.161', 'fts_first_cooldown_pixel_16_28_fts_scan_7.5mm_sec.nc':'92.16', 'fts_first_cooldown_pixel_16_28_fts_scan_7.5mm_sec.nc':'92.28'
    #netcdf_files={'fts_first_cooldown_2_peaked_on_229_W92.nc':'92.229',  'fts_first_cooldown_3_peaked_up_pixel_108_actual_fts_data.nc':'92.108', 'fts_first_cooldown_3_pixel_108.nc':'92.108','fts_first_cooldown_4_peaking_up_pixel_105_106_107.nc':'92.105','fts_first_cooldown_4_peaking_up_pixel_105_106_107.nc':'92.106' ,'fts_first_cooldown_4_peaking_up_pixel_105_106_107.nc':'92.107','fts_first_cooldown_4_peaking_up_pixel_160_161.nc','92.160','fts_first_cooldown_4_peaking_up_pixel_160_161.nc','92.161','fts_first_cooldown_4_peaking_up_pixel_230_242.nc':'92.230','fts_first_cooldown_4_peaking_up_pixel_230_242.nc':'92.242','fts_first_cooldown_4_peaking_up_pixel_231_243.nc':'92.231','fts_first_cooldown_4_peaking_up_pixel_231_243.nc':'92.243','fts_first_cooldown_4_peaking_up_pixel_233_234.nc':'92.233','fts_first_cooldown_4_peaking_up_pixel_233_234.nc':'92.234','fts_first_cooldown_pixel_230_242_fts_scan.nc':'92.230','fts_first_cooldown_pixel_230_242_fts_scan.nc':'92.234' }
    #netcdf_files={'W94_254_20161020_2340.nc':'94.197',  'W94_231px_20161021_0017.nc':'94.197', 'W94_218px_20161021_0035.nc':'94.197','W95_249px_20161022_0510.nc':'95.197', 'W95_250px_20161022_0534.nc':'95.197','W95_245px_20161022_1727.nc':'95.197','W95_219px_20161022_1835.nc':'95.197'}
    #netcdf_files={'W98_px3_20161020_2047.nc':'98.3','W98_px13_20161020_2205.nc':'98.13','W98_px36_20161020_2237.nc':'98.36', 'W94_254_20161020_2340.nc':'94.254',  'W94_231px_20161021_0017.nc':'94.231', 'W94_218px_20161021_0035.nc':'94.218','W99_178px_20161021_0126.nc':'99.178', 'W99_179px_20161021_0149.nc': '99.179', 'W99_161px_20161021_0219.nc':'99.161' }
    # 20161020/20161020_232102_drop_bolos/data/IceBoard_0117.Mezz_2.ReadoutModule_1_OUTPUT.pkl
    #
    #'W95_249px_20161022_0510.nc':'95.249', 'W95_250px_20161022_0534.nc':'95.250'
    ########20161022_200208_drop_bolos
    #
    #'W98_px36_unsaturated_20161020_2247.nc':'98.36','W94_254px_unsaturated_20161020_2347.nc':'94.254', 'W99_179px_unsaturated_20161021_0149.nc':'99.179', 'W99_161px_unsaturated_20161021_0219.nc':'99.161'
    # 20161020/20161020_232102_drop_bolos/data/IceBoard_0117.Mezz_2.ReadoutModule_1_OUTPUT.pkl
    #
    #'W95_245px_20161022_1727.nc':'95.245','W95_219px_20161022_1835.nc':'95.219', 'W98_114px_20161022_2021.nc':'98.114'
    #20161022_200208_drop_bolos
    #
    #'W98_17px_20161023_0527.nc': '98.17' #,'W98_114px_20161023_0442.nc':'98.114'
    #fts/20161023/RELATIVE_OPTICAL_EFFICIENCY_TESTS/PLATE_ON_DROP_BOLOS/20161023_085308_drop_bolos
    '''
    Arguments: netcdf_files, distance_range.
    netcdf_files should be a dictionary in the format of {'location of the netcdf file': centered pixel number(int)}
    distance_range is the distance to the centered pixel you want to get spectrum on.
    '''
    pixel_position=pickle.load(open('pixel_positions.pkl','rb'))
    #print pixel_position[5]
    k=0
    for file_name in sorted(netcdf_files.keys()):
       center_pixel= netcdf_files[file_name].split('.')[1]
       waf=netcdf_files[file_name].split('.')[0]
       #print center_pixel
       #print center_pixel
       center_position= pixel_position[int(center_pixel)]
       #print center_pixel
       neighbor_pixel=[]
       all_pixel=[]
       centered_pixel=[]
       second_pixel=[]
       far_pixel=[]
       for pixel in pixel_position.keys():
          coord= pixel_position[pixel]
          dist= ((coord[0]-center_position[0])**2+(coord[1]-center_position[1])**2)**(0.5)
          if dist<10 and dist>=1:
              neighbor_pixel.append(pixel)
          if dist<12:
              all_pixel.append(pixel)
          if dist<1:
              centered_pixel.append(pixel)
          if dist>10 and dist<12:
              second_pixel.append(pixel)
          if dist>20 and dist<25:
              far_pixel.append(pixel)
       #print centered_pixel
       bololist=['90.x','90.y','150.x', '150.y','220.x', '220.y']
       band=[90,90,150,150,220,220]
       #if you want to plot the first neighbor, change "centered_pixel" below to neighbor_pixel.
       #it can also be changed to be all_pixel, second_pixel, far_pixel to plot different pixel groups.
       #when doing this, it's better to also change data_dir above to be a different folder
       #this way you can store the fts spectra for the centered, neighbor, and second pixels separately 
       for i in range(len(centered_pixel)):
            for j in range(len(bololist)):
                   print file_name
                   bolo_data=glob(file_name)[0]
                   bolo_data=Dataset(bolo_data,mode='r')
                   if 'W'+str(waf)+'_'+str(centered_pixel[i])+'.'+bololist[j]+'_I' in bolo_data.variables.keys():
                    k=k+1
                    fts_result= fts(netcdf=file_name.split('.')[0],boloname='W'+str(waf)+'_'+str(centered_pixel[i])+'.'+bololist[j], center_pixel=center_pixel, band=band[j],speed=.5,chop=False,optics_correction=False,save=True,return_data=True,plot=False, num=k)
                    pl.close("all")
    return 0

#What are the centered pixels?
observed_pixels = ['176.150', '176.93','176.213', '176.97', '176.25', '176.204', '176.45', '176.208', '176.258', '176.141','187.49', '187.26', '187.97', '187.111', '187.214', '187.230', '187.45', '187.88', '187.107', '187.160', '187.249', '174.18', '174.54', '174.119', '174.205', '174.245', '174.193', '174.107','174.26','177.26','177.210', '177.158', '177.206', '177.247', '177.54','177.35', '177.25', '177.78' , '177.59', '177.107', '188.111', '188.63', '188.26', '188.214', '188.150', '188.66', '188.28', '188.258', '188.142', '188.54']


netcdfs=[]
#Find the netcdf files for these pixels
for pixel in observed_pixels:
  netcdf_file= sorted(glob('/poleanalysis/fts/summer2017-18/'+pixel+'/'+pixel+'*.nc'))
  if len(netcdf_file)>0:
    netcdfs.append(netcdf_file[-1])

#netcdf_dict is the input of netcdf_plot
netcdf_dict={}
for netcdf_file in netcdfs:
  netcdf_dict[netcdf_file]=  netcdf_file.split('/')[2]



netcdf_plot(netcdf_dict)   
#pl.show()


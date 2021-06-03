# -*- coding: utf-8 -*-
"""
Created on Sun May 30 14:04:02 2021

@author: danang eko nuryanto (R&D BMKG)
modified by donaldi permana (R&D BMKG)
"""

# from netCDF4 import Dataset,date2num,num2date
import netCDF4 as nc
import matplotlib.pyplot as pl
pl.rcParams.update({'figure.max_open_warning': 0})
# from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pytz
# import glob
import datetime
from PIL import Image
import sys, getopt
import os
# os.environ['PROJ_LIB'] = r'/usr/share/proj'
os.environ['PROJ_LIB'] = 'C:\\Users\\puslitbang\\anaconda3\\Library\\share' 
from mpl_toolkits.basemap import Basemap
# from wrf import getvar, interplevel, to_np, get_basemap, latlon_coords, extract_times
import warnings
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)

def main(argv):
    year = ''
    month = ''
    day = ''
    cycle = ''
    global degree, unit1, path, ds1, ds2, logofile, shpkabkotfname, shpprovfname, times1, str_times1, times2, str_times2
    degree='degC'
    unit1='per 100000 s'
    try:
       opts, args = getopt.getopt(argv,"hy:m:d:c:e",["year=","month=","day=","cycle="])
    except getopt.GetoptError:
       print ('test.py -y <year> -m <month> -d <day> -c <cycle>')
       sys.exit(2)
    for opt, arg in opts:
       if opt in ("-h","--help"):
          print ('test.py -y <year> -m <month> -d <day> -c <cycle>')
          sys.exit()
       elif opt in ("-y", "--year"):
          year = arg
       elif opt in ("-m", "--month"):
          month = arg
       elif opt in ("-d", "--day"):
          day = arg
       elif opt in ("-c", "--cycle"):
          cycle = arg
    print ("Initial time is ", year+month+day+cycle)
    
    # setting for windows
    dirout = 'E:\\DATA\\BMKG\\2021_BMKG\\05\\30_post_processing_InaNWP\\'
    logofile = dirout+'Logo-BMKG-new.png'
    shpkabkotfname = 'E:\\DATA\\BMKG\\2018_BMKG\\12\\backup\\pyscript\\shp\\Indo_Kab_Kot1'
    shpprovfname = 'E:\\DATA\\BMKG\\2018_BMKG\\12\\backup\\pyscript\\shp\\INDONESIA_PROP1'
    
    # setting for linux
    # dirout = '/home/wrfadmin/install-wrf/'
    # logofile = '/home/wrfadmin/install-wrf/WORK/shp/Logo-BMKG-new.png'
    # shpkabkotfname = '/home/wrfadmin/install-wrf/WORK/shp/Indo_Kab_Kot1'
    # shpprovfname = '/home/wrfadmin/install-wrf/WORK/shp/INDONESIA_PROP1'
    
    # Read WRF data via GDS server 
    inittime = year+month+day+cycle
    # inittime = '2021052912'
    url1 = "http://182.16.248.173:8080/dods/INA-NWP/"+inittime+"/"+inittime+"-d01-asim"
    url2 = "http://182.16.248.173:8080/dods/INA-NWP/"+inittime+"/"+inittime+"-d02-asim"
    
    print('Reading forecast data from '+ url1)
    print('Reading forecast data from '+ url2)
    
    ds1 = nc.Dataset(url1)
    # model = '9km InaNWPv0.8'
    
    ds2 = nc.Dataset(url2)
    # model = '3km InaNWPv1.0'
    
    # ds1 = nc.Dataset("http://182.16.248.173:8080/dods/INA-NWP/"+inittime+"/"+inittime+"-d01-no_asim")
    # ds2 = nc.Dataset("http://182.16.248.173:8080/dods/INA-NWP/"+inittime+"/"+inittime+"-d02-no_asim")
    
    # lons = ds2.variables['lon'][:]
    # lats = ds2.variables['lat'][:]
    # levels = ds2.variables['lev'][:]
    times1 = ds1.variables['time'][:]
    t_unit = ds1.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"
    try :
        t_cal = ds1.variables['time'].calendar
    except AttributeError : # Attribute doesn't exist
        t_cal = u"gregorian" # or standard
    times1 = nc.num2date(times1,units = t_unit,calendar = t_cal)

    str_times1 = []
    idx_times = range(len(times1))
    for i in idx_times:
        wkt = str(times1[i])
        local_wkt = datetime.datetime(int(wkt[:4]),int(wkt[5:7]),int(wkt[8:10]),int(wkt[11:13]),int(wkt[14:16]),int(wkt[17:19])).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('Asia/Jakarta')).strftime('%HWIB %d %b %Y')
        str_times1.append(local_wkt)
    
    times2 = ds2.variables['time'][:]
    t_unit = ds2.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"
    try :
        t_cal = ds2.variables['time'].calendar
    except AttributeError : # Attribute doesn't exist
        t_cal = u"gregorian" # or standard
    times2 = nc.num2date(times2,units = t_unit,calendar = t_cal)

    str_times2 = []
    idx_times = range(len(times2))
    for i in idx_times:
        wkt = str(times2[i])
        local_wkt = datetime.datetime(int(wkt[:4]),int(wkt[5:7]),int(wkt[8:10]),int(wkt[11:13]),int(wkt[14:16]),int(wkt[17:19])).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('Asia/Jakarta')).strftime('%HWIB %d %b %Y')
        str_times2.append(local_wkt)
        
    tstart = datetime.datetime.now()
    
    path=dirout

    # if cycle=='00':
    #     nn=24
    # elif cycle=='12':
    #     nn=36

    plot_ch24(ds1,1,str_times1[0],inittime,int(cycle))
    plot_ch24(ds2,2,str_times2[0],inittime,int(cycle))
    plot_ch24(ds2,3,str_times2[0],inittime,int(cycle))
        
    plot_all(ds1,1,inittime)
    plot_all(ds2,2,inittime)
    plot_all(ds2,3,inittime)

    print("Looping took:", datetime.datetime.now() - tstart)
  

def peta_indo(lons,lats,bound,data1,data2,data3,data4,parameters,path,fnam,step,init,pred,ver):     
# 	import matplotlib.pyplot as pl   
# 	from mpl_toolkits.basemap import Basemap
	if os.path.exists(path)==False:
		os.makedirs(path)
	fig = pl.figure(figsize=(10,6))
	m = Basemap(resolution='h',projection='merc',llcrnrlat=bound[0],urcrnrlat=bound[1],llcrnrlon=bound[2],urcrnrlon=bound[3])
	xii, yii = m(lons,lats)
	if parameters[-21:]=='Relative Humidity (%)':
		clevs=[50,55,60,65,70,75,80,85,90,95]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu', extend='both')
	elif parameters=='2-m Temperature ('+degree+')':
		clevs=[18,20,22,24,26,28,30,32,34,36]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='850 mb Temperature ('+degree+')':
		clevs=[14,15,16,17,18,19,20,21,22,23]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='700 mb Temperature ('+degree+')':
		clevs=[8,8.5,9,9.5,10,10.5,11,11.5,12,12.5]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='500 mb Temperature ('+degree+')':
		clevs=[-9,-8,-7,-6,-5,-4,-3,-2,-1,0]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='200 mb Temperature ('+degree+')':
		clevs=[-55,-54,-53,-52,-51,-50,-49,-48,-47,-46]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters[-19:]=='Wind & isotach(m/s)':
		# clevs=[8,10,12,14,16,18,20,22,24]
		clevs=[8,12,16,18,20,30,40,50,60]
		ccols=['white','palegreen','lime','limegreen','forestgreen','orange','orangered','red','firebrick','darkred']
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	elif parameters[-25:]=='divergence ('+unit1+')':
		# clevs=[0,0.5,1,1.5,2,2.5,3,4,5]
		clevs=[0,2,4,6,8,10,15,20,25]
		cs = m.contourf(xii,yii,data1,clevs,cmap='OrRd', extend='both')
	elif parameters[-26:]=='convergence ('+unit1+')':
		# clevs=[-5,-4,-3,-2.5,-2,-1.5,-1,-0.5,0]
		clevs=[-25,-20,-15,-10,-8,-6,-4,-2,0]
		cs = m.contourf(xii,yii,data1,clevs,cmap='Greens_r', extend='both')
	elif parameters[-33:]=='relative vorticity ('+unit1+')':
		clevs=[-20,-15,-10,-5,-2,2,5,10,15,20]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='10-m Wind(vector, m/s), 24hr Prec(shaded, mm), MSLP (contour, mb)':
		clevs=[0.5, 1, 5, 10, 15, 20, 40, 50, 65, 80, 100, 150]
		ccols=['#BEBEBE','#E8E8E7','#BDF2BA','#88F487','#68F422','#A4EE1B','#F2F220','#EFD216','#EBA91C','#ED8E1D','#EA661F','#EE251E','#E719B5'] 
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	else:
		clevs=[0.5, 1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50]
		ccols=['#BEBEBE','#E8E8E7','#BDF2BA','#88F487','#68F422','#A4EE1B','#F2F220','#EFD216','#EBA91C','#ED8E1D','#EA661F','#EE251E','#E719B5'] 
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	if np.logical_and(np.logical_and(data2!='none',data3=='none'),data4=='none'):
		ct=m.contour(xii,yii,data2,colors='slategrey',linewidths=.9)    
		pl.gca().clabel(ct, inline=1, fontsize=8,fmt='%1.0i')
	elif np.logical_and(np.logical_and(data2!='none',data3!='none'),data4=='none'):
		yy = np.arange(0, yii.shape[0], 10)
		xx = np.arange(0, xii.shape[1], 10)
		points = np.meshgrid(yy, xx)
		Q=m.quiver(xii[points], yii[points], data2[points],data3[points], color='grey')
		qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'10 m/s', labelpos='E', coordinates='figure')
	elif data4!='none':
		ct=m.contour(xii,yii,data4,colors='slategrey',linewidths=.9)    
		pl.gca().clabel(ct, inline=1, fontsize=8,fmt='%1.0i')
		yy = np.arange(0, yii.shape[0], 10)
		xx = np.arange(0, xii.shape[1], 10)
		points = np.meshgrid(yy, xx)
		Q=m.quiver(xii[points], yii[points], data2[points],data3[points], color='grey')
		qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'10 m/s', labelpos='E', coordinates='figure')
		# qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'0 \frac{m}{s}$', labelpos='E', coordinates='figure')
	m.readshapefile(shpprovfname,'INDONESIA_PROP1')
	plane = np.array(Image.open(logofile))
	m.drawcoastlines(0.5)
	m.drawcountries()
	m.drawparallels(np.arange(-15,35,5), labels=[1,0,0,0], dashes=[6, 900], fontsize=12)
	m.drawmeridians(np.arange(90,150,5), labels=[0,0,0,1], dashes=[6, 900], fontsize=12)
	ax = pl.axes([0.06,0.2, 0.06, 0.06], frameon=True)  
	ax.imshow(plane)
	ax.axis('off') 
	cbar_ax = fig.add_axes([0.07,0.048,0.9,0.03])
	cbar = pl.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=clevs)
	strings = ["%.1f" % number for number in clevs]
# 	cbar.ax.set_yticklabels(strings)
	pl.gcf().text(0.06, 0.96, parameters, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.06, 0.91, 'Forecast: '+pred+' (T+'+step+')', rotation='horizontal',fontsize=15)
	pl.gcf().text(0.80, 0.86, ver, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.06, 0.86, '$\it{Initial :}$ '+init, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.58, 0.11, '@ $\it{Center}$ $\it{for}$ $\it{Research}$ $\it{and}$ $\it{Development}$ $\it{BMKG}$', rotation='horizontal',fontsize=12)
	fig.tight_layout()
	fig.subplots_adjust(hspace=0,wspace=0,left=0.06,right=0.98,bottom=0.06,top=0.96)
	pl.savefig(path+fnam, format='png', dpi=90, bbox_inches='tight')
	pl.cla()

def peta_jawa(lons,lats,bound,data1,data2,data3,data4,parameters,path,fnam,step,init,pred,ver):
# 	import matplotlib.pyplot as pl           
# 	from mpl_toolkits.basemap import Basemap
	if os.path.exists(path)==False:
		os.makedirs(path)
	fig = pl.figure(figsize=(10,5.4))
	m = Basemap(resolution='h',projection='merc',llcrnrlat=bound[0],urcrnrlat=bound[1],llcrnrlon=bound[2],urcrnrlon=bound[3])
	xii, yii = m(lons,lats)
	if parameters[-21:]=='Relative Humidity (%)':
		clevs=[50,55,60,65,70,75,80,85,90,95]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu', extend='both')
	elif parameters=='2-m Temperature ('+degree+')':
		clevs=[18,20,22,24,26,28,30,32,34,36]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='850 mb Temperature ('+degree+')':
		clevs=[14,15,16,17,18,19,20,21,22,23]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='700 mb Temperature ('+degree+')':
		clevs=[8,8.5,9,9.5,10,10.5,11,11.5,12,12.5]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='500 mb Temperature ('+degree+')':
		clevs=[-9,-8,-7,-6,-5,-4,-3,-2,-1,0]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='200 mb Temperature ('+degree+')':
		clevs=[-55,-54,-53,-52,-51,-50,-49,-48,-47,-46]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters[-19:]=='Wind & isotach(m/s)':
		# clevs=[8,10,12,14,16,18,20,22,24]
		clevs=[8,12,16,18,20,30,40,50,60]
		ccols=['white','palegreen','lime','limegreen','forestgreen','orange','orangered','red','firebrick','darkred']
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	elif parameters[-25:]=='divergence ('+unit1+')':
		# clevs=[0,0.5,1,1.5,2,2.5,3,4,5]
		clevs=[0,2,4,6,8,10,15,20,25]
		cs = m.contourf(xii,yii,data1,clevs,cmap='OrRd', extend='both')
	elif parameters[-26:]=='convergence ('+unit1+')':
		# clevs=[-5,-4,-3,-2.5,-2,-1.5,-1,-0.5,0]
		clevs=[-25,-20,-15,-10,-8,-6,-4,-2,0]
		cs = m.contourf(xii,yii,data1,clevs,cmap='Greens_r', extend='both')
	elif parameters[-33:]=='relative vorticity ('+unit1+')':
		clevs=[-20,-15,-10,-5,-2,2,5,10,15,20]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='10-m Wind(vector, m/s), 24hr Prec(shaded, mm), MSLP (contour, mb)':
		clevs=[0.5, 1, 5, 10, 15, 20, 40, 50, 65, 80, 100, 150]
		ccols=['#BEBEBE','#E8E8E7','#BDF2BA','#88F487','#68F422','#A4EE1B','#F2F220','#EFD216','#EBA91C','#ED8E1D','#EA661F','#EE251E','#E719B5'] 
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	else:
		clevs=[0.5, 1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50]
		ccols=['#BEBEBE','#E8E8E7','#BDF2BA','#88F487','#68F422','#A4EE1B','#F2F220','#EFD216','#EBA91C','#ED8E1D','#EA661F','#EE251E','#E719B5'] 
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	if np.logical_and(np.logical_and(data2!='none',data3=='none'),data4=='none'):
		ct=m.contour(xii,yii,data2,colors='slategrey',linewidths=.9)    
		pl.gca().clabel(ct, inline=1, fontsize=8,fmt='%1.0i')
	elif np.logical_and(np.logical_and(data2!='none',data3!='none'),data4=='none'):
		yy = np.arange(0, yii.shape[0], 20)
		xx = np.arange(0, xii.shape[1], 20)
		points = np.meshgrid(yy, xx)
		Q=m.quiver(xii[points], yii[points], data2[points],data3[points], color='grey')
		qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'10 m/s', labelpos='E', coordinates='figure')
	elif data4!='none':
		ct=m.contour(xii,yii,data4,colors='slategrey',linewidths=.9)    
		pl.gca().clabel(ct, inline=1, fontsize=8,fmt='%1.0i')
		yy = np.arange(0, yii.shape[0], 20)
		xx = np.arange(0, xii.shape[1], 20)
		points = np.meshgrid(yy, xx)
		Q=m.quiver(xii[points], yii[points], data2[points],data3[points], color='grey')
		qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'10 m/s', labelpos='E', coordinates='figure')
		# qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'0 \frac{m}{s}$', labelpos='E', coordinates='figure')
	m.readshapefile(shpprovfname,'INDONESIA_PROP1')
	plane = np.array(Image.open(logofile))
	# m.drawcoastlines(0.5)
	# m.drawcountries()
	m.drawparallels(np.arange(-15,35,1), labels=[1,0,0,0], dashes=[6, 900], fontsize=12)
	m.drawmeridians(np.arange(90,150,1), labels=[0,0,0,1], dashes=[6, 900], fontsize=12)
	ax = pl.axes([0.06,0.2, 0.06, 0.06], frameon=True)  
	ax.imshow(plane)
	ax.axis('off') 
	cbar_ax = fig.add_axes([0.07,0.045,0.9,0.03])
	cbar =pl.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=clevs)
# 	cbar.ax.set_yticklabels(clevs) 
	pl.gcf().text(0.06, 0.96, parameters, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.06, 0.91, 'Forecast: '+pred+' (T+'+step+')', rotation='horizontal',fontsize=15)
	pl.gcf().text(0.80, 0.86, ver, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.06, 0.86, '$\it{Initial :}$ '+init, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.58, 0.1, '@ $\it{Center}$ $\it{for}$ $\it{Research}$ $\it{and}$ $\it{Development}$ $\it{BMKG}$', rotation='horizontal',fontsize=12)
	fig.tight_layout()
	fig.subplots_adjust(hspace=0,wspace=0,left=0.06,right=0.98,bottom=0.06,top=0.96)
	pl.savefig(path+fnam, format='png', dpi=90, bbox_inches='tight')
	pl.cla()

def peta_jkt(lons,lats,bound,data1,data2,data3,data4,parameters,path,fnam,step,init,pred,ver):      
# 	import matplotlib.pyplot as pl     
# 	from mpl_toolkits.basemap import Basemap
	if os.path.exists(path)==False:
		os.makedirs(path)
	fig = pl.figure(figsize=(8.9,10))
	m = Basemap(resolution='h',projection='merc',llcrnrlat=bound[0],urcrnrlat=bound[1],llcrnrlon=bound[2],urcrnrlon=bound[3])
	xii, yii = m(lons,lats)
	if parameters[-21:]=='Relative Humidity (%)':
		clevs=[50,55,60,65,70,75,80,85,90,95]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu', extend='both')
	elif parameters=='2-m Temperature ('+degree+')':
		clevs=[18,20,22,24,26,28,30,32,34,36]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='850 mb Temperature ('+degree+')':
		clevs=[14,15,16,17,18,19,20,21,22,23]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='700 mb Temperature ('+degree+')':
		clevs=[8,8.5,9,9.5,10,10.5,11,11.5,12,12.5]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='500 mb Temperature ('+degree+')':
		clevs=[-9,-8,-7,-6,-5,-4,-3,-2,-1,0]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='200 mb Temperature ('+degree+')':
		clevs=[-55,-54,-53,-52,-51,-50,-49,-48,-47,-46]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters[-19:]=='Wind & isotach(m/s)':
		# clevs=[8,10,12,14,16,18,20,22,24]
		clevs=[8,12,16,18,20,30,40,50,60]
		ccols=['white','palegreen','lime','limegreen','forestgreen','orange','orangered','red','firebrick','darkred']
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	elif parameters[-25:]=='divergence ('+unit1+')':
		# clevs=[0,0.5,1,1.5,2,2.5,3,4,5]
		clevs=[0,2,4,6,8,10,15,20,25]
		cs = m.contourf(xii,yii,data1,clevs,cmap='OrRd', extend='both')
	elif parameters[-26:]=='convergence ('+unit1+')':
		# clevs=[-5,-4,-3,-2.5,-2,-1.5,-1,-0.5,0]
		clevs=[-25,-20,-15,-10,-8,-6,-4,-2,0]
		cs = m.contourf(xii,yii,data1,clevs,cmap='Greens_r', extend='both')
	elif parameters[-33:]=='relative vorticity ('+unit1+')':
		clevs=[-20,-15,-10,-5,-2,2,5,10,15,20]
		cs = m.contourf(xii,yii,data1,clevs,cmap='RdBu_r', extend='both')
	elif parameters=='10-m Wind(vector, m/s), 24hr Prec(shaded, mm), MSLP (contour, mb)':
		clevs=[0.5, 1, 5, 10, 15, 20, 40, 50, 65, 80, 100, 150]
		ccols=['#BEBEBE','#E8E8E7','#BDF2BA','#88F487','#68F422','#A4EE1B','#F2F220','#EFD216','#EBA91C','#ED8E1D','#EA661F','#EE251E','#E719B5'] 
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	else:
		clevs=[0.5, 1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50]
		ccols=['#BEBEBE','#E8E8E7','#BDF2BA','#88F487','#68F422','#A4EE1B','#F2F220','#EFD216','#EBA91C','#ED8E1D','#EA661F','#EE251E','#E719B5'] 
		cs = m.contourf(xii,yii,data1,clevs,colors=ccols, extend='both')
	if np.logical_and(np.logical_and(data2!='none',data3=='none'),data4=='none'):
		ct=m.contour(xii,yii,data2,colors='slategrey',linewidths=.9)    
		pl.gca().clabel(ct, inline=1, fontsize=8,fmt='%1.0i')
	elif np.logical_and(np.logical_and(data2!='none',data3!='none'),data4=='none'):
		yy = np.arange(0, yii.shape[0], 5)
		xx = np.arange(0, xii.shape[1], 5)
		points = np.meshgrid(yy, xx)
		Q=m.quiver(xii[points], yii[points], data2[points],data3[points], color='grey')
		qk = pl.quiverkey(Q, 0.92, 0.92, 10, r'10 m/s', labelpos='E', coordinates='figure')
	elif data4!='none':
		ct=m.contour(xii,yii,data4,colors='slategrey',linewidths=.9)    
		pl.gca().clabel(ct, inline=1, fontsize=8,fmt='%1.0i')
		yy = np.arange(0, yii.shape[0], 5)
		xx = np.arange(0, xii.shape[1], 5)
		points = np.meshgrid(yy, xx)
		Q=m.quiver(xii[points], yii[points], data2[points],data3[points], color='grey')
		qk = pl.quiverkey(Q, 0.92, 0.92, 10, r'10 m/s', labelpos='E', coordinates='figure')
		# qk = pl.quiverkey(Q, 0.1, 0.1, 10, r'0 \frac{m}{s}$', labelpos='E', coordinates='figure')
	m.readshapefile(shpkabkotfname,'Indo_Kab_Kot1')
	plane = np.array(Image.open(logofile))
	# m.drawcoastlines(0.5)
	# m.drawcountries()
	m.drawparallels(np.arange(-15,35,.2), labels=[1,0,0,0], dashes=[6, 900], fontsize=12)
	m.drawmeridians(np.arange(90,150,.2), labels=[0,0,0,1], dashes=[6, 900], fontsize=12)
	ax = pl.axes([0.10,0.12, 0.06, 0.06], frameon=True)  
	ax.imshow(plane)
	ax.axis('off') 
	cbar_ax = fig.add_axes([0.93,0.1,0.03,0.79])
	cbar = pl.colorbar(cs, cax=cbar_ax, orientation='vertical', ticks=clevs)
# 	cbar.ax.set_yticklabels(clevs) 
	pl.gcf().text(0.07, 0.97, parameters, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.07, 0.94, 'Forecast: '+pred+' (T+'+step+')', rotation='horizontal',fontsize=15)
	pl.gcf().text(0.65, 0.91, ver, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.07, 0.91, '$\it{Initial :}$ '+init, rotation='horizontal',fontsize=15)
	pl.gcf().text(0.44, 0.05, '@ $\it{Center}$ $\it{for}$ $\it{Research}$ $\it{and}$ $\it{Development}$ $\it{BMKG}$', rotation='horizontal',fontsize=12)
	fig.tight_layout()
	fig.subplots_adjust(hspace=0,wspace=0,left=0.08,right=0.9,bottom=0.1,top=0.9)
	pl.savefig(path+fnam, format='png', dpi=90, bbox_inches='tight')
	pl.cla()
    
def initpred(wk,npred,step,pil):
    wkt=str(wk[0])[:13]
    xxx=wkt[:4]+wkt[5:7]+wkt[8:10]+wkt[11:]
    init=datetime.datetime(int(wkt[:4]),int(wkt[5:7]),int(wkt[8:10]),int(wkt[11:])).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('Asia/Jakarta')).strftime('%HWIB %d %b %Y')
    wkt=str(wk[npred-8])[:13]
    init1=datetime.datetime(int(wkt[:4]),int(wkt[5:7]),int(wkt[8:10]),int(wkt[11:])).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('Asia/Jakarta')).strftime('%HWIB %d %b %Y')
    wkt=str(wk[npred])[:13]
    pred=datetime.datetime(int(wkt[:4]),int(wkt[5:7]),int(wkt[8:10]),int(wkt[11:])).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('Asia/Jakarta')).strftime('%HWIB %d %b %Y')
    if pil==1:
        pred=init1+' - '+pred
    fnam=xxx+'_'+wkt[:4]+wkt[5:7]+wkt[8:10]+wkt[11:]+'('+step+').png'
    return(init,pred,fnam)

def plot_ch24(ds,pil,init,xxx,nn):
    if pil==1:
        ijj='indo9'
        yyy='Indonesia(9km)'
        ver='9km InaNWPv0.8'
        bound=[-13.,7.,94.,142.5]
        times = times1
        str_times = str_times1
    elif pil==2:
        ijj='jawa3'
        yyy='Jawa(3km)'
        ver='3km InaNWPv1.0'
        bound=[-9.5,-4.,102.,116.]
        times = times2
        str_times = str_times2
    elif pil==3:
        ijj='jkt3'
        yyy='Jabodetabek(3km)'
        ver='3km InaNWPv1.0'
        bound=[-7.,-5.7,106.2,107.4]
        times = times2
        str_times = str_times2
    
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    lons, lats = np.meshgrid(lon, lat)
    levels = ds.variables['lev'][:]
        
    indices = [i for i, s in enumerate(str_times) if '07WIB' in s]
    
    # Extract variables from the data
    # tempisobar = np.transpose(fn.variables['tk'][:,0,:,:]-273.15)
    # rh = np.transpose(fn.variables['rh'][:,0,:,:])
    # ugrid = np.transpose(fn.variables['umet'][:,0,:,:])
    # vgrid = np.transpose(fn.variables['vmet'][:,0,:,:])
    
    # tempht = fn.variables['t2'][:,0,:,:]-273.15
    # dewtempht = fn.variables['td2'][:,0,:,:]
    # rhht = fn.variables['rh2'][:,0,:,:]
    u10m = ds.variables['u10m'][:,0,:,:]
    v10m = ds.variables['v10m'][:,0,:,:]
    slp = ds.variables['slp'][:,0,:,:]
    rainc = ds.variables['rainc'][:,0,:,:]
    rainsh = ds.variables['rainsh'][:,0,:,:]
    rainnc = ds.variables['rainnc'][:,0,:,:]
    rain03 = rainc[1:,:,:] - rainc[:-1,:,:] + rainsh[1:,:,:] - rainsh[:-1,:,:] + rainnc[1:,:,:] - rainnc[:-1,:,:]
    # rain03 = np.ma.append(0, rain03)
    # clflo = fn.variables['clflo'][:,0,:,:]
    # clfmi = fn.variables['clfmi'][:,0,:,:]
    # clfhi = fn.variables['clfhi'][:,0,:,:]
     
    for npred in range(1,len(indices)):
        init1=str_times[indices[npred-1]]
        ch = np.nansum(rain03[indices[npred-1]:indices[npred],:,:],axis=0)
        mslp = slp[indices[npred]]
        u10 = u10m[indices[npred]]
        v10 = v10m[indices[npred]]
        pred=str_times[indices[npred]]
        pred=init1+' - '+pred
        if nn > 0:
            step=str((24-nn)+npred*24)
        else:
            step=str(nn+npred*24)
        if int(step)<10:
            fnam=xxx+'(0'+step+').png'
        else:
            fnam=xxx+'('+step+').png'
        fname='24h-prec_'+fnam
        parameters='10-m Wind(vector, m/s), 24hr Prec(shaded, mm), MSLP (contour, mb)'
        pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/Surface/24hr-Prec-mslp-wind/'
        print('Plotting '+fname + ' for domain '+ yyy + ' ...')
        if pil==1:
            peta_indo(lons,lats,bound,ch,u10.data,v10.data,mslp.data,parameters,pathn,fname,step,init,pred,ver)
        elif pil==2:
            peta_jawa(lons,lats,bound,ch,u10.data,v10.data,mslp.data,parameters,pathn,fname,step,init,pred,ver)
        elif pil==3:
            peta_jkt(lons,lats,bound,ch,u10.data,v10.data,mslp.data,parameters,pathn,fname,step,init,pred,ver)
            
def plot_all(ds,pil,xxx):
    if pil==1:
        ijj='indo9'
        yyy='Indonesia(9km)'
        ver='9km InaNWPv0.8'
        bound=[-13.,7.,94.,142.5]
        times = times1
        str_times = str_times1
    elif pil==2:
        ijj='jawa3'
        yyy='Jawa(3km)'
        ver='3km InaNWPv1.0'
        bound=[-9.5,-4.,102.,116.]
        times = times2
        str_times = str_times2
    elif pil==3:
        ijj='jkt3'
        yyy='Jabodetabek(3km)'
        ver='3km InaNWPv1.0'
        bound=[-7.,-5.7,106.2,107.4]
        times = times2
        str_times = str_times2
        
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    lons, lats = np.meshgrid(lon, lat)
    levels = ds.variables['lev'][:]
    init = str_times[0]
        
    # indices = [i for i, s in enumerate(str_times) if '07WIB' in s]
    
    # Extract variables from the data
    # tempisobar = np.transpose(fn.variables['tk'][:,0,:,:]-273.15)
    # rh = np.transpose(fn.variables['rh'][:,0,:,:])
    # ugrid = np.transpose(fn.variables['umet'][:,0,:,:])
    # vgrid = np.transpose(fn.variables['vmet'][:,0,:,:])
    
    t2 = ds.variables['t2'][:,0,:,:]-273.15
    # dewtempht = fn.variables['td2'][:,0,:,:]
    rh2 = ds.variables['rh2'][:,0,:,:]
    u10m = ds.variables['u10m'][:,0,:,:]
    v10m = ds.variables['v10m'][:,0,:,:]
    slp = ds.variables['slp'][:,0,:,:]
    rainc = ds.variables['rainc'][:,0,:,:]
    rainsh = ds.variables['rainsh'][:,0,:,:]
    rainnc = ds.variables['rainnc'][:,0,:,:]
    rain03 = rainc[1:,:,:] - rainc[:-1,:,:] + rainsh[1:,:,:] - rainsh[:-1,:,:] + rainnc[1:,:,:] - rainnc[:-1,:,:]
    # rain03 = np.ma.append(0, rain03)
    # clflo = fn.variables['clflo'][:,0,:,:]
    # clfmi = fn.variables['clfmi'][:,0,:,:]
    # clfhi = fn.variables['clfhi'][:,0,:,:]

    lv=[1000,850,200]
    for qq in range(len(lv)):
        idx = np.where(levels == lv[qq])[0][0]
        u_ = ds.variables['u'][:,idx,:,:]
        v_ = ds.variables['v'][:,idx,:,:]
        w_ = ds.variables['wspd'][:,idx,:,:]
        r_ = ds.variables['rh'][:,idx,:,:]
        t_ = ds.variables['tc'][:,idx,:,:]
           
        for npred in range(1,len(str_times)):
            # init1=str_times[npred-1]
            ch = np.nansum(rain03[npred-1:npred,:,:],axis=0)
            mslp = slp[npred]
            u10 = u10m[npred]
            v10 = v10m[npred]
            pred=str_times[npred]
            # pred=init1+' - '+pred
            
            step=str(npred*3)
            if int(step)<10:
                fnam=xxx+'(0'+step+').png'
            else:
                fnam=xxx+'('+step+').png'
      
            u = u_[npred,:,:]
            v = v_[npred,:,:]
            w = w_[npred,:,:]
            r = r_[npred,:,:]
            t = t_[npred,:,:]
                
            if lv[qq]==1000:
                t = t2[npred,:,:]
                r = rh2[npred,:,:]
                lev='Surface'
                levv='Surface'
                levt='2-m'
                fname=ijj+'km_wind10m-prec_'+fnam
                print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
                pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/wind10m-Precipitation/'
                parameters='10-m Wind(vector, m/s), Prec(shaded, mm)'
                if pil==1:
                    peta_indo(lons,lats,bound,ch,u10.data,v10.data,'none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==2:
                    peta_jawa(lons,lats,bound,ch,u10.data,v10.data,'none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==3:
                    peta_jkt(lons,lats,bound,ch,u10.data,v10.data,'none',parameters,pathn,fname,step,init,pred,ver)
                
                fname=ijj+'km_mslp-prec_'+fnam
                print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
                pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/mslp-Precipitation/'
                parameters='MSLP(contour, mb), Prec(shaded, mm)'
                if pil==1:
                    peta_indo(lons,lats,bound,ch,mslp.data,'none','none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==2:
                    peta_jawa(lons,lats,bound,ch,mslp.data,'none','none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==3:
                    peta_jkt(lons,lats,bound,ch,mslp.data,'none','none',parameters,pathn,fname,step,init,pred,ver)
            else:                
                lev=str(lv[qq])
                levv=str(lv[qq])+' mb'
                levt=str(lv[qq])+' mb'
            
            fname=ijj+'km_streamline_'+fnam
            print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
            pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/streamline/'
            parameters=levv+' Wind & isotach(m/s)'
            if pil==1:
                peta_indo(lons,lats,bound,w,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
            elif pil==2:
                peta_jawa(lons,lats,bound,w,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
            elif pil==3:
                peta_jkt(lons,lats,bound,w,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
            
            fname=ijj+'km_rh_'+fnam
            print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
            pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/rh/'
            parameters=levt+' Relative Humidity (%)'
            if pil==1:
                peta_indo(lons,lats,bound,r,'none','none','none',parameters,pathn,fname,step,init,pred,ver)
            elif pil==2:
                peta_jawa(lons,lats,bound,r,'none','none','none',parameters,pathn,fname,step,init,pred,ver)
            elif pil==3:
                peta_jkt(lons,lats,bound,r,'none','none','none',parameters,pathn,fname,step,init,pred,ver)
            
            fname=ijj+'km_temp_'+fnam
            print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
            pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/temp/'
            parameters=levt+' Temperature ('+degree+')'
            if pil==1:
                peta_indo(lons,lats,bound,t,'none','none','none',parameters,pathn,fname,step,init,pred,ver)
            elif pil==2:
                peta_jawa(lons,lats,bound,t,'none','none','none',parameters,pathn,fname,step,init,pred,ver)
            elif pil==3:
                peta_jkt(lons,lats,bound,t,'none','none','none',parameters,pathn,fname,step,init,pred,ver)
            
            dX = np.array(np.gradient(lons))*111139. 
            dY = np.array(np.gradient(lats))*111139.
            dV = np.array(np.gradient(v.data))
            Vgradient =  dV[1]/dX[1];Vgrad = dV[0]/dY[0]
            dU = np.array(np.gradient(u.data))
            Ugradient =  dU[0]/dY[0];Ugrad = dU[1]/dX[1]
            div = (Ugradient + Vgradient)*100000
            
            if np.logical_or(lv[qq]==1000,lv[qq]==200):
                fname=ijj+'km_divergence_'+fnam
                print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
                pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/divergence/'
                parameters=levv+' divergence ('+unit1+')'
                if pil==1:
                    peta_indo(lons,lats,bound,div,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==2:
                    peta_jawa(lons,lats,bound,div,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==3:
                    peta_jkt(lons,lats,bound,div,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
            
            if np.logical_or(lv[qq]==1000,lv[qq]==850):
                fname=ijj+'km_convergence_'+fnam
                print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
                pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/convergence/'
                parameters=levv+' convergence ('+unit1+')'
                if pil==1:
                    peta_indo(lons,lats,bound,div,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==2:
                    peta_jawa(lons,lats,bound,div,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
                elif pil==3:
                    peta_jkt(lons,lats,bound,div,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)
            
            vort = (Vgradient - Ugradient)*100000
            fname=ijj+'km_vorticity_'+fnam
            print('Plotting '+fname + ' for domain '+ yyy + ' at ' + levv + ' ...')
            pathn=path+'PRODUCT/GFS0.25deg+Ground+Radar/'+xxx+'/'+yyy+'/'+lev+'/vorticity/'
            parameters=levv+' relative vorticity ('+unit1+')'
            if pil==1:
                peta_indo(lons,lats,bound,vort,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)   
            elif pil==2:     
                peta_jawa(lons,lats,bound,vort,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)   
            elif pil==3:     
                peta_jkt(lons,lats,bound,vort,u.data,v.data,'none',parameters,pathn,fname,step,init,pred,ver)     
        
    
if __name__ == "__main__":
   main(sys.argv[1:])

"""
Modified on Sun May 16 21:19:06 2021

@author: donaldi permana - puslitbang
"""

"""
NCL_meteo_1.py
===============

This script illustrates the following concepts:
   - Drawing a meteogram
   - Creating a color map using hex values
   - Reversing the Y axis
   - Explicitly setting tickmarks and labels on the bottom X axis
   - Increasing the thickness of contour lines
   - Drawing wind barbs
   - Drawing a bar chart
   - Changing the width and height of a plot
   - Overlaying wind barbs and line contours on filled contours
   - Changing the position of individual plots on a page

See following URLs to see the reproduced NCL plot & script:
    - Original NCL script: https://www.ncl.ucar.edu/Applications/Scripts/meteo_1.ncl
    - Original NCL plot: https://www.ncl.ucar.edu/Applications/Images/meteo_1_lg.png
"""

import cartopy.crs as ccrs
import numpy as np
import netCDF4 as nc
import datetime as dt
import pandas as pd
import sys, getopt
import pytz
import datetime
import os
# os.environ['PROJ_LIB'] = r'/usr/share/proj'
os.environ['PROJ_LIB'] = 'C:\\Users\\puslitbang\\anaconda3\\Library\\share' 
import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

###############################################################################
# Import necessary packages
# import xarray as xr
from geocat.viz import util as gvutil
from matplotlib import pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from matplotlib.colors import BoundaryNorm, ListedColormap

###############################################################################

def main(argv):    
    global dirout, ds1, ds2, logofile, shpkabkotfname, shpprovfname, inittime
    global lons1, lats1, levels1, times1, str_times1, lons2, lats2, levels2, times2, str_times2
    arg_lon = 'nan'
    arg_lat = 'nan'
    arg_strpos = 'nan'
    try:
       opts, args = getopt.getopt(argv,"hy:m:d:c:i:j:s:r",["year=","month=","day=","cycle=","strlon=","strlat=","strpos="])
    except getopt.GetoptError:
       print ('test.py -y <year> -m <month> -d <day> -c <cycle> -i <longitude> -j <latitude> -s <strpos>')
       sys.exit(2)
    for opt, arg in opts:
       # print(opt)
       if opt in ("-h","--help"):
          print ('test.py -y <year> -m <month> -d <day> -c <cycle> -lon <longitude> -lat <latitude> -strpos <strpos>')
          sys.exit()
       elif opt in ("-y", "--year"):
          year = arg
       elif opt in ("-m", "--month"):
          month = arg
       elif opt in ("-d", "--day"):
          day = arg
       elif opt in ("-c", "--cycle"):
          cycle = arg
       elif opt in ("-i", "--longitude"):
          arg_lon = arg
       elif opt in ("-j", "--latitude"):
          arg_lat = arg
       elif opt in ("-s", "--strpos"):
          arg_strpos = arg
    print ("Initial time is ", year+month+day+cycle)
    
    # print (arg_lon + arg_lat + arg_strpos)
    
    # setting for windows
    dirout = 'E:\\DATA\\BMKG\\2021_BMKG\\05\\16_NCL_meteogram\\'
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
    
    lons1 = ds1.variables['lon'][:]
    lats1 = ds1.variables['lat'][:]
    levels1 = ds1.variables['lev'][:]
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
    
    lons2 = ds2.variables['lon'][:]
    lats2 = ds2.variables['lat'][:]
    levels2 = ds2.variables['lev'][:]
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
    
    if (arg_strpos != 'nan' and arg_lon != 'nan' and arg_lat != 'nan'):

        check_plot_meteogram(arg_strpos, arg_lon, arg_lat) 
    else:        
        
        # Read coordinate for plotting meteogram
        ifile = dirout + 'titik_meteogram.xlsx'
        
        print('Reading point coordinates from '+ ifile)
        
        metadata = pd.read_excel(ifile,sheet_name=0)
        
        ##############################
        
        # for idx in range(0, 10):#len(metadata)):
        for idx in range(len(metadata)):
        
            str_pos = str(metadata[metadata.columns[0]][idx])
            str_lon = str(metadata[metadata.columns[1]][idx])  
            str_lat = str(metadata[metadata.columns[2]][idx])
        
            # str_pos = 'Kebayoran Lama'
            # str_lon = '106.7453617'
            # str_lat = '-6.2493198'
            
            # str_pos = 'Bekasi Timur'
            # str_lon = '106.6559969'
            # str_lat = '-6.3307958'
    
            if not(str_pos != 'nan' and str_lon != 'nan' and str_lat != 'nan'):
                print(str(idx+1) + ' Plotting meteogram for '+ str_pos + ' is incomplete due to unsufficient information')
                continue
            
            check_plot_meteogram(str_pos,str_lon,str_lat)
            
    print("Looping took:", datetime.datetime.now() - tstart)

def check_plot_meteogram(str_pos,str_lon,str_lat):
    
    # set lon dan lat
    lat = round(float(str_lat), 3)
    lon = round(float(str_lon), 3)

    if (lon>=np.nanmin(lons1) and lon<=np.nanmax(lons1) and lat>=np.nanmin(lats1) and lat<=np.nanmax(lats1)):
        if (lon>=np.nanmin(lons2) and lon<=np.nanmax(lons2) and lat>=np.nanmin(lats2) and lat<=np.nanmax(lats2)):
            plot_meteogram(ds2, 2, str_pos, str_lon, str_lat)
        else:
            plot_meteogram(ds1, 1, str_pos, str_lon, str_lat)
    else:
        print('the coordinate position is out of the domain area')
                    
def plot_meteogram(ds,id,str_pos,str_lon,str_lat):

    print('Plotting meteogram for '+ str_pos + ' (lon = ' + str_lon + ', lat = ' + str_lat + ')' )
    
    lat = round(float(str_lat), 3)
    lon = round(float(str_lon), 3)
    
    if (id==1):
        model = '9km InaNWPv0.8'
        lons = lons1
        lats = lats1
        levels = levels1
        str_times = str_times1
        times = times1
    elif(id==2):
        model = '3km InaNWPv1.0'
        lons = lons2
        lats = lats2
        levels = levels2
        str_times = str_times2
        times = times2
        
    # set maintitle
    area_name = str_pos
    str_title = area_name + '\n 0~3day 3-hourly Forecast Meteogram for (' + str(lon) + '; ' + str(lat) + ') by ' + model + '\n@ Center for Research and Development BMKG' 
    idxlon = np.where(np.abs(lons-lon) == np.nanmin(np.abs(lons-lon)))[0][0]
    idxlat = np.where(np.abs(lats-lat) == np.nanmin(np.abs(lats-lat)))[0][0]

    # Extract variables from the data
    tempisobar = np.transpose(ds.variables['tk'][:,:,idxlat,idxlon]-273.15)
    rh = np.transpose(ds.variables['rh'][:,:,idxlat,idxlon])
    ugrid = np.transpose(ds.variables['umet'][:,:,idxlat,idxlon])
    vgrid = np.transpose(ds.variables['vmet'][:,:,idxlat,idxlon])
    
    tempht = ds.variables['t2'][:,0,idxlat,idxlon]-273.15
    dewtempht = ds.variables['td2'][:,0,idxlat,idxlon]
    rhht = ds.variables['rh2'][:,0,idxlat,idxlon]
    ws10 = ds.variables['ws10'][:,0,idxlat,idxlon]
    u10m = ds.variables['u10m'][:,0,idxlat,idxlon]
    v10m = ds.variables['v10m'][:,0,idxlat,idxlon]
    rainc = ds.variables['rainc'][:,0,idxlat,idxlon]
    rainsh = ds.variables['rainsh'][:,0,idxlat,idxlon]
    rainnc = ds.variables['rainnc'][:,0,idxlat,idxlon]
    rain03 = rainc[1:] - rainc[:-1] + rainsh[1:] - rainsh[:-1] + rainnc[1:] - rainnc[:-1]
    rain03 = np.ma.append(0, rain03)
    clflo = ds.variables['clflo'][:,0,idxlat,idxlon]
    clfmi = ds.variables['clfmi'][:,0,idxlat,idxlon]
    clfhi = ds.variables['clfhi'][:,0,idxlat,idxlon]
    
    idx_times = range(len(times))
    taus = idx_times
    
    ###############################################################################
    # Plot:
    
    # Generate figure (set its size (width, height) in inches)
    # fig = plt.figure(figsize=(15, 20))
    fig = plt.figure(figsize=(8, 14))
    spec = fig.add_gridspec(ncols=1, nrows=6, height_ratios=[4, 1, 1, 1, 1, 1], hspace=0.1)
    
    # Create axis for contour/wind barb plot
    ax1 = fig.add_subplot(spec[0, 0], projection=ccrs.PlateCarree())
    
    # Add coastlines to first axis
    ax1.coastlines(linewidths=0.5, alpha=-1)
    
    # Set aspect ratio of the first axis
    ax1.set_aspect(0.95)
    
    # Create a color map with a combination of matplotlib colors and hex values
    colors = ListedColormap(
        np.array([
            'white', 'white', 'white', 'mintcream', "#DAF6D3",
            "#DAF6D3", "#B2FAB9", 'springgreen', 'lime', "#54A63F"
        ]))
    # colors1 = ListedColormap(
    #     np.array([
    #         'white', 'white', 'white', 'white', 'white', 'mintcream', "#DAF6D3",
    #         "#B2FAB9", "#B2FAB9", 'springgreen', 'lime', "#54A63F"
    #     ]))
    contour_levels = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # contour_levels1 = [-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    normalized_levels = BoundaryNorm(boundaries=contour_levels, ncolors=10)
    # normalized_levels1 = BoundaryNorm(boundaries=contour_levels1, ncolors=12)
    
    level_to_plot = 100
    idx_lvl = np.where(levels==level_to_plot)[0][0]
    
    # Plot filled contours for the rh variable
    contour1 = ax1.contourf(rh[:idx_lvl,:],
                            transform=ccrs.PlateCarree(),
                            cmap=colors,
                            norm=normalized_levels,
                            levels=contour_levels,
                            zorder=2)
    
    # Plot black outlines on top of the filled rh contours
    contour2 = ax1.contour(rh[:idx_lvl,:],
                           transform=ccrs.PlateCarree(),
                           colors='black',
                           levels=contour_levels,
                           linewidths=0.1,
                           zorder=3)
    
    # Plot contours for the tempisobar variable
    contour3 = ax1.contour(tempisobar[:idx_lvl,:],
                           transform=ccrs.PlateCarree(),
                           cmap=plt.get_cmap('gist_rainbow'),
                           levels=np.append(np.arange(-80,-10,10), np.arange(-10,40,5)),
                           linewidths=0.7,
                           linestyles='solid',
                           zorder=4)
    
    contour3.collections[9].set_linestyle('dashed')
    contour3.collections[9].set_linewidth(1)
    contour3.collections[9].set_color('black')
    
    # Create lists of coordinates where the contour labels are going to go
    # Before creating an axes over top of the contour plot, hover your
    # mouse over the locations where you want to plot the contour labels.
    # The coordinate will show up on the bottom right of the figure window.
    # cont2Labels = [(1.71, 3.82), (5.49, 3.23), (9.53, 4.34), (9.27, 3.53),
    #                (14.08, 4.81), (19.21, 2.24), (17.74, 1.00), (22.23, 3.87),
    #                (12.87, 2.54), (10.45, 6.02), (11.51, 4.92)]
    
    # cont3Labels = [(7.5, 6.1), (10.0, 4.58), (19.06, 1.91), (8.68, 0.46),
    #                (19.52, 4.80), (16.7, 6.07), (8.62, 5.41), (18.53, 5.46)]
    
    # Manually plot contour labels
    cont2labels = ax1.clabel(contour2,
                             # manual=cont2Labels,
                             fmt='%d',
                             inline=True,
                             fontsize=7)
    cont3labels = ax1.clabel(contour3,
                             # manual=cont3Labels,
                             fmt='%d',
                             inline=True,
                             fontsize=7)
                             # colors='black')
    
    # Set contour label backgrounds white
    [
        txt.set_bbox(dict(facecolor='none', edgecolor='none', pad=.5))
        for txt in contour2.labelTexts
    ]
    [
        txt.set_bbox(dict(facecolor='none', edgecolor='none', pad=.5))
        for txt in contour3.labelTexts
    ]
    
    # Determine the labels for each tick on the x and y axes
    yticklabels = np.array(levels, dtype=np.int)
    
    # xticklabels = str_times
    
    # xticklabels = [
    #     '12z', '15z', '18z', '21z', 'Apr29', '03z', '06z', '09z', '12z', '15z',
    #     '18z', '21z', 'Apr30', '03z', '06z', '09z', '12z', '15z', '18z', '21z',
    #     'May01', '03z', '06z', '09z', '12z'
    # ]
    
    # Make an axis to overlay on top of the contour plot
    axin = fig.add_subplot(spec[0, 0])
    
    # Use the geocat.viz function to set the main title of the plot
    gvutil.set_titles_and_labels(axin,
                                 maintitle=str_title,
                                 maintitlefontsize=16,
                                 ylabel='RH (shaded,%) Temp (line,\N{DEGREE SIGN}C) Wind (barbs, m/s)'
                                 '\n (millibars)',
                                 labelfontsize=12)
    
    # Add a pad between the y axis label and the axis spine
    axin.yaxis.labelpad = 3
    
    # Use the geocat.viz function to set axes limits and ticks
    gvutil.set_axes_limits_and_ticks(axin,
                                      xlim=[taus[0], taus[-1]],
                                      ylim=[levels[0], levels[idx_lvl]],
                                      xticks=np.array(taus),
                                      yticks=np.array(levels),
                                      xticklabels=[],
                                      yticklabels=yticklabels)
    
    # Make axis invisible
    axin.patch.set_alpha(0.0)
    
    # Make ticks point inwards
    axin.tick_params(axis="x", direction="out", length=5)
    axin.tick_params(axis="y", direction="out", length=5, labelsize=9)
    axin.tick_params(bottom=True, left=True, right=True, top=False)
    # axin.tick_params(which='minor', top=False, bottom=False)
    
    # Rotate the labels on the x axis so they are vertical
    for tick in axin.get_xticklabels():
        tick.set_rotation(90)
    
    # Set aspect ratio of axin so it lines up with axis underneath (ax1)
    axin.set_aspect(0.018)
    
    # Plot wind barbs
    barbs = axin.barbs(idx_times,
                       levels[:idx_lvl],
                       ugrid[:idx_lvl,:],
                       vgrid[:idx_lvl,:],
                       color='black',
                       lw=0.3,
                       length=5)
    
    # Create text box at lower right of contour plot
    # ax1.text(1.0,
    #           -0.28,"",
    #           # "CONTOUR FROM -20 TO 60 BY 10",
    #           horizontalalignment='right',
    #           transform=ax1.transAxes,
    #           bbox=dict(boxstyle='square, pad=0.15',
    #                     facecolor='white',
    #                     edgecolor='black'))
    
    # Create two more axes, one for the bar chart and one for the line graph
    axin1 = fig.add_subplot(spec[1, 0])
    axin2 = fig.add_subplot(spec[2, 0])
    axin3 = fig.add_subplot(spec[3, 0])
    axin4 = fig.add_subplot(spec[4, 0])
    axin5 = fig.add_subplot(spec[5, 0])
    
    # Plot line chart
    
    # Plot lines depicting the wind10m variable
    axin1.plot(taus, ws10, color='orange', marker='o', markerfacecolor='white')
    axin1.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
    Xq,Yq = np.meshgrid(taus, (np.nanmin(ws10)+np.nanmax(ws10))/2)
    axin1.barbs(Xq,Yq,u10m,v10m,color='black',lw=0.3,length=5)
    
    # Use the geocat.viz function to set the y axis label
    gvutil.set_titles_and_labels(axin1, ylabel='10m Wind \nSpeed & Barbs \n(m/s)', labelfontsize=12)
    
    # Determine the labels for each tick on the x and y axes
    w10m_levels = np.arange(np.nanmin(ws10),np.nanmax(ws10),5)
    # yticklabels = ['59.0', '60.0', '61.0', '62.0', '63.0', '64.0']
    # xticklabels = [
    #     '12z', '', '18z', '', 'Apr29', '', '06z', '', '12z', '', '18z', '', 'Apr30',
    #     '', '06z', '', '12z', '', '18z', '', 'May01', '', '06z', '', '12z'
    # ]
    
    # Use the geocat.viz function to set inset axes limits and ticks
    gvutil.set_axes_limits_and_ticks(axin1,
                                     xlim=[taus[0], taus[-1]],
                                     ylim=[w10m_levels[0], np.nanmax(ws10)],
                                     xticks=np.array(taus),
                                     # yticks=t2m_levels,
                                     xticklabels=[],
                                     # yticklabels=t2m_levels
                                     )
    
    # Use the geocat.viz function to add minor ticks
    gvutil.add_major_minor_ticks(axin1, y_minor_per_major=5, labelsize="small")
    
    # Make ticks only show up on bottom, right, and left of inset axis
    axin1.tick_params(bottom=True, left=True, right=True, top=False)
    axin1.tick_params(which='minor', top=False, bottom=False)
    
    # Rotate the labels on the x axis so they are vertical
    for tick in axin1.get_xticklabels():
        tick.set_rotation(90)
    
    # Plot lines depicting the tempht variable
    axin2.plot(taus, tempht, color='red')
    axin2.plot(taus, dewtempht, color='red', linestyle = '--')
    axin2.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
    
    # Use the geocat.viz function to set the y axis label
    gvutil.set_titles_and_labels(axin2, ylabel='2m Temp \n2m DewPt \n(\N{DEGREE SIGN}C)', labelfontsize=12)
    
    # Determine the labels for each tick on the x and y axes
    t2m_levels = np.arange(np.nanmin(dewtempht),np.nanmax(tempht),0.1)
    # yticklabels = ['59.0', '60.0', '61.0', '62.0', '63.0', '64.0']
    # xticklabels = [
    #     '12z', '', '18z', '', 'Apr29', '', '06z', '', '12z', '', '18z', '', 'Apr30',
    #     '', '06z', '', '12z', '', '18z', '', 'May01', '', '06z', '', '12z'
    # ]
    
    # Use the geocat.viz function to set inset axes limits and ticks
    gvutil.set_axes_limits_and_ticks(axin2,
                                     xlim=[taus[0], taus[-1]],
                                     ylim=[np.nanmin(dewtempht), np.nanmax(tempht)],
                                     xticks=np.array(taus),
                                     # yticks=t2m_levels,
                                     xticklabels=[],
                                     # yticklabels=t2m_levels
                                     )
    
    # Use the geocat.viz function to add minor ticks
    gvutil.add_major_minor_ticks(axin2, y_minor_per_major=5, labelsize="small")
    
    # Make ticks only show up on bottom, right, and left of inset axis
    axin2.tick_params(bottom=True, left=True, right=True, top=False)
    axin2.tick_params(which='minor', top=False, bottom=False)
    
    # Rotate the labels on the x axis so they are vertical
    # for tick in axin2.get_xticklabels():
    #     tick.set_rotation(90)
        
    # Plot lines depicting the rhht variable
    axin3.plot(taus, rhht, color='green', marker='o', markerfacecolor='white')
    axin3.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
    
    # Use the geocat.viz function to set the y axis label
    gvutil.set_titles_and_labels(axin3, ylabel='2m RH \n(%)', labelfontsize=12)
    
    # Determine the labels for each tick on the x and y axes
    rh2m_levels = np.arange(np.nanmin(rhht),np.nanmax(rhht),10)
    # yticklabels = ['59.0', '60.0', '61.0', '62.0', '63.0', '64.0']
    # xticklabels = [
    #     '12z', '', '18z', '', 'Apr29', '', '06z', '', '12z', '', '18z', '', 'Apr30',
    #     '', '06z', '', '12z', '', '18z', '', 'May01', '', '06z', '', '12z'
    # ]
    
    # Use the geocat.viz function to set inset axes limits and ticks
    gvutil.set_axes_limits_and_ticks(axin3,
                                     xlim=[taus[0], taus[-1]],
                                     ylim=[np.nanmin(rhht), np.nanmax(rhht)],
                                     xticks=np.array(taus),
                                     # yticks=t2m_levels,
                                     xticklabels=[],
                                     # yticklabels=t2m_levels
                                     )
    
    # Use the geocat.viz function to add minor ticks
    gvutil.add_major_minor_ticks(axin3, y_minor_per_major=5, labelsize="small")
    
    # Make ticks only show up on bottom, right, and left of inset axis
    axin3.tick_params(bottom=True, left=True, right=True, top=False)
    axin3.tick_params(which='minor', top=False, bottom=False)
    
    # Rotate the labels on the x axis so they are vertical
    # for tick in axin2.get_xticklabels():
    #     tick.set_rotation(90)
        
    # Plot bars depicting the cloud fraction variable
    barWidth = 0.15
    br1 = taus
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]
     
    # Make the plot
    axin4.bar(br1, clflo*100, color ='darkblue', width = barWidth,
            edgecolor ='darkblue', label ='low')
    axin4.bar(br2, clfmi*100, color ='blue', width = barWidth,
            edgecolor ='blue', label ='middle')
    axin4.bar(br3, clfhi*100, color ='lightblue', width = barWidth,
            edgecolor ='lightblue', label ='high')

    # axin4.bar(taus,
    #           rain03,
    #           width=1.0,
    #           color='limegreen',
    #           edgecolor='black',
    #           linewidth=.2)
    axin4.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
    axin4.legend(fontsize = 5, ncol = 3, loc = 'upper right')
    
    # Use the geocat.viz function to set the y axis label
    gvutil.set_titles_and_labels(axin4, ylabel='Cloud Cover \n(%)', labelfontsize=12)
    
    # Determine the labels for each tick on the x and y axes
    clf_levels = np.arange(0,101,20)
            
    # yticklabels = ['0.0', '0.10', '0.20', '0.30', '0.40', '0.50']
    # xticklabels = str_times
    # xticklabels = [
    #     '12z', '', '18z', '', 'Apr29', '', '06z', '', '12z', '', '18z', '', 'Apr30',
    #     '', '06z', '', '12z', '', '18z', '', 'May01', '', '06z', '', '12z'
    # ]
    
    # Use the geocat.viz function to set axes limits and ticks
    gvutil.set_axes_limits_and_ticks(axin4,
                                     xlim=[taus[0], taus[-1]],
                                     ylim=[0, 120],
                                     xticks=np.array(taus),
                                     yticks=clf_levels,
                                      xticklabels=[],
                                     # yticklabels=rain3h_levels
                                     )
    
    # Use the geocat.viz function to add minor ticks
    gvutil.add_major_minor_ticks(axin4, y_minor_per_major=5, labelsize="small")
    
    # Make ticks only show up on bottom, right, and left of inset axis
    axin4.tick_params(bottom=True, left=True, right=True, top=False)
    axin4.tick_params(which='minor', top=False, bottom=False)
    
    # Rotate the labels on the x axis so they are vertical
    # for tick in axin4.get_xticklabels():
    #     tick.set_rotation(90)

    # Plot bars depicting the rain03 variable
    axin5.bar(taus,
              rain03,
              width=1.0,
              color='limegreen',
              edgecolor='black',
              linewidth=.2)
    axin5.bar(taus[0], np.nanmax(rain03),width=1,color='lightgrey',edgecolor='lightgrey',linewidth=.2)
    axin5.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
        
    # Use the geocat.viz function to set the y axis label
    gvutil.set_titles_and_labels(axin5, ylabel='3hr rain total \n(mm)', labelfontsize=12)
    
    # Determine the labels for each tick on the x and y axes
    rain3h_levels = np.arange(0,np.nanmax(rain03), 0.1)
    if np.nanmax(rain03) < 0.1:
        ymax = 0.1
    else:
        ymax = rain3h_levels[-1]
    
    axin5.text(len(taus) - 7, 0.8*ymax, r'3-Day Total = ' + str(round(np.nansum(rain03), 3)))
        
    # yticklabels = ['0.0', '0.10', '0.20', '0.30', '0.40', '0.50']
    # xticklabels = str_times
    # xticklabels = [
    #     '12z', '', '18z', '', 'Apr29', '', '06z', '', '12z', '', '18z', '', 'Apr30',
    #     '', '06z', '', '12z', '', '18z', '', 'May01', '', '06z', '', '12z'
    # ]
    
    # Use the geocat.viz function to set axes limits and ticks
    gvutil.set_axes_limits_and_ticks(axin5,
                                     xlim=[taus[0], taus[-1]],
                                     ylim=[0, ymax],
                                     xticks=np.array(taus),
                                     # yticks=rain3h_levels,
                                      xticklabels=str_times,
                                     # yticklabels=rain3h_levels
                                     )
    
    # Use the geocat.viz function to add minor ticks
    gvutil.add_major_minor_ticks(axin5, y_minor_per_major=5, labelsize="small")
    
    # Make ticks only show up on bottom, right, and left of inset axis
    axin5.tick_params(bottom=True, left=True, right=True, top=False)
    axin5.tick_params(which='minor', top=False, bottom=False)
    
    # Rotate the labels on the x axis so they are vertical
    for tick in axin5.get_xticklabels():
        tick.set_rotation(90)
   
    # Adjust space between the first and second axes on the plot
    plt.subplots_adjust(hspace=-0.1)
    
    # plt.show()
    
    try:
        os.mkdir(dirout + inittime + '\\')
    except:
        a = 1    
    fig.savefig(dirout + inittime + '\\meteogram_'+inittime+'_'+area_name.replace('/' , '-')+'.png', format='png', dpi=100, bbox_inches='tight')
    plt.cla()
    
if __name__ == "__main__":
   main(sys.argv[1:])
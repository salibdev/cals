#!/usr/bin/python

# Calculating stellar chromospheric activity parameters and plotting spectrum diagram of Ca II H&K lines with LAMOST LRS spectra.

import math
import os

import numpy as np
from astropy.io import fits
import scipy
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator



class CaIISindex:
    '''
    Calculate stellar chromospheric activity parameters and plot spectrum diagram of Ca II H&K lines with LAMOST LRS spectra.
    
    If you want to display the values of stellar parameters (teff, logg, feh) in the spectrum figure,
    you should provide them via the 'stellar_parameters',
    otherwise they will be displayed as 'unknown' in the figure.
    
    Example:
    data_path = "example.fits"
    cals = CaIISindex(data_path)
    cals.Sindex_savepath = 'S_index_example.csv'                # set csv file save path
    cals.figure_savepath = 'example.png'                        # set figure file save path
    cals.stellar_parameters = ['5534.22','4.423','-0.025']      # set stellar parameters [teff (K), logg (dex), feh (dex)]
    cals.calcSindex()                 # return a dict of stellar activity parameters
    cals.calcError()                  # return a dict of activity parameter errors
    cals.plotSpectrum()               # plot spectrum diagram
    cals.saveRecord()                 # save calculated activity parameters to a csv file and spectrum diagram to a figure file

    '''

    figure_savepath = None
    Sindex_savepath = None
    stellar_parameters = ['unknown','unknown','unknown']
    L_K = 3934.78
    L_H = 3969.59
    L_R = 4002.20
    L_V = 3902.17
    step1 = 0.01
    step2 = 0.001
    FWHM = 1.09
    para_catalog = ['R_mean','R_mean_err','V_mean','V_mean_err',
             	'H_mean_tri','H_mean_tri_err','K_mean_tri','K_mean_tri_err','S_tri','S_tri_err',
         	'S_MWL','S_MWL_err',
         	'H_mean_rec','H_mean_rec_err','K_mean_rec','K_mean_rec_err','S_rec','S_rec_err',
		'condition_tag']
    begin_index = 0
    end_index = None

    def __init__(self,path):
        '''

        Initialize the fits data path.

        '''
        self.path = path
        self.__prework()
        self.__loadData(path)

    def wave_band(self,band_center,band_width,step):
        rec_band = np.arange(band_center-band_width/2+0.5*step,band_center+band_width/2,step)
        return rec_band

    def tri_func(self,step):
        tem_tri = list(np.arange(step/self.FWHM/2,1,step/self.FWHM))
        trig = tem_tri+tem_tri[::-1]
        return trig
        
    def __prework(self):    
        self.band_V = self.wave_band(band_center=self.L_V, band_width=20, step=self.step1)
        self.band_R = self.wave_band(band_center=self.L_R, band_width=20, step=self.step1)
        self.band_K = self.wave_band(band_center=self.L_K, band_width=1 , step=self.step2)
        self.band_H = self.wave_band(band_center=self.L_H, band_width=1 , step=self.step2)
        self.band_tri_K = self.wave_band(band_center=self.L_K, band_width=self.FWHM*2, step=self.step2)
        self.band_tri_H = self.wave_band(band_center=self.L_H, band_width=self.FWHM*2, step=self.step2)

        self.trig = self.tri_func(self.step2)
    
    def __loadData(self,path,printname=False):
        fitsdata = fits.open(path)
        self.specdata = fitsdata[0]
        if self.specdata.header['TELESCOP']!='LAMOST':
            print('warning: not LAMOST data')
        if printname:
            print('fitsname '+self.specdata.header['FILENAME'])
        #read FITS file name
        file_name = os.path.basename(self.path)
        file_name = os.path.splitext(file_name)[0]
        self.figure_savepath = '{}.png'.format(file_name)
        self.Sindex_savepath = 'S_index_{}.csv'.format(file_name)
        #treat spectrum data
        self.spec = self.specdata.data
        self.flux = self.spec[0][self.begin_index:self.end_index]
        self.error = self.spec[1][self.begin_index:self.end_index]
        self.wavelen = self.spec[2][self.begin_index:self.end_index]

        self.wave_z = [i/(1+float(self.specdata.header['Z'])) for i in self.wavelen]
        self.flux_func = scipy.interpolate.interp1d(self.wave_z,self.flux,kind='linear')
        return self.specdata

    def __calc_trig(self,x):
        y = self.flux_func(x)
        newy = [y[i]*self.trig[i] for i in range(len(y))]
        total = sum(newy)*self.step2
        return total/1.09
    
    def calcSindex(self):
        '''
        
        Return a dict of stellar activity parameters:
        para_dict = {'R_mean':R_mean,           # mean flux of 20Å R band
                     'V_mean':V_mean,           # mean flux of 20Å V band
                     'H_mean_tri':H_mean_tri,   # mean flux of H line in 1.09Å FWHM triangular bandpass
                     'K_mean_tri':K_mean_tri,   # mean flux of K line in 1.09Å FWHM triangular bandpass
                     'S_tri':S_tri,             # S index using 1.09Å FWHM triangular bandpass
                     'S_MWL':S_MWL                # Mount Wilson S index
                     'H_mean_rec':H_mean_rec,   # mean flux of H line in 1Å rectangular bandpass
                     'K_mean_rec':K_mean_rec,   # mean flux of K line in 1Å rectangular bandpass
                     'S_rec':S_rec,             # S index using 1Å rectangular bandpass
                    }
        
        '''
        R_mean = sum(self.flux_func(self.band_R))*self.step1/20. 
        V_mean = sum(self.flux_func(self.band_V))*self.step1/20.
        H_mean_rec = sum(self.flux_func(self.band_H))*self.step2/1.
        K_mean_rec = sum(self.flux_func(self.band_K))*self.step2/1.
        H_mean_tri = self.__calc_trig(self.band_tri_H)
        K_mean_tri = self.__calc_trig(self.band_tri_K)
        S_rec = (H_mean_rec+K_mean_rec)/(R_mean+V_mean)
        S_tri = (H_mean_tri+K_mean_tri)/(R_mean+V_mean)
        S_MWL = 8*1.8*1.09/20.*(H_mean_tri+K_mean_tri)/(R_mean+V_mean)
        S_info = [R_mean,V_mean,H_mean_tri,K_mean_tri,S_tri,S_MWL,H_mean_rec,K_mean_rec,S_rec]
        self.para_dict = dict(zip(self.para_catalog[::2][:-1],S_info))
        return self.para_dict

    def __getOrigError(self):
        original_error = []
        for i in range(len(self.error)):
            if self.error[i]!=0:
                original_error.append(math.sqrt(1/self.error[i]))
            else:
                original_error.append(0)
        orig_error_bin = [2/(self.wave_z[i+1]-self.wave_z[i-1]) for i in range(1,len(self.wave_z)-1)]
        orig_error_bin = [1/(self.wave_z[1]-self.wave_z[0])]+orig_error_bin+[1/(self.wave_z[-2]-self.wave_z[-1])]
        self.f_bin = scipy.interpolate.interp1d(self.wave_z,orig_error_bin,kind='linear')
        self.f_err = scipy.interpolate.interp1d(self.wave_z,original_error,kind='linear')
    
    def __calc_RV_err(self):
        R_err_list = self.f_err(self.band_R)
        V_err_list = self.f_err(self.band_V)
        R_bin_list = self.f_bin(self.band_R)
        V_bin_list = self.f_bin(self.band_V)
        if min(R_err_list)<=0:
            R_err = -9999
        else:
            tem_err = [R_err_list[i]**2/self.step1/R_bin_list[i] for i in range(len(R_err_list))]
            R_err = math.sqrt(sum(tem_err))*self.step1/20.
        if min(V_err_list)<=0:
            V_err = -9999
        else:
            tem_err = [V_err_list[i]**2/self.step1/V_bin_list[i] for i in range(len(V_err_list))]
            V_err = math.sqrt(sum(tem_err))*self.step1/20.
        self.S_info_err[0] = R_err
        self.S_info_err[1] = V_err

    def __calc_HK_tri_err(self):
        H_err_list = self.f_err(self.band_tri_H)
        K_err_list = self.f_err(self.band_tri_K)
        H_bin_list = self.f_bin(self.band_tri_H)
        K_bin_list = self.f_bin(self.band_tri_K)
        if min(H_err_list)<=0:
            H_tri_err = -9999
        else:
            tem_err = [H_err_list[i]**2/self.step2/H_bin_list[i]*self.trig[i]**2 for i in range(len(H_err_list))]
            H_tri_err = math.sqrt(sum(tem_err))*self.step2/1.09
        if min(K_err_list)<=0:
            K_tri_err = -9999
        else:
            tem_err = [K_err_list[i]**2/self.step2/K_bin_list[i]*self.trig[i]**2 for i in range(len(K_err_list))]
            K_tri_err = math.sqrt(sum(tem_err))*self.step2/1.09
        self.S_info_err[2] = H_tri_err
        self.S_info_err[3] = K_tri_err
        if -9999 not in [self.S_info_err[0],self.S_info_err[1],self.S_info_err[5],self.S_info_err[6]]:
            RV_err = math.sqrt(self.S_info_err[0]**2+self.S_info_err[1]**2)
            RV_err = RV_err/(self.para_dict['R_mean']+self.para_dict['V_mean'])
            HK_tri_err = math.sqrt(H_tri_err**2+K_tri_err**2)
            HK_tri_err = HK_tri_err/(self.para_dict['H_mean_tri']+self.para_dict['K_mean_tri'])
            S_tri_err = self.para_dict['S_tri']*math.sqrt(RV_err**2+HK_tri_err**2)
            S_MWL_err = S_tri_err*8*1.8*1.09/20.
        else:
            S_tri_err=-9999
            S_MWL_err=-9999
        self.S_info_err[4] = S_tri_err
        self.S_info_err[5] = S_MWL_err

    def __calc_HK_rec_err(self):
        H_err_list = self.f_err(self.band_H)
        K_err_list = self.f_err(self.band_K)
        H_bin_list = self.f_bin(self.band_H)
        K_bin_list = self.f_bin(self.band_K)
        if min(H_err_list)<=0:
            H_rec_err = -9999
        else:
            tem_err = [H_err_list[i]**2/self.step2/H_bin_list[i] for i in range(len(H_err_list))]
            H_rec_err = math.sqrt(sum(tem_err))*self.step2/1.
        if min(K_err_list)<=0:
            K_rec_err = -9999
        else:
            tem_err = [K_err_list[i]**2/self.step2/K_bin_list[i] for i in range(len(K_err_list))]
            K_rec_err = math.sqrt(sum(tem_err))*self.step2/1.
        self.S_info_err[6] = H_rec_err
        self.S_info_err[7] = K_rec_err
        if -9999 not in self.S_info_err[0:4]:
            RV_err = math.sqrt(self.S_info_err[0]**2+self.S_info_err[0]**2)
            RV_err = RV_err/(self.para_dict['R_mean']+self.para_dict['V_mean'])
            HK_rec_err = math.sqrt(H_rec_err**2+K_rec_err**2)
            HK_rec_err = HK_rec_err/(self.para_dict['H_mean_rec']+self.para_dict['K_mean_rec'])
            S_rec_err = self.para_dict['S_rec']*math.sqrt(RV_err**2+HK_rec_err**2)
        else:
            S_rec_err=-9999
        self.S_info_err[8] = S_rec_err

    def calcError(self):
        '''
        
        Return a dictionary of activity parameter errors:
        para_err_dict = {'R_mean_err':R_mean_err, 
                         'V_mean_err':V_mean_err, 
                         'H_mean_tri_err':H_mean_tri_err, 
                         'K_mean_tri_err':K_mean_tri_err,
                         'S_tri_err':S_tri_err, 
                         'S_MWL_err':S_MWL_err,
                         'H_mean_rec_err':H_mean_rec_err, 
                         'K_mean_rec_err':K_mean_rec_err, 
                         'S_rec_err':S_rec_err
                        }
        
        '''
        
        self.calcSindex()
        original_error = self.__getOrigError()
        self.S_info_err = len(self.para_dict)*[0]
        self.__calc_RV_err()    
        self.__calc_HK_tri_err()
        self.__calc_HK_rec_err()
        para_dict_err = dict(zip(self.para_catalog[1::2],self.S_info_err))
        return para_dict_err

    def __getCalcInfo(self):
        self.calcError()
        if -9999 not in self.S_info_err:
            condition_tag = ' '
        else:
            condition_tag = 'uncertainty=-9999'
        S_info = ['{:.6f}'.format(value) for key,value in self.para_dict.items()]
        S_info_err = ['{:.6f}'.format(i) for i in self.S_info_err]
        new_info = []
        for i in range(len(S_info)):
            new_info+=[S_info[i],S_info_err[i]]
        new_info += condition_tag
        return new_info
    
    def __recordSErr(self,new_info,header=True):
        obsid = self.specdata.header['OBSID']
        fitsname = self.specdata.header['FILENAME']
        record_info = [str(obsid),fitsname] + new_info
        record_catalog = ['obsid','fitsname'] + self.para_catalog
        f=open(self.Sindex_savepath,'w')
        if header:
            f.write('|'.join(record_catalog)+'\n')
        f.write('|'.join(record_info)+'\n')
        f.close()
        return [record_catalog,record_info]
    
    def __data_plot(self):
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.figure(figsize=(14,7))
        ax1 = plt.gca()
        all_band = np.arange(self.L_V-10,self.L_R+10,self.step1)
        plot_flux = self.flux_func(all_band)
        fig1 = ax1.plot(all_band,plot_flux/max(plot_flux),color='red',label='shifted',linewidth='1.5')
        ax1.plot([self.L_V-10,self.L_R+10],[1,1],linestyle='--',color='black')
        ax2 = ax1.twinx()
        fig2 = ax2.plot(self.wavelen,self.flux,label='original',linewidth='1.5',linestyle='-.')
        c1 = 'orange'
        ax1.fill_between([self.L_R-10,self.L_R+10],[1,1],[0,0],color=c1)
        ax1.fill_between([self.L_V-10,self.L_V+10],[1,1],[0,0],color=c1)
        ax1.fill_between( [self.L_H-1,self.L_H+1], [1,1],[0,0],color=c1)
        ax1.fill_between( [self.L_K-1,self.L_K+1], [1,1],[0,0],color=c1)

        #ax1.plot([self.L_H-1.09,self.L_H],[0,1],color='lime',linestyle='--')
        #ax1.plot([self.L_H,self.L_H+1.09],[1,0],color='lime',linestyle='--')
        #ax1.plot([self.L_K-1.09,self.L_K],[0,1],color='lime',linestyle='--')
        #ax1.plot([self.L_K,self.L_K+1.09],[1,0],color='lime',linestyle='--')
        #ax1.fill_between(self.band_tri_H,self.trig,[0]*len(self.trig),color='lime')
        #ax1.fill_between(self.band_tri_K,self.trig,[0]*len(self.trig),color='lime')
        ax1.fill_between([self.L_H-1.09,self.L_H,self.L_H+1.09],[0,0,0],[0,1,0],color='lime')
        ax1.fill_between([self.L_K-1.09,self.L_K,self.L_K+1.09],[0,0,0],[0,1,0],color='lime')
        
        ax1.plot([self.L_H,self.L_H],[0,1],color='black',linestyle='--')
        ax1.plot([self.L_K,self.L_K],[0,1],color='black',linestyle='--')

        ax1.annotate("Ca II H "+str(self.L_H)+r'$\,$'+u'\u00C5',xy=(self.L_H,0.1),xytext=(self.L_H+4,0.15),
                color="black",weight="bold",fontsize=11,arrowprops=dict(arrowstyle="->",color="black"))
        ax1.annotate("Ca II K "+str(self.L_K)+r'$\,$'+u'\u00C5',xy=(self.L_K,0.1),xytext=(self.L_K+4,0.15),
                color="black",weight="bold",fontsize=11,arrowprops=dict(arrowstyle="->",color="black"))
        ax2.legend(handles=fig1+fig2,loc='upper left',fontsize='15')
        ax2.set_ylim(0,max(plot_flux)*1.2)
        return ax1,ax2

    def __plotset(self,ax1,ax2):
        xmajorLocator = MultipleLocator(10)
        ax1.xaxis.set_major_locator(xmajorLocator)
        ymajorLocator = MultipleLocator(0.2)
        ax1.yaxis.set_major_locator(ymajorLocator)
        xminorLocator = MultipleLocator(2)
        ax1.xaxis.set_minor_locator(xminorLocator)
        yminorLocator = MultipleLocator(0.05)
        ax1.yaxis.set_minor_locator(yminorLocator)
        ax1.tick_params(which='major',length=8,width=1.5,labelsize=15)
        ax1.tick_params(which='minor',length=4,width=1)
        ax1.set_xlabel(r'$\rm{Vacuum\ Wavelength\ (\AA)}$',fontsize=15)
        ax1.set_ylabel('Relative Flux',fontsize=15)
        ax1.set_xlim(self.L_V-10,self.L_R+10)
        ax1.set_ylim(0,1.2)
        ax2.set_ylabel('Flux',fontsize=15)
        ax2.tick_params(which='major',length=8,width=1.5,labelsize=15)
        return ax1,ax2

    def __plotInfo(self,ax1,ax2,calced_info):
        para_catalog = self.para_catalog
        rv_str = '{:.2f}'.format(float(self.specdata.header['Z'])*299792.458)
        rv_str = '$'+rv_str+'$'
        ax1.text(self.L_V+9,1.03,'Radial Velocity\n='+r'$\,$'+rv_str+r'$\,$'+'km/s\n', fontsize=15)
        ax1.text(self.L_K-4,1.03,'H_mean_rec={:.3f}\nK_mean_rec={:.3f}\nS_rec={:.3f}'.format(float(calced_info[para_catalog.index('H_mean_rec')]),
                        float(calced_info[para_catalog.index('K_mean_rec')]),float(calced_info[para_catalog.index('S_rec')])),fontsize=15)
        ax1.text(self.L_H-10,1.03,'H_mean_tri={:.3f}\nK_mean_tri={:.3f}\nS_tri={:.3f}'.format(float(calced_info[para_catalog.index('H_mean_tri')]),
                        float(calced_info[para_catalog.index('K_mean_tri')]),float(calced_info[para_catalog.index('S_tri')])),fontsize=15)
        ax1.text(self.L_R-15,1.03,'R_mean={:.3f}\nV_mean={:.3f}\nS_MWL={:.3f}'.format(float(calced_info[para_catalog.index('R_mean')]),
                        float(calced_info[para_catalog.index('V_mean')]),float(calced_info[para_catalog.index('S_MWL')])),fontsize=15)
        ax1.text(self.L_V-5,0.2,r'$V$ band',fontsize=15)
        ax1.text(self.L_R-5,0.2,r'$R$ band',fontsize=15)
        teff_str = '{}'.format(self.stellar_parameters[0])
        if (teff_str != 'unknown'):
            teff_str = teff_str.strip()
            if teff_str[-1] == 'K':
                teff_str = teff_str[:-1]
                teff_str = teff_str.strip()
            teff_str = teff_str + r'$\,$K'
        logg_str = '{}'.format(self.stellar_parameters[1])
        logg_str = logg_str.strip()
        feh_str = '{}'.format(self.stellar_parameters[2])
        feh_str = feh_str.strip()
        feh_str = '$'+feh_str+'$'
        snrg_str = '{}'.format(self.specdata.header['SNRG'])
        snrr_str = '{}'.format(self.specdata.header['SNRR'])
        title_info = r'$T_\mathrm{eff}$=' + teff_str + r'   $\log\,g$=' + logg_str + '   [Fe/H]='+ feh_str + r'   $\mathrm{SNR}_g$=' + snrg_str + r'   $\mathrm{SNR}_r$=' + snrr_str
        ax1.text((self.L_R + self.L_V)/2.0, 1.22, title_info, fontsize=15, horizontalalignment='center')
        ax1.set_title('Obsid-{} {} ({})\n'.format(self.specdata.header['OBSID'],self.specdata.header['FILENAME'],self.specdata.header['DATA_V']),fontsize=15)
        return ax1,ax2

    def plotSpectrum(self, figure_savepath=None, stellar_parameters=None, figshow=True):
        '''
        
        Generate spectrum diagram.
        
        Optional arguments:
            figure_savepath:        set figure file save path
            stellar_parameters:     set stellar parameters [teff (K), logg (dex), feh (dex)]
            figshow:                show spectrum diagram on screen (True/False, default True)
        
        '''
        if stellar_parameters:
            self.stellar_parameters = stellar_parameters
        calced_info = self.__getCalcInfo()
        ax1,ax2 = self.__data_plot()
        ax1,ax2 = self.__plotset(ax1,ax2)
        self.ax1,self.ax2 = self.__plotInfo(ax1,ax2,calced_info)
        if figure_savepath:
            self.figure_savepath = figure_savepath
        plt.savefig(self.figure_savepath,bbox_inches='tight')
        if figshow:
            plt.show()
        return calced_info
 
    def saveRecord(self, Sindex_savepath=None, figure_savepath=None, stellar_parameters=None, csv_header=True, figshow=False):
        '''
        
        Save stellar activity parameters to a csv file and spectrum diagram to a figure file.
            
        Optional arguments:
            Sindex_savepath:        set csv file save path
            figure_savepath:        set figure file save path
            stellar_parameters:     set stellar parameters [teff (K), logg (dex), feh (dex)]
            csv_header:             write out column names in csv file (True/False, default True)
            figshow:                show spectrum diagram on screen (True/False, default False)
        
        '''
        new_info = self.plotSpectrum(figure_savepath, stellar_parameters, figshow)
        if Sindex_savepath:
            self.Sindex_savepath = Sindex_savepath
        record_info = self.__recordSErr(new_info,csv_header)
        record_dict = dict(zip(record_info[0],record_info[1]))
        return record_dict


if __name__ == '__main__':
    #pass
    data_path = "example.fits"
    example = CaIISindex(data_path)
    example.stellar_parameters = ['5534.22','4.423','-0.025']
    example.saveRecord()
    


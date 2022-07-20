#!/usr/bin/python

# Calculating stellar chromospheric activity parameters and plotting spectrum diagram of Ca II H&K lines with LAMOST LRS spectra.

import math
import os
import argparse
import time

import numpy as np
from astropy.io import fits
import scipy
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def wave_band(band_center,band_width,step):
    rec_band = np.arange(band_center-band_width/2+0.5*step,band_center+band_width/2,step)
    return rec_band

def tri_func(step,FWHM):
    tem_tri = list(np.arange(step/FWHM/2,1,step/FWHM))
    trig = tem_tri+tem_tri[::-1]
    return trig

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
    
    save_path = ''
    stellar_parameters = ['unknown','unknown','unknown']

    L_K = 3934.78   #Center wavelength of K windows in vacuum
    L_H = 3969.59   #Center wavelength of H windows in vacuum
    L_R = 4002.20   #Center wavelength of R windows in vacuum
    L_V = 3902.17   #Center wavelength of V windows in vacuum

    step1 = 0.01
    step2 = 0.001
    FWHM = 1.09
    alpha = 1.8
    para_catalog = ('R_mean','R_mean_err','V_mean','V_mean_err',
             	    'H_mean_tri','H_mean_tri_err','K_mean_tri','K_mean_tri_err','S_tri','S_tri_err',
         	    'S_MWL','S_MWL_err',
                    'H_mean_rec','H_mean_rec_err','K_mean_rec','K_mean_rec_err','S_rec','S_rec_err')
                    
    begin_index = 0
    end_index = None
    interp_type = 'linear'

    def __init__(self,path=False):
        '''

        Initialize the fits data path.

        '''

        self.all_para_dict = dict(zip(self.para_catalog,['unknown']*len(self.para_catalog)))
        self.__prework()
        if path:
            #self.path = path
            self.__loadData(path)
    
    def __prework(self):    
        self.band_V = wave_band(band_center=self.L_V, band_width=20, step=self.step1)
        self.band_R = wave_band(band_center=self.L_R, band_width=20, step=self.step1)
        self.band_K_rec = wave_band(band_center=self.L_K, band_width=1, step=self.step2)
        self.band_H_rec = wave_band(band_center=self.L_H, band_width=1, step=self.step2)
        self.band_K_tri = wave_band(band_center=self.L_K, band_width=self.FWHM*2, step=self.step2)
        self.band_H_tri = wave_band(band_center=self.L_H, band_width=self.FWHM*2, step=self.step2)
        self.all_band = np.arange(self.L_V-10,self.L_R+10,self.step1)
        
        self.trig = tri_func(step=self.step2, FWHM=self.FWHM)
    
    def __loadData(self,path,printname=False):
        self.path = path
        fitsdata = fits.open(path)
        self.specdata = fitsdata[0]
        if self.specdata.header['TELESCOP']!='LAMOST':
            print('warning: not LAMOST data!')
        if printname:
            print('fitsname '+self.specdata.header['FILENAME'])
        #read FITS file name
        file_name = os.path.basename(self.path)
        file_name = os.path.splitext(file_name)[0]
        self.figure_savepath = self.save_path + '{}.png'.format(file_name)
        self.Sindex_savepath = self.save_path + 'S_index_{}.csv'.format(file_name)
        #treat spectrum data
        self.spec = self.specdata.data
        self.flux = self.spec[0][self.begin_index:self.end_index]
        self.error = self.spec[1][self.begin_index:self.end_index]
        self.wavelen = self.spec[2][self.begin_index:self.end_index]
        self.z = float(self.specdata.header['Z'])

    def __wave_calib(self):
        self.wave_z = [i/(1+self.z) for i in self.wavelen]
        self.flux_func = scipy.interpolate.interp1d(self.wave_z,self.flux,kind=self.interp_type)


    def __calc_tri_mean(self,bin_center,step,band_width):
        y = self.flux_func(bin_center)
        newy = [y[i]*self.trig[i] for i in range(len(y))]
        tri_mean = sum(newy)*step/band_width
        return tri_mean

    def calc_H_tri_mean(self):
        H_tri_mean = self.__calc_tri_mean(bin_center=self.band_H_tri,step=self.step2,band_width=self.FWHM)
        self.all_para_dict['H_mean_tri'] = H_tri_mean
        return H_tri_mean

    def calc_K_tri_mean(self):
        K_tri_mean = self.__calc_tri_mean(bin_center=self.band_K_tri,step=self.step2,band_width=self.FWHM)
        self.all_para_dict['K_mean_tri'] = K_tri_mean
        return K_tri_mean

    def __calc_rec_mean(self,bin_center,step,band_width):
        rec_mean = sum(self.flux_func(bin_center))*step/band_width
        return rec_mean
    
    def calc_R_mean(self):
        R_mean = self.__calc_rec_mean(bin_center=self.band_R,step=self.step1,band_width=20)
        self.all_para_dict['R_mean'] = R_mean
        return R_mean

    def calc_V_mean(self):
        V_mean = self.__calc_rec_mean(bin_center=self.band_V,step=self.step1,band_width=20)
        self.all_para_dict['V_mean'] = V_mean
        return V_mean

    def calc_H_rec_mean(self):
        H_rec_mean = self.__calc_rec_mean(bin_center=self.band_H_rec,step=self.step2,band_width=1)
        self.all_para_dict['H_mean_rec'] = H_rec_mean
        return H_rec_mean

    def calc_K_rec_mean(self):
        K_rec_mean = self.__calc_rec_mean(bin_center=self.band_K_rec,step=self.step2,band_width=1)
        self.all_para_dict['K_mean_rec'] = K_rec_mean
        return K_rec_mean

    def __check_calc_para(self,calc_para):
        for i in calc_para:
            if i not in self.para_catalog:
                error = ('can not calculate '+i+', please check the parameter name！'
                    +' CaIISindex can only be used to calculate: '+', '.join(para_catalog[:-1]))
                raise NameError(error)
        
    def calcSindex(self):
        '''
        
        Return a dict of stellar activity parameters:
        para_dict = {'R_mean':R_mean,           # mean flux of 20Å R band
                     'V_mean':V_mean,           # mean flux of 20Å V band
                     'H_mean_tri':H_mean_tri,   # mean flux of H line in self.FWHMÅ FWHM triangular bandpass
                     'K_mean_tri':K_mean_tri,   # mean flux of K line in self.FWHMÅ FWHM triangular bandpass
                     'S_tri':S_tri,             # S index using self.FWHMÅ FWHM triangular bandpass
                     'S_MWL':S_MWL                # Mount Wilson S index
                     'H_mean_rec':H_mean_rec,   # mean flux of H line in 1Å rectangular bandpass
                     'K_mean_rec':K_mean_rec,   # mean flux of K line in 1Å rectangular bandpass
                     'S_rec':S_rec,             # S index using 1Å rectangular bandpass
                    }
        
        '''
        
        self.__wave_calib()
        
        R_mean = self.calc_R_mean()
        V_mean = self.calc_V_mean()
        H_mean_rec = self.calc_H_rec_mean()
        K_mean_rec = self.calc_K_rec_mean()
        H_mean_tri = self.calc_H_tri_mean()
        K_mean_tri = self.calc_K_tri_mean()
        
        S_rec = (H_mean_rec+K_mean_rec)/(R_mean+V_mean)
        S_tri = (H_mean_tri+K_mean_tri)/(R_mean+V_mean)
        S_MWL = 8*self.alpha*self.FWHM/20.*(H_mean_tri+K_mean_tri)/(R_mean+V_mean)

        self.all_para_dict['S_rec'] = S_rec
        self.all_para_dict['S_tri'] = S_tri
        self.all_para_dict['S_MWL'] = S_MWL
        
        S_info = [self.all_para_dict[i] for i in self.para_catalog[::2]]
        para_dict = dict(zip(self.para_catalog[::2],S_info))
        return para_dict

    def __calcErr_from_zplus(self,all_para_dict):        
        if float(self.specdata.header['Z_ERR'])<0:
            print(self.specdata.header['OBSID'],self.specdata.header['FILENAME'],self.specdata.header['Z'],self.specdata.header['Z_ERR'])
            for i in ['R_mean','V_mean','H_mean_tri','K_mean_tri','S_tri','S_MWL','H_mean_rec','K_mean_rec','S_rec']:
                self.all_para_dict[i+'_err'] = -9999 
        else:
            self.z = float(self.specdata.header['Z'])+float(self.specdata.header['Z_ERR'])
            self.calcSindex()
            for i in ['R_mean','V_mean','H_mean_tri','K_mean_tri','S_tri','S_MWL','H_mean_rec','K_mean_rec','S_rec']:
                self.all_para_dict[i+'_err'] = abs(all_para_dict[i] - self.all_para_dict[i])
        
        S_info_err = [self.all_para_dict[i] for i in self.para_catalog[1::2]]
        para_dict_err = dict(zip(self.para_catalog[1::2],S_info_err))
        return para_dict_err

    
    def __calcErr_from_zminus(self,all_para_dict):
        if float(self.specdata.header['Z_ERR'])<0:
            print(self.specdata.header['OBSID'],self.specdata.header['FILENAME'],self.specdata.header['Z'],self.specdata.header['Z_ERR'])
            for i in ['R_mean','V_mean','H_mean_tri','K_mean_tri','S_tri','S_MWL','H_mean_rec','K_mean_rec','S_rec']:
                self.all_para_dict[i+'_err'] = -9999 
        else:
            self.z = float(self.specdata.header['Z'])-float(self.specdata.header['Z_ERR'])
            self.calcSindex()
            for i in ['R_mean','V_mean','H_mean_tri','K_mean_tri','S_tri','S_MWL','H_mean_rec','K_mean_rec','S_rec']:
                self.all_para_dict[i+'_err'] = abs(all_para_dict[i] - self.all_para_dict[i])
        
        S_info_err = [self.all_para_dict[i] for i in self.para_catalog[1::2]]
        para_dict_err = dict(zip(self.para_catalog[1::2],S_info_err))
        return para_dict_err

    def calcErr_from_z(self):
        for i in self.para_catalog[::2]:
            if self.all_para_dict[i] == 'unknown':
                error = ('Can not calculate error of {} without the value of {}!'.format(i,i)
                    +' You should run calcSindex() before calculating the error.')
                raise KeyError(error)
            
        all_para_dict = self.all_para_dict.copy()
        zplus = self.__calcErr_from_zplus(all_para_dict)
        zminus = self.__calcErr_from_zminus(all_para_dict)
        self.all_para_dict = all_para_dict.copy()
        for i in ['R_mean','V_mean','H_mean_tri','K_mean_tri','S_tri','S_MWL','H_mean_rec','K_mean_rec','S_rec']:
            self.all_para_dict[i+'_err'] = (zplus[i+'_err'] + zminus[i+'_err'])/2
        
        S_info_err = [self.all_para_dict[i] for i in self.para_catalog[1::2]]
        para_dict_err = dict(zip(self.para_catalog[1::2],S_info_err))
        return para_dict_err
        
        
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

    def __calc_band_err(self,bin_center,step,band_width,tri=False):
        err_list = self.f_err(bin_center)
        bin_list = self.f_bin(bin_center)
        if min(err_list)<=0:
            rec_err = -9999
        else:
            if tri:
                tem_err = [err_list[i]**2/step/bin_list[i]*self.trig[i]**2 for i in range(len(err_list))]
            else:
                tem_err = [err_list[i]**2/step/bin_list[i] for i in range(len(err_list))]
            rec_err = math.sqrt(sum(tem_err))*step/band_width
        return rec_err
            
    def calc_R_err(self):
        R_mean_err = self.__calc_band_err(bin_center=self.band_R,step=self.step1,band_width=20)
        self.all_para_dict['R_mean_err'] = R_mean_err
        return R_mean_err

    def calc_V_err(self):
        V_mean_err = self.__calc_band_err(bin_center=self.band_V,step=self.step1,band_width=20)
        self.all_para_dict['V_mean_err'] = V_mean_err
        return V_mean_err

    def calc_H_rec_err(self):
        H_rec_mean_err = self.__calc_band_err(bin_center=self.band_H_rec,step=self.step2,band_width=1)
        self.all_para_dict['H_mean_rec_err'] = H_rec_mean_err
        return H_rec_mean_err

    def calc_K_rec_err(self):
        K_rec_mean_err = self.__calc_band_err(bin_center=self.band_K_rec,step=self.step2,band_width=1)
        self.all_para_dict['K_mean_rec_err'] = K_rec_mean_err
        return K_rec_mean_err

    def calc_H_tri_err(self):
        H_tri_mean_err = self.__calc_band_err(bin_center=self.band_H_tri,step=self.step2,band_width=self.FWHM,tri = True)
        self.all_para_dict['H_mean_tri_err'] = H_tri_mean_err
        return H_tri_mean_err

    def calc_K_tri_err(self):
        K_tri_mean_err = self.__calc_band_err(bin_center=self.band_K_tri,step=self.step2,band_width=self.FWHM,tri = True)
        self.all_para_dict['K_mean_tri_err'] = K_tri_mean_err
        return K_tri_mean_err

    def __calc_S_tri_err(self,R_mean_err,V_mean_err,H_mean_tri_err,K_mean_tri_err):
        if -9999 not in [R_mean_err,V_mean_err,H_mean_tri_err,K_mean_tri_err]:
            RV_err = math.sqrt(R_mean_err**2+V_mean_err**2)
            RV_err = RV_err/(self.all_para_dict['R_mean']+self.all_para_dict['V_mean'])
            HK_tri_err = math.sqrt(H_mean_tri_err**2+K_mean_tri_err**2)
            HK_tri_err = HK_tri_err/(self.all_para_dict['H_mean_tri']+self.all_para_dict['K_mean_tri'])
            S_tri_err = self.all_para_dict['S_tri']*math.sqrt(RV_err**2+HK_tri_err**2)
            S_MWL_err = S_tri_err*8*self.alpha*self.FWHM/20.
        else:
            S_tri_err=-9999
            S_MWL_err=-9999

        return S_tri_err,S_MWL_err

    def __calc_S_rec_err(self,R_mean_err,V_mean_err,H_mean_rec_err,K_mean_rec_err):
        if -9999 not in [R_mean_err,V_mean_err,H_mean_rec_err,K_mean_rec_err]:
            RV_err = math.sqrt(R_mean_err**2+V_mean_err**2)
            RV_err = RV_err/(self.all_para_dict['R_mean']+self.all_para_dict['V_mean'])
            HK_rec_err = math.sqrt(H_mean_rec_err**2+K_mean_rec_err**2)
            HK_rec_err = HK_rec_err/(self.all_para_dict['H_mean_rec']+self.all_para_dict['K_mean_rec'])
            S_rec_err = self.all_para_dict['S_rec']*math.sqrt(RV_err**2+HK_rec_err**2)
        else:
            S_rec_err=-9999
        return S_rec_err

    def calcErr_from_flux(self):
        for i in self.para_catalog[::2]:
            if self.all_para_dict[i] == 'unknown':
                error = ('Can not calculate error of {} without the value of {}!'.format(i,i)
                    +' You should run calcSindex() before calculating the error.')
                raise KeyError(error)
        
        original_error = self.__getOrigError()
        R_mean_err = self.calc_R_err()
        V_mean_err = self.calc_V_err()
        H_mean_rec_err = self.calc_H_rec_err()
        K_mean_rec_err = self.calc_K_rec_err()
        H_mean_tri_err = self.calc_H_tri_err()
        K_mean_tri_err = self.calc_K_tri_err()
        
        S_rec_err = self.__calc_S_rec_err(R_mean_err,V_mean_err,H_mean_rec_err,K_mean_rec_err)
        S_tri_err,S_MWL_err = self.__calc_S_tri_err(R_mean_err,V_mean_err,H_mean_tri_err,K_mean_tri_err)

        self.all_para_dict['S_rec_err'] = S_rec_err
        self.all_para_dict['S_tri_err'] = S_tri_err
        self.all_para_dict['S_MWL_err'] = S_MWL_err
        
        S_info_err = [self.all_para_dict[i] for i in self.para_catalog[1::2]]
        para_dict_err = dict(zip(self.para_catalog[1::2],S_info_err))
        return para_dict_err

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
        z_err = self.calcErr_from_z()
        flux_err = self.calcErr_from_flux()
            
        for i in ['R_mean','V_mean','H_mean_tri','K_mean_tri','S_tri','S_MWL','H_mean_rec','K_mean_rec','S_rec']:
            if z_err[i+'_err']<0 or flux_err[i+'_err']<0:
                tem_err = -9999
            else:
                tem_err = (z_err[i+'_err']**2 + flux_err[i+'_err']**2)
                tem_err = tem_err**0.5
            self.all_para_dict[i+'_err'] = tem_err
        
        S_info_err = [self.all_para_dict[i] for i in self.para_catalog[1::2]]
        para_dict_err = dict(zip(self.para_catalog[1::2],S_info_err))
        return para_dict_err        
        
    def calc(self,calc_para=False,accuracy=6,datasav=False,fitsinfo = ['OBSID','FILENAME'],header=True):
        if calc_para == False:
            self.calcError()
            for key in self.all_para_dict:
                self.all_para_dict[key] = round(self.all_para_dict[key],accuracy)
        if datasav:
            self.recorddata(fitsinfo = fitsinfo,header=header)
        return self.all_para_dict

    def __get_fitsinfo(self,fitsinfo=['OBSID','FILENAME']):
        info = []
        info_data=[]
        for i in fitsinfo:
            try:
                info_data.append(self.specdata.header[i])
                if i=='FILENAME':
                    info.append('fitsname')
                else:
                    info.append(i)
            except:
                print('Warning:',i,'not in .fits file!')
        return dict(zip(info,info_data))
    
    def __recorddata(self,fitsinfo = ['OBSID','FILENAME'],header=True):
        fits_dict = self.__get_fitsinfo(fitsinfo =fitsinfo)
        
        record_dict = dict(**fits_dict,**self.all_para_dict)
        record_catalog = fitsinfo + list(self.para_catalog)
        record_info = [str(record_dict[i]) for i in record_dict]
        f=open(self.Sindex_savepath,'w')
        if header:
            f.write(','.join(record_catalog)+'\n')
        f.write(','.join(record_info)+'\n')
        f.close()
        return record_dict
    
    def __data_plot(self,fig2=True):
        plot_flux = self.flux_func(self.all_band)
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.figure(figsize=(14,7))
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        ax2.set_ylim(0,max(plot_flux)*1.2)
        if self.z == 0:
            fig1 = ax1.plot(self.all_band,plot_flux/max(plot_flux),label='original',linewidth='1.5',color='red')
            ax1.legend(handles=fig1,loc='upper left',fontsize='15')
        elif fig2==False:
            fig1 = ax1.plot(self.all_band,plot_flux/max(plot_flux),color='red',label='shifted',linewidth='1.5')        
            ax1.legend(handles=fig1,loc='upper left',fontsize='15')
        else:
            fig1 = ax1.plot(self.all_band,plot_flux/max(plot_flux),color='red',label='shifted',linewidth='1.5')
            fig2 = ax2.plot(self.wavelen,self.flux,label='original',linewidth='1.5',linestyle='-.',color='blue')
            ax2.legend(handles=fig1+fig2,loc='upper left',fontsize='15')
        
        ax1.plot([self.L_V-10,self.L_R+10],[1,1],linestyle='--',color='black')

        tri_color = 'lime'
        ax1.fill_between([self.L_H-self.FWHM,self.L_H,self.L_H+self.FWHM],[0,0,0],[0,1,0],color=tri_color)
        ax1.fill_between([self.L_K-self.FWHM,self.L_K,self.L_K+self.FWHM],[0,0,0],[0,1,0],color=tri_color)
        rec_color = 'yellow'
        ax1.fill_between([self.L_R-10,self.L_R+10],  [1,1],[0,0],color=rec_color)
        ax1.fill_between([self.L_V-10,self.L_V+10],  [1,1],[0,0],color=rec_color)
        ax1.fill_between([self.L_H-1/2,self.L_H+1/2],[1,1],[0,0],color=rec_color)
        ax1.fill_between([self.L_K-1/2,self.L_K+1/2],[1,1],[0,0],color=rec_color)

        ax1.plot([self.L_H,self.L_H],[0,1],color='black',linestyle='--')
        ax1.plot([self.L_K,self.L_K],[0,1],color='black',linestyle='--')

        self.fig1 = fig1
        self.fig2 = fig2
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
        ax1.tick_params(axis='x',which='major',length=8,width=1.5,labelsize=15,pad=8)
        ax1.tick_params(axis='y',which='major',length=8,width=1.5,labelsize=15)
        ax1.tick_params(which='minor',length=4,width=1)
        ax1.set_xlabel(r'$\rm{Vacuum\ Wavelength\ (\AA)}$',fontsize=15)
        ax1.set_ylabel('Relative Flux',fontsize=15)
        ax1.set_xlim(self.L_V-10,self.L_R+10)
        ax1.set_ylim(0,1.2)
        ax2.set_ylabel('Flux',fontsize=15)
        ax2.tick_params(which='major',length=8,width=1.5,labelsize=15)
        return ax1,ax2

    def __plotS(self,ax1,ax2):
        ax1.text(self.L_K-4,1.03,'H_mean_rec={:.3f}\nK_mean_rec={:.3f}\nS_rec={:.3f}'.format(self.all_para_dict['H_mean_rec'],
                self.all_para_dict['K_mean_rec'],self.all_para_dict['S_rec']),fontsize=15)
        ax1.text(self.L_H-10,1.03,'H_mean_tri={:.3f}\nK_mean_tri={:.3f}\nS_tri={:.3f}'.format(self.all_para_dict['H_mean_tri'],
                self.all_para_dict['K_mean_tri'],self.all_para_dict['S_tri']),fontsize=15)
        ax1.text(self.L_R-15,1.03,'R_mean={:.3f}\nV_mean={:.3f}\nS_MWL={:.3f}'.format(self.all_para_dict['R_mean'],
                self.all_para_dict['V_mean'],self.all_para_dict['S_MWL']),fontsize=15)

        rv_str = '{:.2f}'.format(float(self.z)*299792.458)
        rv_str = '$'+rv_str+'$'
        ax1.text(self.L_V+9,1.03,'Radial Velocity\n='+r'$\,$'+rv_str+r'$\,$'+'km/s\n', fontsize=15)
        return ax1,ax2
        
    def __plotInfo(self,ax1,ax2):
        ax1,ax2 = self.__plotS(ax1,ax2)
        
        ax1.text(self.L_V-5,0.2,r'$V$ band',fontsize=15)
        ax1.text(self.L_R-5,0.2,r'$R$ band',fontsize=15)
        ax1.annotate("Ca II H "+str(self.L_H)+r'$\,$'+u'\u00C5',xy=(self.L_H,0.1),xytext=(self.L_H+4,0.15),
                color="blue",weight="bold",fontsize=11,arrowprops=dict(arrowstyle="->",color="blue"))
        ax1.annotate("Ca II K "+str(self.L_K)+r'$\,$'+u'\u00C5',xy=(self.L_K,0.1),xytext=(self.L_K+4,0.15),
                color="blue",weight="bold",fontsize=11,arrowprops=dict(arrowstyle="->",color="blue"))

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
        if (feh_str != 'unknown'):
            feh_str = '$'+feh_str+'$'
        snrg_str = '{}'.format(self.specdata.header['SNRG'])
        snrr_str = '{}'.format(self.specdata.header['SNRR'])
        title_info = r'$T_\mathrm{eff}$=' + teff_str + r'   $\log\,g$=' + logg_str + '   [Fe/H]='+ feh_str + r'   $\mathrm{SNR}_g$=' + snrg_str + r'   $\mathrm{SNR}_r$=' + snrr_str
        #ax1.text((self.L_R + self.L_V)/2.0, 1.22, title_info, fontsize=15, horizontalalignment='center')
        ax1.set_title('Obsid-{} {} ({})\n{}'.format(self.specdata.header['OBSID'],self.specdata.header['FILENAME'],self.specdata.header['DATA_V'],title_info),fontsize=15)
        return ax1,ax2

    def plotSpectrum(self, figure_savepath=None, stellar_parameters=None, figsav=True, figshow=True):
        '''
        
        Generate spectrum diagram.
        
        Optional arguments:
            figure_savepath:        set figure file save path
            stellar_parameters:     set stellar parameters [teff (K), logg (dex), feh (dex)]
            figshow:                show spectrum diagram on screen (True/False, default True)
        
        '''
        if stellar_parameters:
            self.stellar_parameters = stellar_parameters
        calced_info = self.calc()
        ax1,ax2 = self.__data_plot()
        ax1,ax2 = self.__plotset(ax1,ax2)
        self.ax1,self.ax2 = self.__plotInfo(ax1,ax2)
        if figure_savepath:
            self.figure_savepath = figure_savepath
        if figsav:
            plt.savefig(self.figure_savepath,bbox_inches='tight')
        if figshow:
            plt.show()
        plt.close()
        return calced_info
 
    def saveRecord(self, figure_savepath=None, stellar_parameters=None, figsav=True, figshow=False,
                   Sindex_savepath=None, fitsinfo=['OBSID','FILENAME'], data_header=True, datasav=True):
        '''
        
        Save stellar activity parameters to a csv file and spectrum diagram to a figure file.
            
        Optional arguments:
            Sindex_savepath:        set csv file save path
            figure_savepath:        set figure file save path
            stellar_parameters:     set stellar parameters [teff (K), logg (dex), feh (dex)]
            csv_header:             write out column names in csv file (True/False, default True)
            figshow:                show spectrum diagram on screen (True/False, default False)
        
        '''
        calced_info = self.plotSpectrum(figure_savepath=figure_savepath, stellar_parameters=stellar_parameters, figsav=figsav, figshow=figshow)
        if Sindex_savepath:
            self.Sindex_savepath = Sindex_savepath
        if datasav:
            record_dict = self.__recorddata(fitsinfo=fitsinfo, header=data_header)

        return record_dict

def main():
    print('Start cals\nVesion 1.0.0')
    parser = argparse.ArgumentParser(description='description\n')
    parser.add_argument('--cs','-c',help = 'Calculate the S by setting file path. | e.g. cals -c example.fits',required=True)
    parser.add_argument('--savepath','-sp',help = 'Setting save path.',default='')
    parser.add_argument('--figsav','-fs',help = 'Save/not save figure by setting figsav as 1/0, default=1.',default=1,choices=[0,1],type=int)
    parser.add_argument('--figdisplay','-fd',help = 'Display/not display figure by setting figshow as 1/0, default=0.',default=0,choices=[0,1],type=int)
    parser.add_argument('--datasav','-ds',help = 'Save/not save data by setting datasav as 1/0, default=1.',default=1,choices=[0,1],type=int)
    parser.add_argument('--dataprint','-dp',help = 'Display/not display data by setting dataprint as 1/0, default=1.',default=1,choices=[0,1],type=int)
    args = parser.parse_args()

    if args.cs:
        CaIISindex.save_path = args.savepath
        cs = CaIISindex(args.cs)
        para_dict = cs.saveRecord(figsav=args.figsav,figshow=args.figdisplay,datasav=args.datasav)
    if args.dataprint:
        print('*'*18+'Result'+'*'*18)
        for key in para_dict:
            print('{:^20}: {}'.format(key, para_dict[key]))
        print('*'*19+'End'+'*'*19)

if __name__ == '__main__':

    data_path = "example.fits"
    cs = CaIISindex(data_path)
    #cs.calcSindex()
    
    '''
    cs.calcErr_from_z()
    print(cs.all_para_dict)
    cs.calcErr_from_flux()
    print(cs.all_para_dict)
    cs.calcError()
    print(cs.all_para_dict)
    '''

    cs.stellar_parameters = ['5534.22','4.423','-0.025']
    print(cs.saveRecord())

    


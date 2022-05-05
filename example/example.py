#!/usr/bin/python

from cals import CaIISindex

if __name__ == '__main__':
    #help(cals)
    #help(cals.CaIISindex)
    data_path = "example.fits"
    cs = CaIISindex(data_path)
    #cs.Sindex_savepath = 'S_index_example.csv'               # set csv file save path
    #cs.figure_savepath = 'example.png'                       # set figure file save path
    cs.stellar_parameters = ['5534.22', '4.423', '-0.025']    # set stellar parameters [teff (K), logg (dex), feh (dex)]
    cs.saveRecord()                          # save calculated parameters to a csv file and spectrum diagram to a figure file

    # print(cs.calcSindex())                 # return a dict of stellar activity parameters
    # print(cs.calcError())                  # return a dict of parameter errors
    # cs.plotSpectrum()                      # plot spectrum diagram
    
    

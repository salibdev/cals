#!/usr/bin/python

from cals import CaIISindex

if __name__ == '__main__':
    #help(cals)
    #help(cals.CaIISindex)
    data_path = "example.fits"
    example = CaIISindex(data_path)
    #example.Sindex_savepath = 'S_index_example.csv'               # set csv file save path
    #example.figure_savepath = 'example.png'                       # set figure file save path
    example.stellar_parameters = ['5534.22', '4.423', '-0.025']    # set stellar parameters [teff (K), logg (dex), feh (dex)]
    example.saveRecord()                          # save calculated parameters to a csv file and spectrum diagram to a figure file

    # print(example.calcSindex())                 # return a dict of stellar activity parameters
    # print(example.calcError())                  # return a dict of parameter errors
    # example.plotSpectrum()                      # plot spectrum diagram
    
    

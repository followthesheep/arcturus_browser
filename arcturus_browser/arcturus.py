import pylab as plt
import numpy as np
import pandas as pd
import os
import overlap
from specutils import plotlines
import seaborn as sns

def plot_spectrum(wave_range = [2.15,2.2]):
    # plot the Arcturus spectrum

    datadir = os.path.join(os.path.dirname(__file__),'data')
    # plot all the data into a PDF form


    # get the line list
    atomic_file = os.path.join(datadir,'arcturus_atomic_lines.txt')
    molecular_file = os.path.join(datadir,'arcturus_molecular_lines.txt')    
    atomic_lines, atomic_line_names = np.loadtxt(atomic_file,delimiter=',',unpack=True,dtype=str)
    molecular_lines, molecular_line_names = np.loadtxt(molecular_file,delimiter=',',unpack=True,dtype=str)    
    atomic_lines = np.array(atomic_lines,dtype=float)-.00001
    molecular_lines = np.array(molecular_lines,dtype=float)-0.00001
    molecular_line_names = ['$'+lamb+'$' for lamb in molecular_line_names]
    atomic_line_names = ['$'+lamb+'$' for lamb in atomic_line_names]

    sns.set_context('paper',font_scale=1.0, rc={"lines.linewidth": 0.8})
    sns.set_style('white')
    sns.set_style('ticks')
    
    fontscale = 0.75    
    wave, flux = spectrum(wave_range = wave_range)
    plt.clf()
    plt.plot(wave,flux)
    plt.ylim(0.2,1.1)
    plt.xlim(np.min(wave),np.max(wave))
    plotlines.oplotlines(label=True,lines=molecular_lines,line_names=molecular_line_names,size=12*fontscale,alpha=0.6)
    plotlines.oplotlines(label=True,lines=atomic_lines,line_names=atomic_line_names,size=12*fontscale,color='r',alpha=0.9,linestyle=':')
    
    

def spectrum(wave_range = None, linear = False,winter=True,debug=False):
    # Returns the arcturus spectrum
    #
    # Keywords
    # --------
    # wave_range : wavelength range to return in microns (default: None, returns all available data)
    #
    # linear : By default, the spectrum is returned at the sampling of
    # the original atlas. Set this keyword to returns linearly sampled
    # points using average spacing within wave_range. This will be a
    # problem for parts of the spectrum that is missing. (default:
    # False)
    #
    # winter : By default will return the winter spectrum. Set this to
    # False to return the summer spectrum
    #
    
    datafile = os.path.join(os.path.dirname(__file__),'data/arcturus_combined_unique.npz')
    if os.path.exists(datafile):
        npfiles = np.load(datafile)

        if winter:
            wave_corr = npfiles['wave_num']*.9999129
            flux = npfiles['rat_w']
        else:
            wave_corr = npfiles['wave_num']*1.0000478
            flux = npfiles['rat_s']
    
        # convert to microns
        wave_corr = (1.0/wave_corr)*1e4

        if wave_range is None:
            wave_range = [np.min(wave_corr),np.max(wave_corr)]

        good = np.where((wave_corr >= wave_range[0]) & (wave_corr <= wave_range[1]))[0]

        wave_corr = wave_corr[good]
        flux = flux[good]

        # sort by increasing wavelength
        sort_ind = np.argsort(wave_corr)
        wave_corr = wave_corr[sort_ind]
        flux = flux[sort_ind]

        if debug:
            plt.clf()
            plt.plot(wave_corr,flux)

        # resample linearly
        if linear:
            mean_delt = np.mean(wave_corr[1:]-wave_corr[:-1])
            nsamples = np.fix((np.max(wave_corr) - np.min(wave_corr))/mean_delt)
            wave_samples = np.linspace(np.min(wave_corr),np.max(wave_corr),num=nsamples)

            flux_samples = np.interp(wave_samples,wave_corr,flux)

            if debug:
                plt.plot(wave_samples,flux_samples)
            return wave_samples,flux_samples
        else:
            return wave_corr, flux
    else:
        print 'file not found:'+datafile
        return 0,0


def test_edges():
    datadir = os.path.join(os.path.dirname(__file__),'data/')

    file1 = os.path.join(datadir,'ab4275.')
    file2 = os.path.join(datadir,'ab4250.')

    wave_num, obs_s, tel_s, rat_s, obs_w, tel_w, rat_w = np.loadtxt(file1,unpack=True)
    wave_num2, obs_s2, tel_s2, rat_s2, obs_w2, tel_w2, rat_w2 = np.loadtxt(file2,unpack=True)
    
    summer_wave = wave_num*1.0000478
    summer_wave2 = wave_num2*1.0000478    

    winter_wave = wave_num*.9999129
    winter_wave2 = wave_num2*.9999129    
    
    # convert from wavenumber (1/cm) into microns
    summer_wave = (1.0/summer_wave)*1e4
    summer_wave2 = (1.0/summer_wave2)*1e4        
    winter_wave = (1.0/winter_wave)*1e4
    winter_wave2 = (1.0/winter_wave2)*1e4        

    plt.clf()
    #plt.plot(winter_wave,rat_w)
    #plt.plot(winter_wave2,rat_w2)


    # try to make a spectrum that has unique wavelength region
    non_overlap = np.in1d(wave_num,wave_num2,invert=True)
    ind1, ind2 = overlap.overlap(wave_num,wave_num2)
    print ind1,ind2
    plt.plot(winter_wave[ind1],rat_w[ind1])
    plt.plot(winter_wave2[ind2],rat_w2[ind2])    

    # average the two pieces
    rat_w2[ind2] = (rat_w2[ind2]+rat_w[ind1])/2.0
    
    combine_wave = np.append(winter_wave2,winter_wave[non_overlap])
    combine_flux = np.append(rat_w2,rat_w[non_overlap])
    

    plt.plot(combine_wave,combine_flux-0.2)
    

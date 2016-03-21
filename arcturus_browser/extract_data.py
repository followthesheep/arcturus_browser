import pylab as plt
import numpy as np
import os
from specutils import plotlines
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

def extract_data():
    # convert the data in text file into a different form that is easier to load
    
    pass

def plot_orig_data(datafile = None,micron=True):
    # plot the original text files
    if datafile is None:
        #datafile = os.path.join(os.path.dirname(__file__),'data/ab4475.')
        datafile = os.path.join(os.path.dirname(__file__),'data/ab10900.')
    wave_num, obs_s, tel_s, rat_s, obs_w, tel_w, rat_w = np.loadtxt(datafile,unpack=True)

    summer_wave = wave_num*1.0000478

    winter_wave = wave_num*.9999129
    if micron:
        # convert from wavenumber (1/cm) into microns
        summer_wave = (1.0/summer_wave)*1e4
        winter_wave = (1.0/winter_wave)*1e4        
    plt.clf()
    plt.plot(summer_wave,rat_s)
    plt.plot(winter_wave,rat_w)
    atomic_file = os.path.join(os.path.dirname(__file__),'data/arcturus_atomic_lines.txt')
    molecular_file = os.path.join(os.path.dirname(__file__),'data/arcturus_molecular_lines.txt')    
    atomic_lines, atomic_line_names = np.loadtxt(atomic_file,delimiter=',',unpack=True,dtype=str)
    molecular_lines, molecular_line_names = np.loadtxt(molecular_file,delimiter=',',unpack=True,dtype=str)    
    atomic_lines = np.array(atomic_lines,dtype=float)-.00001
    molecular_lines = np.array(molecular_lines,dtype=float)-0.00001
    molecular_lines = np.array(molecular_lines,dtype=float)
    molecular_line_names = ['$'+lamb+'$' for lamb in molecular_line_names]
    atomic_line_names = ['$'+lamb+'$' for lamb in atomic_line_names]
    plt.ylim(0,1.1)
    plt.xlim(np.min(winter_wave),np.max(winter_wave))    
    #plotlines.oplotlines(spec_wave=winter_wave,spec_flux=rat_w,lines=atomic_lines,line_names=atomic_line_names)
    plotlines.oplotlines(lines=atomic_lines,line_names=atomic_line_names,label=True)
    plotlines.oplotlines(lines=molecular_lines,line_names=molecular_line_names,label=True)        
    #plotlines.oplotlines(spec_wave=winter_wave,spec_flux=rat_w,lines=molecular_lines,line_names=molecular_line_names)

    plt.tight_layout()

def plot_data_pdf(filename = 'Arcturus_IR_spectrum.pdf',per_page = 2,micron=True,plot_both=False):
    datadir = os.path.join(os.path.dirname(__file__),'data')
    # plot all the data into a PDF form
    pdf_pages = PdfPages(filename)

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
    fontscale = 0.75


    # spectra files
    specfiles = os.path.join(datadir,'filenames.txt')
    specnames = np.loadtxt(specfiles,dtype=str,unpack=True)
    

    for i in xrange(len(specnames)):
        datafile = os.path.join(datadir,specnames[i])
        wave_num, obs_s, tel_s, rat_s, obs_w, tel_w, rat_w = np.loadtxt(datafile,unpack=True)        
        
        summer_wave = wave_num*1.0000478
        
        winter_wave = wave_num*.9999129
        if micron:
            # convert from wavenumber (1/cm) into microns
            summer_wave = (1.0/summer_wave)*1e4
            winter_wave = (1.0/winter_wave)*1e4        


        if ((i % per_page ) == 0):
            fig = plt.figure(figsize=(11,8.5),dpi=100)
            plt.clf()
        plt.subplot(per_page,1,(i % per_page)+1)
        plt.ylim(0,1.1)
        plt.xlim(np.min(winter_wave),np.max(winter_wave))
        

        plotlines.oplotlines(label=True,lines=molecular_lines,line_names=molecular_line_names,size=12*fontscale)
        plotlines.oplotlines(label=True,lines=atomic_lines,line_names=atomic_line_names,size=12*fontscale,color='b')

        if plot_both:
            plt.plot(summer_wave,rat_s)
            plt.plot(winter_wave,rat_w)
        else:
            plt.plot(winter_wave,rat_w,'k')
        #plotlines.oplotlines(spec_wave=summer_wave,spec_flux=rat_s,lines=atomic_lines,line_names=atomic_line_names,size=12*fontscale)
        #plotlines.oplotlines(spec_wave=summer_wave,spec_flux=rat_s,lines=molecular_lines,line_names=molecular_line_names,size=12*fontscale)
        plt.xlabel(r'$Wavelength (\mu m)$')
        plt.ylabel(r'$Flux$')

        if ((i != 0) and (((i+1) % per_page) == 0)):
            pdf_pages.savefig(fig,orientation='landscape')
            plt.close(fig)

    pdf_pages.close() 
    plt.close('all')
            

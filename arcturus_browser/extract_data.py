import pylab as plt
import numpy as np
import os
from specutils import plotlines
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import overlap

def extract_data():
    # convert the data in text file into a different form that is easier to load
    
    pass

def clip(arr,n=4):
    # clip some points from the beginning and end of the spectra
    # remove n points from the beginning and end
    return arr[n:-n]
    
def combine_data_unique(debug=False,save=False):
    # combine all the abs* files into one without overlap in wavelength sampling
    
    datadir = os.path.join(os.path.dirname(__file__),'data/')
    savefile = os.path.join(datadir,'arcturus_combined_unique')
    specfiles = os.path.join(datadir,'filenames.txt')
    specnames = np.loadtxt(specfiles,dtype=str,unpack=True)

    file1 = os.path.join(datadir,specnames[0])
    wave_num, obs_s, tel_s, rat_s, obs_w, tel_w, rat_w = np.loadtxt(file1,unpack=True)
    ## wave_num = clip(wave_num)
    ## obs_s = clip(obs_s)
    ## tel_s = clip(tel_s)
    ## rat_s = clip(rat_s)
    ## obs_w = clip(obs_w)
    ## tel_w = clip(tel_w)
    ## rat_w = clip(rat_w)
    wave_micron = (1.0/wave_num)*1e4
    if debug:
        plt.clf()
        plt.plot(wave_micron,rat_s)
        
    for i in xrange(1,len(specnames)):
        file2 = os.path.join(datadir,specnames[i])
        wave_num2, obs_s2, tel_s2, rat_s2, obs_w2, tel_w2, rat_w2 = np.loadtxt(file2,unpack=True)
        ## wave_num2 = clip(wave_num2)
        ## obs_s2 = clip(obs_s2)
        ## tel_s2 = clip(tel_s2)
        ## rat_s2 = clip(rat_s2)
        ## obs_w2 = clip(obs_w2)
        ## tel_w2 = clip(tel_w2)
        ## rat_w2 = clip(rat_w2)
        
        wave_micron2 = (1.0/wave_num2)*1e4

        if debug:
            plt.plot(wave_micron2,rat_s2)            
        non_overlap = np.in1d(wave_num,wave_num2,invert=True)

        # find the indices in common between the two
        ind1, ind2 = overlap.overlap(wave_num,wave_num2)

        # average the two pieces
        if len(ind1) > 0:
            obs_s2[ind2] = (obs_s2[ind2]+obs_s[ind1])/2.0
            tel_s2[ind2] = (tel_s2[ind2]+tel_s[ind1])/2.0
            rat_s2[ind2] = (rat_s2[ind2]+rat_s[ind1])/2.0                
            obs_s2[ind2] = (obs_s2[ind2]+obs_s[ind1])/2.0
            obs_w2[ind2] = (obs_w2[ind2]+obs_w[ind1])/2.0
            tel_w2[ind2] = (tel_w2[ind2]+tel_w[ind1])/2.0
            rat_w2[ind2] = (rat_w2[ind2]+rat_w[ind1])/2.0

        
        wave_num = np.append(wave_num2,wave_num[non_overlap])
        obs_s = np.append(obs_s2,obs_s[non_overlap])
        rat_s = np.append(rat_s2,rat_s[non_overlap])
        tel_s = np.append(tel_s2,tel_s[non_overlap])

        obs_w = np.append(obs_w2,obs_w[non_overlap])        
        tel_w = np.append(tel_w2,tel_w[non_overlap])
        rat_w = np.append(rat_w2,rat_w[non_overlap])


    if save:
        print 'saving: '+savefile+'.txt'
        output = open(savefile+'.txt','w')
        for j in xrange(len(wave_num)):
            output.write('%f %f %f %f %f %f %f\n' % (wave_num[j],obs_s[j], tel_s[j], rat_s[j], obs_w[j], tel_w[j], rat_w[j]))
        output.close()
        print 'saving: '+savefile+'.npz'
        np.savez(savefile,wave_num=wave_num,obs_s=obs_s,tel_s=tel_s,rat_s=rat_s,
                 obs_w= obs_w,tel_w=tel_w,rat_w=rat_w)
    if debug:
        plt.plot((1.0/wave_num)*1e4,rat_s)
        plt.ylim(0,1.8)

    
        
def plot_orig_data(datafile = None,micron=True,alpha=0.8):
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

    #plt.plot(summer_wave,rat_s)
    plt.plot(winter_wave,rat_w,label='Arcturus')
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
    plotlines.oplotlines(lines=atomic_lines,line_names=atomic_line_names,label=True,color='r',linestyle=':')
    plotlines.oplotlines(lines=molecular_lines,line_names=molecular_line_names,label=True,alpha=alpha)        
    #plotlines.oplotlines(spec_wave=winter_wave,spec_flux=rat_w,lines=molecular_lines,line_names=molecular_line_names)



def plot_data_pdf(filename = 'Arcturus_IR_spectrum.pdf',per_page = 2,micron=True,plot_both=False,label_file=False):
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
    sns.set_style('ticks')
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
        

        plotlines.oplotlines(label=True,lines=molecular_lines,line_names=molecular_line_names,size=12*fontscale,alpha=0.6)
        plotlines.oplotlines(label=True,lines=atomic_lines,line_names=atomic_line_names,size=12*fontscale,color='r',alpha=0.9,linestyle=':')

        if plot_both:
            plt.plot(summer_wave,rat_s)
            plt.plot(winter_wave,rat_w)
        else:
            plt.plot(winter_wave,rat_w,'k')
        #plotlines.oplotlines(spec_wave=summer_wave,spec_flux=rat_s,lines=atomic_lines,line_names=atomic_line_names,size=12*fontscale)
        #plotlines.oplotlines(spec_wave=summer_wave,spec_flux=rat_s,lines=molecular_lines,line_names=molecular_line_names,size=12*fontscale)
        plt.xlabel(r'$Wavelength (\mu m)$')
        plt.ylabel(r'$Flux$')
        if label_file:
            plt.text(winter_wave[-2],1.05,specnames[i],size=9)
        if ((i != 0) and (((i+1) % per_page) == 0)):
            pdf_pages.savefig(fig,orientation='landscape')
            plt.close(fig)

    pdf_pages.close() 
    plt.close('all')
            

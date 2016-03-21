import numpy as np


def read_appendix():
    # read the appendix files with the line lists

    file1 = 'data/appendix.a'  # atomic lines
    file2 = 'data/appendix.b'  # molecular lines
    outfile1 = 'data/arcturus_atomic_lines.txt'
    outfile2 = 'data/arcturus_molecular_lines.txt'
    
    infile = open(file1,'r')

    #skip the first three lines
    lines = infile.readlines()

    print(lines[0].strip())
    print(lines[1].strip())
    print(lines[2].strip())

    elname = ''
    linenames = []
    linepos = []
        
    for i in xrange(3,len(lines)):
        inline = lines[i].strip()
        if inline != '':
            try:
                val = float(inline)
                linenames.append(elname)
                linepos.append(val)
            except ValueError:
                elname = inline
                

    output1 = open(outfile1,'w')
    for j in xrange(len(linenames)):
        output1.write('%f,%s\n' % (1.0/linepos[j]*1e4,linenames[j]))

                
    infile.close()


    infile2 = open(file2,'r')
    mol_lines = infile2.readlines()
    infile2.close()

    mol_wave = []
    mol_name_arr = []
    mol_name = ''
    for i in xrange(3,len(mol_lines)):
        inline = mol_lines[i].split()
        if len(inline) > 0:
            try:
                val = float(inline[0])
                mol_name_arr.append(mol_name+' '+' '.join(inline[1:]))
                mol_wave.append(val)
            except ValueError:
                mol_name = ' '.join(inline)

    output2 = open(outfile2,'w')
    for j in xrange(len(mol_name_arr)):
       output2.write('%f,%s\n' % (1.0/mol_wave[j]*1e4,mol_name_arr[j]))


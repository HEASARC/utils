#!/usr/bin/env python3

# Program to create apec-style files by running a spex model at various temperature
# grid points. The lines and continuum are then gathered up, and stored in an
# apec-format table model for xspec

# Variables set in the next section control the temperature grid to use, the energy
# range and number of energy bins for the continuum and pseudo-continuum, and the
# value of epsilon above which lines will be stored individually in the line file
# and below which lines will be accumulated into the pseudo-continuum.

# This requires spex to be installed with the executable in the directory
# $SPEX90/bin/spex. It also requires astropy and numpy. It does not require heasoft.

# Based on Jeremy Sanders spex_to_xspec.py but reworked to use less memory
# kaa 1/20/2026

import os
import subprocess
import os.path
import re

import numpy as np
from astropy.io import fits

###############################################################################
# Adjustable parameters

# output root (using spex version) for apec format filenames
# creates outroot_(line|coco).fits
outroot = 'spex'

##############
# temperature grid: uncomment line and comment other to select

# old APEC temperature grid
#temperatures = np.logspace(np.log10(0.0008617385), np.log10(86.17385), 51)
# new APEC temperature grid from 3.0.9+ (see http://atomdb.org/interpolation/)
temperatures = np.logspace(np.log10(0.0008617385), np.log10(86.17385), 201)

# denser gridding between useful 0.01 and 100 keV
#temperatures = np.logspace(np.log10(0.01), np.log10(100), 201)

# for quick testing
#temperatures = np.array([1,2,4,8])
##############

# energy range and stepping to sample continuum (log spacing used)
contminenergy = 0.05
contmaxenergy = 50.
contenergysteps = 2048

# energy range and stepping to sample pseudo-continuum (log spacing used)
pcontminenergy = 0.05
pcontmaxenergy = 50.
pcontenergysteps = 4096

# Limit for storing lines separately to save space. Lines lower than
# this flux (photon cm^3/s) are put into a pseudo-continuum rather
# than stored separately.  The APEC default is 1e-20, but this
# produces many fewer lines using this for SPEX. Set to None to
# disable putting weak lines into a pseudo-continuum.
minepsilon = 1e-22
#minepsilon = None

# where to put output files (workdir by default)
tmpdir = os.environ.get('WORKDIR', 'workdir')

# end adjustable parameters
###############################################################################

# location of spex installation
try:
    spexroot = os.environ['SPEX90']
except KeyError:
    raise RuntimeError('Please set SPEX90 and initialize the SPEX environment')

# executable to use to run spex
spexexecutable = os.path.join(spexroot, 'bin/spex')

# for checking for numerical data
digits = set('0123456789')

# conversions
keV_K = 11.6048e6
keV_erg = 1.6022e-9

# convert from unit norm to cm3 in spex
norm_factor_cm3 = 1e58

# identify element names with element numbers in spex
elements = (
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O',
    'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
    'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
    'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
)

# create dict to convert element names into numbers
element_nums = {}
for num, element in enumerate(elements):
    element_nums[element] = num+1

# these are the apec elements (non hydrogen)
# (note this distinction is historical, as apec used fewer elements)
apec_elements = elements[1:]

# with hydrogen
all_elements = elements

# roman numerals (which come out of spex)
roman_numerals = (
    'I', 'II', 'III', 'IV', 'V', 'VI', 'VII',
    'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV',
    'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX',
    'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII',
    'XXVIII', 'XXIX', 'XXX'
)

# dict to convert numerals to numbers
roman_to_number = {}
for num, numeral in enumerate(roman_numerals):
    roman_to_number[numeral] = num+1

class Line:
    def __init__(self, element, ion, lowerlevel, upperlevel, wavelength, epsilon, energy):
        self.element = element
        self.ion = ion
        self.lowerlevel = lowerlevel
        self.upperlevel = upperlevel
        self.wavelength = wavelength
        self.epsilon = epsilon
        self.energy = energy

# end global definitions
###############################################################################
# utility routines

def deleteFile(f):
    """For debugging. Comment out to retain temporary files"""
    if os.path.exists(f) : os.unlink(f)

def cnvtNum(s):
    """Sometimes scientific notation values in spex output are expressed
    without an 'e', e.g.  5.62495-100, so we need to convert by hand.
    """
    try:
        v = float(s)
    except ValueError:
        m = re.match('([0-9.]+)([+-][0-9]{3})', s)
        if not m:
            raise ValueError
        else:
            a, b = m.groups()
            v = float(a)*10**float(b)
    return v

###############################################################################
# routines to initialize the output files

def initializeOutputFiles():
    """Create the output files and write their PARAMETERS extensions"""

    # construct HDU describing parameters for the continuum file
    col_kT = fits.Column(
        name='kT', format='1E', unit='keV',
        array=temperatures)
    col_EDensity = fits.Column(
        name='EDensity', format='1E', unit='cm**-3',
        array=np.zeros(len(temperatures))+1)
    col_NElement = fits.Column(
        name='NElement', format='1J',
        array=np.zeros(len(temperatures))+1)
    col_NCont = fits.Column(
        name='NCont', format='1J',
        array=np.zeros(len(temperatures))+1)
    col_NPseudo = fits.Column(
        name='NPseudo', format='1J',
        array=np.zeros(len(temperatures))+1)
    tabhdu = fits.BinTableHDU.from_columns([
        col_kT, col_EDensity, col_NElement, col_NCont, col_NPseudo
    ])
    tabhdu.name = 'PARAMETERS'

    # make output file containing the PARAMETERS extension
    hdulist = fits.HDUList([fits.PrimaryHDU(), tabhdu])
    cocofilename = outroot + '_coco.fits'
    hdulist.writeto(cocofilename, overwrite=True)

    
    # construct HDU describing parameters for the line file
    # note that col_kT and col_EDensity are the same between
    # continuum and line files
    col_Nelement = fits.Column(
        name='Nelement', format='1J',
        array=np.zeros(len(temperatures))+1)
    col_Nline = fits.Column(
        name='Nline', format='1J',
        array=np.zeros(len(temperatures))+1)
    tabhdu = fits.BinTableHDU.from_columns([
        col_kT, col_EDensity, col_Nelement, col_Nline])
    tabhdu.name = 'PARAMETERS'

    # make output line file containing the PARAMETERS extension
    hdulist = fits.HDUList([fits.PrimaryHDU(), tabhdu])
    linefilename = outroot + '_line.fits'
    hdulist.writeto(linefilename, overwrite=True)

    return (cocofilename, linefilename)

###############################################################################
# routines to run spex

def writeScript(fobj, T):
    """Write script sent to spex."""

    # write the header
    print('egrid log %e:%e %i' % (
        contminenergy, contmaxenergy,
        contenergysteps), file=fobj)
    # switch to the latest spex model
    print('var calc new', file=fobj)
    print('var newmekal all true', file=fobj)
    print('comp cie', file=fobj)
    print('abundance ag', file=fobj)

    # switch on all the elements
    for e in elements:
        print('par %02i val 1.0' % (element_nums[e]), file=fobj)

    # set temperature
    print('par t val %e' % T, file=fobj)

    # compute model
    print('calc', file=fobj)

    # dump out lines
    selfile = os.path.join(tmpdir, 'tmp_lines_T%010f.linesel' % T)
    fsel = open(selfile,'w')
    print('for flux 0', file=fsel)
    print('sel range %g:%g unit kev' % (contminenergy, contmaxenergy), file=fsel)
    outfile = os.path.join(tmpdir, 'tmp_lines_T%010f.fits' % T)
    print('ascdump fits %s 1 1 line key %s' % (outfile, selfile), file=fobj)

    # loop over all elements writing the continuum for that element
    for el in all_elements:
        print('ions mute all', file=fobj)
        print('ions unmute z %i' % element_nums[el], file=fobj)
        selfile = os.path.join(tmpdir, 'tmp_conti_T%010f_%s.consel' % (T, el))
        fsel = open(selfile,'w')
        print('sel range %g:%g unit kev' % (contminenergy, contmaxenergy), file=fsel)
        outfile = os.path.join(tmpdir, 'tmp_conti_T%010f_%s.fits' % (T, el))
        print('ascdump fits %s 1 1 tcl key %s' % (outfile, selfile), file=fobj)

    # end spex run
    print('quit', file=fobj)

def runSpex(T):
    """Run spex for temperature T and create temporary output files."""

    fname = os.path.join(tmpdir, 'tmp_spex_T%010f.script' % T)
    with open(fname, 'w') as fout:
        writeScript(fout, T)

    with open(fname) as fin:
        subprocess.call([spexexecutable], stdin=fin)

    return

###############################################################################
# routines to read the spex temporary output and make the line output file

def makeLineHDU(lines, T, totflux):
    """Given lines list, produce line HDU.

    lines is (element, ion, lowerlevel, upperlevel, wavelength, epsilon, energy) list."""

    # sort lines by element and ion and energy
    lines.sort(key=lambda x: (x.element, x.ion, 1/x.wavelength))

    # construct up FITS table to APEC format
    col_lambda = fits.Column(
        name='Lambda', format='1E', unit='A',
        array=[i.wavelength for i in lines])
# no errors at the moment so ignore
#    col_lambda_err = fits.Column(
#        name='Lambda_Err', format='1E', unit='A',
#        array=np.zeros( (len(lines),) ) + np.nan)
    col_epsilon = fits.Column(
        name='Epsilon', format='1E',
        unit='photons cm^3 s^-1',
        array=[v.epsilon for v in lines])
# no errors at the moment so ignore
#    col_epsilon_err = fits.Column(
#        name='Epsilon_Err', format='1E',
#        unit='photons cm^3 s^-1',
#        array=np.zeros((len(lines),)) + np.nan)
    col_element = fits.Column(
        name='Element', format='1J',
        array=[i.element for i in lines])
    col_ion = fits.Column(
        name='Ion', format='1J',
        array=[i.ion for i in lines])

    col_upperlev = fits.Column(
        name='UpperLev', format='1J',
        array=[i.upperlevel for i in lines])
    col_lowerlev = fits.Column(
        name='LowerLev', format='1J',
        array=[i.lowerlevel for i in lines])

# no errors at the moment so leave out those columns
    tabhdu = fits.BinTableHDU.from_columns([
        col_lambda, col_epsilon,
        col_element, col_ion,
        col_upperlev, col_lowerlev])
#    tabhdu = fits.BinTableHDU.from_columns([
#        col_lambda, col_lambda_err, col_epsilon,
#        col_epsilon_err, col_element, col_ion,
#        col_upperlev, col_lowerlev])

    tabhdu.name = 'EMISSIVITY'
    h = tabhdu.header
    h['HIERARCH TEMPERATURE'] = T*keV_K
    h['XTEMP'] = T

    # fixme wrong below (erg not photon)
    h['TOT_LINE'] = totflux
    h['N_LINES'] = len(lines)
    return tabhdu

def processStrongLines(T, linefilename):
    """Process the strong lines in the spex dumped spectra for temperature T."""

    totflux = 0.
    lines = []
    tempfile = os.path.join(tmpdir, 'tmp_lines_T%010f.fits' % T)
    ff = fits.open(tempfile)
    ftable = ff[1].data

    # skip lines out of energy range
    mask = ftable['ener'] >= contminenergy
    temptable = ftable[mask]
    mask = temptable['ener'] <= contmaxenergy
    goodtable = temptable[mask]

    # convert from total photon flux to normalised photon flux
    goodepsilon = goodtable['flux'] / norm_factor_cm3

    # split out the strong lines
    strong_mask = minepsilon is None or goodepsilon > minepsilon

    # make the HDU
    strongtable = goodtable[strong_mask]
    element = strongtable['iz']
    ion = strongtable['jz']
    lowerlevel = strongtable['il']
    upperlevel = strongtable['iu']
    wavelength = strongtable['wav']
    epsilon = goodepsilon[strong_mask]
    energy_kev = strongtable['ener']
    lines = [Line(*attributes) for attributes in zip(element,ion,lowerlevel,upperlevel,wavelength,epsilon,energy_kev)]
    totflux = sum(energy_kev*keV_erg*epsilon)

    print('T=%g, %i strong lines' % (T, len(lines)))
    tabhdu = makeLineHDU(lines, T, totflux)
    ff.close()

    # append the HDU
    fits.append(linefilename, tabhdu.data, tabhdu.header)

    return (len(elements), len(lines))
        
###############################################################################
# routines to read the spex temporary output and make the continuum output file

def readContinuumFile(filename):
    """Take spex dumped model spectrum file, and extract continuum."""

    ff = fits.open(filename)
    ftable = ff[1].data
    outenergy = ftable['ener']
    outval = ftable['scon'] * 1e44 / norm_factor_cm3

    return (outenergy, outval)

def readWeakLines(T):
    """Read the lines file for T and return the weak lines"""

    totflux = 0.
    lines = []
    tempfile = os.path.join(tmpdir, 'tmp_lines_T%010f.fits' % T)
    ff = fits.open(tempfile)
    ftable = ff[1].data

    # skip lines out of energy range
    mask = ftable['ener'] >= contminenergy
    temptable = ftable[mask]
    mask = temptable['ener'] <= contmaxenergy
    goodtable = temptable[mask]

    # convert from total photon flux to normalised photon flux
    goodepsilon = goodtable['flux'] / norm_factor_cm3

    # split out the weak lines
    weak_mask = minepsilon is not None and goodepsilon <= minepsilon

    weaktable = goodtable[weak_mask]
    element = weaktable['iz']
    ion = weaktable['jz']
    lowerlevel = weaktable['il']
    upperlevel = weaktable['iu']
    wavelength = weaktable['wav']
    epsilon = goodepsilon[weak_mask]
    energy_kev = weaktable['ener']
    weaklines = [Line(*attributes) for attributes in zip(element,ion,lowerlevel,upperlevel,wavelength,epsilon,energy_kev)]

    print('T=%g, %i weak lines' % (T, len(weaklines)))
    
    return weaklines


def computePseudoContinuum(T):
    """Calculate the pseudo continuum array for T"""

    # read the line file and get the weak lines
    weaklines = readWeakLines(T)
    
    energyedges = np.logspace(
        np.log10(pcontminenergy), np.log10(pcontmaxenergy),
        pcontenergysteps+1)

    pseudocontinuua = {}
    pseudoenergies = {}

    for element in all_elements:
        elidx = element_nums[element]
        energies = [line.energy for line in weaklines if line.element==elidx]
        epsilons = [line.epsilon for line in weaklines if line.element==elidx]

        summedlines, edgesout = np.histogram(
            energies, weights=epsilons, bins=energyedges)
        # divide by bin width to convert to photon cm^3/s/keV
        pseudocontinuua[element] = summedlines / (energyedges[1:]-energyedges[:-1])
        pseudoenergies[element] = 0.5*(energyedges[1:]+energyedges[:-1])

    return (pseudocontinuua, pseudoenergies)

def makeContinuumHDU(T, contbins, energies, cont, pcontbins, penergies, pcont, totcoco):
    """Make the output HDU for the continuum and pseudocontinuum for T"""

    contformat = 'PE'
    col_element = fits.Column(
        name='Z', format='1J',
        array=[element_nums[i] for i in all_elements])
    col_rmJ = fits.Column(
        name='rmJ', format='1J', array=np.zeros(len(all_elements)))

    col_N_Cont = fits.Column(
        name='N_Cont', format='1J',
        array= [contbins[i] for i in all_elements])
    col_E_Cont = fits.Column(
        name='E_Cont', format=contformat,
        unit='keV', array=np.array([energies[i] for i in all_elements],dtype=np.object_))
    col_Continuum = fits.Column(
        name='Continuum', format=contformat,
        unit='photons cm^3 s^-1 keV^-1',
        array=np.array([cont[i] for i in all_elements],dtype=np.object_))
    # for the moment we have no errors on the continuum so don't write out the column
    #    col_Cont_Err = fits.Column(
    #        name='Cont_Err', format=contformat,
    #        unit='photons cm^3 s^-1 keV^-1',
    #        array=np.array([conterrors[i] for i in all_elements],dtype=np.object_))

    pcontformat = 'PE'
    col_N_Pseudo = fits.Column(
        name='N_Pseudo', format='1J',
        array= [pcontbins[i] for i in all_elements])
    col_E_Pseudo = fits.Column(
        name='E_Pseudo', format=pcontformat,
        unit='keV',
        array=np.array([penergies[i] for i in all_elements],dtype=np.object_))
    col_Pseudo = fits.Column(
        name='Pseudo', format=pcontformat,
        array=np.array([pcont[i] for i in all_elements],dtype=np.object_))
    # for the moment we have no errors on the pseudo continuum so don't write out the column
    #    col_Pseudo_Err = fits.Column(
    #        name='Pseudo_Err', format=pcontformat,
    #        array=np.array([pconterrors[i] for i in all_elements],dtype=np.object_))

    # when we have actual errors available then use this commented out version of tabhdu
    #    tabhdu = fits.BinTableHDU.from_columns([
    #        col_element, col_rmJ,
    #        col_N_Cont, col_E_Cont,
    #        col_Continuum, col_Cont_Err,
    #        col_N_Pseudo, col_E_Pseudo,
    #        col_Pseudo, col_Pseudo_Err])
    tabhdu = fits.BinTableHDU.from_columns([
        col_element, col_rmJ,
        col_N_Cont, col_E_Cont, col_Continuum,
        col_N_Pseudo, col_E_Pseudo, col_Pseudo])

    tabhdu.name = 'EMISSIVITY'
    h = tabhdu.header
    h['HIERARCH TEMPERATURE'] = T*keV_K
    h['XTEMP'] = T
    h['DENSITY'] = 1.0
    h['TOT_COCO'] = totcoco

    return tabhdu

def processContinuum(T, cocofilename):
    """Process the temporary continuum files to make the output table"""

    # read in continum from each element file
    continuua = {}

    for element in all_elements:
        filename = os.path.join(
            tmpdir, 'tmp_conti_T%010f_%s.fits' % (T, element))
        energy, vals = readContinuumFile(filename)
        continuua[element] = vals

    # now calculate the pseudo continuua from the weak lines
    pseudocontinuua, pseudoenergies = computePseudoContinuum(T)

    # remove trailing zeroes in each continuum so we can save space
    # and save resulting number of bins
    contbins = {}
    energies = {}
    cont = {}
    conterrors = {}
    totcontbins = 0
    for element in all_elements:
        cont[element] = np.trim_zeros(continuua[element],'b')
        nbins = len(cont[element])
        contbins[element] = nbins
        totcontbins += nbins
        energies[element] = energy[:nbins]
        conterrors[element] = np.zeros(nbins)

    # do the same for the pseudo continuua
    pcontbins = {}
    penergies = {}
    pcont = {}
    pconterrors = {}
    totpcontbins = 0
    for element in all_elements:
        pcont[element] = np.trim_zeros(pseudocontinuua[element],'b')
        nbins = len(pcont[element])
        pcontbins[element] = nbins
        totpcontbins += nbins
        penergies[element] = (pseudoenergies[element])[:nbins]
        pconterrors[element] = np.zeros(nbins)

    # sum continuum flux
    totcoco = 0.
    for i in all_elements:
        totcoco += continuua[i].sum()

    # make the continuum HDU
    tabhdu = makeContinuumHDU(T, contbins, energies, cont, pcontbins, penergies, pcont, totcoco)
    
    # append the HDU
    fits.append(cocofilename, tabhdu.data, tabhdu.header)

    return (len(elements), totcontbins, totpcontbins)

###############################################################################
# routines to clean up all the temporary files created when running spex for T

def cleanUp(T):

    for element in all_elements:
        contname = os.path.join(
            tmpdir, 'tmp_conti_T%010f_%s.fits' % (T, element))
        conselname = os.path.join(
            tmpdir, 'tmp_conti_T%010f_%s.consel' % (T, element))
        deleteFile(contname)
        deleteFile(conselname)

    linename = os.path.join(tmpdir, 'tmp_lines_T%010f.fits' % T)
    lineselname = os.path.join(tmpdir, 'tmp_lines_T%010f.linesel' % T)
    deleteFile(linename)
    deleteFile(lineselname)

    scriptname = os.path.join(tmpdir, 'tmp_spex_T%010f.script' % T)
    deleteFile(scriptname)
        
    return

###############################################################################
# main routine

def main():
    """Main routine."""

    # set up the output files with their PARAMETERS extensions
    cocofilename, linefilename = initializeOutputFiles()

    # arrays which will be accumulated for temperatures
    NelementLine = []
    Nline = []
    NelementCont = []
    Ncont = []
    Npseudo = []

    # loop over temperatures. for each temperature run spex then
    # process the temporary file to make the output HDUs for the
    # line and continuum files

    for T in temperatures:

        # run spex for this temperature
        runSpex(T);

        # read the temporary file with lines to make the HDU for this
        # temperature
        numelements, numlines = processStrongLines(T, linefilename)
        NelementLine.append(numelements)
        Nline.append(numlines)

        # read the temporary files with continuua and lines to make the HDU for this
        # temperature
        numelements, numcontbins, numpcontbins = processContinuum(T, cocofilename)
        NelementCont.append(numelements)
        Ncont.append(numcontbins)
        Npseudo.append(numpcontbins)

        # tidy up the temporary files created for this temperature
        cleanUp(T)

    # set Nline and Nelement values in the PARAMETERS extension of the line file
    with fits.open(linefilename, mode='update') as ffl:
        params_hdu = ffl[1]
        params_hdu.data['Nelement'] = NelementLine
        params_hdu.data['Nline'] = Nline

    # set NCont, Npseudo and NElement values in the PARAMETERS extension of
    # the continuum file
    with fits.open(cocofilename, mode='update') as ffc:
        params_hdu = ffc[1]
        params_hdu.data['NElement'] = NelementCont
        params_hdu.data['NCont'] = Ncont
        params_hdu.data['NPseudo'] = Npseudo
        
    return    

if __name__ == '__main__':
    main()

__author__ = 'MFA Corcoran'
__version__ = '0.2'

import os
import sys
from astropy.io import fits

if sys.version_info.major <= 2:
    from xspec2 import AllData, Model, Fit, Plot
    import xspec2 as xspec
else:
    from xspec import AllData, Model, Plot

import numpy as np
import pandas as pd


def raw_input(x):
    return input(x).strip()


# class PHA(pha):
#     """
#     based on the heasp pha class with additional set_pha, set_channel methods
#
#     retrieving the counts array
#     In [39]: spec.__swig_setmethods__['Pha']=arange(10)
#     In [42]: spec.__swig_setmethods__['Pha']
#     Out[42]: array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
#
#     Actually this seems not to work...
#
#     """
#     from heasp import pha

#     def __init__(self, phafile=None):
#         PHA = pha()
#         if phafile is not None:
#             PHA.read(phafile)
#     def set_counts(self, counts):
#         self.__swig_setmethods__['Pha'] = counts
#         return
#     def set_channels(self, channels):
#         self.__swig_setmethods__['Channels'] = channels
#         return
#     def get_counts(self):
#         counts = self.__swig_setmethods__['Pha']
#         return counts
#     def get_channels(self):
#         channels = self.__swig_setmethods__['Channels']
#         return channels


# def read_xcm_old(xcmfile, verbose=True):
#     """
#     Reads and xspec command file XCM and returns the AllData and Model objects
#     :param xcmfile:
#     :return: AllData & Model instances
#     """
#     AllData.clear()
#     datastring = ''
#     ddir = os.path.split(xcmfile)[0]
#     with open(xcmfile,'r') as f:
#         xcm=f.readlines()
#     if verbose:
#         for x in xcm:
#             print(x.strip())
#     for i,l in enumerate(xcm):
#         print(l)
#         if 'data' in l:
#             ll = l.strip().replace('data ', '')
#             dnum = ll.split()[0].split(':')[-1].strip()
#             dfile = os.path.join(ddir,ll.split()[-1]).strip()
#             datastring = "{dstr} {dnum} {dfile}".format(dstr=datastring, dnum=dnum, dfile=dfile)
#             AllData(datastring)
#         elif 'resp' in l:
#             ll = l.strip().split()
#             specnum = int(ll[1].split(':')[0])
#             respfile = ll[-1]
#             AllData(specnum).response = respfile
#         elif 'arf' in l:
#             ll = l.strip().replace('arf', '')
#             arf = ll.split(':')[1].split()
#             arfnum = int(arf[0])
#             arffile = arf[1]
#             AllData(specnum).response.arf = arffile
#         elif 'model' in l:
#             mo = Model(l.replace('model ', '').strip())
#             pars = xcm[i+1:]
#             for j,p in enumerate(pars):
#                 mo(j+1).values = p.strip()
#         elif 'ignore' in l:
#             # TODO: allow setting of ignored channels
#             igstring = ignore.split()[1]
#         elif 'method' in l:
#             Fit.method = l.strip().replace('method','')
#         elif 'abund' in l:
#             ab = l.strip().replace('abund','')
#         elif 'xsect' in l:
#             xs = l.strip().replace('xsect','')
#         elif 'cosmo' in l:
#             co = l.strip().replace('cosmo','')
#         elif 'xset' in l:
#             xset = l.strip().replace('xset','')
#         elif 'systematic' in l:
#             syst = l.strip().replace('systematic','')
#     return AllData, mo

def read_xcm(xcmfile, verbose=False, default_model='TBabs*apec'):
    """
    Reads and xspec command file XCM and returns the AllData and Model objects
    :param xcmfile:
    :return: AllData & Model instances
    """
    curdir= os.getcwd()
    ddir=os.path.split(xcmfile)[0]
    # move to xcmfile directory
    if len(ddir) != 0:
        os.chdir(ddir)
    #print("Now in to directory {0}".format(os.getcwd()))
    AllData.clear()
    with open(xcmfile,'r') as f:
        xcm=f.readlines()
    if verbose:
        for x in xcm:
            print(x.strip())
    # get data string
    datastring = [x.strip() for x in xcm if 'data' in x]
    datastring = ' '.join(datastring).replace('data', '').strip()
    AllData(datastring)
    # get & load response
    resp = [x.strip() for x in xcm if "resp" in x]
    for r in resp:
        if ':' in r:
            specnum = int(r.split(':')[1].split(' ')[0])
        else:
            specnum = 1
        respname = r.split()[-1]
        AllData(specnum).response = respname
    # get & load arf
    arf = [x.strip() for x in xcm if "arf" in x]
    for r in arf:
        if ':' in r:
            specnum = int(r.split(':')[1].split(' ')[0])
            arfname =  r.split(' ')[-1]
        else:
            specnum=1
            arfname = r.split()[-1]
        AllData(specnum).response.arf = arfname
    # get & load backgrnd
    bkg = [x.strip() for x in xcm if "back" in x]
    for r in bkg:
        if ':' in r:
            specnum = int(r.split()[1])
            bkgname =  r.split()[2]
            AllData(specnum).background = bkgname
        else:
            specnum=1
            bkgname = r.split()[-1]
            AllData(specnum).background = bkgname
    # load ignored strings
    ig = [x.strip() for x in xcm if "ignore" in x]
    if len(ig)>0:
        igsplit = ig[0].split()[1:]
        for i in igsplit:
            if ':' in i:
                specnum=int(i.split(':')[0])
                igstring = i.split(':')[1]
            else:
                specnum = 1
                igstring = i
            AllData(specnum).ignore(igstring)
    # get model
    model = read_model_xcm(xcmfile, verbose=verbose, default_model='TBabs*apec')
    return AllData, model

def read_model_xcm(xcmo, default_model='TBabs*apec', verbose=False):
    """
    This function reads a model file created with the xspec
    "save model <filename>" command
    and returns a pyxspec model object
    @param xcmo: name of command file
    @return: xspecmodel
    """
    par = []
    modelfound = False
    if os.path.isfile(xcmo):
        with open(xcmo) as file:
            for line in file:
                #print line
                if 'model' in line:
                    xspecmodel = Model(line[5:].strip())
                    modelfound = True
                if modelfound:
                    # append to parameter array if line is one of the acceptable parameter types
                    #  either begins with a digit, or equal sign or slash
                    testforpar = line.strip()[0]
                    if ((testforpar.isdigit()) or (testforpar == "=") or (testforpar == "/") or (testforpar == "-")):
                        par.append(line.strip('\n').strip())
    if modelfound:
        #par = par[1:]
        for j in np.arange(len(par)):
            i=int(j)
            try:
                xspecmodel(i+1).values = par[i]
                if verbose:
                    print(f'Setting {xspecmodel(i+1).name} to {par[i]}')
            except Exception as errmsg:
                # exception can be caused if an integer parameter is specified as a float
                # so remove the decimal points from the parameter values
                print(errmsg)
                #p = par[i].replace('.00000','')
    else:
        print('Could not find {0}; Setting to {1}'.format(xcmo, default_model))
        xspecmodel = Model(default_model)
    return xspecmodel

def write_xcm(xcmfile, spectrum, model=None, clobber=False):
    """
    This function takes a spectrum object and (optionally) a model object and
    writes out an xspec12 command file
    @param xcmfile: output filename (without .xcm extension)
    @param spectrum: source spectrum object from pyxspec
    @param model: pyxspec model object
    @return:
    """
    xcm=['data '+spectrum.fileName]
    try:
        xcm.append('back ' + spectrum.background.fileName)
    except:
        pass
    xcm.append('resp ' + spectrum.response.rmf)
    try:
        xcm.append('arf '+spectrum.response.arf)
    except:
        pass
    xcm.append('ignore '+spectrum.ignoredString())

    xcm.append('ignore bad')

    xcmfileout = "{xcmfile}.xcm".format(xcmfile=xcmfile)

    if model:
        mo_xcm_list = write_xcm_model(xcmfile,model)
        for m in mo_xcm_list:
            xcm.append(m)

    if os.path.isfile(xcmfileout) and not clobber:

        print ("{0} exists".format(xcmfileout))
        ans = raw_input('Overwrite [y/n]? ')
        if ans.strip().lower() == 'n':
            print ("{0} not overwritten; Returning".format(xcmfileout))
            return
    print ("Writing Data + Model to File {0}".format(xcmfileout))
    f = open(xcmfileout, 'w')
    for i in xcm:
        f.write(i + "\n")
    f.close()
    return


def write_xcm_model(savefile, model, pyxspec = False):
    """
    This writes a model instance from pyxspec as an xspec 12 model command file
    or as a pyxspec- formatted command (need to develop a reader for
    pyxspec-formatted model file!)
    :param model: pyxspec model instance
    :param savefile: Name of xcm file to write to (without _mo.xcm - this will be added)
    :param pyxspec: if False, saves to old-style xspec format
    :return:

    20200404 MFC updated to write linked parameters
    """
    savefileout = "{0}_mo.xcm".format(savefile)
    print("Writing model command file {savefileout}".format(savefileout=savefileout))
    with open(savefileout, mode='wt') as mofile:
        mo_xcm_list=["model {0}".format(model.expression)]
        mofile.write("model {0}".format(model.expression))
        mofile.write("\n")
        if pyxspec:
            for c in model.componentNames:
                #print "component = {0}".format(c)
                for p in model.__getattribute__(c).parameterNames:
                    #print "Values of Parameter {0}".format(p)
                    val= model.__getattribute__(c).__getattribute__(p).values
                    mofile.write("model.{0}.{1}.values = {2}".format(c,p,val))
        else:
            for c in model.componentNames:
                #print "component = {0}".format(c)
                for p in model.__getattribute__(c).parameterNames:
                    #print "Values of Parameter {0}".format(p)
                    # check to see if the parameter is linked to another one
                    if len(model.__getattribute__(c).__getattribute__(p).link.strip()) != 0:
                        val = model.__getattribute__(c).__getattribute__(p).link
                        mofile.write(val)
                        mo_xcm_list.append(val)
                        mofile.write("\n")
                    # if not linked, just write the values string
                    else:
                        val= model.__getattribute__(c).__getattribute__(p).values
                        v = ["{0:.3e}".format(x) for x in val]
                        mofile.write(','.join(v))
                        mo_xcm_list.append(','.join(v))
                        mofile.write("\n")
        #mofile.close()
    return mo_xcm_list


def addspec(phafiles, Xenergy=None,
            rmf='/Users/corcoran/Dropbox/nicer_cal/nicer_resp_ver1.02/nicer_v1.02.rmf',
            arf='/Users/corcoran/Dropbox/nicer_cal/nicer_resp_ver1.02/ni_xrcall_onaxis_v1.02.arf',
            ignore="0.0-0.3", verbose=False):
    """
    from a list of phafiles from a given mission, combines the spectra in energy space by converting rates to counts, 
    combining counts, and calculating the combined total exposure
    :param phafiles: list of phafiles
    :param Xenergy: array of energies (monotonic) for the combined spectrum
    :param rmf: response matrix file (default: NICER)
    :param arf: ancillary response file (default: NICER)
    :param ignore: energy band(s) to ignore
    :return: energy array and array of the combined spectrum
    """
    # define a standard X-axis from 0.2 - 10.2 keV with 500 steps

    if not Xenergy:
        Xenergy = np.arange(500) / 50. + 0.2
    xspec.AllData.clear()
    xspec.AllData += phafiles[0]

    xspec.AllData(1).response = rmf
    xspec.AllData(1).response.arf = arf
    xspec.AllData(1).ignore(ignore)

    xspec.Plot.device = '/null'
    xspec.Plot.xAxis = 'keV'
    xspec.Plot('da')

    X = np.asarray(xspec.Plot.x(1))
    Y = np.asarray(xspec.Plot.y(1))
    CountSum = np.interp(Xenergy, X, Y) # rates
    CountSum = CountSum * xspec.AllData(1).exposure # counts
    TotExpo = xspec.AllData(1).exposure

    if verbose:
        print (xspec.AllData(1).fileName, xspec.AllData(1).exposure, TotExpo)

    for i,P in enumerate(phafiles[1:]):
        xspec.AllData.clear()
        phafile = P
        xspec.AllData += P
        xspec.AllData(1).response = rmf
        xspec.AllData(1).response.arf = arf
        xspec.AllData(1).ignore(ignore)
        xspec.Plot.xAxis = 'keV'
        xspec.Plot('da')
        X = np.asarray(xspec.Plot.x(1))
        Y = np.asarray(xspec.Plot.y(1))
        expo = xspec.AllData(1).exposure
        CountSum = CountSum + np.interp(Xenergy, X, Y)*expo
        TotExpo += expo
        if verbose:
            print (xspec.AllData(1).fileName, xspec.AllData(1).exposure, TotExpo)
    return Xenergy, CountSum, TotExpo

def combinespec(pha2sum, rmflist=None, phaout=''):
    """
    from a list of pha files from a single mission and instrument creates a simple summed counts spectrum
    then writes it out to a pha file.
    :param pha2sum: list of pha files to sum
    :param phaout: output pha file containing the combined spectrum to write, if not an empty string
    :TODO param rmflist: Energy boundaries of channels, corresponding to the phafiles in pha2sum; if not None, use EBOUNDS from the rmf files to sum the spectra
    :return: chan, countstot, expotot, gti - channel, combined counts, and combined exposure
    """
    expotot = 0.0
    countstot = 0
    bingti_start = []
    bingti_stop = []
    for p in pha2sum:
        hdu = fits.open(p)
        cts = hdu['SPECTRUM'].data['COUNTS']
        expotot = expotot + hdu['SPECTRUM'].header['EXPOSURE']
        countstot = countstot + cts
        gtistart_i = hdu['GTI'].data['START']
        gtistop_i = hdu['GTI'].data['STOP']
        bingti_start.extend(gtistart_i)
        bingti_stop.extend(gtistop_i)
    chan = fits.open(pha2sum[0])['SPECTRUM'].data['CHANNEL']

    # write out binned phafile to phasout

    hdu = fits.open(pha2sum[0])  # read in a spectrum to use as template
    hdu['SPECTRUM'].data['COUNTS'] = countstot
    hdu['SPECTRUM'].header['EXPOSURE'] = expotot

    col1 = fits.Column(name='START', format='E', array=bingti_start)
    col2 = fits.Column(name='STOP', format='E', array=bingti_stop)

    cols = fits.ColDefs([col1, col2])
    newgti = fits.BinTableHDU.from_columns(cols)
    newgti.header['EXTNAME'] = 'GTI'
    hdu['GTI'] = newgti
    hdu.writeto(phaout, overwrite=True, checksum=True)
    return chan, countstot, expotot

def showinnotebook(xspecobject):
    """
    shows an xspec model or spectrum object (or anything with a show() method)
    in a jupyter notebook
    :param xspecobject: model or spectrum object
    :return:
    """
    from xspec import Xset
    chat = Xset.chatter
    Xset.chatter=10
    from wurlitzer import sys_pipes
    try:
        xspecobject.__getattribute__('show')
    except AttributeError:
        print ("Object has no show() method; returning")
        return
    with sys_pipes():
        xspecobject.show()
    Xset.chatter = chat
    return

def plotinnotebook(data, model, device='null', xAxis = 'keV',
                   figsize=[8,6], fmt='k.', alpha=0.2, title=''):
    """
    Plots spectrum and fit from an xspec session
    in the default jupyter notebook backend ("inline","notebook", etc)

    :param data: xspec data object
    :param model: xspec model objec
    :param device: xspec Plot device
    :return: x, y, yerr, mo
    """
    from matplotlib.pyplot import figure, errorbar, plot
    Plot.xAxis=xAxis
    Plot.device = device
    Plot('da')
    x = Plot.x()
    y = Plot.y()
    yerr = Plot.yErr()
    xerr = Plot.xErr()
    mo = Plot.model()
    fig = figure(figsize=figsize)
    errorbar(x, y, yerr=yerr, xerr=xerr, fmt=fmt, alpha=alpha)
    plot(x, mo, 'r-')
    return {'x':x, 'y':y, 'xerr':xerr, 'yerr':yerr, 'model':mo}




def lineid(wl, temp, mineps= 1e-18, nei=False, tau=1e11, Te_init=1e4, printit=False):
    """
    This code will produce a list of lines in a given wavelength range at a
    given temperature. It also shows the use of an NEI version, where you
    have to additionally specify the initial ionization temperature (or the
    ionization fraction directly) and the elapsed Ne*t.

    The results of the list_lines codes are numpy arrays which can be sorted any
    way you wish. You can, of course, extract the lines easily at this point. There
    is also a print_lines routine for a fixed format output.

    Based on
    https://atomdb.readthedocs.io/en/master/examples.html?highlight=line%20list#make-line-list

    from code written by Adam Foster 2015-12-02
    Version 0.1

    :param wl: 2 element array with start, end wavelengths (angstrom)
    :param temp: electron temperature in K
    :param mineps: minimum emissivity (epsilon) value
    :param nei: if True, calculate line strength for non-equilibrium plasma (slow), else equilibrium
    :param tau: electron density * time (cm^-3 s) for NEI case
    :param Te_init: initial ionization balance temperature (NEI only) in K
    :return: line_list
    """
    import pyatomdb

    # as of version 0.0.0.3, Introduced new way to calculated spectra,
    # using the Session and Spec objects
    # so check if version later than 0.0.0.3 and if so use new method

    Te = temp

    vers = pyatomdb.__version__.split('.')
    if ((int(vers[0])>=0) &  (int(vers[1])>0) & (int(vers[2])>0)):
        if nei:
            res = pyatomdb.spectrum.NEISession().return_linelist(Te, tau, wl)
        else:
            res = pyatomdb.spectrum.CIESession().return_linelist(Te, wl)
    else:
        if not nei:
            # get equilibrium line list
            res = pyatomdb.spectrum.list_lines(wl, Te=Te, teunit='K', minepsilon=mineps)
        else:
            # now do an NEI version. This is slow at the moment, but functional.
            Te_init = temp
            tau = 1e11
            res = pyatomdb.spectrum.list_nei_lines(wl, Te=Te, teunit='K', \
                                                   minepsilon=mineps, \
                                                   Te_init=Te_init, \
                                                   tau=tau)

    # specify wavelength range, in Angstroms
    # wl = [8.0, 9.0]
    # electron temperature in K

    if printit:
        # re-sort by element, ion then emissivity
        res.sort(order=['Element', 'Ion', 'Epsilon'])
        print ("sorted by Element, Ion, Emissivity:")
        pyatomdb.spectrum.print_lines(res)

    # change byte order to avoid "Big-endian buffer not supported on little-endian compiler error"
    line_list = pd.DataFrame(res.byteswap().newbyteorder())
    line_list.drop(['Epsilon_Err', 'Lambda_Err'], axis=1, inplace=True)

    line_list['Symbol'] = ["{0} {1}".format(pyatomdb.atomic.Ztoelsymb(x),
                                            pyatomdb.atomic.int_to_roman(y)) for x, y
                           in zip(line_list['Element'],line_list['Ion'])]
    #line_list.drop(['Epsilon_Err', 'Lambda_Err'], axis=1, inplace=True)
    return line_list

def get_components(model, xAxis='keV', device='null'):
    """
    for a model with one or more additive components, returns the energy and counts for each component
    :param model: XSPEC model instance
    :return: compdict, dictionary of the x, y values of the individual components
    """
    compdict=dict()
    norms=dict()
    for i in range(model.nParameters):
        k=i+1
        p = model(k)
        if p.name == 'norm':
            norms[k]=p.values[0]
    for k in norms.keys():
        model(k).values=0.0
        #print(model(k).values)
    for k in norms.keys():
        ckey = "comp{k}".format(k=k)
        model(k).values=norms[k]
        #print(model(k).values, norms[k])
        Plot.device = device
        Plot.xAxis = xAxis
        Plot('data')
        plotvals = dict()
        plotvals['x'] = Plot.x()
        plotvals['model'] = Plot.model()
        compdict[ckey]= plotvals
        # set minimum norm to 0.0
        model(k).values="0.0,,0.0,0.0"
        #model(k).values=0.0
        #print('\n\n\nAFTER*******************:')
        #print(model(k).values, norms[k])
        #time.sleep(2.5)
    for k in norms.keys():
        model(k).values = norms[k]
    return compdict

def get_abund_table():
    import pandas as pd
    abund = os.path.join(os.environ['HEADAS'],'../Xspec/src/manager','abundances.dat')
    tab = pd.read_csv(abund, sep='\s+', nrows=7)
    tab.set_index('elts:', inplace=True)
    # isel = np.where(tab.index == 'References:')[0][0]
    # tab = tab.iloc[0:isel-1]
    newind = {}
    for i in tab.index:
        newind[i] = i[:-1]
    tab.rename(newind, axis=0, inplace=True)
    tab.index.rename('Abundance Table', inplace=True)
    return tab

def get_mo_params(model, verbose=True):
    """
    This retrieves the parameters for a model instance from pyxspec
    :param model: pyxspec model instance
    :return: returns a dictionary describing the model
    """
    if verbose:
        print('model("{0}")'.format(model.expression))
    modeldict = dict()
    modeldict['Expression'] = model.expression
    for c in model.componentNames:
        modeldict[c]=dict()
        for p in model.__getattribute__(c).parameterNames:
            #print "Values of Parameter {0}".format(p)
            val= model.__getattribute__(c).__getattribute__(p).values
            if verbose:
                print("model.{0}.{1}.values = {2}".format(c,p,val))
            modeldict[c][p] =val
    return modeldict


def get_specparams(obsID, model, phaname, rmffile, arffile,
                   calc_errors=False,
                   backfile=None,
                   fluxband="2.0 10.0",
                   statMethod="cstat",
                   ignore="0.0-0.45, 7.5-**",
                   gtinum=None,
                   writexcm=True,
                   xcmroot = '',
                   use_xset_save = False,
                   allowPrompting = False,
                   clobber=False,
                   dofit=True,
                   verbose=False):
    """Get spectrum parameters from fit
    This function gets the observation and spectrum parameters from a model fit
    for the given obsid (and optionally the give interval for the obsid)

    :param obsID: id (number or string) for the spectrum to be analyzed
    :param phaname: name of phafile (with directory path); will be constructed if not specified
    :param model: xspec model object to compare to/fit to spectrum
    :param get_errors: if True gets the parameter errors (simple method) for non-frozen parameters
    :param calc_errors: if True calculate errors for non-frozen parameters; if False, get parameter sigma as error
    :param rmffile: response file
    :param arffile: effective area file; if None, then the rmffile is a response with the arf folded into it
    :param workdir: user-defined work directory for output
    :param datadir: user-defined directory holding input pha file
    :param fluxband: band over which to calculate fluxes in keV
    :param statMethod: statistic to use in fit ("chi", "cstat")
    :param ignore: energy range in keV to ignore for fit
    :param gtinum: if not None, append a this number to the output xcm files; this is useful to divide an obsid spectrum by time
    :param writexcm: if True, writes the xspec command file after fit to xcmroot
    :param xcmroot: root name to use for xcm file (constructed if not specified)
    :param use_xset_save: if True use the xspec.Xset.save method to write the xcm files
    :param clobber: if True overwrite xcm file when xcm file written
    :param verbose: increase chattiness of output
    :param dofit: if True, fit the model, otherwise apply existing model to the data and return parameters
    :return: pandas DataFrame of best fit spectrum parameters

    Versions:
       0.2: Don't return parameter error if parameter frozen or linked
    """
    __version__ = 0.2

    import xspec
    from astropy.table import Table
    from astropy.time import Time

    if gtinum is None:
        obs = str(obsID)
    else:
        obs = "{0}_{1}".format(obsID, gtinum)

    xspec.AllData.clear()

    xspec.Xset.allowPrompting=allowPrompting

    try:
        if verbose:
            print('Loading pha file {0}'.format(phaname))
        pha = xspec.Spectrum(phaname)
        # don't forget to ignore bad channels
        xspec.AllData.ignore('bad')
        skip_calc = False
        if backfile:
            try:
                pha.background = backfile
            except Exception as e:
                print("Can't find background file {0} ({1})".format(backfile, e))
    except Exception as errmsg:
        print("Problem analyzing {0} ({1})".format(phaname,errmsg))
        skip_calc = True
    if not skip_calc:
        pha.response = rmffile
        if verbose:
            print(f'Setting RMF to {pha.response.rmf}')
        if arffile is not None:
            pha.response.arf = arffile
        if verbose:
            print(f'Setting ARF to {pha.response.arf}')
        pha.ignore(ignore)

        if xspec.Fit.dof < 1:
            print('Number of degrees of freedom < 1: Cannot perform fit')
            status = -1
            return status
        if dofit:
            if verbose:
                print(f'Initial Model is {model.expression}')
                model.show()
            xspec.Fit.statMethod = statMethod
            print('Fitting')
            try:
                xspec.Fit.perform()
            except Exception as e:
                print(f"Can't perform fit for {obs} ({e}); returning")
                status = -1
                return status
            if writexcm:
                if use_xset_save:
                    xs_xcmfile = f'{xcmroot}_xset.xcm'
                    if os.path.exists(xs_xcmfile):
                        if clobber:
                            os.remove(xs_xcmfile)
                        else:
                            print(f'{xs_xcmfile} exists and clobber = False')
                    print(f'Xset Saving {xs_xcmfile}')
                    xspec.Xset.save(xs_xcmfile, info='a')
                    xs_xcmofile = xs_xcmfile.replace('.xcm','_mo.xcm')
                    if os.path.exists(xs_xcmofile):
                        if clobber:
                            os.remove(xs_xcmofile)
                        else:
                            print(f'{xs_xcmofile} exists and clobber = False')
                    print(f"Xset Saving {xs_xcmofile}")
                    xspec.Xset.save(xs_xcmofile,info='m')
                else:
                    if len(xcmroot) == 0:
                        xcmroot = input('enter xcmroot (with directory path > ')
                    write_xcm(xcmroot, pha, model=model, clobber=clobber)
        if verbose:
            print('Calculating fluxes')
        xspec.AllModels.calcFlux(fluxband)

        flux = pha.flux[0]
        print("{0}    flux = {1:.3e} {2} Fit Statistic = {3:.1f}".format(obs, flux, fluxband,
                                                                     xspec.Fit.statistic))

        # update JDSTART, JDEND, JDMID using GTI info from pha file

        hdu = fits.open(phaname)
        gtiname = 'GTI'
        try:
            gti = Table.read(phaname,hdu=gtiname)
        except KeyError:
            gtiname = 'STDGTI'
            gti = Table.read(phaname, hdu=gtiname)
        if gtinum:
            scstart = gti[gtinum]['START']
            scend = gti[gtinum]['STOP']
        else:
            scstart = gti['START'].min()
            scend = gti['STOP'].max()
        gti['Duration'] = gti['STOP'] - gti['START']
        Expos = gti['Duration'].sum()
        #gtiname = 'GTI'
        try:
            mjdoff = hdu[gtiname].header['MJDREFI'] + hdu[gtiname].header['MJDREFF']
        except KeyError:
            mjdoff = hdu[gtiname].header['MJDREF']

        jdstart = Time(scstart/86400 + mjdoff, format='mjd').jd
        jdend = Time(scend/86400 + mjdoff, format='mjd').jd
        jdmid = jdstart+(jdend-jdstart)/2.0

        if gtinum is not None:
            numgtis = 1
        else:
            numgtis = len(gti)

        nobsDF = pd.DataFrame(data={obsID:model.expression}, index=['model'])

        nobsDF.loc['JDSTART'] = jdstart
        nobsDF.loc['JDEND'] = jdend
        nobsDF.loc['JDMID'] = jdmid
        nobsDF.loc['EXPOSURE'] = Expos
        nobsDF.loc['Num_GTI'] = numgtis

        # store the spectral parameters in the data frame

        nobsDF.loc['Flux'] = flux
        nobsDF.loc['FluxBand'] = fluxband
        nobsDF.loc['FluxErr']= flux*pha.rate[1]/pha.rate[0]
        nobsDF.loc['Fit_Statistic'] = xspec.Fit.statistic
        nobsDF.loc['dof'] = xspec.Fit.dof
        nobsDF.loc['Rate'] = pha.rate[2]  # Total rate (without background subtraction)
        nobsDF.loc['RateErr'] = np.sqrt(pha.rate[2]*Expos)/Expos
        nobsDF.loc['NetRate'] = pha.rate[0] # Net rate
        nobsDF.loc['NetRateErr'] = pha.rate[1] # Net rate error
        nobsDF.loc['BkgRate'] = nobsDF[obsID]['Rate'] - nobsDF[obsID]['NetRate']
        nobsDF.loc['BkgRateErr'] = np.sqrt(nobsDF[obsID]['RateErr']**2 + nobsDF[obsID]['NetRateErr']**2) # Net rate error

        for ip in range(model.nParameters):
            p = model(ip+1)
            pname = '{0}_{1}'.format(p.name, p.index)
            nobsDF.loc[pname] = p.values[0]
            if (len(p.link)==0) and (not p.frozen):
                # get error if parameter not linked or frozen
                pn_err="{0}_err".format(pname)
                if not (calc_errors):
                    nobsDF.loc[pn_err] = p.sigma
                else:
                    xspec.Fit.error("2.706 {0}".format(p.index))
                    lobnd = model(p.index).errors[0]
                    hibnd = model(p.index).errors[1]
                    nobsDF.loc[pn_err] = (hibnd-lobnd)/2.0
        nobsDF.loc['NullHyp'] = xspec.Fit.nullhyp
        # Calculate absorption-corrected flux (set all NH = 0)
        unabsFlux = nhcorrectflux(pha, model, fluxband)[0]
        nobsDF.loc['UnabsFlux'] = unabsFlux
        return nobsDF

def nhcorrectflux(pha, model, fluxband):
    """
    Calculates the absorption-corrected flux for a given xspec spectrum and model instance
    :param pha: xspec xspectrum instance
    :param model: xspec model instance
    :param fluxband: xspec flux band ("2.0 10.0" for example)
    :return: absorption-corrected flux (flux with all NH components set to 0)
    """
    import xspec
    for ip in range(model.nParameters):
        p = model(ip + 1)
        if 'nH' in p.name:
            model(ip + 1).values = [0.0, 0.01, 0.0, 0.0]
    xspec.AllModels.calcFlux(fluxband)
    return pha.flux


def ignore_negative_channels(pha):
    """
    determines if there are negative channels in the xspec Spectrum object
    and if so, ignores those channels
    :param pha: xspec Spectrum object
    :return: pha with ignored set to negative channels
    """
    from xspec import Plot
    import numpy as np
    Plot.xAxis = 'channel'
    Plot.device = 'null'
    Plot('da')
    chan = np.asarray(Plot.x())
    net = np.asarray(Plot.y())
    try:
        i = np.where(net < 0)[0]
    except:
        # all channels > 0 so return
        return
    # need to add 1 since xspec starts channel counting at 1 not 0
    negchan = np.asarray(chan)[i]+1
    ig = pha.ignored # get list of currently ignored channels
    igchan = ig.extend(list(negchan)) # add in negative channels
    # sort and remove duplicates
    igchan.sort()
    igchan=list(set(igchan))
    return igchan

if __name__== "__main__":
    #xcmo = '/Users/corcoran/research/ETA_CAR/CHANDRA2/repro/seqid/200810/10787/work/meg-1_mo.xcm'
    # xcmo = '/Users/corcoran/research/ETA_CAR/RXTE/WORKSPACE/DATA/processed_data/reduced/1996020915/1996020915_pcu2_L1_3_fit.xcm'
    # mo = read_model_xcm(xcmo)
    # modeldict = show_mo_params(mo)
    # for k in modeldict.keys():
    #     print (k)
    #     print (modeldict[k])
    #ad = read_xcm('/Users/corcoran/program/missions/NICER:OSWG/wr140/work/1120010115/ni1120010115_0mpu7_cl_bin20.xcm')
    #print(ad(1).response.arf)
    #moxcm = '/Users/corcoran/program/missions/NICER:OSWG/cygx3/work/2142010109/ni2142010109_0mpu7_cl_2_mo.xcm'
    #mo = read_model_xcm(moxcm)
    #xcm='/Users/corcoran/program/missions/NICER:OSWG/seec/work/Eps_Eri/2300030101/test.xcm'
    #xcm = '/Users/corcoran/program/missions/NICER:OSWG/cygx3/work/2142010109/ni2142010109_0mpu7_cl_2.xcm'
    #ad, mo = read_xcm(xcm)
    #ad.show()
    #compdict = get_components(ad, mo, xAxis='keV', device='/xw')
    #ad, mo = read_xcm(
    #    '/Users/corcoran/program/missions/NICER:OSWG/wr140/work/1120010102/ni1120010102_0mpu7_cl_bgsub.xcm')
    #mo.show()
    #compdict = get_components(ad,mo)
    mo = read_model_xcm('/Users/mcorcora/program/missions/NICER:OSWG/wr140/work/Template_spectrum/Pollock/WR140.3bvapec.template.2020-05-22.model_mfc.xcm')


# def addlines(alldata, model, nlines, modeltype="gauss", Erange =[6.0,8.0],dofit=True,
#             numtrials=100, statmethod="chi"):
#     # adds lines to the model at energies of largest residuals
#     # based on the tcl script addline.tcl version of kaa 1/5/00  modified for xspec v11
#     """
#     Usage: addline nlines modeltype fit|nofit"
#     Adds to the model nlines additional lines of type modeltype"
#      "
#     The nlines additional lines are added one at a time."
#     Line energies are set to that of the largest residual"
#     between the data and the model. For each line a fit is"
#     performed with the line width and normalization as the only"
#     free parameters."
#      "
#     The default options are one line and a gaussian. The other"
#     model type that can be used is lorentz. If no third argument"
#     is given then the sigma and normalization of each line are"
#     fit. If nofit is specified then the fit is not performed but"
#     if fit is specified then all free parameters are fit."
#
#     kaa  1/5/00"
#
#     :param alldata: pyxspec AllData object (may include data groups)
#     :param model:  pyxspec model object
#     :param nlines: number of lines to add
#     :param modeltype: xspec model to use for line (default: gaussian)
#     :param Erange: energy range over which to add lines
#     :param dofit: if True, perform fit with entire model; if False, freeze model parameters except for line sigma and norm then fit
#     :param numtrials: number of iterations to perform when fitting
#     :param statmethod: statistics method to use (chi, cstat, etc)
#     :return: model with lines added
#     """
#
#     from xspec import *
#
#     # get the current number of datagroups
#     ngroups = alldata.nGroups
#     Fit.nIterations = numtrials
#     Fit.statMethod = statmethod
#
#     # loop over ngroups
#     for n in range(ngrops):
#         # loop over the lines required
#         for i in range(nlines):
#             # get the current number of model components
#             # TODO:
#             ncomp = model
#             # get the current number of parameters
#             numparam = model.nParameters
#             # work out where the next model component and model parameter must go
#             set
#             nextcomp[expr $ncomp /$ngroups + 1]
#             set
#             nextparam[expr $numparam /$ngroups + 1]
#
#             # and the number of parameters per datagroup
#
#             set
#             parpergroup[expr $numparam /$ngroups]
#
#             # set up an array of the parameter deltas
#             for {set impar 1} {$impar <= $numparam} {incr impar} {
#             tclout param $impar
#             scan $xspec_tclout "%f %f" value delta($impar)
#             }
#
#             # get the peak residual energy and strength
#
#             tclout
#             peak
#             scan $xspec_tclout
#             "%f %f"
#             peake
#             peakn
#
#             puts
#             " "
#             puts[format
#             "New peak at %f keV..." $peake]
#             puts
#             " "
#
#             # add the new component. first we have to set up the model parameters string.
#
#             set
#             pstring
#             " "
#             append
#             pstring
#             "& $peake  1.e-4 "
#             append
#             pstring
#             "& 0.0     1.e-4 "
#             append
#             pstring
#             "& $peakn  1.e-4 "
#
#             for {set igroup 2} {$igroup <= $ngroups} {incr igroup} {
#             append pstring "& = $nextparam "
#             append pstring "& = [expr $nextparam + 1] "
#             append pstring "& = [expr $nextparam + 2] "
#             }
#
#             puts[format
#             "%s addcomp %d %s %s" $prompt $nextcomp $linetype $pstring]
#             addcomp $nextcomp $linetype $pstring
#
#             # if the standard fit option is required (fit=1, nofit=1) then
#             # freeze all parameters except for the line sigma and norm
#
#             if {$nofit == 1 & & $fit == 1} {
#             for {set igroup 1} {$igroup <= $ngroups} {incr igroup} {
#             set istart[expr ($igroup-1) * ($parpergroup+3) + 1]
#             set iend[expr $istart + $parpergroup]
#             puts[format "%s freeze %d-%d" $prompt $istart $iend]
#             freeze $istart-$iend
#             }
#             }
#
#             # if the nofit option has not been set (nofit=0) then
#             # fit for the sigma and normalization
#
#             if dofit:
#                 #         puts "$prompt fit 100"
#                 Fit.perform()
#
#                 # if the standard fit option is required (fit=1, nofit=1) then
#                 # thaw all the parameters that were originally variable
#
#             if {$nofit == 1 & & $fit == 1} {
#             for {set igroup 1} {$igroup <= $ngroups} {incr igroup} {
#             for {set impar 1} {$impar <= $parpergroup} {incr impar} {
#
#             set thispar[expr ($igroup-1) * ($parpergroup+3) + $impar]
#             if {$delta($thispar) > 0.} {
#             puts[format "%s thaw %d" $prompt $thispar]
#             thaw $thispar
#             }
#
#             }
#             }
#
#             # and thaw the energy of the new line
#
#             puts[format "%s thaw %d" $prompt $nextparam]
#             thaw $nextparam
#
#             }
#
#             }
#
#         return
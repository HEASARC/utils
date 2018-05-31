__author__ = 'MFA Corcoran'
__version__ = '0.1'

import xspec
import numpy as np
import os
from heasp import pha

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



def read_model_xcm(xcmo, chatter=True):
    """
    This function reads a model file created with the xspec "save model <filename>" command
    and returns a pyxspec model object
    @param xcmo: name of command file
    @return: xspecmodel
    """
    par = []
    modelfound = False
    with open(xcmo) as file:
        for line in file:
            # print line
            if 'model' in line:
                if chatter:
                    print "model is {0}".format(line[5:].strip())
                xspecmodel = xspec.Model(line[5:].strip())
                modelfound = True
            if modelfound:
                par.append(line.strip('\n').strip())
    if modelfound:
        par = par[1:]
        for i in np.arange(len(par) + 1):
            try:
                xspecmodel(i).values = par[i - 1]
            except:
                # exception can be caused if an integer parameter is specified as a float
                # so remove the decimal points from the parameter values
                p = par[i-1].replace('.00000','')
        xspecmodel.show()
    else:
        xspecmodel = 0
    return xspecmodel


def write_xcm(xcmfile, spectrum, model=None, clobber=False):
    """
    This function takes a spectrum object and (optionally) a model object and writes out a xspec12 command file
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

    xcmfileout = "{xcmfile}.xcm".format(xcmfile=xcmfile)

    if model:
        mo_xcm_list = write_xcm_model(xcmfile,model)
        for m in mo_xcm_list:
            xcm.append(m)

    if os.path.isfile(xcmfileout) and not clobber:

        print "%s exists" % xcmfileout
        ans = raw_input('Overwrite [y/n]? ')
        if ans.strip().lower() == 'n':
            print "{0} not overwritten; Returning".format(xcmfileout)
            return
    print "Writing Data + Model to File {0}".format(xcmfileout)
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

    TODO: needs to be able to handle when one parameter is tied to another)
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
                    val= model.__getattribute__(c).__getattribute__(p).values
                    v = [str(x) for x in val]
                    mofile.write(','.join(v))
                    mo_xcm_list.append(','.join(v))
                    mofile.write("\n")
        #mofile.close()
    return mo_xcm_list

def get_mo_params(model, chatter=True):
    """
    This retrieves the parameters for a model instance from pyxspec
    :param model: pyxspec model instance
    :return: returns a dictionary describing the model
    """
    if chatter:
        print('model("{0}")'.format(model.expression))
    modeldict = dict()
    modeldict['Expression'] = model.expression
    for c in model.componentNames:
        modeldict[c]=dict()
        for p in model.__getattribute__(c).parameterNames:
            #print "Values of Parameter {0}".format(p)
            val= model.__getattribute__(c).__getattribute__(p).values
            if chatter:
                print("model.{0}.{1}.values = {2}".format(c,p,val))
            modeldict[c][p] =val
    return modeldict


def addspec(phafiles, Xenergy=None,
            rmf='/Users/corcoran/Dropbox/nicer_cal/nicer_resp_ver1.00tmp180215a/nicer_v1.00tmp180210b.rmf',
            arf='/Users/corcoran/Dropbox/nicer_cal/nicer_resp_ver1.00tmp180215a/ni_xrcall_onaxis_v1.00tmp180215a.arf',
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
        print xspec.AllData(1).fileName, xspec.AllData(1).exposure, TotExpo

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
            print xspec.AllData(1).fileName, xspec.AllData(1).exposure, TotExpo
    return Xenergy, CountSum, TotExpo

def showinnotebook(xspecobject):
    """
    shows an xspec model or spectrum object (or anything with a show() method)
    in a jupyter notebook
    :param xspecobject: model or spectrum object
    :return:
    """
    from wurlitzer import sys_pipes
    try:
        xspecobject.__getattribute__('show')
    except AttributeError:
        print "Object has no show() method; returning"
        return
    with sys_pipes():
        xspecobject.show()
    return


if __name__== "__main__":
    #xcmo = '/Users/corcoran/research/ETA_CAR/CHANDRA2/repro/seqid/200810/10787/work/meg-1_mo.xcm'
    xcmo = '/Users/corcoran/research/ETA_CAR/RXTE/WORKSPACE/DATA/processed_data/reduced/1996020915/1996020915_pcu2_L1_3_fit.xcm'
    mo = read_model_xcm(xcmo)
    modeldict = show_mo_params(mo)
    for k in modeldict.keys():
        print k
        print modeldict[k]

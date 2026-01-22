def raw_input2(a):
    import sys
    pyver = sys.version.split()[0].split('.')[0]
    if (int(pyver) >2) :
        return input(a).strip()
    else:
        return raw_input(a)


import os
from subprocess import Popen, PIPE
try:
    from wurlitzer import sys_pipes
except:
    print('Could not import wurlitzer')
from astropy.io import fits


def xselect(obsid, datadir='/Volumes/SXDC/Data/NICER/etacar',
            workdir='/Users/corcoran/research/ETA_CAR/NICER/work', evtdir = None, evtfile=None, tbin=100.,
            gtilim=0.0, notebook = False, command=None, pha_cutoff=[40,1000],gti=None, clobber=True):

    """
    uses xselect to extract a time-filtered spectrum and lightcurve from an events file

    :param obsid: NICER observation ID (integer)
    :param datadir: path to directory containing the NICER obsids (not needed if evtdir specified)
    :param workdir: output will be placed in workdir/obsid subdirectory
    :param evtdir: Directory path containing event file; if not specified will be created
    :param evtfile: name of event file in evtdir; if not specified used cleaned events
    :param tbin: time bin size in seconds for lightcurve
    :param gtilim: if gti has less than this duration (in seconds) do not include the gti in the gti_spectra extraction
    :param notebook: set to True if being used in notebook
    :param command: initial xselect command to use; could be for example tcurs or gtispec
    :param pha_cutoff: a 2-element list with the lower and upper pha cutoffs in channels
    :param clobber: if True overwrite output

    the script will read the events and extract a lightcurve binned to tbin seconds (default value = 100)
    and will then plot the spectrum.

    either enter a PLT command or "quit" to quit the PLT environment

    Then:

    a) filter time cursor (remember to exit from PLT> prompt!); special command tcurs will do this automatically
    b) save time cursor to outdir = workdir/obsid/; special command savetcurs will do this automatically
    c) once a time filter has been saved, a spectrum will be extracted and saved, and the program will exit

    convenience commands:
        tcurs: creates a time cursor filter
        savetcurs: save the time filter, extracts a spectrum, saves the spectrum, extracts the lightcurve, saves the lightcurve
        gtispec: extract spectra within gtis and saves the extracted data
        halt: exits the xselect environment
     """

    if not evtdir:
        evtdir = os.path.join(datadir, str(obsid), 'xti/event_cl')

    if not evtfile:
        evtfile = Popen('ls {evtdir}/ni*mpu7*cl.evt'.format(evtdir=evtdir),shell=True, stdout=PIPE).communicate()[0]
        if len(evtfile) == 0: # see if the .gz file exists
            evtfile = Popen('ls {evtdir}/ni*mpu7*cl.evt.gz'.format(evtdir=evtdir), shell=True, stdout=PIPE).communicate()[0]
        if len(evtfile) == 0:
            print ("Could not find uncompressed or gzipped event "
                   "file in {0}; returning".format(evtdir))
            return
        evtfile = os.path.split(evtfile)[1].strip().decode("ascii")

    #outroot = evtfile.rsplit('.evt')[0]
    outroot = evtfile.replace('.evt','')
    outdir = os.path.join(workdir, str(obsid))
    if not os.path.exists(outdir):
        print ("{outd} does not exist; creating directory".format(outd=outdir))
        os.mkdir(outdir)
    if notebook:
        with sys_pipes():
            run_xselect(evtdir, evtfile, outroot, outdir, tbin,
                        command=command, gtilim=gtilim, gti=gti, clobber=True)
    else:
        run_xselect(evtdir, evtfile, outroot, outdir, tbin, command=command, gtilim=gtilim,
                    pha_cutoff=pha_cutoff, gti=gti, clobber=clobber)
    return


def run_xselect(evtdir, evtfile, outroot, outdir, tbin, gtilim=0.0,
                gti=None, command=None, pha_cutoff=[], clobber=True):
    """
    runs an interactive xselectproc session
    :return:
    """
    xselectproc = Popen('xselect', stdin=PIPE, shell=True, universal_newlines=True, bufsize=0)
    xselectproc.stdin.write("xsel_session\n")
    xselectproc.stdin.write("read events\n")
    xselectproc.stdin.write("{0}\n".format(evtdir))
    xselectproc.stdin.write("{0}\n".format(evtfile))
    xselectproc.stdin.write("yes\n") # reset the mission
    xselectproc.stdin.write("set device /xw\n")
    if len(pha_cutoff) == 2:
        cmd = "filter pha_cut {chanlo} {chanhi}\n".format(chanlo=pha_cutoff[0],chanhi=pha_cutoff[1])
        xselectproc.stdin.write(cmd)

    xselectproc.stdin.write("set bin {tbin}\n".format(tbin=tbin))
    xselectproc.stdin.write("extract curve\n")
    ans = 'dummy'
    plt = False
    if command is not None:
        ans = command.lower().strip()
        if "gtispec" in ans:
            gti_spectra(xselectproc, evtdir, evtfile, outdir, outroot, tbin, gtilim=gtilim, gti=gti)
            ans="halt"
        elif "tcurs" in ans:
            ans = 'filter time cursor\n'.format(out=os.path.join(outdir, outroot))
            plt = True
        elif "get_spec" in ans:
            get_spec(xselectproc, evtdir, evtfile, outdir, outroot, tbin, clobber=clobber)
            ans="halt"
        else:
            print("Command must be either 'tcurs', 'get_spec', or 'gtispec'; returning")
            ans='halt'
    while ans.lower().strip() != 'halt':
        if plt:
            prompt = 'PLT> '
        else:
            prompt = 'Enter XSELECT Command (type "halt" to exit XSELECT)> '
        ans = raw_input2(prompt)
        if ans.strip().lower() == 'tcurs':
            ans = 'filter time cursor\n'.format(out=os.path.join(outdir, outroot))
            plt = True
        if "gtispec" in ans.strip().lower():
            gti_spectra(xselectproc, evtdir, evtfile, outdir, outroot, tbin, gtilim=gtilim, gti=gti)
        if 'savetcurs' in ans.strip().lower():
            # check to see if time cursor already exists
            foutroot = os.path.join(outdir, outroot)
            gtifile = "{foutroot}.curs_gti".format(foutroot=foutroot)
            ans = "save time cursor {out}\n".format(out=foutroot)
            xselectproc.stdin.write(ans)
            if os.path.isfile(gtifile):
                # file exists, so check clobber
                if clobber:
                    print ("Clobbering File")
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "save time cursor {out}\n".format(out=foutroot)
                    xselectproc.stdin.write(ans)
                    ans = "no\n"
                    xselectproc.stdin.write(ans)
            ans = "extract events\n"
            xselectproc.stdin.write(ans)
            evtfile = "{foutroot}.evt".format(foutroot=foutroot)
            ans = "save events {out}\n".format(out=evtfile)
            xselectproc.stdin.write(ans)
            if os.path.isfile(evtfile):
                # file exists, so check clobber
                if clobber:
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "no\n"
                    xselectproc.stdin.write(ans)
            # response to "Use filtered events as input data file?" is "yes"
            ans = "yes\n"
            xselectproc.stdin.write(ans)
            ans = "extract spectrum\n"
            xselectproc.stdin.write(ans)
            phafile = "{foutroot}.pha".format(foutroot=foutroot)
            ans = "save spectrum {out}\n".format(out=phafile)
            xselectproc.stdin.write(ans)
            if os.path.isfile(phafile):
                # file exists, so check clobber
                if clobber:
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "no\n"
                    xselectproc.stdin.write(ans)
            ans = "extract curve\n"
            xselectproc.stdin.write(ans)
            lcfile = "{foutroot}.lc".format(foutroot=foutroot)
            ans = "save curve {out}\n".format(out=lcfile)
            xselectproc.stdin.write(ans)
            if os.path.isfile(lcfile):
                # file exists, so check clobber
                if clobber:
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "no\n"
                    xselectproc.stdin.write(ans)
            ans = ""
        if (ans.strip().lower() == "quit") and plt:
            plt = False
            # xselectproc.stdin.write("quit\n")
        if ans.strip().lower()[0:4] == "plot":
            plt = True
            xselectproc.stdin.write("{ans}\n".format(ans=ans))
        else:
            xselectproc.stdin.write("{ans}\n".format(ans=ans))
    xselectproc.stdin.write('exit\n')
    xselectproc.stdin.write('no\n')
    xselectproc.wait()

def get_spec(xselectproc, evtdir, evtfile, outdir, outroot, tbin, gtilim=0.0, clobber=True):
    foutroot = os.path.join(outdir, outroot)
    xselectproc.stdin.write('\n')

    # clear all time cursors
    ans = 'clear time keyboard all\n'
    xselectproc.stdin.write(ans)
    # response to "Proceed? (yes or no) >[yes]" is yes
    # ans = 'yes\n'
    # xselectproc.stdin.write(ans)

    ans = "show filter time\n"
    xselectproc.stdin.write(ans)

    # clear Events
    ans = "clear events\n"
    xselectproc.stdin.write(ans)

    # extract, save spectrum
    ans = "extract spectrum\n"
    xselectproc.stdin.write(ans)
    phafile = "{foutroot}.pha".format(foutroot=foutroot)
    ans = "save spectrum {out}\n".format(out=phafile)
    xselectproc.stdin.write(ans)
    if os.path.isfile(phafile):
        # file exists, so check clobber
        if clobber:
            ans = "yes\n"
            xselectproc.stdin.write(ans)
        else:
            ans = "no\n"
            xselectproc.stdin.write(ans)


def gti_spectra(xselectproc, evtdir, evtfile, outdir, outroot, tbin, gtilim=0.0,
                gti=None, clobber=True):
    """
    extract spectra (and lightcurve and events) within gtis if the duration of the gti is greater than gtilim

    :param xselectproc: xselect process opened by Popen
    :param evtdir: directory holding events file
    :param evtfile: name of events file
    :param outdir: output directory
    :param outroot: root name for output files
    :param tbin: time bin size in seconds for lightcurve
    :param gtilim: use only gtis which have a duration in seconds greater than gtilim
    :param gti: list of [start, stop] times for the data extraction, in mjd
    :param clobber: if True overwrite files
    :return:
    """
    hdu = fits.open(os.path.join(evtdir,evtfile))
    if gti is None:
        gti = list(hdu['GTI'].data)
        mjdoff = hdu['GTI'].header['MJDREFI'] + hdu['GTI'].header['MJDREFF']
        gti = [[x / 86400.0 + mjdoff, y / 86400. + mjdoff, y-x] for x, y in gti]
    # create a spectrum and ligthcurve for each gti
    for i, g in enumerate(gti):
        if g[2] > gtilim:
            outrooti = "{0}_{1}".format(outroot, i)
            foutroot = os.path.join(outdir, outrooti)
            xselectproc.stdin.write('\n')

            # clear all time cursors
            ans = 'clear time keyboard all\n'
            xselectproc.stdin.write(ans)
            # response to "Proceed? (yes or no) >[yes]" is yes
            # ans = 'yes\n'
            # xselectproc.stdin.write(ans)

            ans = "show filter time\n"
            xselectproc.stdin.write(ans)

            # clear Events
            ans = "clear events\n"
            xselectproc.stdin.write(ans)

            # reload events
            # xselectproc.stdin.write("read events\n")
            # ans = "{0}\n".format(evtdir)
            # xselectproc.stdin.write(ans)
            # ans = "{0}\n".format(evtfile)
            # xselectproc.stdin.write(ans)

            xselectproc.stdin.write("show data\n")

            # then extract a lightcurve
            ans="extract curve\n"
            xselectproc.stdin.write(ans)

            # now do the time selection
            ans = 'filter time mjd\n'
            xselectproc.stdin.write(ans)
            # quit the PLT environment
            xselectproc.stdin.write("quit\n")
            # input the MJD filters
            tstart = gti[i][0]
            tstop  = gti[i][1]
            ans ='{0:.5f}, {1:.5f}\n exit\n'.format(tstart,tstop)
            xselectproc.stdin.write(ans)
            # exit time MJD filtering environment
            #xselectproc.stdin.write("x\n")


            # extract, save spectrum
            ans = "extract spectrum\n"
            xselectproc.stdin.write(ans)
            phafile = "{foutroot}.pha".format(foutroot=foutroot)
            ans = "save spectrum {out}\n".format(out=phafile)
            xselectproc.stdin.write(ans)
            if os.path.isfile(phafile):
                # file exists, so check clobber
                if clobber:
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "no\n"
                    xselectproc.stdin.write(ans)
            # extract, save lightcurve
            ans = "extract curve\n"
            xselectproc.stdin.write(ans)
            lcfile = "{foutroot}.lc".format(foutroot=foutroot)
            ans = "save curve {out}\n".format(out=lcfile)
            xselectproc.stdin.write(ans)
            if os.path.isfile(lcfile):
                # file exists, so check clobber
                if clobber:
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "no\n"
                    xselectproc.stdin.write(ans)

            # # now extract, save events
            ans = "extract events\n"
            xselectproc.stdin.write(ans)
            evtf = "{foutroot}.evt".format(foutroot=foutroot)
            ans = "save events {out}\n".format(out=evtf)
            xselectproc.stdin.write(ans)
            if os.path.isfile(evtf):
                # file exists, so check clobber
                if clobber:
                    ans = "yes\n"
                    xselectproc.stdin.write(ans)
                else:
                    ans = "no\n"
                    xselectproc.stdin.write(ans)

            # response to "Use filtered events as input data file?" is "no"
            xselectproc.stdin.write("no\n")

def test_gtispec(obsid = '1142010101',datadir = '/Volumes/SXDC/Data/NICER/cyg_x3/proc',
                 workdir = '/Users/corcoran/program/missions/NICER:OSWG/cygx3/work',
                 gtilim = 5.0, command=None):
    xselect(obsid, datadir=datadir, workdir=workdir, gtilim=gtilim,
            notebook=False, clobber=True, command=command)
    pass

def test_getspec(obsid = '1142010101',datadir = '/Volumes/SXDC/Data/NICER/cyg_x3/proc',
                 workdir = '/Users/corcoran/program/missions/NICER:OSWG/cygx3/work',
                 command=None):
    xselect(obsid, datadir=datadir, workdir=workdir,
            notebook=False, clobber=True, command=command)
    pass

def test_xsel_py3():
    from heasarc.utils import xselect_utils as xs
    o = 1120010181
    workdir = '/Users/mcorcora/research/WR140/NICER/work/'
    datadir = '/Users/mcorcora/SXDC/Data/NICER/wr140'
    xs.xselect(o, datadir=datadir, workdir=workdir, tbin = 10, notebook=False, clobber=True, pha_cutoff=[40, 600])


if __name__ == '__main__':
    #obs = 1110010140
    #datadir = '/Volumes/SXDC/Data/NICER/etacar/repro'
    #workdir = "/Users/corcoran/research/ETA_CAR/NICER/work/repro"
    #xselect(obs, datadir=datadir, workdir=workdir, tbin=10., notebook=False)

    #test_gtispec(command='gtispec')
    #test_getspec(command='gtispec')
    test_xsel_py3()

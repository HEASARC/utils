import datetime
import re
import urllib2
import requests
from bs4 import BeautifulSoup
import re
import pandas as pd

def cdscat2browse_request(bibcode, browse_request_file, browse_table_name='', observatory='', master_catalog='',
                          table_type='', freq_regime='', server='cdsarc.u-strasbg.fr'):
    """
    from a CDS data table, use the CDS ReadMe file to create a request for creation of a Browse table
    Sample CDS ReadMe file ftp://cdsarc.u-strasbg.fr/pub/cats/J/ApJS/208/17/ReadMe
    param: bibcode - string, should be of the form YYYYJournal..volume...page, for example '2013ApJS..208...17A'
    param: browse_request_file: Name of output browse request file
    return: creates browse table readme file
    """
    #
    # find the tables in the CDS ReadMe file
    # they can be found in the lines beginning with "Byte-by-byte Description of file:'
    #
    b = bibcode.split('..')
    yr = b[0][0:4]
    journal = b[0][4:]
    vol = b[1]
    pg = b[2][1:-1]
    pgsuffix = b[2][-1]
    rurl = 'ftp://{3}/pub/cats/J/{0}/{1}/{2}/ReadMe'.format(journal, vol, pg, server)
    print 'Retrieving {0}'.format(rurl)
    #freadme = open(ReadMe_file, mode='r')
    try:
        freadme = urllib2.urlopen(rurl)
    except Exception, errmsg:
        print 'Error reading {0}'.format(rurl)
        print 'Error message: "{0}"'.format(errmsg)
        print 'Returning'
        return

    readme = freadme.readlines()
    freadme.close()
    catindx = []
    catname = []
    # generate table author name
    auth = readme[0].split('   ')[-1].strip().strip('(').strip(')').replace('+',' et al.')
    for i, r in enumerate(readme):
        if 'Byte-by-byte Description of file:' in r:
            catindx.append(i)
            catname.append(r.split('file:')[1].strip())
    catendindx = [x - 1 for x in catindx[1:]]
    catendindx.append(len(readme) - 1)
    #
    # print available catalogs with description and let user select one
    #
    catdes = []
    print "Available catalogs"
    for i, c in enumerate(catname):
        # print "#{0} {1} {}".format(i, c)
        pat = re.compile(c)
        for r in readme:
            m = pat.match(r)
            if m:
                d = r.split('  ')[-1].strip()
                catdes.append(d)
                # print catdes
                print "#{0} {1} {2}".format(i, c.replace('.dat', ''), d)
    catselindx = raw_input('Enter number of selected catalog >')
    catselindx = int(catselindx)

    # 0 or more spaces followed by one or more digits followed by zero or more "-" followed by 1 or more spaces followed by one or more digits
    # p1 = re.compile('\s*\d*-*\s*\d*')

    s = '--------------------------------------------------------------------------------'

    p1 = re.compile(
        '\s*\d*-\s*\d*\s*')  # 0 or more spaces followed by one or more digits followed by a - followed by 1 or more spaces followed by one or more digits

    dbytes = []
    dformat = []
    dunit = []
    dname = []
    dhname = []
    dcomment = []
    err = re.compile('e_', re.I) # re.I means case insensitive matching
    for r in readme[catindx[catselindx] + 4:catendindx[catselindx] - 1]:
        if not s in r:
            if r[0:10].strip():  # if first 10 characters are not spaces process line
                # print r
                comment = ''
                m = p1.match(r)
                if m:
                    # process line to get data
                    dbytes.append(m.group())
                    test = r[m.end():].split(' ')
                    l = [x for x in test if x != '']
                    dformat.append(l[0])
                    dunit.append(l[1])
                    name = l[2].strip()
                    dname.append(name)
                    em = err.match(name)
                    if em:  # if name begins with an "e_" or "E_", assume it's an error value
                        name = "{0}_error".format(name[em.end():])
                    if name=='RAdeg':
                        name ='RA' # assume this is right ascension
                    if name == 'DEdeg':
                        name = 'Dec' # assume this is  Dec
                    dhname.append(name)
                    comment = ' '.join(l[3:])
                    dcomment.append(comment.strip().replace("\n", " "))

            else:
                comment = "{0} {1}".format(comment.strip(),
                                             r.strip())  # first 10 comments are blank so append as continuation comment
                # print comment
                dcomment[-1] = comment.strip()  # replace last comment with appended comment

        else:
            break  # end of data reached
    #
    # Write output
    #
    delim = '--------------------------------------------------------------------------------'
    metadata = [['Request to Create New Database', str(datetime.datetime.now())],
                ['New Database Name', browse_table_name],
                ['db description', catdes[catselindx]],
                ['w3b description', catdes[catselindx]],
                ['Observatory', observatory],
                ['Master Catalog', master_catalog],
                ['Table Type', table_type],
                ['Table Author', auth],
                ['Row Type', ' '],
                ['Frequency Regime', freq_regime],
                ['Default Search Radius', 1],
                ['Line Summary Display', ''],
                ['Unique Key', ''],
                ['Catalog Bibcode', "{0}{1}.{2}.{3}{4}".format(yr,journal,vol,pg, pgsuffix)],
                ['Add Class Field', ' '],
                ['Directory', ' '],
                ['Data File', 'ReadMe'],
                ['This File', 'README']
                ]
    print "\n"
    for m in metadata:
        print "{0:<30s}: {1:<}".format(m[0], m[1])
    print delim
    print "\n"
    print "Byte-by-byte Description of file: {0}".format(catname[catselindx])
    title = "{0:>10s}  {1:>5s} {2:>9s}  {3:>10s} {4:>15s}  {5:<}".format('Bytes', 'Format', 'Units', 'Label',
                                                                         'HEASARC Name', 'Comment')
    print delim
    print title
    print delim

    a = [len(x) for x in dbytes]
    dbml = max(a)

    for db, df, du, dn, dh, dc in zip(dbytes, dformat, dunit, dname, dhname, dcomment):
        print "{0:>10s}  {1:>5s} {2:>9s}  {3:>10s} {4:>15s}  {5:<}".format(db.strip(), df.strip(), du, dn.strip(),
                                                                           dh.strip(), dc)
    #
    # write to output request file
    #
    print "\n Writing to {0}".format(browse_request_file)
    fout = open(browse_request_file, 'w')
    fout.write("\n")
    for m in metadata:
        fout.write("{0:<30s}: {1:<}\n".format(m[0], m[1]))
    fout.write(delim)
    fout.write("\n")
    fout.write("Byte-by-byte Description of file: {0}\n".format(catname[catselindx]))
    fout.write(delim)
    fout.write("\n")
    fout.write("{0}\n".format(title))
    fout.write(delim)
    fout.write("\n")

    a = [len(x) for x in dbytes]
    dbml = max(a)

    for db, df, du, dn, dh, dc in zip(dbytes, dformat, dunit, dname, dhname, dcomment):
        fout.write("{0:>10s}  {1:>5s} {2:>9s}  {3:>10s} {4:>15s}  {5:<}\n".format(db.strip(), df.strip(), du, dn.strip(),
                                                                           dh.strip(), dc))
    fout.close()

    return

def cds_get_ucds(vizier_url):
    """
    Get the UCD values for a vizier table for the parameters in a browse README table creation request file
    :param: browse_readme = name of the browse README table creation request file
    :param: vizier_url = url of the vizier page for the table
    :return:
    """
    #vizier_url = 'http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/588/A103'
    req = requests.get(vizier_url)
    soup = BeautifulSoup(req.text,'lxml')
    tables = soup.findAll('table')
    table = tables[7] # parameter table seems to be table 7 - not sure how general this is
    rows = table.find_all('tr')
    ucd_dict = dict()
    for r in rows:
        col = r.findAll('td')
        try:
            cdspar = col[2].text.strip()
        except:
            cdspar = "Not Found"
        #pat = re.compile(r'\([a-z]*\.\D*\)')
        # this pattern matches an open parenthesis, combinations of letters, periods, semicolons, numbers and dashes and a closing parenthesis
        pat = re.compile('\([a-z,A-Z]*\.{1}[a-z,;,.,A-Z,-]*[0-9]?\)')
        try:
            cdsucd = pat.search(col[4].text).group().split(') (')[0].strip('(').strip(')')
        except:
            cdsucd='Not Found in {0}'.format(cdspar)
        #print " Parameter {0} has UCD {1}".format(cdspar, cdsucd)
        ucd_dict[cdspar]=cdsucd
    return ucd_dict

def README_to_dict(readmefile, header=0,footer=0):
    """
    reads the browse README file then
    returns a dictionary of the cds parameter names and the corresponding brows parameter name, i.e.
    readme_dict={'BROWSE_PARNAME1':'CDS_PARNAME1', ...}
    :param readmefile: BROWSE readme file
    :param header: number of header lines to skip
    :param footer: number of footer lines to skip
    :return:
    """
    f = open(readmefile,'r')
    lines = f.readlines()
    numlines=len(lines)
    readme_dict = dict()
    for l in lines[header:numlines-footer]:
        # first get rid of the integer and dash, which begins some (not all) of the lines
        try:
            x = re.split('\d+\-', l)[1]
        except IndexError:
            x=l
        # then split x on 1 or more spaces - strip whitespace from x first; element 3 and 4 are the parameter names
        cols = re.split(r'\s{1,}', x.strip())
        try:
            cdspar = cols[3].strip()
        except IndexError:
            cdspar = 'cds_par_smissing'
        try:
            browsepar = cols[4].strip()
        except IndexError:
            browsepar = 'browse_par_missing'
        readme_dict[browsepar] = cdspar
    return readme_dict

def strip(text):
    try:
        return text.strip()
    except AttributeError:
        return text

def update_ucdfile(vizier_url, UCDfile, README, header=0, footer=0, UpdateAll = False):
    """
    Updates the initial UCD file produced when ingesting a table into BROWSE based on the
    UCDs assigned by the CDS
    :param vizier_url: the vizier URL of the original table
    :param UCDfile: the name of the initial ucd file
    :param README: the name of the BROWSE README file used to request the table ingest
    :param header: header lines in BROWSE README file to skip
    :param footer: footer lines in BROWSE README file to skip
    :param UpdateAll: if False, updates UCDs only the parameters that are missing UCDs
    :return: creates new ucd file
    """
    # read the initial ucd file as a pandas dataframe
    df_ucd = pd.read_csv(UCDfile, sep="|", names=['browse_table', 'browse_parameter', 'ucd'],
                         converters = {'browse_parameter':strip, 'ucd':strip})
    # get the dictionary of  browse parameters and their cds names from the README file
    rdict = README_to_dict(README, header=header, footer=footer)
    # get ucds from the vizier table
    ucd_dict = cds_get_ucds(vizier_url)
    # update ucds in the dataframe
    for i, bp in enumerate(df_ucd['browse_parameter']):
        # get cds parameter corresponding to browse parameter
        try:
            cdspar = rdict[bp]
        except KeyError:
            print "Browse parameter {0} not found in README dictionary".format(bp)
            cdspar = 'missing'
        if not cdspar == 'missing':
            try:
                ucd = ucd_dict[cdspar]
            except KeyError:
                print "CDS parameter {0} not found in UCD dictionary".format(cdspar)
                ucd = ' '
            if UpdateAll:
                df_ucd.iloc[i]['ucd'] = ucd.strip()
            else:
                # only updating ucds if they are missing
                if len(df_ucd.iloc[i]['ucd'].strip()) == 0:
                    df_ucd.iloc[i]['ucd'] = ucd.strip()
    return df_ucd

def write_ucdfile(ucdfile,ucd_dataframe):
    """
    Writes out the ucd file
    :param ucdfile: path and name of the output ucd file
    :param ucd_dataframe: dataframe containing the ucd for each parameter
    :return:
    """
    fout = open(ucdfile,mode='w')
    for i,tn  in enumerate(ucd_dataframe['browse_table']):
        fout.write("{0}|{1:24s} |{2}\n".format(tn, ucd_dataframe.iloc[i]['browse_parameter'],
                                               ucd_dataframe.iloc[i]['ucd'] ))
    fout.close()
    return

def test_cdscat2browse():
    # cdsr="/Users/corcoran/program/HEASARC/missions/Fermi/Browse_Tables/2PC_auxiliary_files_v04/CDS/ReadMe"
    cdsr = '2013ApJS..208...17A'
    obs = 'Fermi'
    cdscat2browse_request(cdsr,
                          '/Users/corcoran/program/HEASARC/missions/Fermi/Browse_Tables/2PC_auxiliary_files_v04/CDS/browse_request_test.txt',
                          observatory=obs, master_catalog='X-ray', browse_table_name='LAT2PSC',
                          freq_regime='gamma-ray', table_type='Object')
    # cdsr = '2015MNRAS..448...3766K'
    # cdscat2browse_request(cdsr, '/Users/corcoran/browse_request_test.txt')


def test_cds_get_ucds():
    vizier_url = 'http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/588/A103'
    ucd_dict = cds_get_ucds(vizier_url)
    print ucd_dict.keys()


def test_readme_to_dict():
    rfile = 'sample_data/README'
    header = 27
    footer = 12
    rdict = README_to_dict(rfile, header=header, footer = footer)
    print rdict.keys()

def test_update_ucdfile():
    vizier_url = 'http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/588/A103'
    rfile = 'sample_data/Readme'
    ucdfile = 'sample_data/rass2rxs.ucd'
    header = 27
    footer = 12
    df_ucd = update_ucdfile(vizier_url,ucdfile,rfile,header=header, footer=footer)
    # print df_ucd['ucd']
    ucdfile = 'sample_data/rass2rxs.ucd.update'
    write_ucdfile(ucdfile, df_ucd)
    return

if __name__=="__main__":
    #test_cdscat2browse()
    #test_cds_get_ucds()
    #test_readme_to_dict()
    test_update_ucdfile()
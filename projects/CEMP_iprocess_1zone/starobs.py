import h5py
import data_plot as dp
import os
import time
import utils
import numpy as np
import data_support as ds
import matplotlib

'''
GENERAL INFORMATIONS :
Contact : laurent.dardelet@ens.fr

To see how this db can be used immediately, see HOW TO USE THIS MODULE section, later.

This module has been developped since June 2014.
Its main purpose is to fill the lack of long-term efficient way to import and use observational data to compare to the work done by the Nugrid collaboration.
It includes the entire interface with a database file, composed by functions to import data, retrieve this data, and interact with it.
The idea of this module is that one can either build his own database, or group it with several people.

The current world :

Mr. A wants to compare his simulation with data from stars 1 and 2. 
He looks for this data in papers or web databases, imports this data in a - probably - random format, writes a script to use this data for a particular purpose. 
It should work at some point, the plot is done. The data is most likely lost, or forgotten, and probably the same data can be found in Mr. A's HDD several times.

The amazing world with the starobs module.

Mr. A still wants to compare his simulations with data from stars 1 and 2.
He looks in the common database, and finds data about star 1, added by Mr. B not so long ago. He can retrieve immediately this data, interactively with python, maybe the plotting tools Mr. A intended to use are already in place.
Then he looks for star 2. Less lucky, star 2 is not in the database. So he imports the data from the paper, with the requested information for properly feed the database. Then he uses this data, to do the plots, and the data is safe.
Mr. B, two days later, is happy to find star 2 in the database.

================================

This module is at its first stages. 
All the basics are already here. One can add data, retrieve and work with it.

Main ideas :

 - All the data is added the way we find things in the paper, that depends on a solar abundances set. It would be interesting to be able to add data specifying this solar abundances set, and retrieve with an other set. Maybe in that case the data should be stored in absolute abundance, and not in bracket notations...
 - Now, the very basic information needed is Z and [X/Fe]. The data can be stored, but needs to have references, precise and to be safe. Adding data must request proper information about the paper, and no data should be stored without error value.
 - The oxydation levels for each element have to be taken into account. For now, one element is stored, whatever its oxydation level is. It makes some paper tables harder to  handle, and the stored data potentially wrong.
 - Add the possibility to specify only a min cycle for autofit.

Main code improvements to make : (see later the scheme of the file)

 - The groups of the H5 file are currently named by the star name, that is bad, because one cannot then rename a star, incase there was an error when entering data.
Groups should be named by the number of this star in this database when it is added, and the name given as an attribute of that group, as it is done  for the datasets.
 - The row references '!' are not automatic now. They should be, in order to  be able to add all the data we want, no matter what is usualy stored. 
   => for row in data:
        attri.add(row_ID)='!'+str(data.index(row)) for example

File scheme :
Most of this file works as a dictionnary, for groups, datasets and attributes.

------------------
MAIN COMMENTS
------------------
GROUP 'HE0143-0113' : | Attributes :
                      | '12' = '1' These key attributes are the Z for the available elements for this star.
                      | '13' = '2' These value attributes are the default dataset for this element for this star.
                      | '30' = '1'
                      | '60' = '1'
                      | '62' = '2'
                      | ... This allows to set default datasets for each element, and change it easily.
                      | ...
                      | DATASET '1' | Attributes : 
                      |             | 'Author' = 'This guy' First important information on this dataset, to keep track of the data.
                      |             | 'Year' = '2002' Second important information on this dataset. We should know now the paper it refers.
                      |             | 'Added on' = '06/30/2014' More data about when that data has been added in that db.
                      |             | 'Z' = '!0' '!' mark is used to know to what row in the values this attributes refers to. 
                      |             | The data is stored by numpy.ndarrays, but to make sure that we are retrieving the good attribute (incase something changes
                      |             | in the db module) and to have the possibility to add more than these 5 following rows, the db adds itself the row. 
                      |             | Z is obviously the list of charge number of the elements present in this dataset.
                      |             | 'XFe' = '!1' The [X/Fe] abundances.
                      |             | 'Err' = '!2' The sigma adopted error.
                      |             | 'Err_type' = '!3' A try, about if it is a random or systematic error. Not used.
                      |             | 'Upper' = '!4' 0 if it is a regular abundance, 1  if it is an upper limit. We could imagine -1 for lower limits.
                      |             | Values :
                      |             | [[12, 13, 30, 60],
                      |             | [0.5, 4.2, -1.02, 3.7],
                      |             | [0.2, 0.3,  0.12, 0.15],
                      |             | [0.0, 0.0, 0.0, 0.0],
                      |             | [0.0, 0.0, 1.0, 0.0]]
                      | DATASET '2' | Attributes : 
                      |             | 'Author' = 'This other guy'
                      |             | 'Year' = '2002'
                      |             | 'Added on' = '06/30/2014'
                      |             | 'Z' = '!0' These informations are stored here and not in the group attributes or in the main comments because we can imagine 
                      |             | that two datasets about one star don't give the same kind of information or precision. Two datasets won't have the same 
                      |             | attributes.
                      |             | 'XFe' = '!1'
                      |             | 'Err' = '!2'
                      |             | 'Err_type' = '!3'
                      |             | 'Upper' = '!4'
                      |             | Values :
                      |             | [[13, 62],
                      |             | [1.3, 3.1],
                      |             | [0.17, 0.17],
                      |             | [0.0, 0.0],
                      |             | [0.0, 0.0]]


=======================================

HOW TO USE THIS MODULE

The starobs module intends to create a database for observational star data.
It goes along with the database file, a HDF5 file, that can not be modified manually. 
It uses the python module h5py, but is internally designed in such a format that 
modifying it manually would most likely kill the database.

The starobs module is coupled with different modules for plots and study of this data.
The main idea is to create a common database for the NuGrid collaboration, hosted in 
a remote location.

In order to create and maintain this database user-friendly, structured and correct, 
it is highly advised to check all the data entered, and give author and year of 
publication of the paper the data is extracted from.

All the abundances stored until now are in the astronomical way: (X/Fe)/(Xsolar/Fesolar)

Usage 
======
INITIATE THE DATABASE :

Import the module : 

>>> import starobs

And start the instance  :

>>> db=starobs.start()



ADD DATA :

>>> db.add_dataset('CS31062-050', author='Aoki', year='2006', data=...)

Three different ways to enter data :
data is a string : it must be the location of a file, that is an ascii table, with 
the format :
 * & Na & * & * & * & * & 2.52 & * & 0.2
 2.52 in this example is the abundance and 0.2 is sigma.
data is a list :
 [[11, 12], [2.52, 3.23], [0.2, 0.1], [0, 0], [0, 1]]
 The columns are :
 Z
 XFe
 Sigma
 Error type : random or systematic. (Not being used yet)
 Upper limit : 1 for yes, 0 for no.
data is None (default) :
 One can the enter the data manually, type !h for help.

The last way to enter data is to merge two databases :

>>> db.merge('db2.h5', set_default='new')

The imported database will now be included in the current one.



PRIORITIES SYSTEM :

The main idea of this database is to interactively work with plotting and support tools. One 
may want to plot as many elements as possible about one single star. All this informations 
about the elements may not be given in one single paper. We need in that case to merge several 
datasets, to gather information.
Incase two datasets give information about one single element, we need to define priorities 
about this element : which datasets tells the truth (or the truth you want at least).

Each time you add data, you have the possibilitie to define priorities, while adding the data, 
or later with the command :

>>> set_default(dataset='1', star='CS31062-050',  z='10')

Here, I want the Neon information for CS31062-050 to be given by the dataset '1'.



GETTING INFORMATION ABOUT THE DATABASE :

Which stars do we have data about ?
>>> db.starlist

(Variable and not function !)

What datasets are available for a particular star ?
>>> db.datasets(star='CS31062-050')

What attributes can we get from this database ?
>>> db.attri_list(star='CS31062-050', dataset='1')

If no dataset is specified, the attributes of the first dataset of the star will be displayed. 

If you reached this point following step by step the examples, you must be very bored of writing
"CS31062-050' all the time :
>>> db.use_star('CS31062-050')

This command will activate the star as the current default one. If no star is specified in most 
functions, the star data will be taken trom this one.

SEPARATE DATABASE INFO, STUDY ON STAR AND PLOTS :

When the instance is created, 2 subclasses are initialized.
>>> db.plot.* 

will access the different plotting tools. Some come from data_plot, some others are specific to 
the ppn module.
And :
>>> db.data_support.* 

Will access differents tools (that have yet to be added for most of them) to work on this specific
stars data.
'''

def _confirm(default=None):
    '''Very simple private function, to ask confirmation y/n.'''
    n='n'
    y='y'
    if default=='y':
        y='Y'
    elif default=='n':
        n='N'
    print "[{0}/{1}]".format(y,n)
    while(1):
        confirm=raw_input().lower()
        if confirm == 'y':
            return True
        elif confirm=='n':
            return False
        elif confirm=='' and default<>None:
            return default=='y'


def _prepareData(data):
    '''Private function to prepare the data in the good shape for the creation of a new set'''
    XFe=[]
    El=[]
    Err=[]
    Err_Type=[]
    Uppers=[]
    if type(data)==list:
        data[0]=[float(i) for i in data[0]]
        for i in range(len(data[2])):
            data[2][i]=data[2][i]
        return data
    if type(data)==str:
        try:
            datafile=open(data, 'r')
        except:
            return None
        print "File found"
        for line in datafile:
            entries=line.split('&')
            XFe.append(float(entries[6].strip()))
            try:
                nline=float(entries[7].strip())
            except:
                nline=10000
            El.append(int(utils.get_z_from_el(entries[1].strip())))
            Err.append(float(entries[8].strip()))
            if '<' in line:
                Uppers.append(1)
            else:
                Uppers.append(0)
            Err_Type.append(0)
        if not len(El)==len(XFe)==len(Err)==len(Err_Type):
            return None
    if data==None:
        cyclepos=0
        entering = True
        upper=False
        print "Enter your new dataset here. (enter !h for help) :"
        print ""
        while entering:
            if cyclepos%3 == 0:
               print "Element :"
            if cyclepos%3 == 1:
                print "Abundance ["+utils.get_el_from_z(El[-1])+"/Fe] :"
            if cyclepos%3 == 2:
                print "Error :"
            entry=raw_input()
            if entry=='':
                if cyclepos%3 == 0:
                    print "Stop entering data now ?"
                    if _confirm(default='y'):
                        entering=False
                if cyclepos%3 == 1:
                    print "Wrong abundance entry, please resume"
                if cyclepos%3 == 2:
                    Err.append(float(0))
                    Err_Type.append(0)
                    Uppers.append(0)
                    cyclepos+=1
            elif entry[0]=='!':
                if entry[1]=='h':
                    print "Enter successively the symbol, or atomic number of the element, the abundance, and the error value."
                    print "Press ENTER between each entry."
                    print ""
                    print "To signify an upper limit, enter '<' before either the abundance or the error value."
                    print "Enter :"
                    print "'!c' to correct the previous entry,"
                    print "'!n' to change the entire previous element data, "
                    print "'!e' or press ENTER with an empty element entry to end this prompt."
                if entry[1]=='e':
                    entering=False
                    if cyclepos%3 == 1:
                        El=El[:-1]
                    if cyclepos%3 == 2:
                        El=El[:-1]
                        XFe=XFe[:-1]
                if entry[1]=='c' and cyclepos<>0:
                    cyclepos-=1
                    if cyclepos%3 == 0:
                        El=El[:-1]
                    if cyclepos%3 == 1:
                        XFe=XFe[:-1]
                    if cyclepos%3 == 2:
                        Err=Err[:-1]
                        Err_Type=Err_Type[:-1]
                        Uppers=Uppers[:-1]
                if entry[1]=='n' and cyclepos<>0:
                    cyclepos-=1
                    while cyclepos%3 <> 0:
                        if cyclepos%3 == 1:
                            XFe=XFe[:-1]
                        if cyclepos%3 == 2:
                            Err=Err[:-1]
                            Err_Type=Err_Type[:-1]
                            Uppers=Uppers[:-1]
                        cyclepos-=1
                    if cyclepos%3 == 0:
                        El=El[:-1]
            else:
                if cyclepos%3 == 0:
                    missed=False
                    try:
                        z=int(entry)
                        El+=[z]
                        cyclepos+=1
                    except:
                        try:
                            z=int(utils.get_z_from_el(entry[:2]))
                            El_tmp=z
                            entry=entry[2:]
                        except:
                            try:
                                z=int(utils.get_z_from_el(entry[:1]))
                                El_tmp=z
                                entry=entry[1:]
                            except:
                                print "Wrong element ID, please resume"
                                missed=True
                        if not missed:
                            cyclepos+=1
                            El.append(El_tmp)
                elif cyclepos%3 == 1:
                    if '<' in entry:
                        upper=True
                        entry=entry.split('<')[1]
                    try:
                        XFe.append(float(entry))
                        cyclepos+=1
                    except:
                        print "Wrong abundance, please resume"
                elif cyclepos%3 == 2:
                    if '<' in entry:
                        upper=True
                        entry=entry.split('<')[1]
                    try:
                        err=float(entry)
                        cyclepos+=1
                        Err_Type.append(0)
                        Err.append(err)
                        if upper:
                            Uppers.append(1)
                        else:
                            Uppers.append(0)
                    except:
                        print "Wrong error entered, please resume"
            if cyclepos%3 == 0 and cyclepos>0:
                print "Last entry :"
                toprint = "Element : {0}".format(utils.get_el_from_z(El[-1]))
                toprint+=",  {0}".format(El[-1])
                print toprint
                print "Abundance : {0}, error : {1}".format(XFe[-1], Err[-1])
                if Uppers[-1]:
                    print "Upper limit"
                print "---"
                print ""
                upper=False
    data=[El, XFe, Err, Err_Type, Uppers]
    return(data)

def _confirm_dataset(data_prepared, star, author, year):
    '''Private functions that displays the data that will be entered in the database, to make make sure the data seems correct.
    Since it is hard to delete data from the database, it is better to check it...
    '''
    print "The data that is about to be saved is :"
    for i in range(len(data_prepared[0])):
        toprint="[{0}/Fe] = {1}, sigma = {2}".format(utils.get_el_from_z(data_prepared[0][i]),data_prepared[1][i], data_prepared[2][i])
        if data_prepared[4][i]==1:
            toprint+=", upper limit"
        print  toprint
    print "For {0}.".format(star)
    print "Author of publication : {0} in {1}".format(author, year)
    return _confirm(default='y')

def start(filename='starobs_db.h5'):
    '''
    Function that initializes the databese.
    Can be used, as well as starobs directly.
    '''
    return starobs(filename=filename)

class plot(dp.DataPlot, ds.Plots):
    '''
    Class to avoid all data_plot functions to be displayed in starobs, to separate the plotting method from the database methods.
    '''
    def __init__(self, out_class_plot):
        print "Plotting class initialized." 
        self._out_class_plot=out_class_plot

    def get(self, attri_plot, star_plot=None, dataset_plot=None, display_plot=True,solar_norm=[None,None]):
        '''
        Private redirection method from plot to starobs
        '''
        return self._out_class_plot.get(attri=attri_plot, star=star_plot, dataset=dataset_plot, display=display_plot,solar_norm=solar_norm)

    def plotstar(self):
        '''
        Fastest way to plot the element abundance, with errorbar of a star.
        '''
        matplotlib.pyplot.errorbar(self.get('Z'), self.get('XFe'), yerr=self.get('Err'))

class data_support(ds.Tools):
    '''
    Class to avoid all data_support functions to be displayed in starobs, to separate the studying methods and plotting methods proper to ppn from the database methods.
    '''
    def __init__(self, out_class_ds):
        print "Support class initialized." 
        self._out_class_ds=out_class_ds

    def get(self, attri_ds, star_ds=None, dataset_ds=None, display_ds=True,solar_norm=[None,None]):
        '''
        Private redirection method from plot to starobs
        '''
        return self._out_class_ds.get(attri=attri_ds, star=star_ds, dataset=dataset_ds, display=display_ds,solar_norm=solar_norm)


class starobs():
    '''
    Reads the starobs stars database in a hdf5 file.
    Should have several plotting methods in it as well.

    Also allows to edit this database, and ennter new stars or datasets.
    The way to enter these data still have to be precisely defined.
    '''
    def __init__(self, filename='starobs_db.h5'):
        '''
        Initializes the class.
        A default name is given, considering that we shouln't care about the file itself.
        The avoid the file to be lost, We place it in the personnal directory.
        '''
        self.filename=filename
        print "Loading the star database from the file {0}".format(self.filename)
        if not os.path.isfile(filename):
            self._create()
            print "This file doesn't exist, a new database has been created"
        self.starobsfile=h5py.File(self.filename, 'a')
        self.starlist=[]
        self.current_star=None
        
        self.plot=plot(self)
        self.data_support=data_support(self)

        for star_group in self.starobsfile.items():
            self.starlist.append(star_group[0])
        print "{0} star(s) have been found in this database".format(len(self.starlist))

    def close(self):
        '''
        Only way to cleanly close the h5 file.
        If it is not done, some changes MAY remain unsaved (though the module is made to avoid such mistakes.
        '''
        self.starobsfile.close()

    def update(self):
        '''
        Can be useful if several people are working on the same file.
        All changed made to the database will be updated.

        Else, the database if the one at the time the file has been opened.
        '''
        old_star=self.current_star
        self.starobsfile.close()
        self.starobsfile=h5py.File(self.filename, 'a')
        if old_star <> None:
            self.use_star(old_star, display=False)

    def _create(self):
        '''Private function to create the file incase it doesn't exist already'''
        starobsfile=h5py.File(self.filename, 'w')
        starobsfile.attrs.create('Comments', "Database of star data, for the NuGrid collaboration.")
        starobsfile.close()

    def add_dataset(self, star='filename', author='Unknown', year='Unknown' , data='star_data.txt', set_default=True, activate_star=True, autoconfirm=False):
        '''
        Main function to add data.
        If no data is given an interface will be given to enter manually the data.
        The data can be given through a file, just specify the filename.

        Requires either a star name, or the star data file, from which python will try to extract the name.
        In order to keep a well-structured database, it is advised to add the name of the author and the year of publication.
        Parameters :

        star : str
            Name of the star. If no name is given, the name will be the file name given in data, without the extention.
        author : str, optional
            Name of the author of the publication the data comes from.
        year : str, optional
            Year of the publication the data comes from.
        data : str or list, optional
            If a string is given, python will extract the data from the file given in this strig.
            If a list is given, python will try to extract data from it, a star name must be provided.
            If nothing is given, and interface will be given to enter data manually, a star name must be provided.
        activate_star : bool, optional
            Permits to activate this star to work on it later.
        set_default : bool, optional
            Sets the new entry as the default dataset from which value will be taken in priority, incase one element is provided in several datasets.
            default is True
        autoconfirm : bool, optional
            Private argument, to avoid to enter a confirmation.
        '''
        if not activate_star:
            old_star=self.current_star
        if star=='filename':
            if(type(data)==list or data==None):
                print "You must enter a star name. Now returning none."
                return None
            elif(type(data)==str):
                star=data.split('/')[-1].split('.')[0]
        if star in self.starlist:
            print "Star found"
            self.use_star(star, display=False)
            dataset=str(len(self._currentstar.items())+1)
            if not self._compare_othersets(author, year, autoconfirm=autoconfirm):
                print "Dataset canceled, returning None."
                return None
        else:
            print "New star created"
            dataset='1'
        print "This is dataset number {0} for {1}.".format(dataset, star)
        data_prepared=_prepareData(data)
        if(data_prepared == None):
            print "The specified file doesn't exist or has ill defined data."
            return None
        if autoconfirm or _confirm_dataset(data_prepared, star, author, year):
            if dataset=='1':
                try:
                    self.starobsfile.create_group(star)
                    self.starlist.append(star)
                except:
                    print "A previous add_dataset attempt seems to have failed. You can check the integrity of the file via h5dump on Linux."
                self.use_star(star, display=False)
            self._currentstar[dataset]=data_prepared
            self.currentdataset=self._currentstar[dataset]
            self.currentdataset.attrs.create('Added on', time.strftime("%m/%d/%Y"))
            self.currentdataset.attrs.create('Author', author)
            self.currentdataset.attrs.create('Year', year)
            self.currentdataset.attrs.create('Z', '!0')
            self.currentdataset.attrs.create('XFe', '!1')
            self.currentdataset.attrs.create('Err', '!2')
            self.currentdataset.attrs.create('Err_Type', '!3')
            self.currentdataset.attrs.create('Upper', '!4')
            print "New dataset created."
            self.update()
            if set_default or dataset=='1':
                self.set_default(dataset, star=star, display=False)
            else:
                print "filling"
                self.set_default(dataset, star=star, z='fill', display=False)

        else:
            if not activate_star and old_star<> None:
                self.use_star(old_star, display=False)
            print "Dataset canceled, returning None."
            return None
        if not activate_star and old_star<> None:
            self.use_star(old_star, display=False)
        elif activate_star:
            self.use_star(star, display=False)

    def _compare_othersets(self, author, year, autoconfirm=False):
        '''Private method, that compares the entry the the existing ones, to avoid have twice the same dataset.'''
        found_similar=False
        for data_set in self._currentstar.values():
            if(data_set.attrs['Author']==author and data_set.attrs['Year']==year):
                found_similar=True
        if(found_similar):
            if not autoconfirm:
                print "A similar dataset has been found. Do you want to continue anyway ?"
                return _confirm(default='n')
            else:
                return False
        else:
            return True

    def datasets(self, star=None):
        '''
        Very short function, to have an easy way to display the diifferent datasets available for one star.
        Parameters:

        star : str, optional
            If no star is given, a star must have been activated before, through self.use_star.
        '''
        if star <> None:
            if star not in self.starlist:
                print  "Wrong star name (datasets)."
                return None
        if star == None and self.current_star==None:
            print "Currently no star selected. Use starobs.use_star(star) to use one."
            return None
        if star <> None:
            old_star=self.current_star
            self.use_star(star, display=False)
        print "Datasets availables : "
        for data_set in self._currentstar.items():
            print "#{0}".format(data_set[0])
            print "Default for {0}/{1}/{2} elements.".format(self._currentstar.attrs.values().count(data_set[0]), len(self.get('Z', star=self.current_star, dataset=data_set[0])), len(self._currentstar.attrs.values()))
            for key, value in zip(self._currentstar[data_set[0]].attrs.keys(), self._currentstar[data_set[0]].attrs.values()):
                if '!' not in value:
                    print "{0} : {1}".format(key, value)
            print ''
        if star<>None and old_star<> None:
            self.use_star(old_star)

    def attri_list(self, star=None, dataset=None, display=True):
        '''
        Function that displays the different attributes that can be read in self.get()
        Parameters:

        star : str, optional
            If no star is given, a star must have been activated before, through self.use_star.
        dataset : str, optional
            If no dataset is specified, the attributes of the first dataset will be displayed.
        display : bool, optional
            Private argument, for internal use.
        '''
        attributes=[]
        if star<>None:
            if star not in self.starlist:
                if display:
                    print "Wrong star name (attri_list)."
                return None
            else:
                old_star=self.current_star
                self.use_star(star, display=False)
        elif star==None:
            if self.current_star==None:
                if display:
                    print "Currently no star selected. Use starobs.use_star(star) to use one."
                return None
        if dataset==None:
            dataset='1'
        if display:
            print "Getting data from {0}, with dataset {1}".format(self.current_star, dataset)
        dataset=self._currentstar[dataset]
        for attri in dataset.attrs.keys():
            attributes+=[attri]
        if star<>None:
            self.use_star(old_star, display=False)
        return attributes

    def search(self, attri=None, value=None, star='current', dataset='default'):
        '''
        This function allows to search for an attribute, a value or an attribute with a particular value, in the entire database or  a part of it.
        starobs.search () works.

        Parameters:

        attri : str, optional
            A specific attribute you are looking for. The list of attributes for a specific star can be found with starobs.attri_list()
            It can be also an element symbol.
        value : str, float, optional
            A specific value you are looking for. It can be a string (author, year of publication) or a number (z, abundance, ...)
            Yet, it is HIGHLY advised to use the element symbol in attri instead of z here.
        star : str, list, optional
            The star(s) you want to include in that research. Default is 'all'. Use 'current' to use starobs.current_star or a specific star name.
        dataset : str, list, optional
            The dataset(s) you want to include in you research. 
            Default is 'default' to use the default one for every star included, or a specific dataset number or list.
            Use 'all' to search through all the different datasets for this star.
        '''
        old_star=self.current_star
        if attri<>None and type(attri)==str:
            try:
                utils.get_z_from_el(attri)
                isSymbol=True
            except:
                isSymbol=False
        else:
            isSymbol=False
        if star=='current':
            star=self.current_star
            if star==None:
                print "Currently no star selected. Use starobs.use_star(star) to use one."
                return None
        elif star=='all':
            star=self.starlist
        if type(star)==list:
            for starID in star:
                self.search(attri=attri, value=value, star=starID, dataset=dataset)
            if old_star<>None:
                self.use_star(old_star,  display=False)
            return None
        else:
            if star not in self.starlist:
                print "Wrong star name (search)."
                return None
        
        self.use_star(star, display=False)

        if dataset=='default':
            if isSymbol:
                if str(utils.get_z_from_el(attri)) in self._currentstar.attrs.keys():
                    dataset=self._currentstar.attrs[str(utils.get_z_from_el(attri))]
                else:
                    return None
            else:
                dataset=np.unique(self._currentstar.attrs.values()).tolist()
        elif dataset=='all':
            dataset=self._currentstar.keys()
        if type(dataset)==list:
            for datasetID in dataset:
                self.search(attri=attri, value=value, dataset=datasetID)
            if old_star<>None:
                self.use_star(old_star,  display=False)
            return None
        else:
            if dataset not in self._currentstar.keys():
                print "Wrong dataset ID."
                return None
        
        attris=self.attri_list(star=star, dataset=dataset, display=False)
        if attri <> None:
            if attri in attris or (isSymbol and str(utils.get_z_from_el(attri)) in self._currentstar.attrs.keys()):
                attris=[attri]
                out=True
            else:
                out=False
        else:
            out=True
        if out:
            for attriID in attris:
                values=self.get(attriID, dataset=dataset, display=False, oxy=None)
                if value == None:
                    print "{0} found in dataset {1} of {2}, with the value(s) :".format(attriID, dataset, star)
                    if type(values)==list and len(values)>1:
                        el=values[0]
                        z=values[1]
                        xfe=values[2]
                        err=values[3]
                        upper=values[5]
                        if len(values)>6:
                            oxy=values[6]
                        else:
                            oxy=0
                        toprint = "{0}".format(el)
                        if oxy <> 0:
                            toprint +=" {0}".format(oxy)
                        toprint+=" ({0}) : [X/Fe]={1}, sigma={2}".format(int(z), xfe, err)
                        if upper:
                            toprint+=", upper limit"
                        print toprint
                        print ""
                    else:
                        print values
                else:
                    if type(values)<>list and type(values)<> np.ndarray:
                        values=[values]
                    for pos in [a for a, b in enumerate(values) if b == value]:
                        print "{0} found for the attribute {1}, dataset {2} of star {3}".format(value, attriID, dataset, star)
                        if len(values)>1:
                            if not isSymbol:
                                if type(values)==np.ndarray:
                                    values=values.tolist()
                                values=self.get(utils.get_el_from_z(int(self.get('Z', dataset=dataset, display=False)[pos])), dataset=dataset, display=False, oxy=None)
                            el=values[0]
                            z=values[1]
                            xfe=values[2]
                            err=values[3]
                            upper=values[5]
                            if len(values)>6:
                                oxy=values[6]
                            else:
                                oxy=0
                            toprint = "{0}".format(el)
                            if oxy <> 0:
                                toprint +=" {0}".format(oxy)
                            toprint+=" ({0}) : [X/Fe]={1}, sigma={2}".format(int(z), xfe, err)
                            if upper:
                                toprint+=", upper limit"
                            print toprint
                        print ""


    def get(self, attri, star=None, dataset=None, display=True, oxy='sum',solar_norm=[None,None]):
        '''
        This functions permits to get data from a star.
        A star must have been activated before.
        Use starobs.use_star(star) for that.

        Parameter :

        attri : str
            Attribute you want to get or element symbol.
        star : str, optional
            If no star is given, a star must have been activated before, through self.use_star.
        dataset : str, optional
            If no dataset is specified, the attributes of the default elements will be returned.
        display : bool, optional
            Private argument, for internal use.
        oxy : str, int, optional
            Level of oxydation requested incase an element is requested.
            If an int is given, the requested level will be given.
            Default is 'sum' and will give the sum of all the ions, elements... It is the one needed incase of a plot.
            None will
	solar_norm : list, string, optional
	    Allows to convert abundances (in spectroscopic notation) normalized through a
	    particular solar abundance distribution into abundance normalized by another solar
	    distribution of choice. First entry in solar_norm is path to file hosting
            the original solar normalization mass fractions while the second entry points to the file
	    with the new solar normalization abundances of choice.
	    Both files need to be in NuGrids USEEPP format.
 
        '''
        relevant_attri_norm=['XFe','Err']


	#print 'get is executed ',attri,star,dataset
        dataset_ask=dataset
        if star <> None:    
            if star not in self.starlist:
                if display:
                    print "Wrong star name (get)."
                return None
            else:
                old_star=self.current_star
                self.use_star(star, display=False)
        elif star==None:
            if self.current_star <> None:
                if display:
                    print "Getting data from {0}".format(self.current_star)
            else:
                print "Currently no star selected. Use starobs.use_star(star) to use one."
                return None

	#here I select dataset, define dataset variable
	#print 'dataset ask ',dataset_ask
        if dataset_ask <> None:
            dataset=self._currentstar[dataset_ask]
        else:
	    #if not specified, take first entry, _e.g. _currentstar('1')
            dataset=self._currentstar[np.unique(self._currentstar.attrs.values())[0]]
        try:
            z=utils.get_z_from_el(attri)
        except:
            z=0
	# if attri is an element symbol charge Z return..
        if z<>0:
            if str(z) not in self._currentstar.attrs.keys():
                if display:
                    print "Element not available for this star"
                return None
            else:
                data_def = self._currentstar.attrs[str(z)]
            if z in dataset.value[0] and dataset_ask <> None:
                toreturn=[attri]
                toreturn+=[dataset.value[i][dataset.value[0].tolist().index(z)] for i in range(len(dataset.value))]
            else:
                if display:
                    if dataset_ask <> None:
                        print "Element not available in this dataset"
                    print "In the default dataset :"
                toreturn=[attri]
                dataset=self._currentstar[str(data_def)]
                for i in range(len(dataset.value)):
                    toreturn+=[float(dataset.value[i][dataset.value[0].tolist().index(z)])]
            return toreturn
	#if attr in dataset.attrs keys [u'Added on', u'Author', u'Year', u'Z', u'XFe', u'Err', u'Err_Type', u'Upper']
	#attri_value are attrs values e.g.  ['02/18/2016', 'Andreas Koch', '2016', '!0', '!1', '!2', '!3', '!4']
        if attri in dataset.attrs:
            attri_value=dataset.attrs[attri]
        else:
            if display:
                print "This attribute doesn't exist. Use starobs.attrilist() to get the different attributes"
            return None
	## for attribute values with ! and their keys (e.g. 'XFe') which have corresponding values in dataset.value
	# dataset.value values are returned.
        if attri_value[0]=='!':
            if dataset_ask<> None:
		#take value given by attri 
		toreturn=dataset.value[int(attri_value[1])]
		if ((not solar_norm[0]==None) and (not solar_norm[1]==None)):
		   if attri in relevant_attri_norm:
		       print 'renormalize values!'
		       #to get charge numbers from current dataset
		       Z_values=dataset.value[int(dataset.attrs['Z'][1])]
                       converted_data = ds.convert_solar_norm(Z_values,toreturn,solar_norm)
		       print 'original: ',toreturn
		       print 'converted:',converted_data
                       toreturn = converted_data
                return toreturn
	    # if no dataset is specified
            else:
		#data includes all values of  u'Z', u'XFe', u'Err', u'Err_Type', u'Upper'
		data=[[] for i in range(len(dataset.value))]
		# loop over charges
                for Z_def in np.sort(self._currentstar.attrs.keys()):
		    #from my current understanding the dataset does not change through the loop
		    # self._currentstar.attrs[Z_def] is always '1'
                    dataset=self._currentstar[self._currentstar.attrs[Z_def]]
		    #loop over values
                    for i in range(len(dataset.value)):
			data[i]+=[float(dataset.value[i][dataset.value[0].tolist().index(float(Z_def))])]
		return_data=data[int(attri_value[1])]
		if ((not solar_norm[0]==None) and (not solar_norm[1]==None)):
		     if attri in relevant_attri_norm:
			#order of Z as looped over above
			all_Z_def = np.sort(self._currentstar.attrs.keys())	
			converted_data = ds.convert_solar_norm(all_Z_def,return_data,solar_norm)
			#print 'return data', return_data
			#print 'converged data',converted_data
			return_data = converted_data
            return return_data
	# in case attri are values not with !, the corresponding attri values can be returned 
        else:
            data=[]
            if dataset_ask<>None:
                dataset=self._currentstar[dataset_ask]
		toreturn=dataset.attrs[attri] #return value of attri
		'''
		if ((not solar_norm[0]==None) and (not solar_norm[1]==None)):
		    if attri in relevant_attri_norm:
		       attri_valueZ=dataset.attrs['Z']
		       converted_data = ds.convert_solar_norm(dataset.value[int(attri_valueZ)],toreturn,solar_norm)	
		       toreturn = converted_data
		'''
                return toreturn
            else:
		print '######################IAM EXE'
                for datasetID in np.unique(self._currentstar.attrs.values()):
                    dataset=self._currentstar[datasetID]
		    dataID = dataset.attrs[attri]
		    '''
		    if ((not solar_norm[0]==None) and (not solar_norm[1]==None)):
		       print 'attri ',attri,relevant_attri_norm
		       if attri in relevant_attri_norm:
			  attri_valueZ=dataset.attrs['Z']
			  converted_data = ds.convert_solar_norm(dataset.value[int(attri_valueZ)],dataID,solar_norm)
			  converted_data = ds.convert_solar_norm(dataset.attrs['Z'],dataID,solar_norm)
		          print 'old: ',dataset.attrs[attri]  
			  print 'converted: ',converted_data
			  dataID = converted_data
		    '''
                    data+=[dataID]
                return data #CR moved out of foor loop
        self.use_star(old_star, display=False)

    def use_star(self, star, display=True):
        '''
        Function to study a star in particular.
        It will show you  the different datasets available.

        Parameters :

        star : str
            Name of the star you want to activate.
            You can obtain the list of stars with starobs.starlist.
        display : bool, optional
            Private argument, for internal use.
        '''
        if star in self.starlist:
            self.current_star=star
            self._currentstar=self._currentstar=self.starobsfile[star]
            if display:
                self.datasets()
        else:
            print "Wrong star name (use_star)."

    def merge(self, file_to_add, set_default='old'):
        '''
        Function that merges two different database. 
        If the two databases have been built properly with author names and years of publication, it will avoid conflicts.
        This is the only way to share the database, svn (or svn-like) protocoles won't work using the h5 file directly.

        Parameters:
        filename : str
            Filename of the database to include in the current one.
        set_default : str, optional
            Sets the priority to a dataset or to another incase the same star is present in both files. 
            Default is 'old', and keeps the current priority set.
            'new' sets the priority to the datasets added.
            'recent' will set the priority to the most recent publication
        '''
        old_nbr_stars=len(self.starlist)
        old_nbr_datasets=0
        for star in self.starlist:
            self.use_star(star, display=False)
            old_nbr_datasets+=len(self._currentstar.items())

        starobsfile_to_add=starobs(filename=file_to_add)
        stars_to_add=starobsfile_to_add.starlist
        for starname in stars_to_add:
            starobsfile_to_add.use_star(starname, display=False)
            for datasetID in starobsfile_to_add._currentstar.items():
                dataset=datasetID[0]
                author=starobsfile_to_add.get('Author', starname, dataset, display=False)
                year=starobsfile_to_add.get('Year', starname, dataset, display=False)
                z=starobsfile_to_add.get('Z', starname, dataset, display=False)
                xfe=starobsfile_to_add.get('XFe', starname, dataset, display=False)
                err=starobsfile_to_add.get('Err', starname, dataset, display=False)
                err_type=starobsfile_to_add.get('Err_Type', starname, dataset, display=False)
                upper=starobsfile_to_add.get('Upper', starname, dataset, display=False)
                data=[z, xfe, err, err_type, upper]
                if starname not in self.starlist:
                    sd=True
                else:
                    if set_default=='old':
                        sd=False
                    elif set_default=='new':
                        cs=starobsfile_to_add[star]
                        if cs.attrs['default']==dataset:
                            sd=True
                        
                self.add_dataset(starname, author, year, data, set_default=sd, activate_star=False, autoconfirm=True)
        
        new_nbr_stars=len(self.starlist)
        new_nbr_datasets=0
        for star in self.starlist:
            self.use_star(star, display=False)
            new_nbr_datasets+=len(self._currentstar.items())
        if set_default=='recent':
            for star in self.starlist:
                self.set_default('recent')
        elif set_default=='new':
            for star in self.starlist:
                self.set_default('last')
        print "{0} new stars, now containing {1} stars.".format((new_nbr_stars-old_nbr_stars), new_nbr_stars)
        print "{0} datasets added, up to a total of {1} datasets.".format((new_nbr_datasets-old_nbr_datasets), new_nbr_datasets)

    def set_attri(self, dataset, star=None, author=None, year=None):
        '''
        Method to set the author and the year of the publication the star data is extracted from.
        dataset : str
            ID of the dataset to change.
        star : str, optional
            If no star is given, a star must have been activated before, through self.use_star.
        author : str, optional
            If no author is given, and an author is already set, the previous entry will be kept.
        year : str, optional
            If no year is given, and a year is already set, the previous entry will be kept.
        '''
        if star <> None:   
            if star not in self.starlist:
                if display:
                    print "Wrong star name (set_attri)."
                return None
            else:
                old_star=self.current_star
                self.use_star(star, display=False)
        elif star==None:
            if self.current_star <> None:
                star=self.current_star
                print "Changing attribute of {0}.".format(star)
            else:
                print "Currently no star selected. Use starobs.use_star(star) to use one."
                return None
        self.currentdataset=self._currentstar[dataset]
        if author<>None:
            self.currentdataset.attrs.modify('Author', author)
        if year<>None:
            self.currentdataset.attrs.modify('Year', year)
        if old_star<>None:
            self.use_star(old_star)

    def set_default(self, dataset, star=None, z=None, display=True):
        '''
        Method to change the default dataset for a certain star.

        Parameters:
        dataset : str
            It can be either the ID of a dataset, or one of the following.
            'last' sets default to the last dataset added in the Database
            'recent' sets default to the most recent publication.
        star : str, optional
            If no star is given, a star must have been activated before, through self.use_star.
        z : int, str, optional
            Element number or symbol of the specific element one wants to change. If None, all the elements of this dataset will be set as default.
            Private argument 'fill' will just set the default position for new elements added from a dataset not set as default.
        display : bool, optional
            Private argument, for internal use.
        '''
        if z<>None:
            if type(z)==str:
                try:
                    z=int(utils.get_z_from_el(z))
                except:
                    print "Wrong element symbol. Returning None."
                    return None
        if star <> None:   
            if star not in self.starlist:
                if display:
                    print "Wrong star name (set_default)."
                return None
            else:
                old_star=self.current_star
                self.use_star(star, display=False)
        elif star==None:
            if self.current_star <> None:
                if display:
                    print "Getting data from {0}".format(self.current_star)
            else:
                print "Currently no star selected. Use starobs.use_star(star) to use one."
                return None
        if dataset=='recent' or dataset=='last':
                years=[]
                added=[]
                datasets=self._currentstar.items()
                for dataset in datasets:
                    years+=[int(self.get('Year', dataset=dataset[0], display=False))]
                    added+=[int(self.get('Added on', dataset=dataset[0], display=False))]
                if dataset=='recent':
                    self.set_default(datasets[years.index(max(years))][0])
                if dataset=='last':
                    self.set_default(datasets[added.index(max(added))][0])
        else:
            if z==None:
                for z_in in self.get('Z', dataset=dataset, display=False):
                    self.set_default(dataset=dataset, z=z_in, display=False)
            elif z=='fill':
                for z_in in self.get('Z', dataset=dataset, display=False):
                    if str(int(z_in)) not in self._currentstar.attrs.keys():
                        self.set_default(dataset=dataset, z=z_in, display=False)
            else:
                if z in self.get('Z', dataset=dataset, display=False):
                    self._currentstar.attrs[str(int(z))]=dataset
                    self.update()
                else:
                    if display:
                        print "This element doesn't exists for this dataset. Retuning none."
                    if old_star<>None:
                        self.use_star(old_star)
                    return None
        try:
            if old_star<>None:
                self.use_star(old_star)
        except:
            None

    def export(self, filename, star=None, zrange=None, row=None, upper_symbol='$<$', comments='', separator='&',  end_of_line='\\\\', spaces=' ', autoconfirm=False):
        '''
        Simple method that exports star data in an ascii file, in order of example to be included in a tex file.
        
        Parameters :
            filename : str
                Name of the file to export
            columns : list
                List of the attributes to be printed, side by side. 'Z' and 'El' will be added automatically.
            star : str, list, optional
                Star to include in the file.
                If no star is specified, starobs.current_star will be default.
            zrange : list, optional
                Limit if z values of the elements to include in the file
            row : str, list, optional
                Additional row under the star names
            upper_symbol : str, optional
                Symbol displayed in fron of the abundances (if the abundances are requested) to signify an upper limit.
                Default is '<'
            separator : str, optional
                Caracters to be used betwee each entry.
                Default is '&'
            end_of_line : str, optional
                Caracters to be used at the end of each line.
                Default is '\\'
            spaces : str, optional
                spaces one want to put between each separator and entry.
        '''
        form='%1.2F'
        if star==None:
            if self.current_star<>None:
                stars=[self.current_star]
            else:
                print "No star specified or activated as default. Returning None."
                return None
        elif type(star)==str:
            stars=[star]
        elif type(star)==list:
            stars=star
        else:
            print "Wrong star entry. Returning None."
            return None
        if os.path.isfile(filename) and not autoconfirm:
            print "This file already exists. Overwrite ?"
            if not _confirm(default='y'):
                return None
        output=open(filename, "w")
        if comments<>'':
            output.write(comments)
            output.write('\n')
        data=[]

        for star in stars:
            data_tmp= [self.get('Z', star=star)]
            for attri in ['XFe', 'Err', 'Upper']:
                data_tmp+=[self.get(attri, star=star)]
            if zrange==None:
                data+=[data_tmp]
            else:
                data+=[[[] for i in range(len(data_tmp))]]
                for i in range(len(data_tmp[0])):
                    if zrange[0]<=data_tmp[0][i]<=zrange[1]:
                        for j in range(len(data_tmp)):
                            data[stars.index(star)][j]+=[data_tmp[j][i]]
        output.write('\\tablecolumns{'+str(len(stars))+'}\n')
        output.write('\\tablehead{\n')

        line=spaces+spaces
        line+=separator+spaces
        line+=spaces+spaces
        line+=separator+spaces
        for star in stars:
            line+='\multicolumn{2}{c}{'+star+'}'+spaces
            line+=separator+spaces
            line+=separator+spaces
        line=line[:-(len(spaces)+len(separator))]
        line=line[:-(len(spaces)+len(separator))]
        line+=' \\\\'
        output.write(line)
        output.write('\n')

        for star in range(len(stars)):
            indic="{0}-{1}".format(3*(star+1), 3*(star+1)+1)
            output.write('\cline{'+indic+'}\n')
        line='Elem.'+spaces
        line+=separator+spaces
        line+='Z'+spaces
        line+=separator+spaces
        for star in stars:
            line+='[X/Fe]'+spaces
            line+=separator+spaces
            line+='$\sigma$'+spaces
            line+=separator+spaces
            line+=separator+spaces
        line=line[:-(len(spaces)+len(separator))]
        line=line[:-(len(spaces)+len(separator))]
        line+=' }'
        output.write(line)
        output.write('\n')
        
        output.write('\startdata\n')
        
        for z in range(1, 93):
            appears=False
            line=utils.get_el_from_z(str(int(z)))+' '
            line+=separator+spaces
            line+=str(z)+spaces
            line+=separator+spaces
            for star in range(len(stars)):
                for i in range(2):
                    if z in data[star][0]:
                        appears=True
                        if i==0:
                            if data[star][3][data[star][0].index(z)]:
                                line+=upper_symbol
                        line+=form%data[star][i+1][data[star][0].index(z)]+spaces
                    else:
                        line+='\dots'+spaces
                    line+=separator+spaces
                line+=separator+spaces
            line=line[:-(len(spaces)+len(separator))]
            line=line[:-(len(spaces)+len(separator))]
            line+=spaces+end_of_line
            if appears:
                output.write(line)
                output.write('\n')
        output.close()

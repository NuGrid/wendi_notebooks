import data_plot as dp
import ppn
import sys
import os
import astronomy as ast
import utils
import glob
import multiprocessing as mp
from scipy import integrate
from os import chdir
import numpy as np
from matplotlib.pyplot import *
import matplotlib
import re
import math

'''
This module is a super class that works with the starobs database module.
It will most likely not work alone !
'''

import errno

def convert_solar_norm(Z,XFe,solar_norm):

    '''
    Takes abundance distribution in spectroscopic notation
    and converts to a new solar normalization.

    Z : list of charged numbers
    XFe : list of [X_i/Fe] corresponding to Z
    solar_norm : list of two entries, both path to NuGrids USEEPP files
                 first entry represents the original normalization while
		 the latter represents the solar abundance of choice

    '''
    Zint=[]
    for k in range(len(Z)):
	Zint.append(int(Z[k]))
    Z=Zint + [26]
    solar_old=get_solar_abu(el_Z=Z,solar_norm_file=solar_norm[0])
    solar_new=get_solar_abu(el_Z=Z,solar_norm_file=solar_norm[1])
    
    solar_old_Fe=solar_old[-1]
    solar_new_Fe=solar_new[-1]
    #print 'solar Fe',solar_old_Fe,solar_new_Fe	
    XFe_converted=[]
    for k in range(len(XFe)):
	#check for nan values
	if math.isnan(XFe[k]):
		#print 'skip ',k,XFe[k]
		XFe_converted.append(XFe[k])
		continue	
	offset=np.log10(solar_old[k]/solar_old_Fe) - np.log10(solar_new[k]/solar_new_Fe)
	XFe_converted.append( XFe[k] + offset )

    return XFe_converted

#CR added function
def get_solar_abu(el_name=None,el_Z=None,solar_norm_file=None):
	'''
	   Get solar system abundance mass fraction
	   for elements.

	   elements: List of form ['H','He']
	   Xsolarfile: path to NuGrid initial abundance files (USEEPP)

	'''
	if (not el_Z == None) and (el_name == None):
		#get names from el_Z
	   elements=[]
	   for k in range(len(el_Z)):
		elements.append(utils.get_el_from_z(str(el_Z[k]))) 
	else:
		elements=el_name
        X_sol=np.zeros(len(elements))
	solardata=utils.iniabu(solar_norm_file)
        for name in solardata.habu.keys():
	   name1=name.replace(" ","")
	   match = re.match(r"([a-z]+)([0-9]+)",name1, re.I)
	   if match:
    		items = match.groups()
	   else:
		print 'problem with ',name1
	   ele=items[0].capitalize() #form such as Fe
	   for k in range(len(elements)):
		if elements[k] == ele:
	        	X_sol[k] += solardata.habu[name] 

	return X_sol

#CR: added special handling of queue.get to prevent error
def my_queue_get(queue, block=True, timeout=None):
    while True:
        try:
            return queue.get(block, timeout)
        except IOError, e:
            if e.errno != errno.EINTR:
                raise

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

def _getData(p, thread, out_q, advancement_q, mode, clist, zlist, star_data=None, sim_0=None, Xref=None):
    '''
    Private method to get decayed ratios from a ppn simulation, with multiprocessing. Only way to efficiently get data in a reasonable time
    mode defines the purpose of it. If a new function uses it, one has to put the new arguments with default None !
    '''
    if mode=='ppn_compare':
        data=[]
    if mode=='autofit' or mode=='element_chart':
        sys.stdout = open(os.devnull, 'w')
        outdict = {}
        cycledict={}
    for cycle in clist:
        p.get(cycle, decayed=True)
        z_el=np.unique(p.z_iso_to_plot)

        a_el=[]; el_abu=[]; el_abu_hash={}

        if mode=='element_chart':
            out=[]
            for i in range(len(zlist)):
                out+=[p.abunds[np.where(p.el_iso_to_plot==utils.get_el_from_z(str(zlist[i])))[0].tolist()].sum()]
        else:
            for z in zlist:
                el=p.el_iso_to_plot[np.where(p.z_iso_to_plot==z)[0].tolist()[0]]
                X_el=p.abunds[np.where(p.el_iso_to_plot==el)[0].tolist()].sum()
                a_el.append(p.a_iso_to_plot[np.where(p.z_iso_to_plot==z)[0].tolist()[0]])
                el_abu.append(X_el)
                el_abu_hash[el]=X_el
        if mode=='autofit':
            out=0
            el_abu /= np.array(sim_0[0])
            XFe_offset=star_data[1][star_data[0].index(utils.get_z_from_el(Xref))]-np.log10(el_abu_hash[Xref]/sim_0[1][Xref])
            
            sum_w=0

            for z in star_data[0]:
                try:
                    weight_z=star_data[3][z]
                except:
                    weight_z=1
                sum_w+=weight_z

                if(star_data[1][star_data[0].index(z)]-star_data[2][star_data[0].index(z)])< (np.log10(el_abu[zlist.index(z)])+XFe_offset) < (star_data[1][star_data[0].index(z)]+star_data[2][star_data[0].index(z)]) and (not star_data[4][star_data[0].index(z)] or np.log10(el_abu[zlist.index(z)])+XFe_offset < star_data[1][star_data[0].index(z)]):
                    out+=0
                else :
                    if star_data[4][star_data[0].index(z)]:
                        if np.log10(el_abu[zlist.index(z)])+XFe_offset > star_data[1][star_data[0].index(z)]:
                            out+=10000*weight_z
                    out+=weight_z*min(abs(np.log10(el_abu[zlist.index(z)])+XFe_offset-star_data[1][star_data[0].index(z)]-star_data[2][star_data[0].index(z)]), abs(np.log10(el_abu[zlist.index(z)])+XFe_offset-star_data[1][star_data[0].index(z)]+star_data[2][star_data[0].index(z)]))

            out/=sum_w

        if mode=='autofit' or mode=='element_chart':
            outdict[cycle]=out
            cycledict[cycle]=cycle
            advancement_q.put(cycledict)
        
        elif mode=='ppn_compare':
            data+=[[zlist, el_abu, el_abu_hash]]

    if mode=='autofit' or mode=='element_chart':
        out_q.put(outdict)
    if mode=='ppn_compare':
        return data

class Plots():
    def element_charts(self, simu_folder=None, Xsolarfile=None, crange=None, star='current', nbr_of_threads=0):
        '''
        Method that plots data in the following way : several pannels of elements, with [A/B] vs [C/D] for each pannel.
        One can add stars on it, incase the star countains the considered elements.
        One can also plot mean values of several elements : [{A,B}/C] vs [D/{E,F,G}]

        The elements are asked after calling the function, or can be entered while calling.

        Since getting the data is long, all the data is saved temporary. A second plot should then be much faster.

        Parameters :
        simu_folder : str, list, optional
            path or list of paths of the folder the ppn simulation(s) to plot.
        Xsolarfile : str, optional
            The solar abundances are requested here to normalize things properly, incase a simulation is plotted too.
        crange : list of int, list of lists, optional
            [cmin, cmax] are the limit cycles of the simulation for this plot.
            If cmin and cmax are ints, the sames cycles are considered for every simulation.
            If they are lists, they must be lists of the same lengh as the number of simulation, and can have particular cycles.
            Ex : crange = [[250, 300, 310], 935]
            If cmax is not specified, the last cycle will be considered.
        star : str, list, optional.
            Star(s) to try to plot on these charts.
            default is 'current',  it will plot the activated star from the database if these is an activated one.
            None will plot only the simulation.
            A string will specify one star, a list or strings will specify several stars
            'all' will plot all stars from the database.
        nbr_of_threads : int, optional
            Number of threads you want to imply, incase a simulation is plotted. Must fit you cpu number of cores. First check the cpu use with top.
            Default is None, python will use as many threads as cores present in this cpu.
        '''
        if simu_folder==None:
            simu_folder=[]
        elif type(simu_folder)==str:
            simu_folder=[simu_folder]
        if crange<>None and len(simu_folder)<>0:
            if type(crange[0])==list and len(crange[0])<>len(simu_folder):
                print "Not the same lengh of cmin and simulations folders. Returning None."
                return None
            if len(crange)==2 and type(crange[1])==list and len(crange[1])<>len(simu_folder):
                print "Not the same lengh of cmin and simulations folders. Returning None."
                return None
        
        if nbr_of_threads==0:
            nbr_of_threads=mp.cpu_count()
        Cx=[[],[]]
        Cy=[[],[]]
        entering = True
        chart=0
        legendmarkers=['o','s','p','*','h','+','x','D','H']
        corres=[Cx, Cy]
        
        try:
            self._alldatasim<>None
            print "Simulations data found in memory"
        except:
            self._alldatasim={}

        print "Enter your charts parameters here (enter !h for help) :"
        newchart=True
        while entering:
            if (chart%4) == 0:
                if chart<>0 and newchart:
                    print ""
                    print "Chart #{0} : [{1}/{2}] vs [{3}/{4}].".format((chart/4), Cy[0][(chart/4)-1], Cy[1][(chart/4)-1], Cx[0][(chart/4)-1], Cx[1][(chart/4)-1])
                    print "======"
                    print ""
                if newchart:
                    print "Chart #{0}".format(chart/4+1)
                print "X-axis :"
            if chart%4 == 1:
                print "Over :"
            if chart%4 == 2:
                print "Y-axis :"
            if chart%4 == 3:
                print "Over :"
            entry=raw_input()
            newchart=True
            if entry=='':
                print "Stop now ?"
                if _confirm(default='y'):
                    entering=False
            elif entry[0]=='!' and len(entry)>1:
                if entry[1]=='h':
                    print "Enter the symbol of the element,  eventually separeted by ',' for a mean value."
                    print "Press ENTER between each entry."
                    print ""
                    print "Enter :"
                    print "'!c' to correct the previous element,"
                    print "'!n' to change the entire previous chart,"
                    print "'!e' or press ENTER with an empty element entry to end this prompt."
                    print "" 
                if entry[1]=='c':
                    if chart <> 0:
                        chart-=1
                        corres[(chart%4)/2][(chart%4)%2]=corres[(chart%4)/2][(chart%4)%2][:-1]
                elif entry[1]=='n' and chart<>0:
                    chart-=1
                    corres[(chart%4)/2][(chart%4)%2]=corres[(chart%4)/2][(chart%4)%2][:-1]
                    while chart%4<>0:
                        chart-=1
                        corres[(chart%4)/2][(chart%4)%2]=corres[(chart%4)/2][(chart%4)%2][:-1]
                elif entry[1]=='e':
                    if chart/4==0:
                        print "Less than one chart filled, returning none."
                        return None
                    else:
                        entering = False
            else:
                entries=entry.split(',')
                try:
                    for el in entries:
                        z=utils.get_z_from_el(el)
                    corres[(chart%4)/2][(chart%4)%2]+=[entries]
                    chart+=1
                except:
                    print "Wrong entry"
                    newchart=False

        while chart%4<>0:
            chart-=1
            corres[(chart%4)/2][(chart%4)%2]=corres[(chart%4)/2][(chart%4)%2][:-1]
        
        print "Starting the plot..."

        stdout_screen=sys.stdout
        if len(simu_folder)>0:
            sys.stdout = open(os.devnull, 'w')

        el_list=[]
        pannel_nbr=len(Cx[0])

        CxR=[[[] for a in range(2)] for b in range(pannel_nbr)]
        CyR=[[[] for a in range(2)] for b in range(pannel_nbr)]

        for C in range(len(corres)):
            for topbot in range(len(corres[C])):
                 for pannel in range(len(corres[C][topbot])):
                     for el in corres[C][topbot][pannel]:
                         if el not in el_list:
                             el_list+=[el]
                         if C==0:
                             CxR[pannel][topbot]+=[el_list.index(el)]
                         if C==1:
                             CyR[pannel][topbot]+=[el_list.index(el)]
        print CxR
        print CyR
        if len(simu_folder)>0:
            EL_LIST=[]
            for el in el_list:
                if (len(el) == 1):
                    EL_LIST += [[el[0]+' ']]
                else:
                    EL_LIST += [[el[0]+el[1].upper()]]
            column=ppn.xtime(simu_folder[0]).cols
            for colnbr in range(len(column)):
                for nEL in range(len(EL_LIST)):
                    if(EL_LIST[nEL][0] == column[colnbr][:2]):
                        EL_LIST[nEL] += [column[colnbr]]
            for nEL in range(len(EL_LIST)):
                if(EL_LIST[nEL][0] == 'H '):
                    EL_LIST[nEL] += ['PROT']

            Xsolar=utils.iniabu(Xsolarfile)
            elementref=np.zeros((len(EL_LIST)))
            speciesIDchanged=''
            for elementnbr in range(len(EL_LIST)):
                for species in EL_LIST[elementnbr][1:]:
                    if (species == 'PROT'):
                        speciesIDchanged='H-1'
                    else:
                        speciesIDchanged += species[0]
                        if (species[1] != ' '):
                            speciesIDchanged += species[1].lower()
                        speciesIDchanged+='-'
                        if (species[2] != ' '):
                            speciesIDchanged += species[2]
                        if (species[3] != ' '):
                            speciesIDchanged += species[3]
                        speciesIDchanged += species[4]
                    for name in Xsolar.habu:
                        if((name[:2].upper()+name[2:]) == species):
                            elementref[elementnbr] += [Xsolar.iso_abundance(speciesIDchanged)]
                    speciesIDchanged=''
        if type(star)==str:
            if star=='current' and self._out_class_plot.current_star<>None:
                stars=[self._out_class_plot.current_star]
            elif star=='current' and self._out_class_plot.current_star==None:
                stars=[]
            elif star=='all':
                stars=self._out_class_plot.starlist
            else:
                stars=[star]
        elif star==None:
            stars=[]
        else:
            stars=star
        stardata=[[[] for a in range(len(el_list))] for b in range(len(stars))]
        can_plot=[[[] for a in range(pannel_nbr)] for b in range(len(stars))]
        for star in stars:
            for el in el_list:
                if el<>'Fe':
                    tmp=self.get(el, star_plot=star, display_plot=False)
                else:
                    tmp=['Fe', 26, 0, 0, 0]
                if tmp<>None:
                    stardata[stars.index(star)][el_list.index(el)]=tmp[2:4]
            for pannel in range(pannel_nbr):
                can_plot_tmp=[]
                for axis in CxR[pannel]:
                    for el in axis:
                        if len(stardata[stars.index(star)][el])==0:
                            can_plot_tmp+=[el]
                for axis in CyR[pannel]:
                    for el in axis:
                        if len(stardata[stars.index(star)][el])==0:
                            can_plot_tmp+=[el]
                can_plot[stars.index(star)][pannel]=can_plot_tmp
        taus=[]
        tmin=[]
        tmax=[]

        for subfolder in simu_folder:
            X=ppn.xtime(subfolder)
            if crange<>None:
                if type(crange[0])==list:
                    cmin=crange[0][simu_folder.index(subfolder)]
                else:
                    cmin=crange[0]
                if len(crange)==2:
                    if type(crange[1])==list:
                        cmax=crange[1][simu_folder.index(subfolder)]+1
                    else:
                        cmax=crange[1]+1
                else:
                    cmax=X.get('cycle')[-1]+1
            else:
                cmax=X.get('cycle')[-1]+1
                cmin=0
            #Getting the neutron exposure
            xn=X.get("NEUT")
            age=X.get("time")
            rho=X.get("rho")
            v_therm=3.e8
            oneyear=utils.constants.one_year
            #Nn=1.e-19*ast.avogadro_constant*1000.
            #Nn/1.e7
            Nn=xn*ast.avogadro_constant*rho
            age_s=age*oneyear
            tau=integrate.cumtrapz(Nn*v_therm,age_s)
            tau=np.log10(tau*1.e-27)
            tau=tau[cmin:cmax]
            Nn=np.log10(Nn)
            Nn=Nn[cmin:cmax]
            taus+=[tau]
            tmin+=[min(tau)]
            tmax+=[max(tau)]
            sys.stdout = stdout_screen
            
        scalemin=min(tmin)
        scalemax=max(tmax)

        if len(simu_folder)>0:
            print "Neutron exposures loaded."
            print ""
            plotcolorbar=True

        for subfolder in simu_folder:
            print 'Doing : '+subfolder
            sys.stdout = open(os.devnull, 'w')
            P=ppn.abu_vector(subfolder)
            sys.stdout = stdout_screen
            print "ppn instance loaded"
            print ""
            if crange<>None:
                if type(crange[0])==list:
                    cmin=crange[0][simu_folder.index(subfolder)]
                else:
                    cmin=crange[0]
                if len(crange)==2:
                    if type(crange[1])==list:
                        cmax=crange[1][simu_folder.index(subfolder)]+1
                    else:
                        cmax=crange[1]+1
                else:
                    cmax=len(P.files)
            else:
                cmax=len(P.files)
                cmin=0
            
            if subfolder not in self._alldatasim or cmin not in self._alldatasim[subfolder] or (cmax-1) not in self._alldatasim[subfolder]:
                if subfolder in self._alldatasim:
                    if cmin in self._alldatasim[subfolder]:
                        cmin=max(self._alldatasim[subfolder].keys())
                    if cmax in self._alldatasim[subfolder]:
                        cmax=min(self._alldatasim[subfolder].keys())
                threads=[]
                out_q=mp.Queue()
                advancement_q=mp.Queue()
                resultdict={}
                cycledict={}
                zlist=np.arange(1, 92)
                
                for thread in range(nbr_of_threads):
                    clist=range(cmin+((cmax-cmin)*thread)/nbr_of_threads, cmin+((cmax-cmin)*(thread+1))/nbr_of_threads)
                    print "Launching thread {0} for cycles {1} to {2}.".format(thread , cmin+((cmax-cmin)*thread)/nbr_of_threads , cmin+((cmax-cmin)*(thread+1))/nbr_of_threads-1)
                    threads.append(mp.Process( target=_getData, args=(P, thread, out_q, advancement_q, 'element_chart', clist, zlist)))
                    threads[-1].start()
                old=-1
                while(len(cycledict)<(cmax-cmin)):
                    if(len(cycledict)*100/(cmax-cmin) <> old):
                        print '{0}%\r'.format(str(len(cycledict)*100/(cmax-cmin))),
                        sys.stdout.flush()
                        old=len(cycledict)*100/(cmax-cmin)
                    cycledict.update(advancement_q.get())
                print "100%"
                for i in range(nbr_of_threads):
                    resultdict.update(out_q.get())
                print "Dictionary updated"
                
                for i in threads:
                    i.join()
                print "Threads ended."
                
                if subfolder not in self._alldatasim:
                    self._alldatasim[subfolder]=resultdict
                else:
                    self._alldatasim[subfolder].update(resultdict)
            else:
                print "Data already loaded found in memory."

            elementData=np.zeros((len(EL_LIST),(cmax-cmin)))
            for cycle in self._alldatasim[subfolder]:
                for el in range(len(el_list)):
                    elementData[el][cycle-cmin]=self._alldatasim[subfolder][cycle][(int(utils.get_z_from_el(el_list[el]))-1)]

            for subp in range(pannel_nbr):
                if pannel_nbr%2==1 and subp==(pannel_nbr-1):
                    subplotnbr=(100*(subp/2+1))+10+(subp/2+1)
                    subplot(subplotnbr)
                else:
                    subplotnbr=((pannel_nbr-1)/2+1)*100+20+subp+1
                    subplot(subplotnbr)

                up=0
                down=0
                for element in CxR[subp][0]:
                    up+=(elementData[element])/elementref[element]
                for element in CxR[subp][1]:
                    down+=(elementData[element])/elementref[element]
                up/=len(CxR[subp][0])
                down/=len(CxR[subp][1])
                xdata=np.log10(up/down)
                       
                up=0
                down=0
                for element in CyR[subp][0]:
                    up+=(elementData[element])/elementref[element]
                for element in CyR[subp][1]:
                    down+=(elementData[element])/elementref[element]
                up/=len(CyR[subp][0])
                down/=len(CyR[subp][1])
                ydata=np.log10(up/down)

                tau=np.array(taus[simu_folder.index(subfolder)])

                points=np.array([xdata,ydata]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = matplotlib.collections.LineCollection(segments, cmap=get_cmap('jet'), norm=Normalize(scalemin, scalemax))
                lc.set_linewidth(2)
                lc.set_array(tau)
                gca().add_collection(lc)

                #This is a cheat, to normalize the colorbar (python not made fot this...)
                while len(tau)<len(xdata):
                    xdata=xdata[:-1]
                    ydata=ydata[:-1]
                if plotcolorbar:
                    scatter(xdata, ydata, c=tau, s=0, cmap='jet', vmin=scalemin,vmax=scalemax)
                    colorbar()
                plot(xdata,ydata,linestyle='o',marker=legendmarkers[simu_folder.index(subfolder)],markevery=30)
                plotcolorbar=False

                text(xdata[0]-abs(xdata[0]-xdata[1]),ydata[0]+abs(ydata[0]-ydata[1]),str(cmin))
                text(xdata[-1]-abs(xdata[-1]-xdata[-2]),ydata[-1]+abs(ydata[-1]-ydata[-2]),str(cmax))

        for subp in range(pannel_nbr):
            if pannel_nbr%2==1 and subp==(pannel_nbr-1):                    
                subplotnbr=(100*(subp/2+1))+10+(subp/2+1)
                subplot(subplotnbr)
            else:
                subplotnbr=((pannel_nbr-1)/2+1)*100+20+subp+1
                subplot(subplotnbr)
            for star_nbr in range(len(stars)):
                if(len(can_plot[star_nbr][subp])==0):
                    up=0
                    down=0
                    for element in CxR[subp][0]:
                        up+=(stardata[star_nbr][element][0])
                    for element in CxR[subp][1]:
                        down+=(stardata[star_nbr][element][0])
                    up/=len(CxR[subp][0])
                    down/=len(CxR[subp][1])
                    x=(up-down)
                    
                    up=0
                    down=0
                    for element in CyR[subp][0]:
                        up+=(stardata[star_nbr][element][0])
                    for element in CyR[subp][1]:
                        down+=(stardata[star_nbr][element][0])
                    up/=len(CyR[subp][0])
                    down/=len(CyR[subp][1])
                    y=(up-down)
                    
                    up=[]
                    down=[]
                    for elementnbr in CxR[subp][0]:
                        up+=[stardata[star_nbr][elementnbr][1]]
                    up=max(up)
                    for elementnbr in CxR[subp][1]:
                        down+=[stardata[star_nbr][elementnbr][1]]
                    down=max(down)
                    errx=max(up, down)
                    up=[]
                    down=[]
                    for elementnbr in CyR[subp][0]:
                        up+=[stardata[star_nbr][elementnbr][1]]
                    up=max(up)
                    for elementnbr in CyR[subp][1]:
                        down+=[stardata[star_nbr][elementnbr][1]]
                    down=max(down)
                    erry=max(up, down)

                    if('LP' in stars[star_nbr]):
                        color='g'
                    elif('HE' in stars[star_nbr]):
                        color='b'
                    elif('CS' in stars[star_nbr]):
                        color='r'
                    else:
                        color='b'
                    errorbar(x,y,xerr=errx,yerr=erry, ecolor=color)

                else:
                    print "Pannel {0} can't be displayed for star {1}, due to lack of data about :".format((subp+1), stars[star_nbr])
                    for el in can_plot[star_nbr][subp]:
                        print el_list[el]
                    print ""

        for subp in range(pannel_nbr):
            if pannel_nbr%2==1 and subp==(pannel_nbr-1):
                subplotnbr=(100*(subp/2+1))+10+(subp/2+1)
                subplot(subplotnbr)
            else:
                subplotnbr=((pannel_nbr-1)/2+1)*100+20+subp+1
                subplot(subplotnbr)
            xname='['
            for el in CxR[(subp)][0]:
                xname+=el_list[el]
                xname+=','
            xname=xname[:-1]+'/'
            for el in CxR[(subp)][1]:
                xname+=el_list[el]
                xname+=','
            xname=xname[:-1]+']'
            yname='['
            for el in CyR[(subp)][0]:
                yname+=el_list[el]
                yname+=','
            yname=yname[:-1]+'/'
            for el in CyR[(subp)][1]:
                yname+=el_list[el]
                yname+=','
            yname=yname[:-1]+']'

            xlabel(xname)
            ylabel(yname)

    def ppn_compare(self, simu_folder, cycle, star=None, Xref='Ba', zrange=None, linestyle=None, text_to_add='', legend_simu=None, grey_zone=False, plot_all_simu=True,solar_norm_file=None,convert_solar_norm_star=[None,None]):
        '''
        Method that adjusts a simulation with star data and plots boths, abundances vs z.
        If only one folder is specified, the curve of the simulation will appear.
        Else, the zone reached by the different folders will be colored in grey.

        Parameters :
        simu_folder : str, list
            path or list of paths of the folder the simulation(s) to plot.
        cycle : int, list
            Cycle of the simulation you want to plot.
            If a list is given, it must be of the same lenght as simu_folder.
        star : str, optional
            The star you want to plot. Default is None, uses starobs.current_star or enter a specific star name.
        Xref : str, list, optional
            Element that will be used as reference for the fit. Used to simulation the mixing factor for the star. Can be any of the element in the star data.
            Default is 'Ba'
        zrange : list, optional
            [zmin, zmax], that define the range of elements to plot.
        linestyle : str, list, optional
            Linestyle(s) one wants to give to the simulation lines. If only one is given, it will be applied to all lines.
            If a list is given, it must be the same lenght as simu_folder, but needs plot_all_simu=True.
        text_to_add : str, optional
            Text to add to the legendbox, with the starname and the neutron information.
        legend_simu : str, list, optional
            Legend one want to give to the simulation line. If a single simu_folder is given and legend=None (default), Nn and  tau will be given.
            If several folders are given, legend (not None) must be of the same lenght as simu_folder
        grey_zone : bool, optional
            Incase several simulation folders are given, allows to grey the zone between the different lines.
        plot_all_simu : bool, optional
            Defines if all simulations curves have to be plotted anyway.
            Default is True.
            Use False only with several folders, or to display only the star data.

        solar_norm_file : string, optional
            To normalize the simulation data with the same solar system abundances as done for the star data provide a path
            to the solar system abundance data. NuGrid initial abundance data type (USEEPP) used as input.
            If None, uses initial abundance as it was originally implemented.

        '''
        if type(simu_folder)==str:
            simu_folder=[simu_folder]
        
        if type(cycle)==list and type(simu_folder)==list and len(cycle)<>len(simu_folder):
            print "Cycles ill defined. Returning None"
            return None
        if type(linestyle)==list and type(simu_folder)==list and len(linestyle)<>len(simu_folder):
            print "Linestyles ill defined. Returning None"
            return None
        if type(legend_simu)==list and type(simu_folder)==list and len(legend_simu)<>len(simu_folder):
            print "Legends ill defined. Returning None"
            return None

        if type(legend_simu)==str or legend_simu==None:
            legend_simu=[legend_simu for i in range(len(simu_folder))]
        if type(cycle)<>list:
            cycle=[cycle for a in range(len(simu_folder))]
        if type(simu_folder)==list and type(linestyle)<>list and linestyle<>None:
            linestyle=[linestyle for i in range(len(simu_folder))]
        if star=='None' and self.__class__==starobs.plot:
            star=self._out_class_ds.current_star
        fig, ax=subplots(1)
        Z_star_tmp=self.get('Z', star)
        XFe_star_tmp=self.get('XFe', star,solar_norm=convert_solar_norm_star)
        Err_star_tmp=self.get('Err', star,solar_norm=convert_solar_norm_star)
        Upper_star_tmp=self.get('Upper', star,solar_norm=convert_solar_norm_star)
        Z_star=[]
        XFe_star=[]
        Err_star=[]
        Upper_star=[]
        if zrange==None:
            zrange=[min(Z_star_tmp), max(Z_star_tmp)]
        for z, x, e, u in zip(Z_star_tmp,XFe_star_tmp,Err_star_tmp,Upper_star_tmp):
            if zrange[0]<=z<=zrange[1]:
                Z_star.append(z)
                XFe_star.append(x)
                Err_star.append(e)
                Upper_star.append(u)

        form_str='%5.2F'
        el_abu_max=[]
        el_abu_min=[]
        toreturn=[]
        Xref_fol=Xref
        if type(Xref_fol)==str:
           Xref_fol=[Xref_fol for i in range(len(simu_folder))] 

        for subfolder_n in range(len(simu_folder)):
            Xref=Xref_fol[subfolder_n]
            print 'Doing : ',simu_folder[subfolder_n]
            p=ppn.abu_vector(simu_folder[subfolder_n])
            cycle_tmp=cycle[subfolder_n]
            p.get(cycle_tmp, decayed=True)
            zlist_tmp=np.arange(zrange[0], zrange[1]+1)
            zlist=[]
            for el in zlist_tmp:
                if len(np.where(p.z_iso_to_plot==el)[0].tolist())<>0:
                    zlist+=[el]
            cycles_data=_getData(p, thread=None, out_q=None, advancement_q=None, mode='ppn_compare', clist=[0,cycle_tmp], zlist=zlist, star_data=None, sim_0=None, Xref=None)

            #CR: get the solar mass fractions for zlist elements
            #replaces first cycle as initial abundance
	    if not solar_norm_file==None:
	        #cycles_data[0] needs to be replaced	   
		#[zlist, el_abu, el_abu_hash]
		temp_lists=cycles_data[0]

                print 'Use solar normalization from ',solar_norm_file
		print cycles_data[0][2].keys()
                solar_norm=get_solar_abu(el_name=cycles_data[0][2].keys(),solar_norm_file=solar_norm_file)
                el_abu=[];el_abu_hash={} 
                for k in range(len(zlist)):
                        el_abu.append(solar_norm[k])
                        el_abu_hash[cycles_data[0][2].keys()[k]]=solar_norm[k]

		cycles_data[0][1]=el_abu
		cycles_data[0][2]=el_abu_hash


	    cycles_data[1][1]/=np.array(cycles_data[0][1])
	    XFe_offset=XFe_star[Z_star.index(utils.get_z_from_el(Xref))]-np.log10(cycles_data[1][2][Xref]/cycles_data[0][2][Xref])
            
            if len(simu_folder)>1:
                if(el_abu_max==[]):
                    el_abu_max[:]=(np.log10(cycles_data[1][1][:])+XFe_offset)
                    el_abu_min[:]=(np.log10(cycles_data[1][1][:])+XFe_offset)
                else:
                    for i in range(len(cycles_data[1][0])):
                        if((np.log10(cycles_data[1][1][i])+XFe_offset)>el_abu_max[i]):
                            el_abu_max[i]=(np.log10(cycles_data[1][1][i])+XFe_offset)
                        if((np.log10(cycles_data[1][1][i])+XFe_offset)<el_abu_min[i]):
                            el_abu_min[i]=(np.log10(cycles_data[1][1][i])+XFe_offset)
            
            x=ppn.xtime(simu_folder[subfolder_n])
            age=x.get("time")[:cycle_tmp]
            rho=x.get("rho")[:cycle_tmp]
            v_therm=3.e8
            xn=x.get('NEUT')[:cycle_tmp]
            Nn=xn*ast.avogadro_constant*rho
            oneyear=utils.constants.one_year
            age_s=age*oneyear
            tau=integrate.cumtrapz(Nn*v_therm,age_s)
            tau=np.log10(tau*1.e-27)

            if subfolder_n==0:
                labelname="{0}\nlog tau={1}".format(star, form_str%tau[-1])
            if text_to_add<>'' and subfolder_n==0:
                labelname+='\n'+text_to_add
            if (plot_all_simu==True or len(simu_folder)==1) and linestyle<>None:
                plot(zlist,(np.log10(cycles_data[1][1])+XFe_offset), linestyle[subfolder_n], label=legend_simu[subfolder_n])
            elif (plot_all_simu==True or len(simu_folder)==1) and linestyle==None:
                plot(zlist,(np.log10(cycles_data[1][1])+XFe_offset), label=legend_simu[subfolder_n])
            legend(loc='best')

            for z in range(len(Z_star)):
                text(Z_star[z]-0.5,XFe_star[z]+0.05,utils.get_el_from_z(int(Z_star[z])), fontsize=15)

            xlabel('Z', fontsize=17)
            ylabel('[X/Fe]', fontsize=17)
            
            toreturn+=[(np.log10(cycles_data[1][1])+XFe_offset)]

        if len(simu_folder)>1 and grey_zone:
            fill_between(Z_star,el_abu_min,el_abu_max,color='0.6')
        errorbar(Z_star, XFe_star, yerr=Err_star, fmt='r*', linestyle='None', markersize=12)
        props=dict(boxstyle='round', fc='white', alpha=0.5)
        text(0.65, 0.18, labelname, transform=ax.transAxes, fontsize=15, va='top',  bbox=props)


        arrows=[]
        for el_nbr in range(len(Z_star)):
            if(Upper_star[el_nbr]):
                if len(el_abu_min)==len(Z_star):
                    arrows+=[[Z_star[el_nbr],XFe_star[el_nbr],0,-abs(XFe_star[el_nbr]-el_abu_min[zlist.index(Z_star[el_nbr])])*3/10]]
                else:
                    arrows+=[[Z_star[el_nbr],XFe_star[el_nbr],0,-0.5-abs(XFe_star[el_nbr]-(np.log10(cycles_data[1][1][zlist.index(Z_star[el_nbr])])+XFe_offset))*3/10]]
        if len(arrows)>0:
            arrowarray=np.array(arrows)
            X,Y,U,V = zip(*arrowarray)

            quiver(X,Y,U,V, headwidth=4, headlength=4, width=0.002)

        return zlist, toreturn

class Tools():

    def autofit(self, simu_folder, star=None, rsch_rg=None, step=1, Xref='Ba', weight=None, zrange=None, nbr_of_threads=0, sol_nbr=1,solar_norm_file_sim=None,convert_solar_norm_star=[None,None]):
        '''
        Method to automatically fit star data with one or several ppn simulations.
        It takes quite a while, since the get method from ppn.abu_vector is quite slow.
        It is highly advised to use a multi-core computer, and to restrict the research range !
    
        Parameters :
    
        simu_folder : str, list
            path or list of paths of the simulation folder to try to fit.
        star : str, optional
            The star you want to autofit. Default is None ans uses CEMP.current_star or enter a specific star name.
        rsch_rg : list, optional
            List with min and max cycle to explore the folder. Very important to redute the computing time. Usual ppn simulations are only interesting after cycle 300
            Default is None, e.g all the simulation.
        step : int, optional
            The simulation will only be explored every (step) cycles. Best way to reduce the time. Default is 1.
        Xref : str, optional
            Element that will be used as reference for the fit. Used to simulation the mixing factor for the star. Can be any of the element in the star data.
            Default is 'Ba'
        zrange : list,  optional
            [zmin, zmax], that define the range of rearch for the elements.
        weight : dict, optional
            Parameter to change the weight applied to each element. Some element can be not much relevant to fit, and so can be ignored putting a weight of 0.
            Use as weight['X']=x, X is the symbol of the isotope.
            Default is None, all elements from zrange will be equals.
        nbr_of_threads : int, optional
            Number of threads you want to imply. Must fit you cpu number of cores. First check the cpu use with top. 
            Default is None, python will use as many threads as cores present in this cpu.
        sol_nbr : int, optional.
            The autofit method tries to find different best fits. Since the simulation in continuous, it will search for best fit zones, described by the best cycle in it, to avoid to have, as best fits, cycles 800, 801, 799...
            sol_nbr defines the number of these zones you want to know.
	solar_norm_file_sim : string, optional
	    To normalize the simulation data with the same solar system abundances as done for the star data provide a path
	    to the solar system abundance data. NuGrid initial abundance data type (USEEPP) used as input.
            If None, uses initial abundance as it was originally implemented.
	convert_solar_norm_star : list, optional
	    Takes abundance distribution in spectroscopic notation
	    and converts to a new solar normalization.
	    list of two entries, both path to NuGrids USEEPP files
	    first entry represents the original normalization while
	    the latter represents the solar abundance of choice
 
        '''
        result={}
       
        if star=='None' and self.__class__==starobs.ppn_study:
            star=self._out_class_ds.current_star

        if nbr_of_threads==0:
            nbr_of_threads=mp.cpu_count()
        if weight==None:
            weight={}
        weight_tmp={}
        for el in weight:
            try:
                z=utils.get_z_from_el(el)
                weight_tmp[z]=weight[el]
            except: 
                print "Wrong element defined in weight. Returning None"
                return None
        weight=weight_tmp
        star_data_tmp=[]
        star_data_tmp+=[self.get('Z', star)]
        star_data_tmp+=[self.get('XFe', star,solar_norm=convert_solar_norm_star)]
        star_data_tmp+=[self.get('Err', star,solar_norm=convert_solar_norm_star)]
        star_data_tmp+=[self.get('Upper', star)]
        star_data=[[],[],[],{}, []]

	elements_for_fit=[]
        for z, x, e, u in zip(star_data_tmp[0], star_data_tmp[1], star_data_tmp[2], star_data_tmp[3]):
            if zrange==None or zrange[0]<=z<=zrange[1]:
	    #CR: added checks to exclude 'nan' values
		elem_name=utils.get_el_from_z(str(int(z)))
		if np.isnan([z,x,e,u]).any():
		    print 'Nan value found for element ',elem_name,'(Z=',z,')!!!'
		    print 'XFe= ',x,', Err=',e,', Upper=',u
		    print 'Fit not possible. Skip element!'
		    continue
                star_data[0]+=[z]
                star_data[1]+=[x]
                star_data[2]+=[e]
                star_data[4]+=[u]
		elements_for_fit.append(elem_name)
                if z in weight:
                    star_data[3][z]=weight[z]
	print '### Elements used for fit : ',elements_for_fit
        zlist=star_data[0]
        try:
            zref=utils.get_z_from_el(Xref)
        except:
            print "Wrong refence element. Returning None"
            return None
        if zref not in star_data[0]:
            print "Element reference not in the elements to be searched or not in star data. Retuning None"
            return None

        if type(simu_folder)==str:
            simu_folder=[simu_folder]
        for subfolder in simu_folder:
            print 'Studying : ',subfolder
            resultfolder={}
            p=ppn.abu_vector(subfolder)
            p.get(0,decayed=True)
            if rsch_rg<>None:
                if rsch_rg[0]<rsch_rg[1]<=len(p.files):
                    cycles=range(rsch_rg[0], rsch_rg[1], step)
                else:
                    print "Wrong research range provided. Returning None."
                    return None
            else:
                cycles=range(0, len(p.files), step)

            a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
            zlist=star_data[0]
            sim_0=[]
            for z in zlist:
                el=p.el_iso_to_plot[np.where(p.z_iso_to_plot==z)[0].tolist()[0]]
                X_el=p.abunds[np.where(p.el_iso_to_plot==el)[0].tolist()].sum()
                a_el.append(p.a_iso_to_plot[np.where(p.z_iso_to_plot==z)[0].tolist()[0]])
                el_abu.append(X_el)
                el_name.append(el)
                el_abu_hash[el]=X_el

            sim_0+=[el_abu]
            sim_0+=[el_abu_hash]
	    print 'test 1#########'
	    print sim_0
	    #CR: get the solar mass fractions for zlist elements
	    #replaces first cycle as initial abundance
            if not solar_norm_file_sim==None:
		print 'Use solar normalization from ',solar_norm_file_sim
                solar_norm=get_solar_abu(el_name,solar_norm_file_sim)
		el_abu=[];el_abu_hash={}	
		for k in range(len(el_name)):
			el_abu.append(solar_norm[k])
			el_abu_hash[el_name[k]]=solar_norm[k]			
		sim_0=[]	
                sim_0+=[el_abu]
                sim_0+=[el_abu_hash]
		print 'test 2########'
		print sim_0

            threads=[]
            out_q=mp.Queue()
            advancement_q=mp.Queue()
            resultdict={}
            cycledict={}

            for thread in range(nbr_of_threads):

                print "Launching thread {0} for cycles {1} to {2}.".format(thread , cycles[len(cycles)*thread/nbr_of_threads], cycles[len(cycles)*(1+thread)/nbr_of_threads-1])
                threads.append(mp.Process( target=_getData, args=(p, thread, out_q, advancement_q, 'autofit', cycles[len(cycles)*thread/nbr_of_threads:len(cycles)*(1+thread)/nbr_of_threads], zlist, star_data, sim_0, Xref)))
                threads[-1].start()
            old=-1

            while len(cycledict)<len(cycles):
                if(len(cycledict)*100/len(cycles) <> old):
                    print '{0}%\r'.format(str(len(cycledict)*100/len(cycles))),
                    sys.stdout.flush()
                    old=len(cycledict)*100/len(cycles)
                #cycledict.update(advancement_q.get())
		cycledict.update(my_queue_get(advancement_q)) #CR
            print "100%"
            for i in range(nbr_of_threads):
                resultdict.update(out_q.get())
            print "Dictionary updated"
            for i in threads:
                i.join()

            print "Threads ended."

            for cycle in resultdict:
                result[subfolder+'&'+str(cycle)]=resultdict[cycle]
        
        print 'Simulations loaded'
        print ''
        bestdist=[]
        bestcycle=[]
        bestfolder=[]
        sol=[]
        i=0
        while i < sol_nbr:
            bestdist+=[result[min(result, key=result.get)]]
            bestcycle+=[int(min(result, key=result.get).split('&')[1])]
            bestfolder+=[min(result, key=result.get).split('&')[0]]
            del result[min(result, key=result.get)]

            if len(bestfolder)<=1 or ((bestfolder[-1]==bestfolder[-2] and (bestcycle[-1]+step not in bestcycle[:-1] and bestcycle[-1]-step not in  bestcycle[:-1])) or bestfolder[-1] not in bestfolder[:-1]):
                print '--------'
                print "The best fit #{0} for {1} is obtained for :" .format((i+1), star)
                print "Cycle : "+str(bestcycle[-1])
                print "Folder : "+bestfolder[-1]
                print "Mean distance to fit : "+str(bestdist[-1])
                i+=1
                sol+=[[bestfolder[-1], bestcycle[-1], bestdist[-1]]]
            if len(result)==0:
                i=sol_nbr
    
        return sol 

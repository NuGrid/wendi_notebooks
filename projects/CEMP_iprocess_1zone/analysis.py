import starobs

import ascii_table
import utils

#read from fits file filename, created by fit.py
#filename="fits_all.txt"
cadc_data_path='/home/nugrid/CADC/nugrid/'

def results(db,star):

	filename="fits_all.txt"
	solfile=open(filename, "r")

	#solar systmen normalization consistent with star data
	#solar_norm_file='/nfs/rpod2/critter/PPN/forum.astro.keele.ac.uk_update/frames/mppnp/USEEPP/iniab1.3E-02As09.ppn'
	solar_norm_file=None

	#write results (figures and tables) in directory outputdir
	import os
	import shutil
	print 'test ',os.path.isdir('analysis')
	if os.path.isdir('analysis'):
		shutil.rmtree('analysis')
	os.mkdir('analysis')
	print 'mkdir analysis'
	outputdir='analysis'

	import matplotlib.pyplot as plt
	cycle=0
	for line in solfile:
		#star entry
	    if cycle%6==0:
		star=line.strip()
		#run entry
	    elif cycle%6==3 or cycle%6==4 or cycle%6==5:
		sol=line.split('&')
		plt.figure(1)
		#
		sol_elems, sol_abu=db.plot.ppn_compare(simu_folder=[sol[1].strip()], star=star, cycle=[int(sol[2].strip())], Xref=['Ba'], zrange=[31, 80], linestyle=['b-'], text_to_add='',solar_norm_file=solar_norm_file)

		#directory I am writing figure and table to
		#name='./Fits/Hvar_final/'+star+'_'+'sol'+sol[0].strip()+'_'+sol[2].strip()
		name = outputdir+'/'+star+'_'+'sol'+sol[0].strip()+'_'+sol[2].strip()
		plt.savefig(name+'.png', dpi=100)
		plt.close(1)
		el=[]
		for z in sol_elems:
		    el+=[utils.get_el_from_z(str(int(z)))]
		ascii_table.write(name+'.txt', headers=[star, sol[1].strip(), 'distance = '+sol[3].strip()], dcols=['El', 'Z', 'XFe'], data=[el, sol_elems, sol_abu[0]])
	    cycle+=1


	#download stuff
	#analysis
	print 'tar data'
	import os
	os.chdir('analysis')
	os.system('zip data *')
	os.chdir('../')
	print 'Data ready for download'
	

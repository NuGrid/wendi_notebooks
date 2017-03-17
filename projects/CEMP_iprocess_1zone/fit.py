import starobs
import ascii_table

#file to write data to
#create new fits_new.txt file
file="fits_all.txt"

solfile=open(file, "w")

#solar systmen normalization consistent with star data
#solar_norm_file='/nfs/rpod2/critter/PPN/forum.astro.keele.ac.uk_update/frames/mppnp/USEEPP/iniab1.3E-02As09.ppn'
#modify star data to match GN93
#convert_solar_norm_star=['/nfs/rpod2/critter/PPN/forum.astro.keele.ac.uk_update/frames/mppnp/USEEPP/iniab1.3E-02As09.ppn','/nfs/rpod2/critter/PPN/forum.astro.keele.ac.uk_update/frames/mppnp/USEEPP/iniab1.0E-03GN93.ppn']
convert_solar_norm_star=[None,None]
cadc_data_path='/home/nugrid/CADC/NuGrid/'
def start_fit(db,star):
	convert_solar_norm_star=[None,None]	
	starlist=[star]
	runsfolder=[cadc_data_path+'/projects/CEMP_iprocess_1zone/ppn_run_0.2_1.0_0.2/', cadc_data_path+'projects/CEMP_iprocess_1zone/ppn_run_0.2_1.0_0.1/', cadc_data_path+'/projects/CEMP_iprocess_1zone/ppn_run_0.2_1.0_0.05/']

	for star in starlist:
	    print 'Fit star ',star
	    #do the fit for star, test cycles in interval of 2, take 3 best values of all runs in runsfolder
	    # for elements between Z=31 (Ga) and Z=80 (Hg)
	    sols=db.data_support.autofit(simu_folder=runsfolder, star=star, zrange=[31,80], Xref='Ba', step=1, sol_nbr=3,solar_norm_file_sim=None,convert_solar_norm_star=convert_solar_norm_star)#,nbr_of_threads=10
	    #write the results into file
	    solfile.write(star+"\n")
	    if sols<>None:
		author=db.get('Author', star=star)
		year=db.get('Year', star=star)
		solfile.write(author[0]+"\n")
		solfile.write(year[0]+"\n")
		for sol in range(3):
		    towrite="{0} & {1} & {2} & {3}\n".format(sol+1, sols[sol][0], sols[sol][1], sols[sol][2])
		    solfile.write(towrite)

	solfile.close()


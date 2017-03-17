#add data to database
from shutil import copyfile
def add_star(star, author, year):
	copyfile('starobs_db_default.h5','starobs_db.h5')
	import starobs
	db=starobs.start('starobs_db.h5')
	db.add_dataset(star=star, author=author, year=year)
	return db

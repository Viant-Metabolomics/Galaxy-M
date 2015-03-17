#!/usr/bin/python

import os
import sqlite3

def Classes(db_name):

	fn_MIDB = os.path.join(os.getcwd(), 'tool-data', 'MI_Pack','MIDB_v0_9b.sqlite')

	conn = sqlite3.connect(fn_MIDB)
	cursor = conn.cursor()

	l = []

	if db_name != None and db_name != "None":
		
		cursor.execute("select class from database_class where db_name = '%s'" % (db_name))		
		records = cursor.fetchall()
		for record in records:
			if record[0] == "*":
				l.append([record[0], record[0], ""])
			else:
				l.append([record[0], record[0], ""])
	else:
		l.append(["*", "*", ""])

	#print l
	return l

#Classes("KEGG_COMPOUND")
#Classes("LIPIDMAPS")
#Classes("HMDB")
#Classes("BIOCYC")


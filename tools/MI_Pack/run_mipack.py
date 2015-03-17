#!/usr/bin/python

import sys
import os
import sqlite3

from copy import deepcopy

import MI_Pack
from MI_Pack import Ions
from MI_Pack import Atoms
from MI_Pack import Isotopes
from MI_Pack import Identification
from MI_Pack import Formula
from MI_Pack import PeakList


def main():

	##########################################################################
	##########################################################################

	fn_MIDB = os.path.join(os.path.abspath(os.path.dirname(__file__)).replace("tools","tool-data"), 'MIDB_v0_9b.sqlite')

	if sys.argv[1] == "SPS":
		
		peaklist = sys.argv[2]
		ppm = float(sys.argv[3])
		ion_mode = sys.argv[4]
		ions_usr = sys.argv[5]
		database = sys.argv[6]
		classes = sys.argv[7]
		#DB_KEGG_COMPOUND = sys.argv[6].replace("__ob____sq__","").replace("__sq____cb__","")
		#DB_KEGG_DRUG = sys.argv[7].replace("__ob____sq__","").replace("__sq____cb__","")
		#DB_LIPIDMAPS = sys.argv[8].replace("__ob____sq__","").replace("__sq____cb__","")
		#DB_BIOCYC = sys.argv[9].replace("__ob____sq__","").replace("__sq____cb__","")
		#DB_DRUG_BANK = sys.argv[10].replace("__ob____sq__","").replace("__sq____cb__","")
		#DB_HMDB = sys.argv[11].replace("__ob____sq__","").replace("__sq____cb__","")
		outfile_sql = sys.argv[8]

		#define databases for searching
		DBsToSearch = {}

		if classes != "None" and "*" not in classes:
			DBsToSearch[database] = classes.split(",")
		else:
			DBsToSearch[database] = "*"

		atoms = Atoms()
		ions = Ions(ion_mode)
		isotopes = Isotopes(ion_mode)

		if os.path.isfile(ions_usr) == False: #if input is not a file
			ions_usr = ions_usr.replace("__ob__","[").replace("__cb__","]")
			ions_usr = ions_usr.split(",")
			for ion in ions.lib.keys():
				if ion not in ions_usr:
					ions.delete(ion)
		else:
			ions.delete("*") 
			inp = open(ions_usr, "r")
			for line in inp:
				line = line.replace("\n", "").split("\t")
				ions.add(line[0],float(line[1]))
			inp.close()
		
		ID = Identification(peaklist, outfile_sql, ion_mode, ppm, atoms, ions, isotopes, fn_MIDB)
		
		ID.MIDB("MIDB", DBsToSearch)

		ID.SPS()


	elif sys.argv[1] == "TM":

		peaklist = sys.argv[2]
		ppm = float(sys.argv[3])
		ion_mode = sys.argv[4]
		ions_usr = sys.argv[5]
		DB_KEGG_COMPOUND = sys.argv[6]
		outfile_sql = sys.argv[7]

		DBsToSearch = {}
		if DB_KEGG_COMPOUND != "None" and "*" not in DB_KEGG_COMPOUND:
			DBsToSearch["KEGG_COMPOUND"] = DB_KEGG_COMPOUND.split(",")
		elif DB_KEGG_COMPOUND == "*":
			DBsToSearch["KEGG_COMPOUND"] = DB_KEGG_COMPOUND

		atoms = Atoms()
		ions = Ions(ion_mode)
		isotopes = Isotopes(ion_mode)

		if os.path.isfile(ions_usr) == False:
			ions_usr = ions_usr.replace("__ob__","[").replace("__cb__","]")
			ions_usr = ions_usr.split(",")
			for ion in ions.lib.keys():
				if ion not in ions_usr:
					ions.delete(ion)
		else:
			ions.delete("*") 
			inp = open(ions_usr, "r")
			for line in inp:
				line = line.replace("\n", "").split("\t")
				ions.add(line[0],float(line[1]))
			inp.close()

		ID = Identification(peaklist, outfile_sql, ion_mode, ppm, atoms, ions, isotopes, fn_MIDB)

		ID.MIDB("MIDB", DBsToSearch)
		ID.TM()


	elif sys.argv[1] == "EFS":

		peaklist = sys.argv[2]
		ppm = float(sys.argv[3])
		ion_mode = sys.argv[4]
		ions_usr = sys.argv[5]
		DBCal = sys.argv[6]
		atoms_usr = sys.argv[7]
		outfile_sql = sys.argv[8]

		atoms = Atoms()
		ions = Ions(ion_mode)
		isotopes = Isotopes(ion_mode)

		if os.path.isfile(ions_usr) == False:
			ions_usr = ions_usr.replace("__ob__","[").replace("__cb__","]")
			ions_usr = ions_usr.split(",")
			for ion in ions.lib.keys():
				if ion not in ions_usr:
					ions.delete(ion)
		else:
			ions.delete("*")
			inp = open(ions_usr, "r")
			for line in inp:
				line = line.replace("\n", "").split("\t")
				ions.add(line[0],float(line[1]))
			inp.close()

		if os.path.isfile(atoms_usr) == False:
			atoms_usr = atoms_usr.split(",")
			for atom in atoms.lib:
				if atom not in atoms_usr:
					atoms.delete(atom)
			
		else:
			atoms.delete("*")
			inp = open(atoms_usr, "r")
			lines = inp.readlines()
			for line in lines:
				line = line.replace("\n", "").split("\t")
				atoms.add(line[0], float(line[1]))
			for line in lines:    
				line = line.replace("\n", "").split("\t")
				atoms.add_limit(line[0], int(line[2]), int(line[3]))
			inp.close()

		ID = Identification(peaklist, outfile_sql, ion_mode, ppm, atoms, ions, isotopes, fn_MIDB)

		ID.EFS()

	elif sys.argv[1] == "PPS":

		peaklist = sys.argv[2]
		ppm = float(sys.argv[3])
		ion_mode = sys.argv[4]
		ions_usr = sys.argv[5]
		isotope_mode = sys.argv[6]
		isotopes_usr = sys.argv[7]
		outfile_sql = sys.argv[8]
		
		atoms = Atoms()
		ions = Ions(ion_mode)
		isotopes = Isotopes(ion_mode)

		if os.path.isfile(ions_usr) == False: #if input is not a file
			ions_usr = ions_usr.replace("__ob__","[").replace("__cb__","]")
			ions_usr = ions_usr.split(",")
			for ion in ions.lib.keys():
				if ion not in ions_usr:
					ions.delete(ion)
		else:
			ions.delete("*") 
			inp = open(ions_usr, "r")
			for line in inp:
				line = line.replace("\n", "").split("\t")
				ions.add(line[0],float(line[1]))
			inp.close()

		if os.path.isfile(isotopes_usr) == False: #if input is not a file
			isotopes_usr = isotopes_usr.replace("__ob__","").replace("__cb__","").split(",")
			for isotope in deepcopy(isotopes.lib):
				if " ".join(map(str,isotope)) not in isotopes_usr:
					isotopes.delete(isotope)
		else:
			isotopes.delete("*") 
			inp = open(isotopes_usr, "r")
			for line in inp:
				line = line.replace("\n","").split("\t")
				isotopes.add([line[0], line[1], float(line[2]), float(line[3]), float(line[4])])
			inp.close()
	
		ID = Identification(peaklist, outfile_sql, ion_mode, ppm, atoms, ions, isotopes, fn_MIDB)
		ID.PPS(ions)
		ID.PPS(isotopes)


	elif sys.argv[1] == "Output":
		
		peaklist = sys.argv[2]
		SPS_TM_sql = sys.argv[3]
		EFS_sql = sys.argv[4]
		PPS_sql = sys.argv[5]
		outfile_txt = sys.argv[6]
		outfile_sql = sys.argv[7]

		OutputType = None
		TypeDefined = False
		ion_mode = "POS" # DEFAULT

		for db in [SPS_TM_sql, EFS_sql, PPS_sql]:
			if db != "None" and TypeDefined == False:

				con = sqlite3.connect(db)
				cursor = con.cursor()
				cursor.execute("select name from sqlite_master where type = 'table';")
				tables = cursor.fetchall()
				con.close()

				for table in tables:
					
					if "_POS_" in table[0]:
						ion_mode = "POS"
					elif "_NEG_" in table[0]:
						ion_mode = "NEG"

					if "SPS" in table[0]:
						OutputType = "SPS"
						TypeDefined = True
					elif "TM" in table[0]:
						OutputType = "TM"
						TypeDefined = True
					elif "EFS" in table[0]:
						OutputType = "EFS"
						TypeDefined = True
					elif "PPS" in table[0]:
						OutputType = "PPS"
						TypeDefined = True

		con = sqlite3.connect(outfile_sql)
		cursor = con.cursor()

		for db in [SPS_TM_sql, EFS_sql, PPS_sql]:
			if db != "None":
				
				cursor.execute("""attach \"%s\" as toMerge""" % (db))
				cursor.execute("select name from toMerge.sqlite_master where type = 'table';")
				tables = cursor.fetchall()
				
				for table in tables:
					if "SPS" in table[0] or "TM" in table[0] or "EFS" in table[0] or "PPS" in table[0]:
						cursor.execute("create table %s as select * from toMerge.%s;" % (table[0], table[0]))
				con.commit()
				cursor.execute("detach database toMerge;")

		con.close()

		atoms = Atoms()
		ions = Ions(ion_mode)
		isotopes = Isotopes(ion_mode)

		ID = Identification(peaklist, outfile_sql, ion_mode, 1.0, atoms, ions, isotopes, fn_MIDB)		
		ID.output(OutputType, 1000000, outfile_txt)
		
	
if __name__ == "__main__": 
	main()



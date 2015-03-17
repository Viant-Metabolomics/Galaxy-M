#!/usr/bin/python

import os

def paths(path = "/home/weberrj/Desktop/data/"):

	l = []

	for d in os.listdir(path):
		if os.path.isdir(os.path.join(path, d)):
			if "@" in d:
				for dd in os.listdir(os.path.join(path, d)):
					if os.path.isdir(os.path.join(path, d, dd)):
						l.append((os.path.join("/", d, dd), os.path.join(path, d, dd), False))
			else:
				l.append(("/%s" % (d), "/%s" % (d), False))
	#print l
	return l


def files(dataTypes = [], path = "/home/weberrj/Desktop/data/"):

	l = []

	for d in os.listdir(path):
		if os.path.isdir(os.path.join(path, d)):
			if "@" in d:
				for dd in os.listdir(os.path.join(path, d)):
					if os.path.isdir(os.path.join(path, d, dd)):
						for fn in os.listdir(os.path.join(path, d, dd)):
						    for dataType in dataTypes:
						        if dataType in fn and fn[-1] != "~":
								    l.append((os.path.join("/", d, dd, fn), os.path.join(path, d, dd, fn), False))

			else:
				for fn in os.listdir(os.path.join(path, d)):
				    for dataType in dataTypes:
				        if dataType in fn and fn[-1] != "~":
				    		l.append((os.path.join("/", d, fn), os.path.join(path, d, fn), False))
	#print l
	return l

#studyPaths()
#files()


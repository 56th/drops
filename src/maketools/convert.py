#! /usr/bin/python

import sys
import re

target = None
target_list = [];

for line in sys.stdin:
	res = re.search("(?P<target>[a-zA-Z0-9_]+):",line)
	if res != None:
		target = (res.group('target'),[])
		continue

	if target == None:
		continue
	
	if re.search("LFLAGS",line):
		target_list.append(target);
		continue
	
	#since target != None && not last line of a target, this must be an interesting "object" line
	
	for token in line.split(' '):
		res = re.search("(?P<package>[a-zA-Z0-9_]+)/(?P<library>[a-zA-Z0-9_]+)\.o",token)
		if res != None:
			target[1].append("%s-%s" % (res.group('package'), res.group('library')))

for target in target_list:
	print "add_executable({name}.exe {name})".format(name=target[0])
	if target[1] != []:
		libs = "target_link_libraries({name}.exe".format(name=target[0])
		for req in target[1]:
			libs += " " + req
		libs += ")\n"
		print libs

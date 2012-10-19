import sys
import re

target_re = re.compile("(?P<name>[a-zA-Z0-9]+):")
end_re = re.compile("LFLAGS")
white_re = re.compile("\s+")
target = None

for line in sys.stdin:
	if target == None:
		res = target_re.match(line)
		if res != None:
			target = res.group('name')
	elif end_re.match(line):
		target = None
	else:
		print line.split(' ')
		for token in line.split(' '):
			if token == "" or token == '\\\n' or white_re.match(token):
				continue
			print '\'' + token + '\''

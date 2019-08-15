from __future__ import print_function

from re import sub,search,finditer,findall
from yaml import load,FullLoader
from numpy import unique

from colors import colors
from utils import isnumber


def crn_parameters(file_path) :
	'''returns a dictionary of crn parameters for use with PyCont'''

	# open and read file as string
	kwargs = { 'parameters':{} }
	with open(file_path, 'r') as file:

		# read and remove comments
		file_string = file.read()
		file_string = strip_comments(file_string)
		file_string = strip_brackets(file_string)

		# parse directives and parameters/rates
		kwargs = parse_directives(file_string,kwargs,ignore=['sweeps','inference','data','simulator','simulation','rates','spatial'])
		kwargs['parameters'] = { key:value for param in kwargs['parameters']
								 for key,value in param.items() if len(param)==1 }

	# hardcoded requirements for PyCont
	kwargs['parameters']['eps76'] = kwargs['parameters']['e76']; del kwargs['parameters']['e76']
	kwargs['parameters']['eps81'] = kwargs['parameters']['e81']; del kwargs['parameters']['e81']

	return kwargs['parameters']


def strip_comments(file_string,comment_symbol='//'):
	'''removes commented lines from string'''
	return sub(r'({}).*'.format(comment_symbol),' ',file_string)


def strip_brackets(file_string):
	'''removes excess brackets [] surrounding state variables'''

	# regex pattern for matching state variables
	pattern = r'\[\b[\w ]+\b\]'

	for state in unique(findall(pattern,file_string)) :
		state = state[1:-1]

		# strip brackets
		file_string = sub(r'\[{}\]'.format(state),state,file_string)

	return file_string


def parse_delimited(file_string,delimiter=r'\|',contains=''):
	'''find all strings containing symbols between delimiters'''

	# regex pattern for matching reaction syntax
	pattern = r'(?:[^{}\n]*{}[^{}\n]*)'.format(delimiter,contains,delimiter)
	reactions = []

	for match in finditer(pattern,file_string) :
			reaction = match.group()

			# parse with python compatible syntax
			reactions += [ reaction.replace('^','**') ]
	return reactions


def to_yaml(file_string) :
	'''convert file string to yaml compatible format'''

	pattern = r'\w+' # regex for mathcing words
	for word in unique(findall(pattern,file_string)) :

		# convert all words to strings
		if not isnumber(word[0]) :
			file_string = sub(r'\b{}\b(?!\()'.format(word),'"{}"'.format(word),file_string)

	# cast to yaml format
	file_string = file_string.replace(';',',')
	file_string = file_string.replace('=',':')

	# python compatible syntax
	file_string = file_string.replace('^','**')

	return file_string


def parse_directives(file_string,kwargs,ignore=['sweeps','inference','data','simulator']) :
	'''parse directives into kwargs'''

	# parse each directive to dictionary
	pattern = r'(?<=directive ).*?(?=(directive|init|\|))'
	for match in finditer(pattern,file_string.replace('\n',' ')) :

		# convert to yaml format
		name,directive = match.group().split(" ",1)
		directive = to_yaml(directive)

		if name.strip() in ignore :
			print(colors.lightgrey,'[parser]','ignoring directive {}'.format(name),colors.reset)

		########################################### TODO(gszep) ugly way to handle rates
		elif name == 'rates' :
			try :
				for rate_string in parse_delimited(directive,delimiter=r'\[,\]'):
					if rate_string.strip() :

						key,rate = rate_string.replace('"','').split(':')
						kwargs['rates'][key.strip()] = rate.replace(' ','')

			except : # parse manually
				raise Exception('failed to parse directive rates')
		###########################################

		else : # attempt to parse automagically
			try :
				kwargs[name.strip()] = load(directive,Loader=FullLoader)

			except : # parse manually
				raise Exception('failed to parse directive {}'.format(name))

	return kwargs


def parse_initials(file_string,kwargs) :
	'''parse initial conditions to kwargs; defaults to zero if not set'''

	for var in unique(findall(r'\b\w+\b',','.join(kwargs['reactions']))) :
		if var not in kwargs['rates'] and var not in kwargs['parameters'] :
			kwargs['init'][var] = 0.0

	########################################### TODO(gszep) ugly way to handle syntax
	for initial in parse_delimited(file_string,contains=r'\w+ \w+') :
		init_string = initial.strip().split(' ')

		if init_string[0] == 'init' :
			name,value = init_string[1],init_string[2]
			kwargs['init'][name.strip()] = float(value) if isnumber(value) else value

		elif init_string[0] != 'directive' :
			name,value = init_string[1],init_string[0]
			kwargs['init'][name.strip()] = float(value) if isnumber(value) else value
	###########################################

	return kwargs

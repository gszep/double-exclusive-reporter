from argparse import ArgumentTypeError


def isnumber(value):
	'''retrun boolean depending on whether the
	current value is parsable as a number.'''

	try:
		float(value)
		return True

	except (ValueError,TypeError):
		return False


def str2bool(v):
	'''converts strings to bools for commandline argparser'''
	if isinstance(v, bool):
		return v

	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True

	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False

	else:
		raise ArgumentTypeError('Boolean value expected.')
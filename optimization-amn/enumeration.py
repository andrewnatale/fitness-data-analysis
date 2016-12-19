
def enum(*args):
	"""
	Hack-y way for having enums in python.
	Declare with e = enum("val1", "val2", "val3", ...)
	Each value will me mapped to an integer on [0, nValues) based on order
	There is an additional field .size that says how many values are in this enum
	"""
	enums = dict(zip(args, range(len(args))));
	enums['size'] = len(args);	# attach a "size" field to the enum
	return type('Enumeration', (), enums);
import math as magic
import ast
import numpy
import warnings
from io import *
from enumeration import enum
from copy import *

# should the macrostates be hard-coded? probably not if this ends up being actually used for tuning other models...
MACROSTATES_T = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

class Model:
	"""
	A multistate design model
	"""
	# NOTE: Java documentation style. mainly because idk what's a standard one for python and Java's is most readable

	# I need these statements for intellisense.
	# vars used by all instances
	# TODO: is the enum even useful?
	MACROSTATES = enum();							# enum of the macrostates
	nMacrostates = 0;								# int
	ensembleSize = 0;								# int
	nPositions = 0;									# int, number of positions on sequence examined
	contiguousPositions = True;						# are the positions contiguous?
	positionMap = {};								# if the positions are not contiguous, a map from non-contiguous to [0, nPositions)
	positionOffset = 0;								# offset for positions when the positions examined are contiguous
	backrubTemp = 0.0;								# float
	boltzmannTemp = 0.0;							# float
	weights = numpy.array(0);						# float[] relative weights of the macrostates
	steepness = 0.0;								# float steepness of siggmoid (s)
	fitnesses = numpy.array(0);						# double[position][residue fitness] calculated fitnesses of residues at each location
	frequencies = numpy.array(0);					# double[position][residue frequency] calculated frequencies of residues at each location
	macrostateResidueEnergies = numpy.array(0);		# double[position][residue energy][macrostate]
	recovery = -1.0;								# float assigned by the outside similiartyMeasure to how well this mode recovers the sequence
	# TODO: check if macrostatesUsed is actually useful
	macrostatesUsed = numpy.array(False);			# bool[], macrostates examined during optimization
	
	# vars used when microstate data is involved
	useAltAveragingMethod = False;					# use accepted avg method or (the alt) xinjie's expression?
	isFrequenciesCalculated = False;				# prevent unnecessary calculations
	useMicrostateData = False;						# do we have data from individual microstates?
	areMicrostatesPicked = False;					# have microstates been selected to be used in the ensemble?
	microstateResidueEnergies = numpy.array(0);		# double[position][residue energy][macrostate][microstate]
	selectedMicrostateEnergies = numpy.array(0);	# double[position][residue energy][macrostate][microstate], subset of microstateREsidueEnergies
	microstateCounts = numpy.array(0);				# double[position][macrostate] number of microstates
	microstatesUsed = numpy.array(0);				# int[position][macrostate][microstate index], the microstates used to calculate the macrostates

	def __init__(self, macrostates:enum, ensembleSize:int, backrubTemp:float, boltzmannTemp:float, weights:numpy.array, steepness:float, positions:int, positionOffset:int, useMicrostateData:bool = False, posMap:dict = None, useAltAverageMethod:bool = False):
		"""
		Default constructor

		@param macrostates				enumeration of the different macrostates in this model
		@param ensembleSize				int of macrostate ensemble sizes
		@param backrubTemp				float of the backrub temperature
		@param boltzmannTemp			float, Boltzmann averaging temeprature
		@param weights					float[], wights for the macrostates
		@param steepness				float, steepness in fitness funciton
		@param positions				int, number of positions examined in this midel
		@param positionOffset			int, the index of the lowest index to be examined
		@param useMicrostateData		bool, are we to actually average microstate data?
		@param posMap					dict<int, int> a remapping of position values if the positions are not contiguous. ONLY pass an object if the positions are not contiguous
		@param useAltAveragingMethod	bool, use the other Boltzmann averaging calculation method?
		"""
		self.MACROSTATES = macrostates;
		self.nMacrostates = macrostates.size;
		self.ensembleSize = ensembleSize;
		self.backrubTemp = backrubTemp;
		self.boltzmannTemp = boltzmannTemp;
		self.weights = weights;
		self.steepness = steepness;
		self.nPositions = positions;
		self.positionOffset = positionOffset;
		self.isFrequenciesCalculated = False;
		self.useMicrostateData = useMicrostateData;
		self.useAltAveragingMethod = useAltAverageMethod;
		self.macrostatesUsed = numpy.array([True] * self.nMacrostates);
		self.positionMap = deepcopy(posMap);
		if posMap is not None:
			self.contiguousPositions = False;
		else:
			self.contiguousPositions = True;

		# allocate arrays
		self.macrostateResidueEnergies = numpy.zeros([self.nPositions, 20, self.nMacrostates], dtype = numpy.float64);
		self.fitnesses = numpy.ones([self.nPositions, 20], dtype = numpy.float64);
		self.frequencies = numpy.zeros([self.nPositions, 20], dtype = numpy.float64);

		# a little checker to prevent two identical entries, used in function addMacrostateData()
		for i in range(self.nMacrostates):
			for j in range(self.nPositions):
				self.macrostateResidueEnergies[j][0][i] = 65536.0; 

		if self.useMicrostateData:
			self.areMicrostatesPicked = False;
			self.microstatesUsed = numpy.zeros([0]);
			self.microstateCounts = numpy.zeros([self.nPositions, self.nMacrostates], dtype = int);
			self.microstateResidueEnergies = numpy.zeros([self.nPositions, 20, self.nMacrostates, 700], dtype = numpy.float64); # magic number 700 - max expected number of microstates		

	def constructFromExisting(existing, ensembleSize:int, backrubTemp:float, boltzmannTemp:float, weights:numpy.array, steepness:float):
		"""
		"Overloaded" "constructor" that uses a pre-existing Model as a template

		@param existing			pre-exisiting Model object
		@param ensembleSize		int, new ensemble size
		@param backrubTemp		float, new backrub temperate
		@param boltzmannTemp	float, new Boltzmann temerature
		@param weights			float[], new weights
		@return Model
		"""
		
		new = Model(existing.MACROSTATES, ensembleSize, backrubTemp, boltzmannTemp, weights, steepness, existing.nPositions, existing.positionOffset, existing.useMicrostateData);
		new.macrostatesUsed = existing.macrostatesUsed;
		new.microstatesUsed = existing.microstatesUsed;
		new.contiguousPositions = existing.contiguousPositions;

		# TODO: is deepy copy necessary for all these values?
		#new.positionMap = deepcopy(existing.positionMap);	# deep copy dict, keep instances completely separate
		#if (not existing.useMicrostateData):				# not using microstate data
		#	new.macrostateResidueEnergies = numpy.array(existing.macrostateResidueEnergies, dtype = numpy.float64);
		#elif ensembleSize == existing.ensembleSize:				# using microstate data, already collapsed
		#	new.selectedMicrostateEnergies = numpy.array(existing.selectedMicrostateEnergies);
		#	new.areMicrostatesPicked = True;
		#else:													# using microstate data, not collapsed
		#	#print("!", end='')
		#	new.microstateResidueEnergies = numpy.array(existing.microstateResidueEnergies);
		#	new.microstateCounts = numpy.array(existing.microstateCounts);

		# SHALLOW COPY of raw data for faster since they should not need modification during any runs
		new.positionMap = existing.positionMap;
		if (not existing.useMicrostateData):				# not using microstate data
			new.macrostateResidueEnergies = existing.macrostateResidueEnergies;
		else:
			if ensembleSize == existing.ensembleSize:				# using microstate data, already collapsed
				new.selectedMicrostateEnergies = existing.selectedMicrostateEnergies;
				new.areMicrostatesPicked = True;
		#else:													# using microstate data, not collapsed
			#print("!", end='')
			new.microstateResidueEnergies = existing.microstateResidueEnergies;
			new.microstateCounts = existing.microstateCounts;
		
		return new;

	def setPositionMap(self, posMap:dict) -> None:
		"""
		Deprecated. There should not be a reason to EVER call this function

		If this model is using non-contiguous position numbers, sets the Map<int, int>
		used to convert them to contiguous numbers on [0, nPositions]

		@param posMap		dict, Map<int, int> of how to change the numbers
		@return None
		"""
		warnings.warn("Obsolete function. its function was taken care of during construction. Now results are not guaranteed", DeprecationWarning);
		self.positionMap = deepcopy(posMap);
		return None;

	def addMacrostateData(self, macrostate:int, position:int, energies:"float[]") -> None:
		"""
		Inserts a macrostate_position set of fitness values to this model

		@param macrostate		int, the macrostate this corresponds to
		@param position			int, the position the energies corresponds to
		@param energies			float[] of length-20 of the energies
		@return void
		"""
		# convert raw position to an internal index
		if self.contiguousPositions:
			pos = position - self.positionOffset;
		else:
			pos = self.positionMap[position];
			
		if self.macrostateResidueEnergies[pos][0][macrostate] != 65536.0:
			raise Exception("Something something this entry already full");

		for i in range(20):
			self.macrostateResidueEnergies[pos][i][macrostate] = energies[i];

	def addMicrostateData(self, macrostate:int, position:int, energies:"float") -> None:
		"""
		Inserts a microstate_position fitness into this model

		@param macrostate		int, the macrostate this microstate belongs to
		@param position			int, the position the energy corresponds to
		@param energy			float, the energy of this mutation
		@retun void
		"""
		if self.contiguousPositions:
			position -= self.positionOffset;
		else:
			position = self.positionMap[position];

		# TODO: do I need a overwrite check as in adding macrostate data?
		for i in range(20):
			self.microstateResidueEnergies[position][i][macrostate][self.microstateCounts[position][macrostate]] = energies[i];

		self.microstateCounts[position][macrostate] += 1;

		return None;

	def useAltAverageMethod(self, yes:bool) -> None:
		"""
		Changes whether to use the other averaging method

		@param yes		bool whether to use it or not
		@return void
		"""
		self.useAltAveragingMethod = yes;
		self.isFrequenciesCalculated = False;

	# change weights sets
	# TODO ascartain that this is actually necessary
	def setWeights(self, newWeights:numpy.array) -> None:
		warning.warn("Why are you changing the weights directly in a model?", UserWarning);
		self.weights = newWeights;

	# PRIVATE
	# TODO: add flag to only compute once. Then we should be able to remove the deep copy
	def averageMicrostates(self) -> None:
		"""
		Boltzmann averages the microstates to calculate the energy for the macrostate

		@param void
		@return void
		"""
		#print(self.microstateResidueEnergies[0][0]);
		#print();
		if not self.areMicrostatesPicked:
			# pick backbones to use for the ensemble.
			self.microstatesUsed = numpy.zeros([self.nPositions, self.nMacrostates, self.ensembleSize], dtype = int);
			for i in range(self.nPositions):
				for j in range(self.nMacrostates):
					self.microstatesUsed[i][j] = numpy.random.randint(0, self.microstateCounts[i][j], [self.ensembleSize]);
			# cherry-pick out the selected microstates
			self.selectedMicrostateEnergies = numpy.zeros([self.nPositions, 20, self.nMacrostates, self.ensembleSize]);
			for i in range(self.nPositions):
				for j in range(20):
					for k in range(self.nMacrostates):
						for l in range(self.ensembleSize):
							self.selectedMicrostateEnergies[i][j][k][l] = self.microstateResidueEnergies[i][j][k][self.microstatesUsed[i][k][l]]
			
			self.areMicrostatesPicked = True;

		#print(self.selectedMicrostateEnergies[0][0]);
		#print();
		if not self.useAltAveragingMethod:
			if (self.boltzmannTemp == 0.0):
				self.macrostateResidueEnergies = numpy.amin(self.selectedMicrostateEnergies, axis = 3);
			elif (self.boltzmannTemp == -1.0):
				self.macrostateResidueEnergies = numpy.mean(self.selectedMicrostateEnergies, axis = 3); 
			else:
				self.macrostateResidueEnergies = numpy.sum(self.selectedMicrostateEnergies * numpy.exp(self.selectedMicrostateEnergies / -self.boltzmannTemp), axis = 3) / numpy.sum(numpy.exp(self.selectedMicrostateEnergies / -self.boltzmannTemp), axis = 3);
		else:
			if (self.boltzmannTemp == 0.0):
				self.macrostateResidueEnergies = numpy.amin(self.selectedMicrostateEnergies, axis = 3);
			elif (self.boltzmannTemp == -1.0):
				self.macrostateResidueEnergies = numpy.mean(self.selectedMicrostateEnergies, axis = 3);
			else:
				self.macrostateResidueEnergies = -numpy.log(sum(numpy.exp(self.selectedMicrostateEnergies / -self.boltzmannTemp), axis = 3));

		#print(self.macrostateResidueEnergies[0]);
		# After averaging, delete the 4D array to save space and flip the microstate flag
		self.microstateResidueEnergies = numpy.array(0);
		return None;

	# PRIVATE
	def calcFitness(self) -> None:
		"""
		Calculates the fitnesses of the each residue at each position.
		There is no need for this function to be externally called

		@param void
		@return void
		"""

		# collapse microstates into macrostates
		if self.useMicrostateData:
			self.averageMicrostates();

		minEnergies = numpy.amin(self.macrostateResidueEnergies, axis = 1);	# for each position and macrostate, which residue had min energy?
		offsets = minEnergies + numpy.divide(numpy.log(99), self.steepness);		# calculate offset is double[position][macrostate]
		self.fitnesses = numpy.ones([self.nPositions, 20], dtype = numpy.float64);
		for i in range(self.nPositions):
			for j in range(20):
				for k in range(self.nMacrostates):
					f = 1.0 / (1.0 + numpy.exp(self.steepness * (self.macrostateResidueEnergies[i][j][k] - offsets[i][k])));
					self.fitnesses[i][j] *= (1 - self.weights[k] + self.weights[k] * f);

		# PRIVATE
	def calcFrequencies(self) -> None:
		"""
		Calculates the fitnesses and the frequencies of residues at each location

		@param void
		@return void
		"""
		if not self.isFrequenciesCalculated: # only do it once!
			self.isFrequenciesCalculated = True;
			self.calcFitness();

			self.frequencies = numpy.divide(self.fitnesses, numpy.subtract(1.0, self.fitnesses));	# non-normalized frequencies
			sums = numpy.sum(self.frequencies, axis = 1);
			for i in range(self.nPositions):	# normalize
				self.frequencies[i] = numpy.divide(self.frequencies[i], sums[i]);

	# get functions
	# member fields should not be directly accessed; use these get funtions instead
	def getEnsembleSize(self) -> int:
		"""
		Self-explanatory name

		@return int
		"""
		return self.ensembleSize;

	def getBackrubTemp(self) -> float:
		"""
		Self-explanatory name

		@return double
		"""
		return self.backrubTemp;

	def getBoltzmannTemp(self) -> float:
		"""
		Self-explanatory name

		@return double
		"""
		return self.boltzmannTemp;

	def getWeights(self) -> numpy.array:
		"""
		Self-explanatory name. Return a deep copy of the array so it's
		safe to directly do math on the return value

		@return float[]
		"""
		return numpy.array(self.weights);

	def getSteepness(self) -> float:
		"""
		Self-explanatory name

		@return float
		"""
		return self.steepness;

	def getFrequencies(self) -> numpy.array:
		"""
		Self-explanatory name. Returns a deep copy of the array so it's
		safe to directly do math on the return value

		@return float[][]
		"""
		
		# this is a special instance for storing data
		if self.ensembleSize == 0 and self.useMicrostateData:
			raise PermissionError("This object is a raw data storage instance and this call should not have been made");
			return None;

		if not self.isFrequenciesCalculated:
			self.calcFrequencies();
		return numpy.array(self.frequencies);

	def equalTo(self, other) -> bool:
		"""
		Are the data stored in this Model correct? i.e. does everything actually correspond
		to the input file's data? Used for debugging

		@param other		Model object to compare this to
		@return bool	if everything compared correctly
		"""
		retval = True;
		if not isinstance(other, Model):
			print("is not same class")
			retval = False;
		if self.ensembleSize != other.ensembleSize:
			print("ensemble sizes: {:d}, {:d}".format(self.ensembleSize, other.ensembleSize))
			retval = False;
		if self.boltzmannTemp != other.boltzmannTemp:
			print("boltzmann temps: {:.2f}, {:.2f}".format(self.boltzmannTemp, other.boltzmannTemp));
			retval = False;
		if self.backrubTemp != other.backrubTemp:
			print("backrub temps: {:.2f}, {:.2f}".format(self.backrubTemp, other.backrubTemp));
			retval = False;
		if numpy.sum(numpy.abs(self.macrostateResidueEnergies - other.macrostateResidueEnergies)) > 1e-9:
			print(numpy.sum(numpy.abs(self.macrostateResidueEnergies - other.macrostateResidueEnergies)))
			retval = False;
		if not retval:
			print("{:d}\t{:.2f}\t{:.2f}".format(self.ensembleSize, self.boltzmannTemp, self.backrubTemp));
		return retval;

	# comparison operators based on similarity measure
	def __eq__(self, other):
		assert(isinstance(other, Model));
		return self.recovery == other.recovery;		# see, this is where static typing comes in useful. it lets autocomplete see the fields of other
	
	def __le__(self, other):
		assert(isinstance(other, Model));
		return self.recovery <= other.recovery;
	
	def __lt__(self, other):
		assert(isinstance(other, Model));
		return self.recovery < other.recovery;

	def __ge__(self, other):
		assert(isinstance(other, Model));
		return self.recovery >= other.recovery;

	def __gt__(self, other):
		assert(isinstance(other, Model));
		return self.recovery < other.recovery;

	def __ne__(self, other):
		assert(isinstance(other, Model));
		return self.recovery != other.recovery;

# see, my thing with semicolons is that they're the quivalent of the period/full stop of the english language.
# sure, you can technically write a paper without using punctuation and I'll still understand it.
# but, hell, does it look wrong and rambling.
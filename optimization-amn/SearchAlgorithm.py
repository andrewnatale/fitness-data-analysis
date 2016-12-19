from SimilarityMeasure import SimilarityMeasure
from model import Model
from enumeration import enum
from datetime import *
import numpy

class SearchAlgorithm:
	"""
	An abstract class/interface defining the expected behaviour of an optimization algorithm
	"""
	#optimizer = Optimizer();					# the optimizer with this this instance is associated
	#similarityMeasure = SimilarityMeasure();	# a SimilarityMeasure object used to calculate the fitness of models
	models = {};								# Map<hyperparams, models>
	maxIterations = 1024;						# hard cap on number of iterations before termination
	
	elapsedTime = datetime.now();	# how long the search took

	# TODO: determine if we need a macrostate enum as a safety measure

	# boundaries on search
	# TODO: resolve hard-coded things
	ensembleSizes = numpy.array([]);					# set of possible ensemble sizes
	backrubTemps = numpy.array([]);						# set of possible backrub temperatures
	boltzmannTemps = numpy.array([]);					# possible boltzmann temperatures or range of possible temps
	continuousBoltzmann = False;						# are we searching a range, or do we only have discrete macrostate data?
	steepnessRange = numpy.array([]);					# {minSteepness, maxSteepness}
	weightMins = numpy.array([]);						# minimum weights
	weightMaxs = numpy.array([]);						# maximum weights

	# which parameters to search, default is to search all of them.
	searchEnsemble = True;
	searchBackrub = True;
	searchBoltzmann = True;
	searchSteepness = True;
	searchWeights = numpy.array([]);

	# record of best parameters found so far
	bestEnsembleSize = 1;
	bestBackrubTemp = 1.0;
	bestBoltzmannTemp = 1.0;
	bestSteepness = 1.0;
	bestWeights = [];
	bestMatchVal = 0;
	bestFrequencies = numpy.array([]);

	# print things to console?
	suppressOutputs = False;

	def __init__(self, models = None, similarityMeasure:SimilarityMeasure = None, continuousBoltzmann:bool = False):
		"""
		Default constructor

		@param models				Map<string, model> of a set of Model objects holding the data read in
		@param similarityMeasure	the similarity measure to use
		@param continuousBoltzmann	whether to the boltzmann averagin search is contiuous or discrete
		"""
		self.similarityMeasure = similarityMeasure;
		self.models = models;
		self.continuousBoltzmann = continuousBoltzmann;
		self.searchEnsemble = True;
		self.searchBackrub = True;
		self.searchBoltzmann = True;
		self.searchSteepness = True;
		self.suppressOutputs = False;

		#self.optimizer = optimizer;

	def setSimilarityMeasure(self, similarityMeasure:SimilarityMeasure) -> None:
		"""
		Sets the similarity measure to be used

		@param similarityMeasure	new similarity measure to be used
		@return void
		"""
		self.similarityMeasure = similarityMeasure;

	def setMaxIterations(self, iter:int) -> None:
		"""
		Sets the max number of iterations of the algorithm. Typically the algorithm *will*
		run for this number of iterations....

		@param iter		int, max limit
		@return void
		"""
		self.maxIterations = iter;

	def setContinuousBoltzmann(self, yes:bool) -> None:
		"""
		Sets the algorithm whether to search through a continuous range of boltzmann averaging
		temperatures or not

		@param yes		bool, true for yes
		@return void
		"""
		if yes is not bool:
			raise TypeError("Input is not a bool");
		self.continuousBoltzmann = yes;
		return None;

	def getContinuousBoltzmann(self) -> bool:
		return self.continuousBoltzmann;

	def setParamBounds(self, ensembleSizes:"int[]", backrubTemps:"float[]", boltzmannTemps:"float[]", steepnessRange:"float[]", weightMins:"float[]", weightMaxs:"float[]") -> None:
		"""
		Sets the bounds on the parameter space to search through

		@param ensembleSizes		int[] of the ensemble sizes
		@param backrubTemps			float[] of discrete backrub temperatures
		@param boltzmannTemps		float[] of discrete Boltzmann temperatures
		@param steepnessRange		float[] of length 2: {lower bound, upper bound}
		@param weightMins			float[] of minimum weights
		@param weightMaxs			float[] of maximum weights
		@return void
		"""
		self.ensembleSizes = ensembleSizes;
		self.backrubTemps = backrubTemps;
		self.boltzmannTemps = boltzmannTemps;
		self.steepnessRange = steepnessRange;
		self.weightMins = weightMins;
		self.weightMaxs = weightMaxs;
		#self.searchWeights = numpy.array([True] * self.weightMins.shape[0]); # why was this line here anyways?

	def setSearchParameters(self, ensemble:bool, backrub:bool, boltzmann:bool, steepness:bool, weights:"bool[]") -> None:
		"""
		Sets which parameters to search through

		@param ensemble		bool
		@param backrub		bool
		@param boltzmann	bool
		@param steepness	bool
		@param weights		bool[]
		@return void
		"""
		self.searchEnsemble = ensemble;
		self.searchBackrub = backrub;
		self.searchBoltzmann = boltzmann;
		self.searchSteepness = steepness;
		self.searchWeights = weights;

	def setAllSearchToTrue(self) -> None:
		"""
		All parameters are to be searched

		@params void
		@return void
		"""
		self.searchEnsemble = True;
		self.searchBackrub = True;
		self.searchBoltzmann = True;
		self.searchSteepness = True;
		self.searchWeights = numpy.array([True] * self.weightMins.shape[0]);

	def setAllSearchToFalse(self) -> None:
		"""
		No parameters are to be searched, this is verifying some found parameters

		@params void
		@return void
		"""
		self.searchEnsemble = False;
		self.searchBackrub = False;
		self.searchBoltzmann = False;
		self.searchSteepness = False;
		self.searchWeights = numpy.array([False] * self.weightMins.shape[0]);

	# VIRTUAL ABSTRACT
	def iterate(self) -> None:
		"""
		Start searching
		
		@param void
		@return void
		"""
		raise NotImplementedError;

	def boundCheckBoltzmann(self, newBoltzmann:float) -> float:
		"""
		Check whether a new Boltzmann averaging temp is in the bounds
		If it is out of bounds, trim to the bound
		Should only be used when searching on a continuous range

		@param newBoltzmann		value to check against the bounds
		@return float	the input value if it's within bounds or a corrected value if not
		"""
		# validate that this function should be used
		if not self.continuousBoltzmann:
			raise AssertionError("This search is not on a continuous range of Boltzmann averaging temps");
			return None; # in the unlikely event that the call was in a try block, return a null to make sure that things still screw up down the line

		elif not self.searchBoltzmann:
			return self.boltzmannTemps[0];
		elif newBoltzmann < -1:						# less than inf
			return self.boltzmannTemps[1];
		elif newBoltzmann == 0:						# force move away from 0
			return self.boltzmannTemps[0];
		elif newBoltzmann < self.boltzmannTemps[0]: # lower than lower bound set to 0 - min boltzmann
			return 0;
		elif newBoltzmann > self.boltzmannTemps[1]:	# higher than higher bound set to -1 - inf
			return -1;
		else:
			return newBoltzmann;

	def boundCheckSteepness(self, newSteep:float) -> float:
		"""
		Check whether a new steepness value is in the bounds
		If it is out of bounds, trim to the bound

		@param newSteep		float, the value to check
		@return float		the corrected steepness value
		"""
		if not self.searchSteepness:	# check if we're searching steepness
			return self.steepnessRange[0];
		elif newSteep < self.steepnessRange[0]:
			return self.steepnessRange[0];
		elif newSteep > self.steepnessRange[1]:
			return self.steepnessRange[1];
		else:
			return newSteep;

	def boundCheckWeights(self, newWeights:numpy.array) -> numpy.array:
		"""
		Check whether a new set of weights in the specified bounds.
		Values out of bounds will be trimmed to the bound

		@param newWeights:		float[], new weights to check
		@return float[]			the checked set of weights
		"""
		newWeights = numpy.array(newWeights);	# useful for intellisense, otherwise the computer doesn't know type >.<
		for i in range(newWeights.size):
			if not self.searchWeights[i]:	# first check if this weight is to be searched
				newWeights[i] = self.weightMins[i];
			elif newWeights[i] < self.weightMins[i]:
				newWeights[i] = self.weightMins[i];
			elif newWeights[i] > self.weightMaxs[i]:
				newWeights[i] = self.weightMaxs[i];
			else:
				pass;
		return newWeights;

	def getModelByParams(self, param1, param2, param3) -> Model:
		"""
		Gets a model by the specified pre-determined parameters.
		Return object is a not deep copy that should not be modified

		@param ensembleS		int, ensemble size
		@param backrubT			float, backrub temperature
		@param boltzmannT		float, boltzmann averaging temp
		@return Model with specified params
		"""
		# convert to string
		param1 = str(param1);
		param2 = str(param2);
		param3 = str(param3);

		return self.models[param1 + " " + param2 + " " + param3];

	def getBestParameters(self) -> {}:
		"""
		Returns the currently best parameters.
		Keys:
			'ensembleSize'
			'backrubTemp'
			'boltzmannTemp'
			'steepness'
			'weights'
			'match'

		@param void
		@return		a dictionary with keys ensembleSize, backrubTemp, boltzmannTemp, steepness, weights
		"""
		bestParams = {};
		bestParams['ensembleSize'] = self.bestEnsembleSize;
		bestParams['backrubTemp'] = self.bestBackrubTemp;
		bestParams['boltzmannTemp'] = self.bestBoltzmannTemp;
		bestParams['steepness'] = self.bestSteepness;
		bestParams['weights'] = self.bestWeights;
		bestParams['match'] = self.bestMatchVal;
		return bestParams;

	def getBestFrequencies(self) -> numpy.array:
		"""
		Returns the set of frequencies corresponding to the best match

		@param void
		@return float[][] of frequencies
		"""
		return numpy.array(self.bestFrequencies);
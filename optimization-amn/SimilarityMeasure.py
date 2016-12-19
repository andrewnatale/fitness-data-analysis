import numpy

class SimilarityMeasure:
	"""
	An abstract class/interface defining the interface of a method for measuring
	the similarity between two distributions, in this case specifically for
	animo acid residue frequencies
	"""
	targetFrequencies = numpy.array(0);		# double[position][reside] of target frequencies

	# TODO: change this to where target freq is automatically assigned by parent Optimizer obj
	def __init__(self, targetFrequencies = None):
		"""
		Default constructor
		@param targetFrequencies	float[position][residue] of target frequencies to examine
		"""

		if targetFrequencies == None: # handle default constructor case
			self.targetFrequencies = None;
			return;

		self.targetFrequencies = numpy.array(targetFrequencies);

	def setTargetFreqs(self, targetFrequencies:numpy.array) -> None:
		"""
		Sets what target frequencies to use. SHould be double[position][residue]

		@param targetFrequencies	array of frequencies
		@return void
		"""
		self.targetFrequencies = numpy.array(targetFrequencies);

	def clone(self):
		"""
		Makes a copy of this similarity measure

		@param void
		@return SimilarityMeasure
		"""
		raise NotImplementedError;
		return SimilarityMeasure();

	def copy(SM):
		"""
		Static method to copy an existing SimilarityMeasure.

		@param SM	exist similarity measure
		@return SimilarityMeasure
		"""
		return SM.clone();

	# VIRTURAL ABSTRACT
	def getSimilarityMeasure(self, expFrequencies) -> float:
		"""
		Sees how similar an experimental frequency set is to the target/natural set

		@param expFrequencies		float[position][residue] of experimental frequencies
		@return						float in [0, 1], where 0 is no similarity and 1 is perfect match
		"""
		raise NotImplementedError;

	# VIRTUAL ABSTRACT
	def __str__(self, **kwargs):
		"""
		Return the name of the current similarity measure

		@param void
		@return	string of the name of the measure
		"""
		raise NotImplementedError;
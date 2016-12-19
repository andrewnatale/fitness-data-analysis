from SimilarityMeasure import SimilarityMeasure
import numpy

class MutualInformation(SimilarityMeasure):
	"""
	Mutual information. Uses I(X, Y) as a similarity measure.
	Minimizes H(X|Y). Unfortunately needs the joint distribution of X, Y
	and I don't think that that can be calculated...
	"""
	entropies = numpy.array(0);
	nPositions = 0;
	LOG20 = numpy.log(20);

	def __init__(self, targetFrequencies = None):
		super().__init__(targetFrequencies);
		self.nPositions = targetFrequencies.shape[0];
		self.entropies = numpy.zeros([self.nPositions]);
		for i in range(self.nPositions):
			# convert frequencies to probabilities at each position and calc entropy
			self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));
			self.entropies[i] = numpy.sum(numpy.nan_to_num(-1 * numpy.multiply(self.targetFrequencies[i], numpy.log(self.targetFrequencies[i]) / self.LOG20)));

	def clone(self):
		return MutualInformation(self.targetFrequencies);

	def getSimilarityMeasure(self, expFrequencies):
		for i in range(self.nPositions):
			expFrequencies[i] = numpy.divide(expFrequencies[i], numpy.sum(expFrequencies[i]));
		return

	def __str__(self, **kwargs):
		return "Mutual information";
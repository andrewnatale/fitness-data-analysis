from SimilarityMeasure import SimilarityMeasure
import numpy

class EntropyWeightedSimilarity(SimilarityMeasure):
	"""
	Similar to EntropyWeightsMixedSimilarity, but with only a single measure.
	Each position contributes differently to the overall similarity depending on
	the entropy at each position. Low entropy positions are given high weighting.
	"""
	entropies = numpy.array(0);
	weights = numpy.array(0);
	totalWeights = 0;
	nPositions = 0;
	LOG20 = numpy.log(20);
	similarityMeasures = [];

	def __init__(self, measure:SimilarityMeasure, targetFrequencies:numpy.array):
		super().__init__(targetFrequencies);
		self.nPositions = targetFrequencies.shape[0];
		self.entropies = numpy.zeros([self.nPositions]);
		self.weights = numpy.zeros_like(self.entropies);
		similarityMeasures = [];

		for i in range(self.nPositions):
			# convert frequencies to probabilities at each position and calc entropy
			self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));
			self.entropies[i] = numpy.sum(numpy.nan_to_num(-1 * numpy.multiply(self.targetFrequencies[i], numpy.log(self.targetFrequencies[i]) / self.LOG20)));
			sm = measure.clone();
			sm.setTargetFreqs(self.targetFrequencies[i]);
			self.similarityMeasures.append(sm);

		# weights is square of 1- entropy to give even less weight to high entropy positions
		self.weights = numpy.power(1 - self.entropies, 2);
		self.totalWeights = numpy.sum(self.weights);

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		self.nPositions = targetFrequencies.shape[0];
		entropies = numpy.zeros([self.nPositions]);
		self.weights = numpy.zeros_like(self.entropies);

		for i in range(self.nPositions):
			self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));
			self.entropies[i] = numpy.sum(numpy.nan_to_num(-1 * numpy.multiply(self.targetFrequencies[i], numpy.log(self.targetFrequencies[i]) / self.LOG20)));
		
		self.weights = numpy.power(1 - self.entropies, 2);
		self.totalWeights = numpy.sum(self.weights);

	def getSimilarityMeasure(self, expFrequencies):
		out = 0;		# composite similarity score
		for i in range(self.nPositions):
			expFrequencies[i] = expFrequencies[i] / numpy.sum(expFrequencies[i]);
			out += self.weights[i] * self.similarityMeasures[i].getSimilarityMeasure(expFrequencies[i]);
		return out / self.totalWeights;

	def __str__(self, **kwargs):
		return self.similarityMeasures[0].__str__() + " weighted by position entropies";
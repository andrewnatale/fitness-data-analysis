from SimilarityMeasure import SimilarityMeasure
import numpy

class EntropyWeightsMixedSimilarity(SimilarityMeasure):
	"""
	Mixes the scores of two measures, weighting them by entropy of
	the natural sequence at each position
	"""
	similarityMeasures1 = [];
	similarityMeasures2 = [];
	entropies = numpy.array(0);
	nPositions = 0;
	LOG20 = numpy.log(20);

	def __init__(self, SM1:SimilarityMeasure, SM2:SimilarityMeasure, targetFrequencies:numpy.array):
		"""
		A mix between two similarity measures, linearly weighted by the entropy at each position.

		@param SM1		similarity measure to be given more weight with increasing entropy
		@param SM2		similarity measure to be given less weight with increasing entropy
		@param targetFrequencies
		"""
		super().__init__(targetFrequencies);
		self.nPositions = targetFrequencies.shape[0];
		self.entropies = numpy.zeros([self.nPositions]);

		for i in range(self.nPositions):
			# convert frequencies to probabilities at each position and calc entropy
			self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));
			self.entropies[i] = numpy.sum(numpy.nan_to_num(-1 * numpy.multiply(self.targetFrequencies[i], numpy.log(self.targetFrequencies[i]) / self.LOG20)));
			# separate objects for comparison at each position
			SM1.setTargetFreqs(self.targetFrequencies[i]);
			self.similarityMeasures1.append(SM1.clone());
			SM2.setTargetFreqs(self.targetFrequencies[i]);
			self.similarityMeasures2.append(SM2.clone());
	
	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);

		similarityMeasures1 = [];
		similarityMeasures2 = [];
		self.nPositions = targetFrequencies.shape[0];
		entropies = numpy.zeros([self.nPositions]);

		for i in range(self.nPositions):
			# convert frequencies to probabilities at each position and calc entropy
			self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));
			self.entropies[i] = numpy.sum(numpy.nan_to_num(-1 * numpy.multiply(self.targetFrequencies[i], numpy.log(self.targetFrequencies[i]) / self.LOG20)));
			# separate objects for comparison at each position
			newSM1 = SM1.clone();
			newSM1.setTargetFreqs(self.targetFrequencies[i]);
			self.similarityMeasures1.append(newSM1);
			newSM2 = SM2.clone();
			newSM2.setTargetFreqs(self.targetFrequencies[i]);
			self.similarityMeasures2.append(newSM2);

	def getSimilarityMeasure(self, expFrequencies):
		out = 0;		# composite similarity score
		for i in range(self.nPositions):
			s = 0;		# similarity at this position
			expFrequencies[i] = numpy.divide(expFrequencies[i], numpy.sum(expFrequencies[i]));
			s += self.entropies[i] * self.similarityMeasures1[i].getSimilarityMeasure(expFrequencies[i]);
			s += (1 - self.entropies[i]) * self.similarityMeasures2[i].getSimilarityMeasure(expFrequencies[i]);
			out += s;
		return (out / self.nPositions); # normalize to 1

	def __str__(self, **kwargs):
		return "Entropy-weighted mixed between " + self.similarityMeasures1[0].__str__() + " and " + self.similarityMeasures2[0].__str__();
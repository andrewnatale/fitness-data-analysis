from SimilarityMeasure import *
import numpy

class JSDistByPos(SimilarityMeasure):
	"""
	Jensen Shannon distance, calculated for each position independently and then averaged
	i.e. for N positions, the measure is calculated as an average of N distanced in 20-space
	as opposed to a a single distance in N x 20 space
	"""
	nPositions = 0;

	def __init__(self, targetFrequencies = None):
		super().__init__(targetFrequencies);
		if self.targetFrequencies != None:
			self.nPositions = targetFrequencies.shape[0];
			# normalizing over each position and not over all vals
			for i in range(self.nPositions):
				self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		self.nPositions = targetFrequencies.shape[0];
		for i in range(self.nPositions):
			self.targetFrequencies[i] = numpy.divide(self.targetFrequencies[i], numpy.sum(self.targetFrequencies[i]));

	def getSimilarityMeasure(self, expFrequencies):
		sum = 0;
		h = lambda x : -1 * numpy.multiply(x, numpy.log2(x));	
		for i in range(self.nPositions):
			expFrequencies[i] = numpy.divide(expFrequencies[i], numpy.sum(expFrequencies[i]));
			JSDiv = numpy.nan_to_num(h(self.targetFrequencies[i]) + h(expFrequencies[i]) - h(self.targetFrequencies[i] + expFrequencies[i]));
			JSDiv = 0.5 * float(numpy.sum(JSDiv));
			if JSDiv < 0 or JSDiv > 1:
				print(JSDiv);
				exit(1);
			sum += numpy.sqrt(JSDiv);
		return sum / self.nPositions;

	def clone(self):
		return JSDistByPos(self.targetFrequencies);

	def __str__(self, **kwargs):
		return "Jensen-Shannon distance by position"
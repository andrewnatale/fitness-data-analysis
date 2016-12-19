from SimilarityMeasure import *
import numpy

class CosineSimilarityByPosition(SimilarityMeasure):
	"""
	Cosine similarity measure of frequencies. The frequencies at each location
	are normalized independently, rather than all frequencies at all positions
	being normalized altogether
	"""
	def __init__(self, targetFrequencies = None):
		super().__init__(targetFrequencies);

		# normalize to unit vect in (nPosition * 20)-space
		if self.targetFrequencies != None:
			for i in range(self.targetFrequencies.shape[0]):
				self.targetFrequencies[i] = numpy.nan_to_num(self.targetFrequencies[i] / numpy.linalg.norm(self.targetFrequencies[i]));

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		for i in range(self.targetFrequencies.shape[0]):
			self.targetFrequencies[i] = numpy.nan_to_num(self.targetFrequencies[i] / numpy.linalg.norm(self.targetFrequencies[i]));
	
	# override
	def getSimilarityMeasure(self, expFrequencies):
		for i in range(self.targetFrequencies.shape[0]):
			expFrequencies[i] = numpy.nan_to_num(expFrequencies[i] / numpy.linalg.norm(expFrequencies[i]));
		similarity = 0;
		for i in range(self.targetFrequencies.shape[0]):
			similarity += numpy.dot(self.targetFrequencies[i], expFrequencies[i]);
		return similarity / self.targetFrequencies.shape[0];	# since all vals are positive, they'll definitely be >= 0

	def clone(self) -> SimilarityMeasure:
		return CosineSimilarityByPosition(self.targetFrequencies);

	def __str__(self, **kwargs):
		return "Cosine similarity normalized by position"
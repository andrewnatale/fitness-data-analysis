from SimilarityMeasure import *
import numpy

class CosineSimilarity(SimilarityMeasure):
	"""
	Cosine similarity measure of frequencies. The frequency of each
	residue at each location is treated as an independent dimension.
	"""
	def __init__(self, targetFrequencies = None):
		super().__init__(targetFrequencies);

		# normalize to unit vect in (nPosition * 20)-space
		if self.targetFrequencies != None:
			self.targetFrequencies = numpy.nan_to_num(self.targetFrequencies / numpy.linalg.norm(self.targetFrequencies));

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		self.targetFrequencies = numpy.nan_to_num(self.targetFrequencies / numpy.linalg.norm(self.targetFrequencies));
	
	# override
	def getSimilarityMeasure(self, expFrequencies):
		expFrequencies = numpy.nan_to_num(expFrequencies / numpy.linalg.norm(expFrequencies));
		similarity = 0;
		for i in range(self.targetFrequencies.shape[0]):
			similarity += numpy.dot(self.targetFrequencies[i], expFrequencies[i]);
		return similarity;	# since all vals are positive, they'll definitely be >= 0

	def clone(self) -> SimilarityMeasure:
		return CosineSimilarity(self.targetFrequencies);

	def __str__(self, **kwargs):
		return "Cosine similarity"
from SimilarityMeasure import SimilarityMeasure
import numpy

class KLDivergence(SimilarityMeasure):
	"""
	Kullback-Leibler divergence, the general divergence upon which Jensen-Shannon is based
	"""
	def __init__(self, targetFrequencies):
		super().__init__(targetFrequencies);
		# normalize
		self.targetFrequencies = numpy.divide(self.targetFrequencies, numpy.sum(self.targetFrequencies));

	def getSimilarityMeasure(self, expFrequencies):
		# normalize
		expFrequencies = numpy.divide(expFrequencies, numpy.sum(expFrequencies));
		similarity = numpy.nan_to_num(numpy.log10(numpy.divide(self.targetFrequencies, expFrequencies))); # log base 2 since information but it really doesn't matter since it's just essentially scaling
		similarity = numpy.multiply(similarity, self.targetFrequencies);
		similarity = numpy.nan_to_num(similarity);
		similarity = numpy.sum(similarity);
		# K-L divergence is on range of [0, +inf), use exp(-KLD) to translate it to [0, 1]
		similarity = numpy.exp(-1 * similarity);
		return similarity;

	def __str__(self, **kwargs):
		return "Kullback-Leibler Divergence"
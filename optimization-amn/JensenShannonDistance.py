from SimilarityMeasure import *
import numpy

class JensenShannonDistance(SimilarityMeasure):
	"""
	Jensen-Shannon Distance, a symmetric version of Kullback-Leibler divergence.
	Specifically, this returns the J-S *distance*, which is a metric, and not
	the J-S *divergence*, which is not a metric. The J-S distance is simply the
	square root of the J-S divergence

	This class calculated the distance by treating all probabilities as a signle distribution.
	i.e. for N positions the distance is calculated in N x 20 space
	"""
	NOT_ZERO_BUT_CLOSE_ENOUGH = 0.0000000001;

	def __init__(self, targetFrequencies = None):
		super().__init__(targetFrequencies);
		# JSD is technically for probability distributions, so everything nees to sum to 1
		if self.targetFrequencies != None:
			self.targetFrequencies = self.targetFrequencies / numpy.sum(self.targetFrequencies);

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		self.targetFrequencies = self.targetFrequencies / numpy.sum(self.targetFrequencies);

	def getSimilarityMeasure(self, expFrequencies):
		expFrequencies = expFrequencies / numpy.sum(expFrequencies);

		# kernalized implementation - there's a paper on this somewhere
		h = lambda x : -1 * numpy.multiply(x, numpy.log2(x));	# in base 2 JS div is on range [0, 1]
		JSDiv = numpy.nan_to_num(h(self.targetFrequencies) + h(expFrequencies) - h(self.targetFrequencies + expFrequencies));
		JSDiv = 0.5 * float(numpy.sum(JSDiv)); # technically JSDiv is 1 - sum(things), but we're flipping the directions
		if float(JSDiv) < 0 or float(JSDiv) > 1:
			n = h(self.targetFrequencies) + h(expFrequencies) - h(self.targetFrequencies + expFrequencies);
			n = numpy.nan_to_num(n);
			print(n);
			print(numpy.sum(n));
			print(JSDiv);
			exit(1);

		#avgFreq = numpy.divide(numpy.add(expFrequencies, self.targetFrequencies), 2.0);
		#expFreqKLDiv = numpy.multiply(expFrequencies, numpy.divide(numpy.log(numpy.divide(expFrequencies, avgFreq)), numpy.log(20)));
		#expFreqKLDiv = numpy.sum(expFreqKLDiv);
		#tgtFreqKLDiv = numpy.multiply(self.targetFrequencies, numpy.divide(numpy.log(numpy.divide(self.targetFrequencies, avgFreq)), numpy.log(20)));
		#tgtFreqKLDiv = numpy.sum(tgtFreqKLDiv);
		#JSDiv = 0.5 * expFreqKLDiv + 0.5 * tgtFreqKLDiv;
		# normalize to the [0, 1] scale specified by superclass since J-S div is bounded on [0, log_n(2)]
		# then flip directionality to that specified by the superclass
		# JSDiv = (1 - JSDiv / (numpy.log(2) / numpy.log(20)));
		# the sqrt of JS divergence is JS distance
		return numpy.sqrt(JSDiv);
	
	def clone(self) -> SimilarityMeasure:
		return JensenShannonDistance(self.targetFrequencies);
	
	def __str__(self, **kwargs):
		return "Jensen-Shannon distance";
from SimilarityMeasure import SimilarityMeasure
import numpy

class Chi2Kernel(SimilarityMeasure):
	"""
	Implements the Chi-2 kernel between two frequency sets. The Chi-2 distance
	is often used against histograms and bag-o-words in computer vision, which
	is conceptually similar to what we're doing here.

	The kernel is implemented since it convinently transforms the [0, inf)
	Chi-2 distance to the [0, 1] scale used by similarity measure
	"""
	coeff = -1;

	def __init__(self, targetFrequencies = None, coeff:int = 1):
		super().__init__(targetFrequencies);
		self.coeff = -coeff;
		if self.targetFrequencies != None:
			self.targetFrequencies = numpy.nan_to_num(self.targetFrequencies / numpy.linalg.norm(self.targetFrequencies));

	def clone(self):
		return Chi2Kernel(self.getSimilarityMeasure());

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		if self.targetFrequencies != None:
			self.targetFrequencies = numpy.nan_to_num(self.targetFrequencies / numpy.linalg.norm(self.targetFrequencies));

	def getSimilarityMeasure(self, expFrequencies):
		expFrequencies = numpy.nan_to_num(expFrequencies / numpy.linalg.norm(expFrequencies));
		diffs = self.targetFrequencies - expFrequencies;
		diffs = numpy.power(diffs, 2);
		sums = self.targetFrequencies + expFrequencies;
		val = numpy.divide(diffs, sums);
		val = self.coeff * numpy.sum(val);
		return numpy.exp(val);

	def __str__(self, **kwargs):
		return "Chi-2 kernel with parameter " + str(self.coeff);

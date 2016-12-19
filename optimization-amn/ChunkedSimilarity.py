from SimilarityMeasure import SimilarityMeasure
import numpy

class ChunkedSimilarity(SimilarityMeasure):
	"""
	For a large number of positions, ChunkedSimilarity evaluates similarities
	of chunks of positions. A rudimentary method for dealing with high-dimensional
	problems. While JSD runs fine on 8 positios, everything goes funky is the 79 positions
	of the PDZ dataset, and the predictions go way off
	"""
	chunkSize = 0;
	nChunks = 0;
	nPositions = 0;
	similarityMeasures = [];

	def __init__(self, measure:SimilarityMeasure, chunkSize:int, targetFrequencies:numpy.array):
		super().__init__(targetFrequencies);
		self.nPositions = targetFrequencies.shape[0];
		self.chunkSize = chunkSize;
		self.nChunks = int(numpy.ceil(self.targetFrequencies.shape[0] / chunkSize));
		self.similarityMeasures = [];
		for i in range(self.nChunks):
			start = i * self.chunkSize;
			end = (i + 1) * self.chunkSize;
			end = end if end < self.nPositions else self.nPositions;
			sm = measure.clone();
			sm.setTargetFreqs(self.targetFrequencies[start:end]);
			self.similarityMeasures.append(sm);

	def setTargetFreqs(self, targetFrequencies):
		super().setTargetFreqs(targetFrequencies);
		measure = self.similarityMeasures[0].clone();
		self.nPositions = targetFrequencies.shape[0];
		self.chunkSize = chunkSize;
		self.nChunks = int(numpy.ceil(self.targetFrequencies.shape[0] / chunkSize));
		self.similarityMeasures = [];
		for i in range(self.nChunks):
			start = i * self.chunkSize;
			end = (i + 1) * self.chunkSize;
			end = end if end < self.nPositions else self.nPositions;
			sm = measure.clone();
			sm.setTargetFreqs(self.targetFrequencies[start:end]);
			self.similarityMeasures.append(sm);

	def getSimilarityMeasure(self, expFrequencies):
		out = 0;
		for i in range(self.nChunks):
			start = i * self.chunkSize;
			end = (i + 1) * self.chunkSize;
			end = end if end < self.nPositions else self.nPositions;
			# each chunk is given a weight based on the number of positions it contains
			# to avoid overweighting the last chunk, which most likely contains fewer positions
			out += (end - start) * self.similarityMeasures[i].getSimilarityMeasure(expFrequencies[start:end]);
		return out / self.nPositions;

	def clone(self):
		return ChunkedSimilarity(self.similarityMeasures[0], self.targetFrequencies);

	def __str__(self, **kwargs):
		return self.similarityMeasures[0].__str__() + " over " + str(self.nPositions) + " chunked with size " + str(self.chunkSize);
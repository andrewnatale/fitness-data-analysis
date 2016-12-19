from SearchAlgorithm import SearchAlgorithm
from SimilarityMeasure import SimilarityMeasure
import numpy

class GridSearch(SearchAlgorithm):
	"""
	Grid search over hyperparameter space
	"""

	steepnessIncrement = 0;				# finess of search grid for steepness
	weightIncrement = 0;				# finess of search grid for weights

	def __init__(self, similarityMeasure, models):
		super().__init__(similarityMeasure, models);

	def iterate(self):
		# TODO: implement this
		pass;
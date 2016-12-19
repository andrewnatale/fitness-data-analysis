from SearchAlgorithm import SearchAlgorithm
from SimilarityMeasure import SimilarityMeasure
from model import Model
from datetime import *
import numpy
import scipy.stats

class CuckooSearch(SearchAlgorithm):
	"""
	Cuckoo search for parameter optimization.
	See https://en.wikipedia.org/wiki/Cuckoo_search
	"""
	
	scaleParam = 1.0;			# scale parameter c used by the Levy distribution
	elimination = 0.20;			# fraction of individuals elimiated on each generation, i.e. discovery rate by parent birds
	populationSize = 512;		# number of eggs
	population = [];			# model[populationSize]

	# TODO: update the models parameter
	def __init__(self, models, similarityMeasure:SimilarityMeasure, continuousBoltzmann:bool, populationSize:int, scaleParam:float, elimination:float):
		"""
		Default constructor

		@param similarityMeasure	a SimiliartyMeasure object
		@param models				a Map<hyperParams, models>
		@param populationSize		int
		@param scaleParam			a float on (0, +inf) for the parameter c for the Levy distribution
										typical good values range from ~0.5 to ~2. Higher values increase
										tail weight of the distribution
		@param elimination			float on (0, 1), fraction of eggs to be eliminated every generation
		"""
		super().__init__(models, similarityMeasure, continuousBoltzmann);
		self.scaleParam = scaleParam;
		self.populationSize = populationSize;
		self.elimination = elimination;
		self.population = [];

	# PRIVATE
	def initPopulation(self):
		"""
		Initializes the population before the optimization begins

		@params void
		@return void
		"""
		self.bestMatchVal = 0;
		self.population = [];
		for i in range(self.populationSize):
			# rand continuous params
			thisSteepness = (numpy.random.rand() * (self.steepnessRange[1] - self.steepnessRange[0])) + self.steepnessRange[0] if self.searchSteepness else	self.steepnessRange[0];

			thisWeights = numpy.random.rand(self.weightMins.size);	# rand init first
			thisWeights *= (numpy.subtract(self.weightMaxs, self.weightMins));
			thisWeights += self.weightMins;
			thisWeights /= numpy.max(thisWeights); # TODO: normalize max to 1 or sum to 1?
			for j in range(thisWeights.shape[0]):					# then force non-search weights to preset val
				if not self.searchWeights[j]:
					thisWeights[j] = self.weightMins[j];

			# rand discrete params
			thisEnsembleSize = self.ensembleSizes[numpy.random.randint(0, self.ensembleSizes.size)] if self.searchEnsemble else	self.ensembleSizes[0];
			thisBackrubTemp = self.backrubTemps[numpy.random.randint(0, self.backrubTemps.size)] if self.searchBackrub else self.backrubTemps[0];
			thisBoltzmannTemp = self.boltzmannTemps[numpy.random.randint(0, self.boltzmannTemps.size)] if self.searchBoltzmann else self.boltzmannTemps[0];
			
			# different IDs depending on whether boltzmann is continuous or not
			if not self.continuousBoltzmann:
				m = Model.constructFromExisting(self.getModelByParams(thisBackrubTemp, thisEnsembleSize, thisBoltzmannTemp), thisEnsembleSize, thisBackrubTemp, thisBoltzmannTemp, thisWeights, thisSteepness);
			else:
				m = Model.constructFromExisting(self.getModelByParams(thisBackrubTemp, None, None), thisEnsembleSize, thisBackrubTemp, thisBoltzmannTemp, thisWeights, thisSteepness);
			m.macrostatesUsed = self.searchWeights;
			m.recovery = self.similarityMeasure.getSimilarityMeasure(m.getFrequencies());
			self.population.append(m);
		self.population.sort(reverse = True);
		self.recordBestParams();

	def iterate(self):
		self.initPopulation();
		# structure from the cuckoo search implementation in MATLAB by the authors
		# cuckoo = individual solution
		# on each iteration, every cuckoo does a Levy flight from its nest and lays a new egg
		# if the egg is better than its parent, the egg replaces the parent
		# else the egg dies
		# a certain percentage of birds are replaced by new birds
		#	the authors say some of the worse nests are replaced, but their code implies
		#	the all nests have an equal chance of being replaced....

		start = datetime.now();	# track runtime

		# a little progress bar to make the wait bearable
		if not self.suppressOutputs:
			updateStep = int(numpy.ceil(self.maxIterations / 70));		# 70-char width outputs
			if self.maxIterations != 1:
				print("going for {:d} generations".format(self.maxIterations));
				for i in range(int(numpy.floor(self.maxIterations / updateStep))):
					print('_', end='');
			else:
				print('_', end='');
			print();

		for i in range(self.maxIterations):
			for j in range(self.populationSize):
				# Levy flight away from this bird's parameters
				# TODO: should the Levy steps be scales differently for each parameter? probably yes.
				newSteepness = self.boundCheckSteepness(self.nextLevyStep() + self.population[j].getSteepness());
				newWeights = self.boundCheckWeights(self.nextLevySteps(self.bestWeights.size) + self.population[j].getWeights());

				# TODO: implement something that is a single hop in the discrete vars instead of random
				newEnsembleSize = self.ensembleSizes[numpy.random.randint(0, self.ensembleSizes.size)] if self.searchEnsemble else self.ensembleSizes[0];
				newBackrubTemp = self.backrubTemps[numpy.random.randint(0, self.backrubTemps.size)] if self.searchBackrub else self.backrubTemps[0];
				# 2 differents conditions for Boltzmann temperatures
				if not self.continuousBoltzmann:
					newBoltzmannTemp = self.boltzmannTemps[numpy.random.randint(0, self.boltzmannTemps.size)] if self.searchBoltzmann else self.boltzmannTemps[0];
					newModel = Model.constructFromExisting(self.getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteepness);
				else:
					newBoltzmannTemp = self.boundCheckBoltzmann(self.nextLevyStep() + self.population[j].getBoltzmannTemp());
					newModel = Model.constructFromExisting(self.getModelByParams(newBackrubTemp, None, None), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteepness);
					
				newModel.recovery = self.similarityMeasure.getSimilarityMeasure(newModel.getFrequencies());
				# replace parent if better
				if newModel.recovery > self.population[j].recovery:
					self.population[j] = newModel;

				# is this egg to be replaced?
				if (numpy.random.rand() < self.elimination):
					# if so, the replacement egg should be somewhat similar to the original egg
					# see auther's MATLAB implementation. I don't *quite* understand why this method. Yet.
					randParent1 = numpy.random.randint(0, self.populationSize);
					randParent2 = numpy.random.randint(0, self.populationSize);
					multiplier = numpy.random.rand();
					steepnessStep = multiplier * (self.population[randParent1].getSteepness() - self.population[randParent2].getSteepness());
					weightsStep = multiplier * (self.population[randParent1].getWeights() - self.population[randParent2].getWeights());
					newSteepness = self.boundCheckSteepness(self.population[j].getSteepness() + steepnessStep);
					newWeights = self.boundCheckWeights(self.population[j].getWeights() + weightsStep);

					# TODO implement something to move between discrete vars instead of random
					newEnsembleSize = self.ensembleSizes[numpy.random.randint(0, self.ensembleSizes.size)] if self.searchEnsemble else self.ensembleSizes[0];
					newBackrubTemp = self.backrubTemps[numpy.random.randint(0, self.backrubTemps.size)] if self.searchBackrub else self.backrubTemps[0];
					if not self.continuousBoltzmann:
						newBoltzmannTemp = self.boltzmannTemps[numpy.random.randint(0, self.boltzmannTemps.size)] if self.searchBoltzmann else self.boltzmannTemps[0];	
						newModel = Model.constructFromExisting(self.getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteepness);
					else:
						boltzmannStep = multiplier * (self.population[randParent1].getBoltzmannTemp() - self.population[randParent2].getBoltzmannTemp());
						newBoltzmannTemp = self.boundCheckBoltzmann(self.population[j].getBoltzmannTemp() + boltzmannStep);
						newModel = Model.constructFromExisting(self.getModelByParams(newBackrubTemp, None, None), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteepness);
						
					newModel.recovery = self.similarityMeasure.getSimilarityMeasure(newModel.getFrequencies());

					self.population[j] = newModel;

			self.population.sort(reverse = True);
			if (self.population[0].recovery > self.bestMatchVal):
				self.recordBestParams();
			if not self.suppressOutputs:
				if int(numpy.mod(i, updateStep)) == 0:
					print(">", end='');
		if not self.suppressOutputs:
			print();
		self.elapsedTime = datetime.now() - start;

	def nextLevyStep(self) -> float:
		"""
		Generates a random number from this model's Levy distribution.
		The magnitude is drawn from a Levy distribution f(x; 0, scaleParam)
		While the direction is uniform random

		@param void
		@return float
		"""
		# generate a random Levy by transforming from a rand uniform
		# see https://en.wikipedia.org/wiki/L%C3%A9vy_distribution#Random_sample_generation
		r1 = numpy.random.rand();
		randLevy = self.scaleParam * numpy.power(scipy.stats.norm.ppf(1.0 - r1 / 2.0), -2); # ppf is the inverse normal cdf
		# random direction
		r2 = numpy.sign((numpy.random.rand() - 0.5));
		return randLevy * r2 * 0.01; # Cuckoo search authors says to use 1/100 of the scale length

	def nextLevySteps(self, steps:int) -> numpy.array:
		"""
		Generates an array of Levy steps

		@param steps		size of array
		@return float[]		an array of independent Levy steps
		"""
		out = numpy.zeros(steps);
		for i in range(steps):
			out[i] = self.nextLevyStep();
		return out;

	def recordBestParams(self) -> None:
		self.bestBackrubTemp = self.population[0].getBackrubTemp();
		self.bestBoltzmannTemp = self.population[0].getBoltzmannTemp();
		self.bestEnsembleSize = self.population[0].getEnsembleSize();
		self.bestSteepness = self.population[0].getSteepness();
		self.bestWeights = self.population[0].getWeights();
		self.bestFrequencies = self.population[0].getFrequencies();
		self.bestMatchVal = self.population[0].recovery;

	def __str__(self, **kwargs):
		return "Cuckoo search, scale: {:.4f}, elimination: {:.4f}, population: {:d}, generations: {:d}".format(self.scaleParam, self.elimination, self.populationSize, self.maxIterations);
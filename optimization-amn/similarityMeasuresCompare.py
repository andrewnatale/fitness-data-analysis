from SimilarityMeasure import SimilarityMeasure
from enumeration import enum
from JensenShannonDistance import JensenShannonDistance
from CosineSimilarity import CosineSimilarity
from KLDivergence import KLDivergence
from EntropyWeightedSimilarity import EntropyWeightedSimilarity
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity
from Chi2Kernel import Chi2Kernel
from Optimizer import Optimizer
from CuckooSearch import CuckooSearch
from datetime import *
import sys
import numpy

similarityMeasure = int(sys.argv[1]);

#def convMicro():
MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

# only looking at MACROSTATE.TS
# only optimizing backrub temperature and steepness
ensembleSizes = numpy.array([20, 50, 100, 150, 200]);
backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
boltzmannTemps = numpy.array([-1, 5.0]);
steepnessRange = numpy.array([0.5, 5]);
minWeights = numpy.array([0, 0, 0, 0, 0]);
maxWeights = numpy.array([1, 1, 1, 1, 1]);

data = "/netapp/home/tianjiao.zhang/data/microstates.dat";
targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

dataAlt = "C:\\users\\candy\\skydrive\\documents\\rotation 2\\DHFR microstates\\microstates.dat";
targetFreqsAlt = "C:\\users\\candy\\skydrive\\documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

optimizer = Optimizer(MACROSTATES, True);

# slightly different paths on my two computers

optimizer.readTargetFrequencies(targetFreqs);	
optimizer.readFormattedMicrostateData(data);

if similarityMeasure == 0:
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 32, 1, 0.25);
elif similarityMeasure == 1:
	search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), True, 32, 1, 0.25);
elif similarityMeasure == 2:
	search = CuckooSearch(optimizer.models, KLDivergence(optimizer.targetFrequencies), True, 32, 1, 0.25);
elif similarityMeasure == 3:
	search = CuckooSearch(optimizer.models, Chi2Kernel(optimizer.targetFrequencies), True, 32, 1, 0.25);
elif similarityMeasure == 4:
	search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(JensenShannonDistance(), optimizer.targetFrequencies), True, 32, 1, 0.25);
elif similarityMeasure == 5:
	search = CuckooSearch(optimizer.models, EntropyWeightsMixedSimilarity(JensenShannonDistance(), CosineSimilarity(), optimizer.targetFrequencies), True, 32, 1, 0.25);
else:
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 32, 1, 0.25);

print(search.similarityMeasure.__str__());

search.setMaxIterations(2048);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setAllSearchToTrue();
optimizer.useAlgorithm(search);
optimizer.optimize();
now = datetime.now();
optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "Similarity measure compare " + str(similarityMeasure) + now.strftime('%Y%m%d%H%M%S') + ".fasta");
optimizer.writeBestParamsToText("Similarity measure compare " + str(similarityMeasure) + now.strftime('%Y%m%d%H%M%S'));
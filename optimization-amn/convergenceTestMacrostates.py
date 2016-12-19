from SimilarityMeasure import SimilarityMeasure
from JensenShannonDistance import JensenShannonDistance
from CuckooSearch import CuckooSearch
from enumeration import enum
from Optimizer import Optimizer
import numpy

MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

# only looking at MACROSTATE.TS
# only optimizing backrub temperature and steepness
ensembleSizes = numpy.array([20, 50]);
backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
steepnessRange = numpy.array([0.5, 5]);
minWeights = numpy.array([0, 0, 0, 0, 0, 0]);
maxWeights = numpy.array([1, 1, 0, 1, 1, 1]);

data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat1.tsv";
targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

optimizer = Optimizer(MACROSTATES);

# slightly different paths on my two computers
try:
	optimizer.readTargetFrequencies(targetFreqs);	
	optimizer.readData(data);
except:
	optimizer.readTargetFrequencies(targetFreqsAlt);	
	optimizer.readData(dataAlt);

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(16);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("16 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(32);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("32 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(64);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("64 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(128);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("128 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(256);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("256 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(512);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("512 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(1024);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("1024 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(2048);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("2048 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(4096);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("4096 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 128, 1, 0.25);
search.setMaxIterations(16);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
optimizer.useAlgorithm(search);
optimizer.optimize();
print("8192 iterations | runtime {:s} | match {:.6f}".format(str(search.elapsedTime), optimizer.getBestParameters()['match']));
from Optimizer import Optimizer;
from model import Model;
from SimilarityMeasure import SimilarityMeasure;
from SearchAlgorithm import SearchAlgorithm;
from CuckooSearch import CuckooSearch;
from KLDivergence import KLDivergence;
from JensenShannonDistance import JensenShannonDistance;
from CosineSimilarity import CosineSimilarity;
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity;
from EntropyWeightedSimilarity import EntropyWeightedSimilarity
from Chi2Kernel import Chi2Kernel;
from enumeration import enum;
from datetime import *
import numpy;
import threading
from io import *
# tests with small subsets of data

def smallTest(iterations = 1024):
	print("Hello!\n");
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

	print("Initializing objects\n");

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataMicro = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";
	dataMicroAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";

	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");

	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 32, 1, 0.25);
	search.setMaxIterations(iterations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
	optimizer.useAlgorithm(search);

	print("\nJS Dist");
	#search.setSimilarityMeasure(JensenShannonDistance(optimizer.targetFrequencies));
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "fixed ensembles " + now.strftime('%Y%m%d%H%M') + ".fasta");
	optimizer.writeBestParamsToText("fixed ensembles " + now.strftime('%Y%m%d%H%M'));
	print(optimizer.getBestParameters()['match']);
	
	return None;

def smallTestBoltz():
	print("Hello!\n");
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([128]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5]);
	steepnessRange = numpy.array([1, 7]);
	minWeights = numpy.array([0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 1, 1, 1]);

	print("Initializing objects\n");

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	dataMicro = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\microstates.dat";
	dataMicroAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\microstates.dat";

	optimizer = Optimizer(MACROSTATES, True);

	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readFormattedMicrostateData(dataMicro);
	except FileNotFoundError:
		optimizer.readtargetfrequencies(targetfreqsalt);	
		optimizer.readFormattedMicrostatedata(datamicroalt);
		
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 64, 1, 0.25);
	search.setMaxIterations(2048);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True, True, True, True, True]));
	optimizer.useAlgorithm(search);
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "var ensembles " + now.strftime('%Y%m%d%H%M') + ".fasta");
	optimizer.writeBestParamsToText("var ensembles " + now.strftime('%Y%m%d%H%M'));
	
	#for i in range(8):
	#	thread = optimizerThread();
	#	thread.copyOptimizer(optimizer);
	#	thread.run();

	return None;

def testRandUniformInput():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	# only looking at MACROSTATE.TS
	# only optimizing backrub temperature and steepness
	ensembleSizes = numpy.array([50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 0, 1, 1, 1]);

	print("Initializing objects\n");

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	
	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	# make energies uniform
	for model in optimizer.models:
		optimizer.models[model].macrostateResidueEnergies = numpy.ones_like(optimizer.models[model].macrostateResidueEnergies);

	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 1, 1, 0.25);
	search.setMaxIterations(1);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, False, False, False, numpy.array([False, False, False, False, False, False]));
	optimizer.useAlgorithm(search);

	outfile = open("uniform energy similarities.txt", 'w');
	optimizer.optimize();
	outfile.write("JSD: {:.4f}\n".format(optimizer.getBestParameters()['match']));

	search.setSimilarityMeasure(CosineSimilarity(optimizer.targetFrequencies));
	optimizer.optimize();
	outfile.write("Cosine similarity: {:.4f}\n".format(optimizer.getBestParameters()['match']));

	search.setSimilarityMeasure(KLDivergence(optimizer.targetFrequencies));
	optimizer.optimize();
	outfile.write("K-L divergence: {:.4f}\n".format(optimizer.getBestParameters()['match']));

	search.setSimilarityMeasure(EntropyWeightsMixedSimilarity(CosineSimilarity(), JensenShannonDistance(), optimizer.targetFrequencies));
	optimizer.optimize();
	outfile.write("Weighted mixed similarity: {:.4f}\n".format(optimizer.getBestParameters()['match']));
	outfile.close();
	return None;

def smalltestPrevOptimalVals():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	# only looking at MACROSTATE.TS
	# only optimizing backrub temperature and steepness
	ensembleSizes = numpy.array([50]);
	backrubTemps = numpy.array([1.8]);
	boltzmannTemps = numpy.array([0.0]);
	steepnessRange = numpy.array([3.0]);
	minWeights = numpy.array([0.80, 0.55, 0, 0.90, 0.35, 1.00]);
	maxWeights = numpy.array([0.80, 0.55, 0, 0.90, 0.35, 1.00]);

	print("Initializing objects\n");

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat6.tsv";
	dataAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat5.tsv";
	
	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");

	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 1, 1, 0.25);
	search.setMaxIterations(1);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, False, False, False, numpy.array([False, False, False, False, False, False]));
	optimizer.useAlgorithm(search);

	#print("Cos similiarity");
	#optimizer.optimize();	
	#optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "testOutCos.fasta");
	#print(optimizer.getBestParameters()['match']);

	print("\nJS Dist");
	#search.setSimilarityMeasure(JensenShannonDistance(optimizer.targetFrequencies));
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "prev opt vals " + now.strftime('%Y%m%d%H%M') + ".fasta");
	optimizer.writeBestParamsToText("prev opt vals " + now.strftime('%Y%m%d%H%M'))
	print(optimizer.getBestParameters()['match']);
	
	return None;

def TSonlyMacro(i:int, similarity:int):
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([20, 50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 1]);
	maxWeights = numpy.array([0, 0, 0, 0, 0, 1]);

	data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat" + str(i + 1) + ".tsv";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	optimizer = Optimizer(MACROSTATES);
	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readData(data);

	measure = "";
	if similarity == 0:
		search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " JSD";
	elif similarity == 1:
		search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " Cos";
	elif similarity == 2:
		search = CuckooSearch(optimizer.models, KLDivergence(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " KLD";
	else:
		search = CuckooSearch(optimizer.models, EntropyWeightsMixedSimilarity(CosineSimilarity(), JensenShannonDistance(), optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " Mix"
	search.setMaxIterations(2048);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(True, True, True, True, numpy.array([False, False, False, False, False, True]));
	optimizer.useAlgorithm(search);
	optimizer.optimize();

	name = "Macrostates " + str(i + 1) + measure + datetime.now().strftime('%Y%m%d%H%M');
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), name + ".fasta", 3);
	optimizer.writeBestParamsToText(name + ".txt");
	
	return None;

def TSonlyMicro():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([32, 64, 128, 256]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5]);
	steepnessRange = numpy.array([1, 7]);
	minWeights = numpy.array([0, 0, 0, 0, 0]);
	maxWeights = numpy.array([0, 0, 0, 0, 1]);

	data = "/netapp/home/tianjiao.zhang/data/microstates.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	
	optimizer = Optimizer(MACROSTATES, True);

	optimizer.readTargetFrequencies(targetFreqs);	
	optimizer.readFormattedMicrostateData(data);
		
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 128, 1, 0.25);
	search.setMaxIterations(4096);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(True, True, True, True, numpy.array([False, False, False, False, True]));
	optimizer.useAlgorithm(search);
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "TS microstates " + now.strftime('%Y%m%d%H%M') + ".fasta");
	optimizer.writeBestParamsToText("TS Microstates " + now.strftime('%Y%m%d%H%M'));
	

def long(generations:int):
	data = "/netapp/home/tianjiao.zhang/data/microstates.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([16, 24, 32, 64, 128]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5]);
	steepnessRange = numpy.array([1, 7]);
	minWeights = numpy.array([0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 1, 1, 1]);
	
	optimizer = Optimizer(MACROSTATES, True);
	optimizer.readTargetFrequencies(targetFreqs);	
	optimizer.readFormattedMicrostateData(data);
		
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 128, 1.25, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setAllSearchToTrue();
	search.suppressOutputs = True;
	optimizer.useAlgorithm(search);
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "var ensembles " + now.strftime('%Y%m%d%H%M') + ".fasta");
	optimizer.writeBestParamsToText("var ensembles " + now.strftime('%Y%m%d%H%M'));

def onlyMacro():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	i = 0;

	# only looking at MACROSTATE.TS
	# only optimizing backrub temperature and steepness
	ensembleSizes = numpy.array([20, 50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 0, 1, 1, 1]);

	data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat" + str(i + 1) + ".tsv";
	#data =  "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat" + str(i + 1) + ".tsv";
	#targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	optimizer = Optimizer(MACROSTATES);
	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readData(data);

	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 64, 1, 0.25);
	measure = " JSD";
	search.setMaxIterations(2048);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(True, True, True, True, numpy.array([True, True, False, True, True, True]));
	optimizer.useAlgorithm(search);
	optimizer.optimize();

	name = "Macrostates " + str(i + 1) + measure + datetime.now().strftime('%Y%m%d%H%M');
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), name + ".fasta", 3);
	optimizer.writeBestParamsToText(name + ".txt");
	
	return None;

def repeatTest():
	print("Hello!\n");
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 0, 1, 1, 1]);

	print("Initializing objects\n");

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataMicro = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";
	dataMicroAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";

	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");

	for i in range(32):
		search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 8, 1, 0.25);
		search.setMaxIterations(16);
		search.suppressOutputs = True;
		search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
		search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
		optimizer.useAlgorithm(search);

		print("\nJS Dist");
		#search.setSimilarityMeasure(JensenShannonDistance(optimizer.targetFrequencies));
		optimizer.optimize();
	
		params = optimizer.getBestParameters();
		print(params['match']);
		print(optimizer.verifyFoundParams(params['ensembleSize'], params['backrubTemp'], params['boltzmannTemp'], params['steepness'], params['weights']));

		search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), False, 8, 1, 0.25);
		search.setMaxIterations(16);
		search.suppressOutputs = True;
		search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
		search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
		optimizer.useAlgorithm(search);

		print("\nCosine")
		optimizer.optimize();
		params = optimizer.getBestParameters();
		print(params['match']);
		print(optimizer.verifyFoundParams(params['ensembleSize'], params['backrubTemp'], params['boltzmannTemp'], params['steepness'], params['weights']));

	return None;

def simpleRepeatTest():
	print("Hello!\n");
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 1]);
	maxWeights = numpy.array([0, 0, 0, 0, 0, 1]);

	print("Initializing objects\n");

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	dataMicro = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";
	dataMicroAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";

	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");

	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 1, 1, 0.25);
	search.setMaxIterations(1);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	#search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
	search.setAllSearchToFalse();
	search.suppressOutputs = True;
	optimizer.useAlgorithm(search);

	print("\nJS Dist");

	for i in range(64):
		optimizer.optimize();
		params = optimizer.getBestParameters();
		m = search.population[0];
		#print(search.similarityMeasure.getSimilarityMeasure(m.getFrequencies()));
		# TODO: getModelByParams doesn't always return the same object.
		m1 = Model.constructFromExisting(optimizer.getModelByParams(m.backrubTemp, m.ensembleSize, m.boltzmannTemp), m.ensembleSize, m.backrubTemp, m.boltzmannTemp, m.getWeights(), m.steepness);
		#print(search.similarityMeasure.getSimilarityMeasure(m1.getFrequencies()));
		if not m.equalTo(m1):
			print("\t{:s}".format(Optimizer.calcParamsID(m.backrubTemp, m.ensembleSize, m.boltzmannTemp)));
			print("\t{:s}".format(Optimizer.calcParamsID(m1.backrubTemp, m1.ensembleSize, m1.boltzmannTemp)));
		#m2 = Model.constructFromExisting(m, m.ensembleSize, m.backrubTemp, m.boltzmannTemp, m.getWeights(), m.steepness);
		#print(search.similarityMeasure.getSimilarityMeasure(m2.getFrequencies()));

		#print(m.equalTo(m2));

	#print(m2.backrubTemp);
	#print(m2.boltzmannTemp);
	#print(m2.ensembleSize);
	#print(m2.steepness);
	#print(m2.weights);
	#print(search.similarityMeasure.getSimilarityMeasure(m2.getFrequencies()));

	#print(numpy.sum(numpy.sum(numpy.abs( m.getFrequencies() - m2.getFrequencies()))));

	return None;

def testNewKLD(iterations = 64):
	print("Hello!\n");
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

	print("Initializing objects\n");

	#data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat1.tsv";
	data =  "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat1.tsv";
	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	#targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");

	search = CuckooSearch(optimizer.models, KLDivergence(optimizer.targetFrequencies), False, 16, 1, 0.25);
	search.setMaxIterations(iterations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
	optimizer.useAlgorithm(search);

	print("\nKL Dist");
	#search.setSimilarityMeasure(JensenShannonDistance(optimizer.targetFrequencies));
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "KLD new test " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("KLD new test " + now.strftime('%Y%m%d%H%M%S'));
	print(optimizer.getBestParameters()['match']);
	
	return None;

def testChi2(iterations = 64):
	print("Hello!\n");
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

	print("Initializing objects\n");

	data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat1.tsv";
	#data =  "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat" + str(i + 1) + ".tsv";
	#targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");

	search = CuckooSearch(optimizer.models, Chi2Kernel(optimizer.targetFrequencies), False, 64, 1, 0.25);
	search.setMaxIterations(iterations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True, False, True, True, True]));
	optimizer.useAlgorithm(search);

	print("\nChi2 kernel");
	#search.setSimilarityMeasure(JensenShannonDistance(optimizer.targetFrequencies));
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "Chi2 test " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("Chi2 test " + now.strftime('%Y%m%d%H%M%S'));
	print(optimizer.getBestParameters()['match']);
	
	return None;

def DHFRcomparemeasures(similarity:int):
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([20, 50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 0, 1, 1, 1]);

	data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat1.tsv";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	optimizer = Optimizer(MACROSTATES);
	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readData(data);

	measure = "";
	if similarity == 0:
		search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " JSD";
	elif similarity == 1:
		search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " Cos";
	elif similarity == 2:
		search = CuckooSearch(optimizer.models, KLDivergence(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " KLD";
	elif similarity == 3:
		search = CuckooSearch(optimizer.models, EntropyWeightsMixedSimilarity(CosineSimilarity(), JensenShannonDistance(), optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = " Mix"
	elif similarity == 4:
		search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(JensenShannonDistance(), optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = "Weighted JSD";
	else:
		search = CuckooSearch(optimizer.models, Chi2Kernel(optimizer.targetFrequencies), False, 64, 1, 0.25);
		measure = "Chi2 kernel";
	search.setMaxIterations(2048);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(True, True, True, True, numpy.array([True, True, False, True, True, True]));
	optimizer.useAlgorithm(search);
	optimizer.optimize();

	name = "DHFR compare measures " + measure + " " + datetime.now().strftime('%Y%m%d%H%M');
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), name + ".fasta", 3);
	optimizer.writeBestParamsToText(name + ".txt");
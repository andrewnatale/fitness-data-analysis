#/usr/bin/env python3
from Optimizer import Optimizer
from model import Model
from SimilarityMeasure import SimilarityMeasure
from SearchAlgorithm import SearchAlgorithm
from CuckooSearch import CuckooSearch
from KLDivergence import KLDivergence
from JensenShannonDistance import JensenShannonDistance
from CosineSimilarity import CosineSimilarity
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity
from EntropyWeightedSimilarity import EntropyWeightedSimilarity
from Chi2Kernel import Chi2Kernel
from enumeration import enum
#from datetime import *
import numpy
import threading
from io import *
import os
import sys

verbose = True

infile = sys.argv[1]
task_id = sys.argv[2]

#output_path = '/kortemmelab/home/anatale/opt3/output'
output_path = '/Users/anatale/school/UCSF/Kortemme_lab/code/testing'
try:
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise

#input_path = '/kortemmelab/home/anatale/opt3/opt_rnd3'
input_path = '/Users/anatale/school/UCSF/Kortemme_lab/code/multi-state-design/opt_rnd3'

# options to load from jobsfile based on task id:
# 1) targetFreqs_filename
# 2) data_filename
# 3) similarityMeasure (as code, possible: JS, CS, EW, EWM, KL, C2)
# 4) iterations
# 5) bool of states to search

#infile = os.path.join(input_path, 'optimizer_jobs_122_127.lst')
#infile = os.path.join(input_path, 'optimizer_jobs_101_109.lst')

with open(infile, 'r') as jobsfile:
    jobs = jobsfile.readlines()
jobs = [i.strip('\n') for i in jobs]

for job in jobs:
    job_params = job.split()
    job_id = job_params[0]
    job_tag = job_params[1]
    minPosition = int(job_params[2])
    targetFreqs_filename = job_params[3]
    data_filename = job_params[4]
    simMeas_id = job_params[5]
    iterations = int(job_params[6])
    #iterations = 4
    usedstates = numpy.array(job_params[7:]).astype(dtype=bool)
    #print(usedstates)
    if job_id == str(task_id):
        print('using options for job %s\n' % job_id)
        break

#targetFreqs = "/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/fitness-data-analysis/highscale_trim.fasta"
targetFreqs = os.path.join(input_path, targetFreqs_filename)
data = os.path.join(input_path, data_filename)
#data = "/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/fitness-data-analysis/testing_microstates.tsv"

MACROSTATES = enum('1i2m','1a2k','1k5d','3gj0','importin','composite')
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# parameter reanges to optimize
ensembleSizes = numpy.array([60,70,80,90,100])
backrubTemps = numpy.array([0.9])
#boltzmannTemps = numpy.array([-1.0]) # set below
steepnessRange = numpy.array([0.5, 5])
minWeights = numpy.array([0, 0, 0, 0, 0, 0])
maxWeights = numpy.array([1, 1, 1, 1, 1, 1])

optimizer = Optimizer(MACROSTATES)
optimizer.readTargetFrequencies(targetFreqs)
print('pos before loading data: ', optimizer.nPositions)
#optimizer.readData(data)
optimizer.readMicrostateData(data, minPosition=minPosition)
print('pos after loading data: ', optimizer.nPositions)
# examine model
for model_id in optimizer.models:
    if verbose:
        print('id:',model_id)
        print('obj:',optimizer.models[model_id])
        print(optimizer.models[model_id].nMacrostates)
        print(optimizer.models[model_id].ensembleSize)
        print(optimizer.models[model_id].backrubTemp)
        print(optimizer.models[model_id].boltzmannTemp)
        print(optimizer.models[model_id].weights)
        print(optimizer.models[model_id].steepness)
        print(optimizer.models[model_id].nPositions)
        print(optimizer.models[model_id].macrostatesUsed)
        print(optimizer.models[model_id].useMicrostateData)
        print(optimizer.models[model_id].positionOffset)
        print(optimizer.models[model_id].macrostateResidueEnergies.shape)
        print(optimizer.models[model_id].microstateResidueEnergies.shape)

    if optimizer.models[model_id].useMicrostateData:
        print('using microstates')
        microstate_optimize=True
        boltzmannTemps = numpy.array([-1.0, 5.0])
    else:
        print('using macrostates')
        microstate_optimize=False
        boltzmannTemps = numpy.array([-1.0])
print("Done reading data files");

if verbose:
    print(optimizer.targetFrequencies.shape)
    print(optimizer.nPositions)
    print(optimizer.minPosition)

def optimize():
    if microstate_optimize == True:
        # init search algorithm
        #search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 64, 1, 0.25)
        if simMeas_id == 'JS':
            search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 64, 1, 0.25)
        elif simMeas_id == 'CS':
            search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), True, 64, 1, 0.25)
        # elif simMeas_id == 'EW':
        #     search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(CosineSimilarity, optimizer.targetFrequencies), True, 64, 1, 0.25)
        # elif simMeas_id == 'EWM':
        #     search = CuckooSearch(optimizer.models, EntropyWeightsMixedSimilarity(optimizer.targetFrequencies), True, 64, 1, 0.25)
        elif simMeas_id == 'KL':
            search = CuckooSearch(optimizer.models, KLDivergence(optimizer.targetFrequencies), True, 64, 1, 0.25)
        elif simMeas_id == 'C2':
            search = CuckooSearch(optimizer.models, Chi2Kernel(optimizer.targetFrequencies), True, 64, 1, 0.25)

        search.setMaxIterations(iterations)
        # set parameters
        search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights)
        search.setSearchParameters(True, False, True, True, usedstates)
        # load search algorithm
        optimizer.useAlgorithm(search)
        # optimize
        optimizer.optimize()
        #now = datetime.now()
        optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), os.path.join(output_path, "var_ensembles_"+job_tag+".fasta"))
        optimizer.writeBestParamsToText(os.path.join(output_path, "var_ensembles_"+job_tag))
    else:
        # init search algorithm
        if simMeas_id == 'JS':
            search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 64, 1, 0.25)
        elif simMeas_id == 'CS':
            search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), False, 64, 1, 0.25)
        # elif simMeas_id == 'EW':
        #     search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(optimizer.targetFrequencies), False, 64, 1, 0.25)
        # elif simMeas_id == 'EWM':
        #     search = CuckooSearch(optimizer.models, EntropyWeightsMixedSimilarity(optimizer.targetFrequencies), False, 64, 1, 0.25)
        elif simMeas_id == 'KL':
            search = CuckooSearch(optimizer.models, KLDivergence(optimizer.targetFrequencies), False, 64, 1, 0.25)
        elif simMeas_id == 'C2':
            search = CuckooSearch(optimizer.models, Chi2Kernel(optimizer.targetFrequencies), False, 64, 1, 0.25)

        search.setMaxIterations(iterations)
        # set parameters
        search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights)
        search.setSearchParameters(False, False, False, True, usedstates)
        # load search algorithm
        optimizer.useAlgorithm(search)
        # optimize
        optimizer.optimize()
        #now = datetime.now()
        optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), os.path.join(output_path, "fixed_ensembles_"+job_tag+".fasta"))
        optimizer.writeBestParamsToText(os.path.join(output_path, "fixed_ensembles_"+job_tag))
    #print(optimizer.getBestParameters()['match'])

optimize()

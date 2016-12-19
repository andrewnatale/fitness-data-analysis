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

# sge_task_id = 1
# if os.environ.has_key("SGE_TASK_ID"):
#sge_task_id = os.getenv("SGE_TASK_ID", default=1)
sge_task_id = sys.argv[1]

#output_path = '/kortemmelab/home/anatale/multistate/output'
output_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/testing'
try:
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise

#input_path = '/kortemmelab/home/anatale/multistate'
input_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/multi-state-design'

# options to load from jobsfile based on task id:
# 1) targetFreqs_filename
# 2) data_filename
# 3) similarityMeasure (as code, possible: JS, CS, EW, EWM, KL, C2)
# 4) iterations

infile = os.path.join(input_path, 'optimizer_jobs.lst')

with open(infile, 'r') as jobsfile:
    jobs = jobsfile.readlines()
jobs = [i.strip('\n') for i in jobs]

for job in jobs:
    job_id, targetFreqs_filename, data_filename, simMeas_id, iterations = job.split()
    if job_id == str(sge_task_id):
        print('using options for job %s\n' % job_id)
        #iterations = int(iterations)
        iterations = 4
        break

#targetFreqs = "/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/fitness-data-analysis/highscale_trim.fasta"
targetFreqs = os.path.join(input_path, targetFreqs_filename)
data = os.path.join(input_path, data_filename)
#data = "/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/fitness-data-analysis/testing_microstates.tsv"

print("Hello!\n")
MACROSTATES = enum('apo_GEF','GTP_GAP','GDP','GTP_importin')
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# only optimizing backrub temperature and steepness
ensembleSizes = numpy.array([40,50,60])
backrubTemps = numpy.array([0.9])
#boltzmannTemps = numpy.array([-1.0]) # set below
steepnessRange = numpy.array([0.5, 5])
minWeights = numpy.array([0, 0, 0, 0])
maxWeights = numpy.array([1, 1, 1, 1])

optimizer = Optimizer(MACROSTATES)
optimizer.readTargetFrequencies(targetFreqs)
print('pos before loading data: ', optimizer.nPositions)
#optimizer.readData(data)
optimizer.readMicrostateData(data, minPosition=122)
print('pos after loading data: ', optimizer.nPositions)
# examine model
for model_id in optimizer.models:
    #print('id:',model_id)
    #print('obj:',optimizer.models[model_id])
    # print(optimizer.models[model_id].nMacrostates)
    # print(optimizer.models[model_id].ensembleSize)
    # print(optimizer.models[model_id].backrubTemp)
    # print(optimizer.models[model_id].boltzmannTemp)
    # print(optimizer.models[model_id].weights)
    # print(optimizer.models[model_id].steepness)
    # print(optimizer.models[model_id].nPositions)
    # print(optimizer.models[model_id].macrostatesUsed)
    # print(optimizer.models[model_id].useMicrostateData)
    # print(optimizer.models[model_id].positionOffset)
    # print(optimizer.models[model_id].macrostateResidueEnergies.shape)
    # print(optimizer.models[model_id].macrostateResidueEnergies.shape)

    if optimizer.models[model_id].useMicrostateData:
        print('using microstates')
        microstate_optimize=True
        boltzmannTemps = numpy.array([-1.0, 5.0])
    else:
        print('using macrostates')
        microstate_optimize=False
        boltzmannTemps = numpy.array([-1.0])
print("Done reading data files");

print(optimizer.targetFrequencies.shape)
print(optimizer.nPositions)
print(optimizer.minPosition)
#print(optimizer.maxPos)

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
        search.setSearchParameters(True, False, True, True, numpy.array([True, True, True, True]))
        optimizer.useAlgorithm(search)
        # load search algorithm
        optimizer.optimize()
        #now = datetime.now()
        optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), os.path.join(output_path, "var_ensembles_"+str(sge_task_id)+".fasta"))
        optimizer.writeBestParamsToText(os.path.join(output_path, "var_ensembles_"+str(sge_task_id)))
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
        search.setSearchParameters(False, False, False, True, numpy.array([True, True, True, True]))
        optimizer.useAlgorithm(search)
        # load search algorithm
        optimizer.optimize()
        #now = datetime.now()
        optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), os.path.join(output_path, "fixed_ensembles_"+str(sge_task_id)+".fasta"))
        optimizer.writeBestParamsToText(os.path.join(output_path, "fixed_ensembles_"+str(sge_task_id)))
    #print(optimizer.getBestParameters()['match'])

optimize()

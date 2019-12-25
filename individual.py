import random 
from math import exp, sqrt
import copy
import numpy as np
from util import smart_initialization

class Chromosome(object):

	def __init__(self, nodes): 
		self.length = len(nodes)  # length of route
		self.nodes = nodes  # list of nodes
		self.sequence = [] # list contining sequnece of nodes - the route
		self.fitness = 0 # fitness value of chromosome
		self.distance = 0 # total distance of route. Will be same as fitness if no diversity is used

	def create_individual(self, seed_val):
		"""
		Create route randomly
		"""
		random.seed(seed_val)
		self.sequence = random.sample(self.nodes, k = self.length)

	def create_individual_heur(self, seed_val, distanceMatrix):
		"""
		Create routes using nearest neighbour
		"""
		start_point = random.randint(1, self.length)
		self.sequence = smart_initialization(self.length+1, start_point, distanceMatrix)
		
		
	def fitness_func(self, iter, distanceMatrix , distances = 0):
		# alpha2 = 1  # required for using diversity
		
		dist = 0
		seq = copy.deepcopy(self.sequence)
		seq.append(seq[0])
		for i in range(len(seq)-1):
			if (seq[i] == seq[i+1]):
				print(i)
				print(seq)
			dist += distanceMatrix[seq[i], seq[i+1]]

		self.distance = dist

		# below part is used for including diversity in the fitness function
		distances = np.asarray(distances)

		self.fitness = dist
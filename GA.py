from individual import Chromosome
from nodes import Node
import numpy as np
import random
import copy
import datetime
import os

class GeneticAlgorithm(object):
	def __init__(self, mut_prob, cross_prob, pop_size, nodes, num_generation, elite_num, 
	distanceMatrix, seed, tsp_name, cutoff_time = 10, verbose = 0):
		self.mut_prob = mut_prob # mutation prob
		self.cross_prob = cross_prob # crossover prob
		self.pop_size = pop_size # population size
		# creating list of chromosomes objects
		self.Population = [Chromosome(nodes = nodes) for i in range(pop_size)]
		self.nodes = nodes # nodes
		self.length = len(nodes) # length of sequence
		self.distanceMatrix = distanceMatrix # distance matrix with distance between 2 nodes
		self.mutation_type = 'swap' # mutation type - options : 'swap'
		self.crossover_type = 'PMX' # crossover type - options: 'PMX', 'OX'
		self.crossover_points = 2 # crossover points - options: 1 , 2
		self.new_population = [] # list to store fitness of new population generated
		self.new_population_fitness = [] # list to store fitness of new population
		self.elite_num = elite_num # elite chromosomes count
		self.num_generation = num_generation # number of generations allowed to run
		self.fitness = [] # list to store fitness of current population
		self.verbose = verbose # print iteration result level - options: 0 , 1, 2
		self.seed = seed # seed value
		self.seed_vals = [] # list to store multiple seeds for creating different sequences
		self.distance = [] # list to distance. will be same as fitness if no diversity factor is included in fitness
		self.smart_pop = 15 # number of routes initialized using heuristics
		self.cutoff_time = cutoff_time # cutoff time - from input
		self.tsp_name = tsp_name # instance name - from input
		self.trace = {} # dicitonary to store the trace of optimization
		self.prev_best = 0 # stores previous best ovserved result
	
	def get_new_pop_count(self, zz):
		"""
		At each iteration few new randomly routes are generated.
		This count decreases as iteration increases (simulated annealing)
		This function returns that required count of new random routes to be generated
		"""
		return(round(self.pop_size*0.20/((zz+1)**0.2)))

	def init_seed_vals(self):
		"""
		create list of seeds required to crate new different initial sequence
		"""
		self.seed_vals = random.sample(range(100000), k = 350*6)

	def PopInitialize(self):
		"""
		Initial population generation. Two types of chromosomes are generated:
		1. randomly generated
		2. based on nearest beighbour
		"""

		# creating initial routes using heuristic
		smart_pop_count = self.smart_pop
		for i in range(smart_pop_count):
			self.Population[i].create_individual_heur(self.seed_vals[i], self.distanceMatrix)
			self.Population[i].fitness_func(0, self.distanceMatrix, list(np.zeros(self.pop_size)))

		# creating initial random routes
		for i in range(smart_pop_count, self.pop_size):
			self.Population[i].create_individual(self.seed_vals[i])
			self.Population[i].fitness_func(0, self.distanceMatrix, list(np.zeros(self.pop_size)))
		
	def get_individuals(self, type, zz, selection = 'roulette'):
		"""
		Selection step: creates the mating pool for crossover
		Input: type: mutation/crosover
		zz: iteration number
		selection: selection criteria: options: 'roulette', 'rank'
		"""
		fitness = []
		for i in range(self.pop_size):
			fitness.append(self.Population[i].fitness)
			
		if type == 'mutation':
			pass

		elif type == 'crossover':

			fitness_prob = np.asarray(fitness)
			if selection == 'roulette':

				# calculating fitness from the distance
				fitness_prob = 1/ fitness_prob
				fitness_prob = fitness_prob/sum(fitness_prob)

				idx = np.argsort(-fitness_prob)
				fitness_prob = -np.sort(-fitness_prob, kind = 'mergesort')

			else:
				# finding ranks of fitness values
				idx = np.argsort(fitness_prob)
				fitness_prob = list(reversed(range(1, len(fitness_prob)+1)))
				fitness_prob = np.asarray(fitness_prob)/sum(fitness_prob)

			# selecting new sequence based on probability
			fitness_prob = np.cumsum(fitness_prob)
			del_prob = np.append(np.asarray([0]),fitness_prob[0:len(fitness_prob)-1])

			# storing index of mating pool
			new_ind = []
			for i in range(self.pop_size):
				random_n = random.random()
				ind = (fitness_prob >= random_n) & (del_prob < random_n)
				ind = np.where(ind)[0]
				new_ind.append(idx[ind][0])

			new_ind = np.split(np.asarray(new_ind),2)

			return(list(zip(new_ind[0], new_ind[1])))
					
	def crossover(self, zz):
		"""
		Does crossover in the mating pool
		"""

		# this will store the new population
		self.new_population_fitness =[]

		# getting mating pool
		msk = self.get_individuals(type = 'crossover', zz = zz, selection = 'rank')
		cnt = 0

		# iterating over all the mating pool
		for i in msk:
			if random.random() <= self.cross_prob:
				ind1 = copy.deepcopy(np.asarray(self.Population[i[0]].sequence))
				ind2 = copy.deepcopy(np.asarray(self.Population[i[1]].sequence))
				
				if self.crossover_points == 2:
					cut = [random.randint(0,self.length), random.randint(0,self.length)]
				else:
					cut = [random.randint(0,self.length), self.length]

				cut.sort()
				if cut[0] == cut[1]:
					if cut[0] == len(self.nodes):
						cut[0]-=1
					else:
						cut[1]+=1
				cross_ind = list(range(cut[0], cut[1]))
				# these elements are taken from ind1 and ind2 to crossover
				from1 = set(ind1[cross_ind])
				from2 = set(ind2[cross_ind])

				# these indices are not cross-overed
				rem_ind = np.asarray(list(set(list(range(0, len(self.nodes)))) - set(list(range(cut[0], cut[1])))))


				if self.crossover_type == 'PMX':
					# these will be repeated in ind1 after crossover
					rep1 = list(from2-from1)
					# these will be repeated in ind2 after crossover
					rep2 = list(from1-from2)
					tmp = ind1[cross_ind]
					ind1[cross_ind] = ind2[cross_ind]
					ind2[cross_ind] = tmp

					# updating the new chromosomes in case they have repeated nodes
					for j in range(len(rep1)):
					    idx = np.intersect1d(np.where(ind1 == rep1[j])[0], rem_ind)
					    ind1[idx] = rep2[j]
					    idx = np.intersect1d(np.where(ind2 == rep2[j])[0], rem_ind)
					    ind2[idx] = rep1[j]

					new_c1 = Chromosome(nodes = self.nodes) 
					new_c1.sequence = list(ind1)
					new_c1.fitness_func(zz, self.distanceMatrix, self.distance)
					new_c2 = Chromosome(nodes = self.nodes) 
					new_c2.sequence = list(ind2)
					new_c2.fitness_func(zz, self.distanceMatrix, self.distance)

				elif self.crossover_type == 'OX':
					child1 = np.asarray([0]*self.length)
					child2 = np.asarray([0]*self.length)
					child1[cross_ind] = ind1[cross_ind]
					child2[cross_ind] = ind2[cross_ind]

					in1 = 0
					in2 = 0
					for k in range(self.length):
					    if ind2[k] not in child1:
					        child1[rem_ind[in1]] = ind2[k]
					        in1+=1
					    if ind1[k] not in child2:
					        child2[rem_ind[in2]] = ind1[k]
					        in2+=1

					new_c1 = Chromosome(nodes = self.nodes) 
					new_c1.sequence = list(child1)
					new_c1.fitness_func(zz, self.distanceMatrix, self.distance)
					new_c2 = Chromosome(nodes = self.nodes) 
					new_c2.sequence = list(child2)
					new_c2.fitness_func(zz, self.distanceMatrix, self.distance)

				self.new_population.append(new_c1)
				self.new_population.append(new_c2)
				self.new_population_fitness.append(new_c1.fitness)
				self.new_population_fitness.append(new_c2.fitness)

			else:
				self.new_population.append(self.Population[i[0]])
				self.new_population.append(self.Population[i[1]])
				self.new_population_fitness.append(self.Population[i[0]].fitness)
				self.new_population_fitness.append(self.Population[i[1]].fitness)

			cnt+=1

	def mutation(self, zz):
		"""
		Does mutation on the crossovered routes
		"""
		if self.mutation_type == 'swap':
			
			for i in range(len(self.new_population)):
			
				indi = copy.deepcopy(self.new_population[i].sequence)

				# using mutation prob to decide if mutation is required
				msk = np.random.rand(len(self.nodes)) < self.mut_prob
				msk = np.where(msk)[0]

				if(len(msk) > 0):
					for j in msk:
						swap1 = j
						swap2 = random.randint(0,self.length-1)

						# swapping the genes
						tmp = indi[swap1]
						indi[swap1] = indi[swap2]
						indi[swap2] = tmp

					new_m = Chromosome(nodes = self.nodes)
					new_m.sequence = indi
					new_m.fitness_func(zz, self.distanceMatrix, self.distance)

					self.new_population[i] = new_m

	def create_next_gen(self, zz, elite):
		"""
		Creates population from next generation:
		New population is made by three ways:
		1. Elite chromosomes from current gen are copied without doing crossover/mutation
		2. Top chormosomes after crossover/mutation are copied
		3. some random chromosomes are generated
		"""

		# deleting current populations to store next generation
		self.fitness = []
		distance = []
		copy_Population = copy.deepcopy(self.Population)
		self.Population = []

		# calculating fitness value of all the chromosomes
		fitness_values = []
		for i in range(self.pop_size):
			fitness_values.append(copy_Population[i].fitness)

		# storing elite chromosomes for next generation
		if elite > 0:
			elite_ind = np.argsort(fitness_values)[0:elite]
			for i in elite_ind:
				self.Population.append(copy.deepcopy(copy_Population[i]))
				self.fitness.append(copy_Population[i].fitness)
				distance.append(copy_Population[i].distance)

		# selecting top chromosomes from crossovered population
		new_pop_size = self.get_new_pop_count(zz)
		new_fitness_values = []
		for i in range(len(self.new_population)):
			new_fitness_values.append(self.new_population[i].fitness)

		new_ind = np.argsort(new_fitness_values)[0:(self.pop_size - elite - new_pop_size)]
		for i in new_ind:
			self.Population.append(copy.deepcopy(self.new_population[i]))
			self.fitness.append(self.new_population[i].fitness)
			distance.append(self.new_population[i].fitness)

		# creating random sequences for next generation and killing the bad ones from this one
		if new_pop_size > 0:
			for i in range(new_pop_size):
				new_rand = Chromosome(nodes = self.nodes) #, alpha = self.alpha, lamb = self.lamb)
				new_rand.create_individual(self.seed_vals[i])
				new_rand.fitness_func(0, self.distanceMatrix, self.distance)
				self.Population.append(copy.deepcopy(new_rand))
				self.fitness.append(self.new_population[i].fitness)
				distance.append(self.new_population[i].fitness)

		self.distance = distance

		self.new_population = []

	def fit(self):
		"""
		This function fits the GA by running iterations
		"""
		start_time = datetime.datetime.now()

		# creating list of seed values
		self.init_seed_vals()
		# Initialize the population
		self.PopInitialize()
		
		# running the generations
		for i in range(self.num_generation):
			
			# selection & crossover
			self.crossover(zz = i)
			# Mutation
			self.mutation(zz = i)
			# scoring and new generation cretion
			self.create_next_gen(zz = i, elite = self.elite_num)
			# printing generation result
			self.print(i, self.verbose)

			scores = np.asarray(self.fitness)
			best_score = min(scores)

			time_diff = (datetime.datetime.now()-start_time).total_seconds()
			# saving in trace
			if i == 0:
				self.prev_best = best_score
				self.trace[time_diff] = best_score
			else:
				if (self.prev_best > best_score):
					self.prev_best = best_score
					self.trace[time_diff] = best_score

			if (time_diff >= self.cutoff_time):
				break

	def print(self, zz, verbose):
		"""
		This function prints the result of each generation
		Input:
		zz: iteration number
		verbose: printing level (0: No print, 1: iteration number and best solution, 
								 2: Iteration number, best and average solution)
		"""
		if verbose == 1:
			scores = np.asarray(self.fitness)
			best_score = min(scores)
			print("Iteration: %d, Best score: %d" %(zz, best_score))
		elif verbose == 2:
			scores = np.asarray(self.fitness)
			best_score = min(scores)
			avg_score = np.mean(scores)

			print("Iteration: %d, Best score: %.d, Average score: %.2f" %(zz, best_score, avg_score))


	def write_results(self):
		"""
		Writes the solution and trace file
		"""
		best_score = self.Population[0].fitness
		best_seq = self.Population[0].sequence
		for i in range(1,len(self.Population)):
			if self.Population[i].fitness < best_score:
				best_score = self.Population[i].fitness
				best_seq = self.Population[i].sequence

		f1 = open(self.tsp_name + "_LS2_" + str(self.cutoff_time)+ "_" + str(self.seed) + ".sol", 'w+')
		f1.write(str(int(round(best_score))) + '\n')
		f1.write(str([t - 1 for t in best_seq]).replace(' ', '')[1:-1])
		f1.close()

		f2 = open(self.tsp_name + "_LS2_" + str(self.cutoff_time) + "_" + str(self.seed) + ".trace", 'w+')
		for soln in self.trace:
			f2.write(str(round(soln,2)) + ', ' + str(int(round(self.trace[soln]))) + '\n')
		f2.close()

def GA(f_path, cutoff_time, seed):
	"""
	Creates the class objects and fits the genetic algorithm
	"""

	# Reading the file and creating the distance matrix
	sd = datetime.datetime.now()
	read_nodes = 0
	node_list = []
	
	with open(f_path, 'r') as openfileobject:
		for line in openfileobject:
			if 'COMMENT' in line:
				node_count = int(line.split(' ')[1])
			if read_nodes > 0 and read_nodes <= node_count:
				temp_node = [float(x) for x in line.rstrip().split(' ')[0:]]
				node_list.append(Node(temp_node[0], temp_node[1], temp_node[2]))
				read_nodes += 1
			if 'NODE_COORD_SECTION' in line:
				read_nodes = 1
	
	distanceMatrix = {}
	for i in range(len(node_list)-1):
		for j in range(i+1, len(node_list)):
			row = node_list[i].id
			col = node_list[j].id
			distanceMatrix[row,col] = node_list[i].distance(node_list[j])
			distanceMatrix[col,row] = distanceMatrix[row,col]

	tsp_name = os.path.basename(f_path)[0:-4]

	random.seed(seed)
	np.random.seed(seed)

	# creating the generic algorithm object
	GA = GeneticAlgorithm(0.1, 0.8, 100, list(range(1,len(node_list)+1)), 3000, 20, distanceMatrix, seed,
		tsp_name, cutoff_time)
	# Running GA
	GA.fit()
	# write result
	GA.write_results()
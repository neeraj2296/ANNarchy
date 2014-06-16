#
# 	A network comparable to Dinkelbach et al. 2012
#
from ANNarchy import *

setup(num_threads=1)
# setup( paradigm = "cuda", num_threads=1 )

# Defining the neurons
InputNeuron = RateNeuron(
    equations="""
        rate = 0.2 : init = 0.1 
    """
)

OutputNeuron = RateNeuron(
    equations="""
        rate = sum(exc)
    """
)

NEURON = 25
CONN = 10

input_pop = Population(geometry=(NEURON), neuron=InputNeuron)
output_pop = Population(geometry=(2), neuron=OutputNeuron)

proj = Projection(input_pop, output_pop, 'exc').connect_fixed_number_pre(CONN, 1.0)
#proj = Projection(input_pop, output_pop, 'exc').connect_all_to_all(1.0)

compile()  # needed to save connectivity matrix



from ANNarchy.core.Global import _projections, _populations

# computation of the adjacency list of the neural network graph
num_vertices = sum(pop.size for pop in _populations)

# contains for each population the neurons and their pre- and post-synaptic neurons (adjacent neurons)
projections_edges = {}

# create an empty array of adjacencies for each neuron of each population
for pop in _populations:
    projections_edges[pop] = [[] for x in xrange(pop.size)]

# iterate over each projection and add the adjacencies
for proj in _projections:
    for pre in xrange(proj.pre.size):
        projections_edges[proj.pre][pre] += [(proj.post, post) for post in xrange(proj.size()) if pre in proj.dendrite(post).rank]
    for post in xrange(proj.size()):
        projections_edges[proj.post][post] += [(proj.pre, pre) for pre in proj.dendrite(post).rank]

# number of edges in total in the neural network
num_edges = sum([sum([len(e) for e in edges]) for edges in projections_edges.values()]) / 2

# calculate an dictionary of offsets, for translating the population-local ranks of each neuron to a global rank of all neurons
population_offsets = {}
for i, population in enumerate(_populations):
    population_offsets[population] = sum([pop.size for pop in _populations[:i]])

# debug output
print num_vertices, num_edges
print population_offsets

# write the adjacency list into a file for METIS
filename = 'metis_input.csv'
with open(filename, mode='w') as w_file:
    # first two numbers are the vertex count and the edge count
    w_file.write(str(num_vertices) + ' ' + 
                 str(num_edges) + '\n')
    # write for every neuron the adjacent neurons in a new line
    for projection in projections_edges.values():
        for neuron in projection:
            for edge in neuron:
                # add the offset of the corresponding population to the local rank
                # add 1, because METIS uses indices starting at 1
                w_file.write(str(population_offsets[edge[0]] + edge[1] + 1) + ' ')
            w_file.write('\n')

# call METIS with the created input file and the number of partitions
numPartitions = 3
from subprocess import call
call(["gpmetis", filename, str(numPartitions)])

# read the partitioning output from METIS
# and create a csv file containing the population-number, the x, y and z coordinates in the population and the partition of the neuron
# the csv file can be used by drawPartitioning.m to generate a visualization
with open(filename + '.part.' + str(numPartitions) + '.csv', mode='w') as w_file:
    for i, part in enumerate(open(filename + '.part.' + str(numPartitions))):
        popNum = max([j for j, p in enumerate(_populations) if i >= population_offsets[p] and (i + 1 >= len(_populations) or i < population_offsets[_populations[i + 1]])])
        pop = _populations[popNum]
        coordinates = list(pop.coordinates_from_rank(i - population_offsets[pop]))
        for cl in xrange(len(coordinates), 3):
            coordinates += [0]
        w_file.write(', '.join([str(popNum)] + [str(c) for c in coordinates] + [str(int(part))]) + '\n')

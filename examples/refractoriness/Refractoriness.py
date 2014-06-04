#
#   ANNarchy - SimpleSTDP
#
#   A simple model showing the usage of refractoriness in ANNarchy
# 
#   See e.g. Mainen & Sejnowski (1995) for experimental results in vitro.
#
#   Code adapted from the Brian example: https://brian2.readthedocs.org/en/latest/examples/non_reliability.html
#
#   authors: Helge Uelo Dinkelbach, Julien Vitay
#
from ANNarchy import *

Neuron = SpikeNeuron(
parameters = """
    tau = 20.0 : population
    sigma = 0.015 : population
    target = 1.1 : population
""",
equations = """
    noise =  sqrt( 2.0 * tau ) * Normal(0.0, sigma) 
    tau * dx/dt + x = target + noise
""",
spike = """
    x > 1
""",
reset = """
    x = 0
""",
refractory = Uniform(1,5)
)

pop = Population( geometry=25, neuron = Neuron )

compile()

pop.start_record('spike')
simulate ( 500.0 )
data = pop.get_record()

spikes = raster_plot(data['spike'])

# Plot the results
import pylab as plt
plt.plot(spikes[:, 0], spikes[:, 1], '.')
plt.show()
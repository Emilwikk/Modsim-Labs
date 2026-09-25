import config
import utils
from matplotlib import pyplot as plt

# Create and Populate Network 
sfn = utils.SFN(config.V_dan,config.gamma_dan)
sfn.populate(S=config.S_dan, I=config.I_dan)

for _ in range(config.steps):
    sfn.move_individuals()

plt.figure()
plt.title("Number of Susceptible Individuals Per Node Over Time")
for n, node in enumerate(sfn.nodes):
    plt.plot(range(config.steps+1), node.S_history, label=f"Node {n}")
plt.legend()
plt.axis([0,config.steps+1,0,400])
plt.xlabel("Time []")
plt.ylabel("Quantity")
plt.show()



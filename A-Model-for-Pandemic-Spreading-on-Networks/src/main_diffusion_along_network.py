import config
import utils

# Create and Populate Network 
sfn = utils.SFN(config.V_dan,config.gamma_dan)
sfn.populate(S=config.S_dan, I=config.I_dan)

print(sfn.nodes[0].S)
print(sfn.nodes[0].I)
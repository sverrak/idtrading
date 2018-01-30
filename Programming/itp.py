# Status: Incomplete
# Description: This file is an implementation of the MSILP

### ----------- External Packages -----------
from gurobipy import *

### ----------- Model Initialization -----------
model = Model('Intraday_Trading_Problem')


### ----------- System Parameters -----------
print_output = True
default_parameters = True
parameter_file = "itp_parameters.txt"


### ----------- Model Parameters -----------
if (default_parameters == True):
  # Set model parameters
  params = 0
  number_of_trading_stages = 12


else:
  # Retrieve parameters from file
  params = read_parameters(parameter_file)




### ----------- Datastructures -----------

# Inspiration (The quantities of the Refinery Example Problem)
if(False):  
  crude_numbers = range(1,2+1)
  petrols = ["Premium_fuel", "Regular_fuel"]

  crude_bounds = {1:20000, 2:30000}
  lb_lube_oil = 500
  ub_lube_oil = 1000

  max_crude = 45000
  max_reform = 10000
  max_cracking = 8000

  distillation_splitting_coefficients = {"Light_naphtha": (0.1, 0.15), "Medium_naphtha": (0.2, 0.25), "Residuum": (0.13, 0.12)}
  reforming_splitting_coefficients = {"Light_naphtha": 0.6, "Medium_naphtha":0.52, "Heavy_naphtha":0.45}

  end_product_profit = {"Premium_fuel":7, "Regular_fuel":6, "Jet_fuel":4, "Fuel_oil":3.5, "Lube_oil":1.5}
  blending_coefficients = {"Light_oil": 0.55, "Heavy_oil": 0.17, "Cracked_oil": 0.22, "Residuum": 0.055}

  lube_oil_factor = 0.5
  pmf_rmf_ratio = 0.4

  vapor_pressure_constants = [0.6, 1.5, 0.05]


# Intraday Trading Problem 
trading_stages = range(number_of_trading_stages)                                        # Set of trading timeslots
trading_stages_remaining =  [range(number_of_trading_stages) - trading_stages[i] 
                            for i in range(number_of_trading_stages)]                   # Set of trading timeslots
volume_options = [5]#[0,1,2,3,5,10,15,20,25]                                            # Set of volume options (from linearization of cost convexity)
max_bid_price = 9999.00
min_bid_price = -9999.00
p_bm_pos = 9999.00
p_bm_neg = 9999.00
v_da = 35
max_q = 100
min_q = 30
c_c = 0.01
c_p = 40
M_2 = 80




### ----------- Variables -----------

# Inspiration (The quantities of the Refinery Example Problem)
if(False):
  crudes = model.addVars(crude_numbers, ub=crude_bounds, name="cr")    
  end_products = model.addVars(end_product_names, name="end_prod")
  end_products["Lube_oil"].lb= lb_lube_oil
  end_products["Lube_oil"].ub= ub_lube_oil
  distillation_products = model.addVars(distillation_products_names, name="dist_prod")
  reform_usage = model.addVars(naphthas, name="napthas_to_reformed_gasoline")
  reformed_gasoline = model.addVar(name="reformed_gasoline")
  cracking_usage = model.addVars(intermediate_oils,name="intermediate_oils_to_cracked_gasoline")
  cracking_products = model.addVars(cracking_products_names,  name="cracking_prods")
  used_for_regular_motor_fuel = model.addVars(used_for_motor_fuel_names, name="motor_fuel_to_regular_motor_fuel")
  used_for_premium_motor_fuel = model.addVars(used_for_motor_fuel_names, name="motot_fuel_to_premium_motor_fuel")
  used_for_jet_fuel = model.addVars(used_for_jet_fuel_names, name="jet_fuel")
  used_for_lube_oil = model.addVar(vtype=GRB.CONTINUOUS,name="residuum_used_for_lube_oil")


# Intraday Trading Problem

# Production
q = model.addVar(vtype=GRB.CONTINUOUS,lb=min_q,ub=max_q, name="production_quantity")    # Production quantity

# Bid features         
b_p = model.addVars(trading_stages, name="bid_price")                                   # Bid price
b_v = model.addVars(trading_stages, lb=0, name="bid_volume")                            # Bid volume
b_vr = [[0 for k in range(trading_stages[i]+i)] for i in range(trading_stages)]         # Residual bid volume

# Note: here, the idea is to get the residual volume variables on an upper triangular form. Index tuning is probably needed for this to work properly. 

# Some extra care is needed to initialize the residual volume
for i in range(len(trading_stages)):
    for k in range(trading_stages[i]):
        b_vr[i+1][k+i] = model.addVar(vtype=GRB.CONTINUOUS,                             # Residual volume
                                      name="residual_volume_%s" % str(str(i+1)+str(i+1+k)))

# Transaction features
for i in range(len(trading_stages)):
    for k in range(trading_stages[i]):
        p_c[i+1][k+i] = model.addVar(vtype=GRB.CONTINUOUS,                              # Transaction price
                                    name="p_c%s" % str(str(i+1)+str(i+1+k)))

for i in range(len(trading_stages)):
    for k in range(trading_stages[i]):
        v_c[i+1][k+i] = model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=b_vr[i+1][k+i],       # Transaction volume
                                    name="v_c%s" % str(str(i+1)+str(i+1+k)))

v_bm_pos = model.addVar(vtype=GRB.CONTINUOUS,name="v_bm_pos")                           # Balancing market+ volume
v_bm_neg = model.addVar(vtype=GRB.CONTINUOUS,name="v_bm_neg")                           # Balancing market- volume


# Support variables
is_bid_placed = [[0 for k in range(trading_stages[i]+i)] for i in range(trading_stages)]                                                     # Bid is placed
is_bid_cleared = [[0 for k in range(trading_stages[i]+i)] for i in range(trading_stages)]                                                     # Bid is cleared
is_bid_killed = [[0 for k in range(trading_stages[i]+i)] for i in range(trading_stages)]                                                   # Bid is killed

for i in range(len(trading_stages)):
    for k in range(trading_stages[i]):
        is_bid_placed[i+1][k+i] = model.addVar(vtype=GRB.BINARY,                       # Is bid placed (d^p)
                                    name="is_bid_placed%s" % str(str(i+1)+str(i+1+k)))

for i in range(len(trading_stages)):
    for k in range(trading_stages[i]):
        is_bid_cleared[i+1][k+i] = model.addVar(vtype=GRB.BINARY,                      # Is bid cleared (d^c)
                                    name="is_bid_cleared_%s" % str(str(i+1)+str(i+1+k)))

for i in range(len(trading_stages)):
    for k in range(trading_stages[i]):
        is_bid_killed[i+1][k+i] = model.addVar(vtype=GRB.BINARY,                       # Is bid killed (d^k)
                                    name="is_bid_killed_%s" % str(str(i+1)+str(i+1+k)))

### ----------- Constraints -----------
# Inspiration (The quantities of the Refinery Example Problem)
if(False):
  # Max Crude
  model.addConstr(crudes.sum() <= max_crude, "max_crude")

  # Splitting
  model.addConstrs((quicksum(distillation_splitting_coefficients[dpn][crude-1]*crudes[crude] 
    for crude in crudes) == distillation_products[dpn]  
    for dpn in distillation_products_names), "splitting_distillation")

  # Reforming
  model.addConstr(reform_usage.prod(reforming_splitting_coefficients) == reformed_gasoline,
                      "splitting_reforming")

  # Cracking
  model.addConstrs((quicksum(cracking_splitting_coefficients[oil, crack_prod]*cracking_usage[oil]
                             for oil in intermediate_oils) == cracking_products[crack_prod]
                    for crack_prod in cracking_products_names),
                   name="splitting_cracking")

  # Blending
  model.addConstr(used_for_regular_motor_fuel["Reformed_gasoline"] +
                  used_for_premium_motor_fuel["Reformed_gasoline"] ==
                  reformed_gasoline, "continuity_reformed_gasoline")

  # Vapour pressure
  model.addConstr(used_for_jet_fuel["Light_oil"] +
                  vapor_pressure_constants[0]*used_for_jet_fuel["Heavy_oil"] +
                  vapor_pressure_constants[1]*used_for_jet_fuel["Cracked_oil"] +
                  vapor_pressure_constants[2]*used_for_jet_fuel["Residuum"] <= end_products["Jet_fuel"],
                  "vapour_pressure")

# Intraday Trading Problem
# Physical Constraints

# Financial Constraints
model.addConstr(sum(sum(x) for x in v_c) 
                  + v_bm_neg - v_bm_pos - q == -v_da, "sold_equals_generated")

# Question: are the indices right? When should the residual volume be set to zero?
mode.addConstrs(M_2*(is_bid_killed[t_b][t]-1)+b_vr[t_b][t:].sum() <= 0
                for t in range(len(trading_stages_remaining[t_b])) 
                for t_b in range(len(trading_stages)),
                "no_transaction_after_kill")

model.addConstrs(is_bid_killed[t_b][t+1] - is_bid_killed[t_b][t] >= 0 
                for t in range(len(is_bid_killed[t_b])) 
                for t_b in range(len(is_bid_killed)-1),
                "propagate_killing")

model.addConstrs(v_c[t_b][:t].sum() + b_vr[t_b][t] - b_v[t_b] == 0
                for t in range(len(v_c[t_b])) for t_b in range(len(v_c)),
                "set_residual_volume")




### ----------- Objective Function -----------
# Inspiration:

model.setObjective(end_products.prod(end_product_profit), GRB.MAXIMIZE)

# Question: Is it possible to write the objective function in two lines?
z = p_c.prod(v_c) - c_c*v_c - c_p*q - (p_bm_pos.prod(v_bm_pos) - p_bm_neg(v_bm_neg))
model.setObjective(z, GRB.MAXIMIZE)


### ----------- Optimization -----------
model.optimize()



### ----------- Output Results -----------
for v in model.getVars():
    if v.X != 0:
        print("%s %f" % (v.Varname, v.X))

model.write("refinery-output.sol")

### ----------- Support Functions -----------

# To be implemented
def read_parameters(filename):
  return 0

using JuMP
using Ipopt
#using Compat
#using Clp

#READING THE SYSTEM DATA
system_name = "IEEE14"

#SHOW INPUT DATA
IData = 0 # 1 == show; 0 == No show

include("ACPBF_dat.jl")

#SETTING THE SOLVER
solver = IpoptSolver()
m = Model(solver = solver)
#max_iter=3

#INCLUDYNG LIBRARIES
include("Filias_Library.jl")

#READING THE MODEL
include("ACPBF_model.jl")

#PRINTING THE RESULTS
include("ACPBF_inc.jl")
using JuMP
using Mosek
#READING THE SYSTEM DATA
system_name = "IEEE3BusSystem"
#SHOW INPUT DATA
IData = 1 # 1 == show; 0 == No show
#INCLUDING DATA
include("Flow_General_data.jl")

#SETTING THE SOLVER
solver = MosekSolver()
# #INCLUDE MODEL
include("Flow_Bai_model.jl")
#INCLUDE PRINT
# EXPORTDATA = 1  # 1 == EXPORT; 0 == NO EXPORT (.mat)
# include("Flow_General_inc.jl")

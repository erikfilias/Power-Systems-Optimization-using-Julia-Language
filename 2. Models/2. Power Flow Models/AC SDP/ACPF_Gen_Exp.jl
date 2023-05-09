# =============================================================================#
#Description: Export Results .mat extension
#Date: 08.14.2019
#Version: 0.0.1.081419
#Model: Non Lineal Rectangular
# =============================================================================#
#                             Library Import                                   #
# =============================================================================#

results=Dict{String, Array{Float64}}()
results["ResultFO"]=ResultFinal
results["BUS"]=BusOut
results["BRANCH"]=BranchOut
results["PowTo"]=PowerTotal
results["LossTo"]=LossTotal
matwrite("ACFP_$ACModel"*"_$file_name.mat", results)

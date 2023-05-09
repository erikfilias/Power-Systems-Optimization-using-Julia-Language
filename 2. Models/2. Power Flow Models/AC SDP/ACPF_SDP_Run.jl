# =============================================================================#
#Description: File to Run || Main Program
#Date: 08.14.2019
#Version: 0.0.1.081419
#Model: SDP
# =============================================================================#
#                             Change Directory                                 #
# =============================================================================#
#pwd():  you could know the current directory
#clearconsole()
# path_mp  = "$(homedir())/Documents/JuliaProjects/01. PowerFlow/01. N-Lineal"
# path_dat = "$(homedir())/Documents/JuliaProjects/01. PowerFlow/00. Data"
# path_res = "$(homedir())/Documents/JuliaProjects/01. PowerFlow/04. Results"
path_mp  = "$(homedir())/Dropbox/TutorialSDP/01. PowerFlow/01. N-Lineal"
path_dat = "$(homedir())/Dropbox/TutorialSDP/01. PowerFlow/00. Data"
path_res = "$(homedir())/Dropbox/TutorialSDP/01. PowerFlow/04. Results"
# =============================================================================#
#                              Library Import
# =============================================================================#
using JuMP
using Mosek
using MosekTools
using Ipopt
# using Cbc
# using CPLEX
using AmplNLWriter
using Printf
using CSV
# using MathOptInterface
# using MathOptInterfaceMosek
using LinearAlgebra
using Match
using MAT

# =============================================================================#
#                                File name                                     #
# =============================================================================#
file_name   = "IEEE118"
ACModel     = "SDP"
FO          = 0
# =============================================================================#
#                                Data import                                   #
# =============================================================================#
cd(path_dat)
include("ACPF_Gen_Dat.jl")
# =============================================================================#
#                                Model import                                  #
# =============================================================================#
cd(path_mp)
model= Model(with_optimizer(Mosek.Optimizer))
include("ACPF_SDP_Mod.jl")
status=optimize!(model)
# =============================================================================#
#                              Print Results                                   #
# =============================================================================#
include("ACPF_SDP_Out.jl")
include("ACPF_Gen_Inc.jl")
# =============================================================================#
#                              Export Results                                  #
# =============================================================================#
cd(path_res)
include("ACPF_Gen_Exp.jl")

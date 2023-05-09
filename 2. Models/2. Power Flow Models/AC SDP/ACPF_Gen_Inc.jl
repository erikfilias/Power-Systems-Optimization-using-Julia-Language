
println(" *** SUMMARY OF RESULTS  ***")
println("――――――――――――――――――――――――――――")
# Printing results of Buses
FOstring_0 = "Objective value (Power Generation in Slack Bus)"
FOstring_1 = "Objective value (Total Power Generation)"
FOstring_2 = "Objective value (Losses)"
if FO ==0;FOstring = FOstring_0
elseif FO ==1;FOstring = FOstring_1
else FOstring = FOstring_2
end
print(FOstring)
@printf ": (→): %5.2f \n" float(ResultFinal[1,1])
println("――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
println("BUS RESULTS (↓):")
println("――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
println("  Bus  Bus_Id   Type   V0(pu)    Th0(°)    Pg(MW)    Qg(MVar)   Pd(MW)    Qd(MVar) gshb(MW)  bshb(MVar)")
println("――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
for i in 1:nbus
    @printf "%5d"      float(BusOut[i,1])
    @printf "%6d"      float(BusOut[i,2])
    @printf "%7.d"     float(BusOut[i,3])
    @printf "%10.2f"   float(BusOut[i,4])
    @printf "%10.2f"   float(BusOut[i,5])
    @printf "%10.2f"   float(BusOut[i,6])
    @printf "%10.2f"   float(BusOut[i,7])
    @printf "%10.2f"   float(BusOut[i,8])
    @printf "%10.2f"   float(BusOut[i,9])
    @printf "%10.2f"   float(BusOut[i,10])
    @printf "%10.2f\n" float(BusOut[i,11])
end
println("―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
@printf "Total %42.2f" float(PowerTotal[1,1])
@printf "%10.2f"       float(PowerTotal[1,2])
@printf "%10.2f"       float(PowerTotal[1,3])
@printf "%10.2f"       float(PowerTotal[1,4])
@printf "%10.2f"       float(PowerTotal[1,5])
@printf "%10.2f \n"    float(PowerTotal[1,6])
println("―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n")

# Printing results of Branch
println("―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
println("BRANCH RESULTS (↓):")
println("―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
println("Branch From  To    Pde(MW)   Ppa(MW)   Qde(MW)   Qpa(MW)    Pls(MW) Qls(MVAr)")
println("―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")

for i in 1:nbranch
    @printf "%5d"      float(BranchOut[i,1])
    @printf "%5d"      float(BranchOut[i,2])
    @printf "%5d"      float(BranchOut[i,3])
    @printf "%10.2f"   float(BranchOut[i,4])
    @printf "%10.2f"   float(BranchOut[i,5])
    @printf "%10.2f"   float(BranchOut[i,6])
    @printf "%10.2f"   float(BranchOut[i,7])
    @printf "%10.2f"   float(BranchOut[i,8])
    @printf "%10.2f\n" float(BranchOut[i,9])
end
println("―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――")
@printf "Total %59.2f" float(LossTotal[1,1])
@printf "%10.2f \n"    float(LossTotal[1,2])

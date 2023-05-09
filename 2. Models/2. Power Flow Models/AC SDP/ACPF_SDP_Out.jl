#==========================================================#
#Description: Print Results
#Model: Polar Flow Power
#Date: 09.05.2018
#Version: beginning v0
#Name: Hever Alcahuaman
#==========================================================#
# Reults Function Objetive
ResultFinal = [JuMP.objective_value(model)*Sbase]
# Matrix results Bus - Rectangular method
BusOut = zeros(nbus,11)
BusOut[:,1] = Bus.busLoc
BusOut[:,2] = Bus.busId
BusOut[:,3] = Bus.bustype
BusOut[:,8] = Bus.Pd*Sbase
BusOut[:,9] = Bus.Qd*Sbase
for i in 1:nbus
	if value(x5[i+i-1,i+i])>=0
		BusOut[i,5] = atan((value(x5[i+i,i+i]))^0.5/value(x5[i+i-1,i+i-1])^0.5)
		BusOut[i,4] = (value(x5[i+i-1,i+i-1])+value(x5[i+i,i+i]))^0.5
	elseif value(x5[i+i-1,i+i])<=0
		BusOut[i,5] = atan(-(value(x5[i+i,i+i]))^0.5/value(x5[i+i-1,i+i-1])^0.5)
		BusOut[i,4] = (value(x5[i+i-1,i+i-1])+value(x5[i+i,i+i]))^0.5
	end
    BusOut[i,6] = value(Pg[i])*Sbase
    BusOut[i,7] = value(Qg[i])*Sbase
    BusOut[i,10] = Sbase*Bus.gshb[i]* BusOut[i,4]^2
    BusOut[i,11] = Sbase*Bus.bshb[i]* BusOut[i,4]^2
end

# Matrix results Branch - Rectangular method
BranchOut = zeros(nbranch,9)
BranchOut[:,1] = Branch.branchLoc
BranchOut[:,2] = Branch.from
BranchOut[:,3] = Branch.to
for i in 1:nbranch
    BranchOut[i,4] = Sbase*value(Pde[i])
    BranchOut[i,5] = Sbase*value(Ppa[i])
    BranchOut[i,6] = Sbase*value(Qde[i])
    BranchOut[i,7] = Sbase*value(Qpa[i])
    BranchOut[i,8] = Sbase*(value(Pde[i])+value(Ppa[i]))
    BranchOut[i,9] = Sbase*(value(Qde[i])+value(Qpa[i]))
end


# Calculating total power
PowerTotal = zeros(1,6)
PowerTotal[1,1]  = sum(BusOut[i,6] for i in 1:nbus)   #PgT
PowerTotal[1,2]  = sum(BusOut[i,7] for i in 1:nbus)   #QgT
PowerTotal[1,3]  = sum(BusOut[i,8] for i in 1:nbus)   #PdT
PowerTotal[1,4]  = sum(BusOut[i,9] for i in 1:nbus)   #QdT
PowerTotal[1,5]  = sum(BusOut[i,10] for i in 1:nbus)  #gshbT
PowerTotal[1,6]  = sum(BusOut[i,11] for i in 1:nbus)  #bshbT

# Calculating total power Loss
LossTotal = zeros(1,2)
LossTotal[1,1] = sum(BranchOut[i,8] for i in 1:nbranch)   #PlsT
LossTotal[1,2]  = sum(BranchOut[i,9] for i in 1:nbranch)  #QlsT

@printf "---------------------------------------------------------------------------------------------\n"
@printf "-------------------------------------------RESULTS-------------------------------------------\n"
@printf "---------------------------------------------------------------------------------------------\n"
println("Objective value (Losses): ", getobjectivevalue(m)*Sbase)
@printf "---------------------------------------------------------------------------------------------\n"
@printf "                                         BUS__RESULTS                                        \n"
@printf "---------------------------------------------------------------------------------------------\n"
@printf "   Bus  Type    V0       Th0       Pg        Qg        Pd        Qd        gshb      bshb    \n"
@printf "---------------------------------------------------------------------------------------------\n"
for i in 1:nbus
    @printf "%5d" float(i)
    @printf "%5d" float(Bus.bustype[i])
    @printf "%10.4f" float(V1[i])
    @printf "%10.4f" float(th1[i]*180/3.14159265359)
    @printf "%10.4f" float(Pg1[i]*Sbase)
    @printf "%10.4f" float(Qg1[i]*Sbase)
    @printf "%10.4f" float(Bus.Pd[i]*Sbase)
    @printf "%10.4f" float(Bus.Qd[i]*Sbase)
    @printf "%10.4f" float(Sbase*Bus.gshb[i]*V1[i])
    @printf "%10.4f" float(Sbase*Bus.bshb[i]*V1[i])
    @printf "\n"
end
@printf "---------------------------------------------------------------------------------------------\n"
@printf "TOTAL"
a=0
b=0
c=0
d=0
e=0
f=0
for i in 1:nbus
    a = a + Pg1[i]
    b = b + Qg1[i]
    c = c + Bus.Pd[i]
    d = d + Bus.Qd[i]
    e = e + Bus.gshb[i]*V1[i]^2
    f = f + Bus.bshb[i]*V1[i]^2
end
@printf "%35.4f" float(Sbase*a)
@printf "%10.4f" float(Sbase*b)
@printf "%10.4f" float(Sbase*c)
@printf "%10.4f" float(Sbase*d)
@printf "%10.4f" float(Sbase*e)
@printf "%10.4f" float(Sbase*f)
@printf "\n"
@printf "---------------------------------------------------------------------------------------------\n"

@printf "-----------------------------------------------------------------------------\n"
@printf "                               BRANCH__RESULTS                               \n"
@printf "-----------------------------------------------------------------------------\n"
@printf " Branch From  To  Pij[MW]   Pji[MW]   Qij[MW]   Qji[MW]    Pls[MW] Qls[MVAr] \n"
@printf "-----------------------------------------------------------------------------\n"
for i in 1:nbranch
    @printf "%5d" float(i)
    @printf "%5d" float(Branch.from[i])
    @printf "%5d" float(Branch.to[i])
    @printf "%10.4f" float(Sbase*Pde1[i])
    @printf "%10.4f" float(Sbase*Ppa1[i])
    @printf "%10.4f" float(Sbase*Qde1[i])
    @printf "%10.4f" float(Sbase*Qpa1[i])
    @printf "%10.4f" float(Sbase*(Pde1[i]+Ppa1[i]))
    @printf "%10.4f" float(Sbase*(Qde1[i]+Qpa1[i]))
    @printf "\n"
end

@printf "-----------------------------------------------------------------------------\n"
@printf "TOTAL"
# using MAT
# matwrite("results5_$system_name.mat", results)

@printf "---------------------------------------------------------------------------------------------\n"
@printf "-------------------------------------------RESULTS-------------------------------------------\n"
@printf "---------------------------------------------------------------------------------------------\n"
println("Objective value (Loss): ", getobjectivevalue(m)*Sbase)
@printf "---------------------------------------------------------------------------------------------\n"
@printf "                                         BUS__RESULTS                                        \n"
@printf "---------------------------------------------------------------------------------------------\n"
@printf "   Bus  Type    V0       Th0       Pg        Qg        Pd        Qd        gshb      bshb    \n" 
@printf "---------------------------------------------------------------------------------------------\n"
for i in 1:nbus
    @printf "%5d" float(i)
    @printf "%5d" float(Bus.bustype[i])
    @printf "%10.4f" float(sqrt(getvalue(Vsqr[i])))
    @printf "%10.4f" float(getvalue(th[i])*180/3.14159265359)
    @printf "%10.4f" float(getvalue(Pg[i])*Sbase)
    @printf "%10.4f" float(getvalue(Qg[i])*Sbase)
    @printf "%10.4f" float(Bus.Pd[i]*Sbase)
    @printf "%10.4f" float(Bus.Qd[i]*Sbase)
    @printf "%10.4f" float(Sbase*Bus.gshb[i]*getvalue(Vsqr[i]))
    @printf "%10.4f" float(Sbase*Bus.bshb[i]*getvalue(Vsqr[i]))
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
    a = a + getvalue(Pg[i])
    b = b + getvalue(Qg[i])
    c = c + Bus.Pd[i]
    d = d + Bus.Qd[i]
    e = e + Bus.gshb[i]*getvalue(Vsqr[i])
    f = f + Bus.bshb[i]*getvalue(Vsqr[i])
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
#for k in 1:nbus
#    for l in 1:nbus
#        for i in 1:nbranch
#            if k == Branch.from[i]
#                if l == Branch.to[i]
#                    @printf "%5d" float(i)
#                    @printf "%5d" float(Branch.from[i])
#                    @printf "%5d" float(Branch.to[i])
#                    @printf "%10.4f" float(Sbase*(getvalue(P[i])+Branch.r[i]*getvalue(Isqr[i])))
#                    @printf "%10.4f" float(-Sbase*getvalue(P[i]))
#                    @printf "%10.4f" float(Sbase*(getvalue(Q[i])+Branch.x[i]*getvalue(Isqr[i])-Branch.c[i]*getvalue(Vsqr[k])))
#                    @printf "%10.4f" float(-Sbase*(getvalue(Q[i])+Branch.c[i]*getvalue(Vsqr[l])))
#                    @printf "%10.4f" float(Sbase*Branch.r[i]*getvalue(Isqr[i]))
#                    @printf "%10.4f" float(Sbase*(Branch.x[i]*getvalue(Isqr[i])-Branch.c[i]*getvalue(Vsqr[k])-Branch.c[i]*getvalue(Vsqr[l])))
#                    @printf "\n"
#                end
#            end
#        end 
#    end
#end
for i in 1:nbranch
    @printf "%5d" float(i)
    @printf "%5d" float(Branch.from[i])
    @printf "%5d" float(Branch.to[i])
    @printf "%10.4f" float(Sbase*(getvalue(P[i])+Branch.r[i]*getvalue(Isqr[i])))
    @printf "%10.4f" float(-Sbase*getvalue(P[i]))
    @printf "%10.4f" float(Sbase*(getvalue(Q[i])+Branch.x[i]*getvalue(Isqr[i])-Branch.c[i]*getvalue(Vsqr[Branch.from[i]])))
    @printf "%10.4f" float(-Sbase*(getvalue(Q[i])+Branch.c[i]*getvalue(Vsqr[Branch.to[i]])))
    @printf "%10.4f" float(Sbase*Branch.r[i]*getvalue(Isqr[i]))
    @printf "%10.4f" float(Sbase*(Branch.x[i]*getvalue(Isqr[i])-Branch.c[i]*getvalue(Vsqr[Branch.from[i]])-Branch.c[i]*getvalue(Vsqr[Branch.to[i]])))
    @printf "\n"
end

@printf "-----------------------------------------------------------------------------\n"
@printf "TOTAL"
a = 0
b = 0
for k in 1:nbus
    for l in 1:nbus
        for i in 1:nbranch
            if k == Branch.from[i]
                if l == Branch.to[i]
                    a = a + Branch.r[i]*getvalue(Isqr[i])
                    b = b + Branch.x[i]*getvalue(Isqr[i]) - Branch.c[i]*getvalue(Vsqr[k]) - Branch.c[i]*getvalue(Vsqr[l])
                end
            end
        end
    end
end
@printf "%60.4f" float(Sbase*a)
@printf "%10.4f" float(Sbase*b)
@printf "\n"
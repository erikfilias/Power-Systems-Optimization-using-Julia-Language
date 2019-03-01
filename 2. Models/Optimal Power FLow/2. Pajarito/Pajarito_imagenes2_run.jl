using JuMP
using Pajarito
using Mosek
using CPLEX
using DelimitedFiles
#READING THE SYSTEM DATA
system_name = "IEEE3BusSystem"
#SHOW INPUT DATA
IData = 0 # 1 == show; 0 == No show
#INCLUDING DATA
include("Flujo_General_data.jl")

nmax=[1 0 0]
Reg=[0.1 0 0]
nstep=[0 0 5]
###Parameters-----------------------------------------------------------------
    nopt= 2*findnz(nmax)[3]+1
    ncomb=prod(nopt)
    nij=Dict()
    aij=Dict()
    for i=1:nbranch
        nij[i]=zeros(1,2*nmax[i]+1)
        aij[i]=zeros(1,2*nmax[i]+1)
        contador=0
        if TrafoLocation[i]!=0
            for j=1:size(nij[i])[2]
                nij[i][j]=nij[i][j]-findnz(nmax)[3][i]+contador

                aij[i][j]=1+findnz(Reg)[3][i]*(nij[i][j]/findnz(nmax)[3][i])
                contador=contador+1
            end
        else
            aij[i]=1
        end
    end
    Branchz=Branch.r+Branch.x*im
    Branchy=1 ./Branchz
    Y=Dict()
    for i=1:ncomb
        Y[i]=zeros(nbus,nbus)
        Y[i]=complex(Y[i])
    end
    branchaij=collect(Base.product(aij[1], aij[2],aij[3])) #Falta mejorar esta parte
    for i=1:ncomb
        for k=1:nbranch
        Y[i][convert(Int, Branchfrom[k]),convert(Int, Branchto[k])] = Y[i][convert(Int, Branchfrom[k]),convert(Int, Branchto[k])] + Branchy[k]*branchaij[i][k]

        Y[i][convert(Int, Branchto[k]),convert(Int, Branchfrom[k])] = Y[i][convert(Int, Branchfrom[k]),convert(Int, Branchto[k])]
        end
    end
    for k=1:ncomb

        for i=1:nbus
            for j=1:nbranch
                if Branchfrom[j]==i
                    Y[k][i,i] = Y[k][i,i] + Branchy[j] * branchaij[k][j]^2 + Branch.c[j]*im
                elseif Branchto[j]==i
                    Y[k][i,i] = Y[k][i,i] + Branchy[j] + Branch.c[j]*im
                end
             end
        end
    end

X=Dict()

function flowbai(l,nbus,Busnumber,Bus,Branch,Branchto,Branchfrom,out_lines,in_lines,Y)
        #SETTING THE SOLVER
        solver = MosekSolver(LOG=0)

        # #model
        mod=Model(solver=solver)

        # Group of buses:
        @variable(mod, x5[1:2*nbus,1:2*nbus],SDP)

        @SDconstraint(mod, x5>= 0)

        # Simple variable of injected reactive power
        @variable(mod, Qshunt[Busnumber])
        for i in 1:nbus
            setvalue(Qshunt[i], 0)
        end

        # Formulacion Compensation Shunt
        for k=1:nbus
            if Bus.bshb[k]!=0
                auxshunt=zeros(2*nbus,2*nbus)
                auxshunt[k+k-1,k+k-1]   =Bus.bshb[k]
                auxshunt[k+k,k+k]       =Bus.bshb[k]
                @constraint(mod,Qshunt[k]==vecdot(auxshunt,x5))
            else
                @constraint(mod,Qshunt[k]==0)
            end
        end
        # Matricial variable of Active Power
        @variable(mod, Pg[Busnumber])

        mtemporal=  zeros(2*nbus,2*nbus)
        mtemporal1= zeros(2*nbus,2*nbus)

        mtemporal2= zeros(2*nbus,2*nbus)
        mtemporal3= zeros(2*nbus,2*nbus)

        WGp= zeros(2*ngen,2*ngen)
        WGq= zeros(2*ngen,2*ngen)

        WBp= zeros(2*nbus,2*nbus)
        WBq= zeros(2*nbus,2*nbus)

        for k=1:nbus
            if Bus.bustype[k] !=0
                #Encerar matrices

                mtemporal=  zeros(2*nbus,2*nbus)
                mtemporal1= zeros(2*nbus,2*nbus)

                mtemporal2= zeros(2*nbus,2*nbus)
                mtemporal3= zeros(2*nbus,2*nbus)

                WGp= zeros(2*ngen,2*ngen)
                WGq= zeros(2*ngen,2*ngen)

                WBp= zeros(2*nbus,2*nbus)
                WBq= zeros(2*nbus,2*nbus)

                mtemporal[k+k-1,k+k-1]= -2*real(Y[l][k,k])
                mtemporal[k+k,k+k]=     -2*real(Y[l][k,k])

                mtemporal2[k+k-1,k+k-1]= 2*imag(Y[l][k,k])
                mtemporal2[k+k,k+k]=     2*imag(Y[l][k,k])

                a=convert(Int,(size(out_lines[k]))[1])

                for i=1:a
                    j=convert(Int,(Branchto[out_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end

                b=convert(Int,(size(in_lines[k]))[1])

                for i=1:b
                    j=convert(Int,(Branchfrom[in_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]= real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=     real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]= real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=     real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]=  -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=   imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]=  -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=   imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end
                WBp=mtemporal1+mtemporal
                WBq=mtemporal2+mtemporal3

                    @constraint(mod,Pg[k]==-vecdot((0.5)*WBp,x5)+Bus.Pd[k])

                    @constraint(mod,vecdot((-0.5)*WBp,x5)+Bus.Pd[k]<=Bus.Pmax[k])
                    @constraint(mod,vecdot((-0.5)*WBp,x5)+Bus.Pd[k]>=Bus.Pmin[k])
                    @constraint(mod,-Qshunt[k]+vecdot((-0.5)*WBq,x5)+Bus.Qd[k]<=Bus.Qmax[k])
                    @constraint(mod,-Qshunt[k]+vecdot((-0.5)*WBq,x5)+Bus.Qd[k]>=Bus.Qmin[k])

            elseif Bus.bustype[k] ==0

                mtemporal=  zeros(2*nbus,2*nbus)
                mtemporal1= zeros(2*nbus,2*nbus)

                mtemporal2= zeros(2*nbus,2*nbus)
                mtemporal3= zeros(2*nbus,2*nbus)

                WBp= zeros(2*nbus,2*nbus)
                WBq= zeros(2*nbus,2*nbus)

                mtemporal[k+k-1,k+k-1]= -2*real(Y[l][k,k])
                mtemporal[k+k,k+k]=     -2*real(Y[l][k,k])

                mtemporal2[k+k-1,k+k-1]= 2*imag(Y[l][k,k])
                mtemporal2[k+k,k+k]=     2*imag(Y[l][k,k])

                a=convert(Int,(size(out_lines[k]))[1])

                for i=1:a
                    j=convert(Int,(Branchto[out_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end

                b=convert(Int,(size(in_lines[k]))[1])

                for i=1:b
                    j=convert(Int,(Branchfrom[in_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end
                WBp=mtemporal1+mtemporal
                WBq=mtemporal2+mtemporal3


                @constraint(mod,Pg[k]==0)

                @constraint(mod,vecdot((0.5)*WBp,x5)==Bus.Pd[k])
                @constraint(mod,Qshunt[k]+vecdot((0.5)*WBq,x5)==Bus.Qd[k])


            end
        end

         @objective(mod,Min,sum(Pg[k] for k=1:nbus if Bus.bustype[k]!=0 ))

        # 2.Constraint of reference bus:
        for i=1:nbus
            if Bus.bustype[i] !=0
                @constraint(mod,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vg[i]^2)
            end
            if Bus.bustype[i] ==3
                @constraint(mod,x5[i+i,i+i]==0)
            end
        end
        # 4.Limits of amplitude at each bus:
        WV=zeros(2*nbus,2*nbus)
        Wf=zeros(2*nbus,2*nbus)
        Winf=zeros(2*nbus,2*nbus)
        Wsup=zeros(2*nbus,2*nbus)

        for k=1:nbus

            WV=zeros(2*nbus,2*nbus)
            Winf=zeros(2*nbus,2*nbus)
            Wsup=zeros(2*nbus,2*nbus)

            WV[k+k-1,k+k-1]=1
            WV[k+k,k+k]=1

            @constraint(mod,vecdot(WV,x5)<=1.1^2)
            @constraint(mod,vecdot(WV,x5)>=0.95^2)
        end

        status=solve(mod)
        # print(mod)
        x=getvalue(x5)

            if status==:Optimal
                A=getvalue(x5)
                for i=1:2*nbus
                    for j=1:2*nbus
                        if A[i,j]<=0.0001
                            A[i,j]=0
                        else
                        end
                    end
                end
                autovalores=eigvals(A)
                autovalormax=maximum(autovalores)
                Matrixrank=0
                indicador=0
                for i=1:length(autovalores)
                    indicador=autovalores[i]/autovalormax
                    if indicador>=0.0001
                        Matrixrank=Matrixrank+1
                    else
                    end
                end

            else
                Matrixrank=1000
            end

            return(x,Matrixrank)
end

function flowbai2(l,nbus,Busnumber,Bus,Branch,Branchto,Branchfrom,out_lines,in_lines,Y)
        # SETTING THE SOLVER
        mip_solver_drives = true
        rel_gap = 1e-5
        # using CPLEX
        mip_solver = CplexSolver(
            CPX_PARAM_SCRIND=(mip_solver_drives ? 1 : 0),
            # CPX_PARAM_SCRIND=1,
            CPX_PARAM_EPINT=1e-8,
            CPX_PARAM_EPRHS=1e-7,
            CPX_PARAM_EPGAP=(mip_solver_drives ? 1e-5 : 1e-9)
        )
        # using SCS
        # conic_solver = SCSSolver(eps=1e-6, max_iters=1000000, verbose=0)
        # using Mosek
        conic_solver = MosekSolver(LOG=0)
        #MODEL
        mod=Model(solver=PajaritoSolver(
            mip_solver_drives=mip_solver_drives,
            log_level=3,
            rel_gap=rel_gap,
        	mip_solver=mip_solver,
        	cont_solver=conic_solver,
        ))


        # Group of buses:
        @variable(mod, x5[1:2*nbus,1:2*nbus],SDP)

        @SDconstraint(mod, x5>= 0)

        # Simple variable of injected reactive power
        @variable(mod, Qshunt[Busnumber])
        for i in 1:nbus
            setvalue(Qshunt[i], 0)
        end
        epsilon=Dict()
        for i=1:nbus
            if Bus.bshb[i]!=0
                epsilon[i]=@variable(mod,[1],Bin)
            end
        end

        # Formulacion Compensation Shunt
        for k=1:nbus
            if Bus.bshb[k]!=0
                auxshunt=zeros(2*nbus,2*nbus)
                auxshunt[k+k-1,k+k-1]   =Bus.bshb[k]
                auxshunt[k+k,k+k]       =Bus.bshb[k]

                @constraint(mod,Qshunt[k]>=-epsilon[k][1]*Bus.bshb[k]*1.1^2)
                @constraint(mod,Qshunt[k]<=epsilon[k][1]*Bus.bshb[k]*1.1^2)

                @constraint(mod,Qshunt[k]-vecdot(auxshunt,x5)>=-(1-epsilon[k][1])*Bus.bshb[k]*1.1^2)
                @constraint(mod,Qshunt[k]-vecdot(auxshunt,x5)<=(1-epsilon[k][1])*Bus.bshb[k]*1.1^2)

            else
                @constraint(mod,Qshunt[k]==0)

            end
        end
        # Matricial variable of Active Power
        @variable(mod, Pg[Busnumber])

        mtemporal=  zeros(2*nbus,2*nbus)
        mtemporal1= zeros(2*nbus,2*nbus)

        mtemporal2= zeros(2*nbus,2*nbus)
        mtemporal3= zeros(2*nbus,2*nbus)

        WGp= zeros(2*ngen,2*ngen)
        WGq= zeros(2*ngen,2*ngen)

        WBp= zeros(2*nbus,2*nbus)
        WBq= zeros(2*nbus,2*nbus)

        for k=1:nbus
            if Bus.bustype[k] !=0
                #Encerar matrices

                mtemporal=  zeros(2*nbus,2*nbus)
                mtemporal1= zeros(2*nbus,2*nbus)

                mtemporal2= zeros(2*nbus,2*nbus)
                mtemporal3= zeros(2*nbus,2*nbus)

                WGp= zeros(2*ngen,2*ngen)
                WGq= zeros(2*ngen,2*ngen)

                WBp= zeros(2*nbus,2*nbus)
                WBq= zeros(2*nbus,2*nbus)

                mtemporal[k+k-1,k+k-1]= -2*real(Y[l][k,k])
                mtemporal[k+k,k+k]=     -2*real(Y[l][k,k])

                mtemporal2[k+k-1,k+k-1]= 2*imag(Y[l][k,k])
                mtemporal2[k+k,k+k]=     2*imag(Y[l][k,k])

                a=convert(Int,(size(out_lines[k]))[1])

                for i=1:a
                    j=convert(Int,(Branchto[out_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end

                b=convert(Int,(size(in_lines[k]))[1])

                for i=1:b
                    j=convert(Int,(Branchfrom[in_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]= real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=     real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]= real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=     real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]=  -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=   imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]=  -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=   imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end
                WBp=mtemporal1+mtemporal
                WBq=mtemporal2+mtemporal3

                    @constraint(mod,Pg[k]==-vecdot((0.5)*WBp,x5)+Bus.Pd[k])

                    @constraint(mod,vecdot((-0.5)*WBp,x5)+Bus.Pd[k]<=Bus.Pmax[k])
                    @constraint(mod,vecdot((-0.5)*WBp,x5)+Bus.Pd[k]>=Bus.Pmin[k])
                    @constraint(mod,-Qshunt[k]+vecdot((-0.5)*WBq,x5)+Bus.Qd[k]<=Bus.Qmax[k])
                    @constraint(mod,-Qshunt[k]+vecdot((-0.5)*WBq,x5)+Bus.Qd[k]>=Bus.Qmin[k])

            elseif Bus.bustype[k] ==0

                mtemporal=  zeros(2*nbus,2*nbus)
                mtemporal1= zeros(2*nbus,2*nbus)

                mtemporal2= zeros(2*nbus,2*nbus)
                mtemporal3= zeros(2*nbus,2*nbus)

                WBp= zeros(2*nbus,2*nbus)
                WBq= zeros(2*nbus,2*nbus)

                mtemporal[k+k-1,k+k-1]= -2*real(Y[l][k,k])
                mtemporal[k+k,k+k]=     -2*real(Y[l][k,k])

                mtemporal2[k+k-1,k+k-1]= 2*imag(Y[l][k,k])
                mtemporal2[k+k,k+k]=     2*imag(Y[l][k,k])

                a=convert(Int,(size(out_lines[k]))[1])

                for i=1:a
                    j=convert(Int,(Branchto[out_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end

                b=convert(Int,(size(in_lines[k]))[1])

                for i=1:b
                    j=convert(Int,(Branchfrom[in_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end
                WBp=mtemporal1+mtemporal
                WBq=mtemporal2+mtemporal3


                @constraint(mod,Pg[k]==0)

                @constraint(mod,vecdot((0.5)*WBp,x5)==Bus.Pd[k])
                @constraint(mod,Qshunt[k]+vecdot((0.5)*WBq,x5)==Bus.Qd[k])


            end
        end

         @objective(mod,Min,sum(Pg[k] for k=1:nbus if Bus.bustype[k]!=0 ))

        # 2.Constraint of reference bus:
        for i=1:nbus
            if Bus.bustype[i] !=0
                @constraint(mod,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vg[i]^2)
            end
            if Bus.bustype[i] ==3
                @constraint(mod,x5[i+i,i+i]==0)
            end
        end
        # 4.Limits of amplitude at each bus:
        WV=zeros(2*nbus,2*nbus)
        Wf=zeros(2*nbus,2*nbus)
        Winf=zeros(2*nbus,2*nbus)
        Wsup=zeros(2*nbus,2*nbus)

        for k=1:nbus

            WV=zeros(2*nbus,2*nbus)
            Winf=zeros(2*nbus,2*nbus)
            Wsup=zeros(2*nbus,2*nbus)

            WV[k+k-1,k+k-1]=1
            WV[k+k,k+k]=1

            @constraint(mod,vecdot(WV,x5)<=1.1^2)
            @constraint(mod,vecdot(WV,x5)>=0.95^2)
        end

        status=solve(mod)

    x=getvalue(x5)
    LD=Dict()
    #eX
    for i=1:nbus
        if Bus.bshb[i]!=0
            LD[i]=getvalue(epsilon[i][1])
        else
            LD[i]="WITHOUT"
        end
    end

        if status==:Optimal
            A=getvalue(x5)
            for i=1:2*nbus
                for j=1:2*nbus
                    if A[i,j]<=0.0001
                        A[i,j]=0
                    else
                    end
                end
            end
            autovalores=eigvals(A)
            autovalormax=maximum(autovalores)
            Matrixrank=0
            indicador=0
            for i=1:length(autovalores)
                indicador=autovalores[i]/autovalormax
                if indicador>=0.0001
                    Matrixrank=Matrixrank+1
                else
                end
            end

        else
            Matrixrank=1000
        end

    return(x,LD,Matrixrank)
end

function flowbai3(l,nbus,Busnumber,Bus,Branch,Branchto,Branchfrom,out_lines,in_lines,Y,nstep)
        # SETTING THE SOLVER
        mip_solver_drives = true
        rel_gap = 1e-5
        # using CPLEX
        mip_solver = CplexSolver(
            CPX_PARAM_SCRIND=(mip_solver_drives ? 1 : 0),
            # CPX_PARAM_SCRIND=1,
            CPX_PARAM_EPINT=1e-8,
            CPX_PARAM_EPRHS=1e-7,
            CPX_PARAM_EPGAP=(mip_solver_drives ? 1e-5 : 1e-9)
        )
        # using SCS
        # conic_solver = SCSSolver(eps=1e-6, max_iters=1000000, verbose=0)
        # using Mosek
        conic_solver = MosekSolver(LOG=0)
        #MODEL
        mod=Model(solver=PajaritoSolver(
            mip_solver_drives=mip_solver_drives,
            log_level=3,
            rel_gap=rel_gap,
        	mip_solver=mip_solver,
        	cont_solver=conic_solver,
        ))


        # Group of buses:
        @variable(mod, x5[1:2*nbus,1:2*nbus],SDP)

        @SDconstraint(mod, x5>= 0)

        #Auxiliar variables for Reactiva Compensation
        @variable(mod, Qshunt[Busnumber])
        for i in 1:nbus
        	setvalue(Qshunt[i], 0)
        end
        Qp=Dict()
        for i=1:nbus
            if QShuntLocation[i]!=0
                Qp[i]=@variable(mod, [1:nstep[i]])
            else
                Qp[i]=0
            end
        end
        nb=Dict()
        for i=1:nbus
            if QShuntLocation[i]!=0
                nb[i]=@variable(mod, [1:nstep[i]],Bin)
            else
                nb[i]=0
            end
        end

        for i=1:nbus
            if QShuntLocation[i]!=0
                @constraint(mod, sum(nb[i][k] for k=1:nstep[i])>=0)
                @constraint(mod, sum(nb[i][k] for k=1:nstep[i])<=nstep[i])
            end
        end

        for i=1:nbus
            if QShuntLocation[i]!=0
                for k=2:nstep[i]
                @constraint(mod, nb[i][k] <= nb[i][k-1])
                end
            end
        end

        # Formulacion Compensation Shunt
        for k=1:nbus
            if Bus.bshb[k]!=0
            auxshunt=zeros(2*nbus,2*nbus)
            auxshunt[k+k-1,k+k-1]   =Bus.bshb[k]
            auxshunt[k+k,k+k]       =Bus.bshb[k]

                for j=1:nstep[k]
                @constraint(mod, Qp[k][j]-(x5[k+k-1,k+k-1]+x5[k+k,k+k])>=-1.1^2*(1 - nb[k][j]))
                @constraint(mod, Qp[k][j]-(x5[k+k-1,k+k-1]+x5[k+k,k+k])<=0.9^2 *(1 - nb[k][j]))

                @constraint(mod, nb[k][j] * 0.95^2 <= Qp[k][j])
                @constraint(mod, Qp[k][j] <= 1.1^2 * nb[k][j])

                end
                @constraint(mod,Qshunt[k]==(Bus.bshb[k]/nstep[k])*sum(Qp[k][j] for j=1:nstep[k]))

            else
            @constraint(mod,Qshunt[k]==0)

            end
        end

        # Matricial variable of Active Power
        @variable(mod, Pg[Busnumber])

        mtemporal=  zeros(2*nbus,2*nbus)
        mtemporal1= zeros(2*nbus,2*nbus)

        mtemporal2= zeros(2*nbus,2*nbus)
        mtemporal3= zeros(2*nbus,2*nbus)

        WGp= zeros(2*ngen,2*ngen)
        WGq= zeros(2*ngen,2*ngen)

        WBp= zeros(2*nbus,2*nbus)
        WBq= zeros(2*nbus,2*nbus)

        for k=1:nbus
            if Bus.bustype[k] !=0
                #Encerar matrices

                mtemporal=  zeros(2*nbus,2*nbus)
                mtemporal1= zeros(2*nbus,2*nbus)

                mtemporal2= zeros(2*nbus,2*nbus)
                mtemporal3= zeros(2*nbus,2*nbus)

                WGp= zeros(2*ngen,2*ngen)
                WGq= zeros(2*ngen,2*ngen)

                WBp= zeros(2*nbus,2*nbus)
                WBq= zeros(2*nbus,2*nbus)

                mtemporal[k+k-1,k+k-1]= -2*real(Y[l][k,k])
                mtemporal[k+k,k+k]=     -2*real(Y[l][k,k])

                mtemporal2[k+k-1,k+k-1]= 2*imag(Y[l][k,k])
                mtemporal2[k+k,k+k]=     2*imag(Y[l][k,k])

                a=convert(Int,(size(out_lines[k]))[1])

                for i=1:a
                    j=convert(Int,(Branchto[out_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end

                b=convert(Int,(size(in_lines[k]))[1])

                for i=1:b
                    j=convert(Int,(Branchfrom[in_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]= real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=     real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]= real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=     real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]=  -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=   imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]=  -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=   imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end
                WBp=mtemporal1+mtemporal
                WBq=mtemporal2+mtemporal3

                    @constraint(mod,Pg[k]==-vecdot((0.5)*WBp,x5)+Bus.Pd[k])

                    @constraint(mod,vecdot((-0.5)*WBp,x5)+Bus.Pd[k]<=Bus.Pmax[k])
                    @constraint(mod,vecdot((-0.5)*WBp,x5)+Bus.Pd[k]>=Bus.Pmin[k])
                    @constraint(mod,-Qshunt[k]+vecdot((-0.5)*WBq,x5)+Bus.Qd[k]<=Bus.Qmax[k])
                    @constraint(mod,-Qshunt[k]+vecdot((-0.5)*WBq,x5)+Bus.Qd[k]>=Bus.Qmin[k])

            elseif Bus.bustype[k] ==0

                mtemporal=  zeros(2*nbus,2*nbus)
                mtemporal1= zeros(2*nbus,2*nbus)

                mtemporal2= zeros(2*nbus,2*nbus)
                mtemporal3= zeros(2*nbus,2*nbus)

                WBp= zeros(2*nbus,2*nbus)
                WBq= zeros(2*nbus,2*nbus)

                mtemporal[k+k-1,k+k-1]= -2*real(Y[l][k,k])
                mtemporal[k+k,k+k]=     -2*real(Y[l][k,k])

                mtemporal2[k+k-1,k+k-1]= 2*imag(Y[l][k,k])
                mtemporal2[k+k,k+k]=     2*imag(Y[l][k,k])

                a=convert(Int,(size(out_lines[k]))[1])

                for i=1:a
                    j=convert(Int,(Branchto[out_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end

                b=convert(Int,(size(in_lines[k]))[1])

                for i=1:b
                    j=convert(Int,(Branchfrom[in_lines[k]])[i])

                    mtemporal1[k+k-1,j+j-1]=real(Y[l][k,j])
                    mtemporal1[k+k,j+j]=    real(Y[l][k,j])
                    mtemporal1[j+j-1,k+k-1]=real(Y[l][k,j])
                    mtemporal1[j+j,k+k]=    real(Y[l][k,j])

                    mtemporal1[k+k-1,j+j]= -imag(Y[l][k,j])
                    mtemporal1[k+k,j+j-1]=  imag(Y[l][k,j])
                    mtemporal1[j+j,k+k-1]= -imag(Y[l][k,j])
                    mtemporal1[j+j-1,k+k]=  imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j-1]=-imag(Y[l][k,j])
                    mtemporal3[k+k,j+j]=    -imag(Y[l][k,j])
                    mtemporal3[j+j-1,k+k-1]=-imag(Y[l][k,j])
                    mtemporal3[j+j,k+k]=    -imag(Y[l][k,j])

                    mtemporal3[k+k-1,j+j]=  -real(Y[l][k,j])
                    mtemporal3[k+k,j+j-1]=   real(Y[l][k,j])
                    mtemporal3[j+j,k+k-1]=  -real(Y[l][k,j])
                    mtemporal3[j+j-1,k+k]=   real(Y[l][k,j])

                end
                WBp=mtemporal1+mtemporal
                WBq=mtemporal2+mtemporal3


                @constraint(mod,Pg[k]==0)

                @constraint(mod,vecdot((0.5)*WBp,x5)==Bus.Pd[k])
                @constraint(mod,Qshunt[k]+vecdot((0.5)*WBq,x5)==Bus.Qd[k])


            end
        end

         @objective(mod,Min,sum(Pg[k] for k=1:nbus if Bus.bustype[k]!=0 ))

        # 2.Constraint of reference bus:
        for i=1:nbus
            if Bus.bustype[i] !=0
                @constraint(mod,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vg[i]^2)
            end
            if Bus.bustype[i] ==3
                @constraint(mod,x5[i+i,i+i]==0)
            end
        end
        # 4.Limits of amplitude at each bus:
        WV=zeros(2*nbus,2*nbus)
        Wf=zeros(2*nbus,2*nbus)
        Winf=zeros(2*nbus,2*nbus)
        Wsup=zeros(2*nbus,2*nbus)

        for k=1:nbus

            WV=zeros(2*nbus,2*nbus)
            Winf=zeros(2*nbus,2*nbus)
            Wsup=zeros(2*nbus,2*nbus)

            WV[k+k-1,k+k-1]=1
            WV[k+k,k+k]=1

            @constraint(mod,vecdot(WV,x5)<=1.1^2)
            @constraint(mod,vecdot(WV,x5)>=0.95^2)
        end

        status=solve(mod)

    x=getvalue(x5)
    LD=Dict()
    #eX
    for i=1:nbus
        if Bus.bshb[i]!=0
            LD[i]=getvalue(nb[i])
        else
            LD[i]="WITHOUT"
        end
    end

        if status==:Optimal
            A=getvalue(x5)
            for i=1:2*nbus
                for j=1:2*nbus
                    if A[i,j]<=0.0001
                        A[i,j]=0
                    else
                    end
                end
            end
            autovalores=eigvals(A)
            autovalormax=maximum(autovalores)
            Matrixrank=0
            indicador=0
            for i=1:length(autovalores)
                indicador=autovalores[i]/autovalormax
                if indicador>=0.0001
                    Matrixrank=Matrixrank+1
                else
                end
            end

        else
            Matrixrank=1000
            A="NOT"
            autovalores="NOT"
            autovalormax="NOT"
        end
        FO=getobjectivevalue(mod)

        return(x,LD,Matrixrank,A,autovalores,autovalormax,FO)
end

for l=1:ncomb
X[l] = flowbai3(l,nbus,Busnumber,Bus,Branch,Branchto,Branchfrom,out_lines,in_lines,Y,nstep)
println("------------------------------------------------ Ok ------------------------------------------------\n")
end

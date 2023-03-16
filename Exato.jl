using Gaston, JuMP, Gurobi, Combinatorics, LightGraphs, Random, LinearAlgebra, Printf,MathOptInterface, FileIO, JLD2

#--------------------------------------------------------------------------
function maxima_ordem(Y,ind,nmin,nmax)
    Y = Int.(round.(Y,digits=0))
    v = [sum(k*Y[k,i] for i in ind) for k in nmin:nmax]
    Kmax = maximum(v) - 1
    return Kmax
end

function Visitados(n,ind,Phi,dia)
    visitados=[]
    for i in ind
        ind = findall(Phi[i,dia] .> 0.5)
        if length(ind) > 0
            push!(visitados,i)
        end
    end
    return visitados
end

function Candidatos(c,Q,dia,n_max,visitados,ultima_visita,atracoes)
    if dia == 1
        ind=union(1:atracoes+1,R[1],H[1])
    else
        #
        cand = setdiff(1:atracoes+1,visitados)
        #
        if dia < maximum(T)
            ind = union(cand,R[dia],H[dia],ultima_visita)
        else
            ind = union(cand,R[dia],ultima_visita)
        end
    end
    return ind
end

function Plotasolucao(n,L,X,atracoes,R,H)
    scatter([L[2,1]], [L[2,2]],size="400,400",linewidth="3", pointtype="fsquare", plotstyle="linespoints",
    linecolor = "'white'",
    Axes(key = :off,xrange=(0,13), yrange=(0,15)))
    p=[]
    for i in 1:n, j in 1:n
        if X[i,j] == 1
            push!(p,plot!([L[i,1],L[j,1]],[L[i,2],L[j,2]],linewidth="0.5",linecolor = "'blue'",plotstyle = "linespoints", pointtype=""))
        end
    end
    plot!(L[2:atracoes+1,1], L[2:atracoes+1,2],linewidth="3", pointtype="fcircle",linecolor = "'cyan'")
    plot!(L[R[1],1], L[R[1],2],linewidth="3", pointtype="ftrianup",linecolor = "'green'")
    plot!(L[H[1],1], L[H[1],2],linewidth="3", pointtype="fdmd",linecolor = "'coral'")
    plot!([L[1,1]], [L[1,2]],linewidth="3", pointtype="fsquare",linecolor = "'red'")
    for i in 2:atracoes+1
        push!(p,plot!([L[i,1]].+0.3,[L[i,2]],supp=["$i"],w = "labels"))
    end
    for i in H[1]
        push!(p,plot!([L[i,1]].+0.3,[L[i,2]],supp=["H"],w = "labels"))
    end
    for i in R[1]
        push!(p,plot!([L[i,1]].+0.3,[L[i,2]],supp=["R"],w = "labels"))
    end
    push!(p,plot!([L[1,1]].+0.3,[L[1,2]],supp=["A"],w = "labels"))
    return p
end

function No_por_periodo_e_ordem(Visita,Ordem,n,T)
    No_por_periodo=[]
    for t in T
        push!(No_por_periodo,setdiff(findall(Visita[1:n,t] .> 0.5),1))
    end
    No_por_ordem=[]
    for k in 1:n
        push!(No_por_ordem,findall(Ordem[k,:] .> 0.5))
    end
    Ordem_por_periodo=[]
    for t in T
        push!(Ordem_por_periodo,length(No_por_periodo[t]))
    end
    return (No_por_periodo,No_por_ordem,Ordem_por_periodo)
end

function Dados(T,atracoes,Hoteis,Restaurantes)
    auxR = 1 .+ 9*rand(Restaurantes,2) #[1+4*rand() 1+4*rand();2+5*rand() 1+4*rand();5+9*rand() 5+9*rand() ]
    auxH =  1 .+ 9*rand(Hoteis,2) #[1+4*rand() 1+4*rand();2+5*rand() 1+4*rand();5+9*rand() 5+9*rand() ]
    AuxR=zeros(0,2);AuxH=zeros(0,2)
    auxTR=zeros(0);auxTH=zeros(0)
    auxPR=zeros(0);auxPH=zeros(0)
    for t in 1:maximum(T)
        AuxR=[AuxR;auxR]
        auxTR=[auxTR;2*ones(Restaurantes)]
        auxPR=[auxPR;[10;5]]
        if t<maximum(T)
            AuxH=[AuxH;auxH]
            auxTH=[auxTH;8*ones(Restaurantes)]
            auxPH=[auxPH;[5;10]]
        end
    end

    auxXA = 10*round.(randperm(1+atracoes)/(1+atracoes),digits=2)
    auxYA = 10*round.(randperm(1+atracoes)/(1+atracoes),digits=2)

    indx = sortperm(auxXA)
    auxXA=auxXA[indx]
    L = [
         [auxXA auxYA];
         AuxR;
         AuxH
    ]

    n = atracoes + 1 + maximum(T)*Restaurantes + (maximum(T)-1)*Hoteis

    c = round.([norm(L[i,:] - L[j,:]) for i in 1:n, j in 1:n],digits=2) #Custo das rotas do arco i para j

    aux3=zeros(atracoes)
    for i in 1:atracoes
        r=rand()
        if r<0.4
            aux3[i]=0.5
        elseif r>=0.4 && r<0.6
            aux3[i]=1
        elseif r>=0.6 && r<0.85
            aux3[i]=1.5
        elseif r>=0.85 && r<1
            aux3[i]=2
        end
    end

    d = round.([0; aux3;[auxTR;auxTH]],digits=2)

    #e = round.(c.*(0.3 .+ 0.7*rand(n,n)),digits=2)
    #for i in 1:n, j in 1:n
    #    if e[i,j] > 1.5
    #        e[i,j] = 0.5 + rand()
    #    end
    #end
    e = c/20 + 0.5*rand(n,n)
    e=round.(e,digits=2)

    P = [zeros(n) rand(0:10,n,atracoes) zeros(n,n-atracoes-1)] #Premio das prioridades (que leva em conta a ordem)
    Q = round.([0; 10*rand(atracoes);auxPR;auxPH],digits=2) #Premio das visitas (sem levar em conta a ordem)

    aux = 8*ones(atracoes)
    aux2= 19*ones(atracoes)
    #cedo=Int(ceil(0.05*atracoes))
    #tarde=Int(ceil(0.05*atracoes))
    #aux[1:cedo] .= 13
    #aux2[cedo+1:tarde+cedo] .= 12
    #perm=randperm(atracoes)
    #aux=aux[perm]
    #aux2=aux2[perm]

    auxAi=zeros(atracoes,0)
    auxBi=zeros(atracoes,0)
    for t in 1:maximum(T)
        auxAi=[auxAi aux .+ 24*(t-1)]
        auxBi=[auxBi aux2 .+ 24*(t-1)]
    end

    #############################
    cedo=Int(ceil(0.15*atracoes))
    tarde=Int(ceil(0.15*atracoes))

    for t in 1:maximum(T)
        pem=randperm(atracoes)
        ind_cedo=pem[1:cedo]
        ind_tarde=pem[cedo+1:cedo+tarde]
        auxBi[ind_cedo,t] .= auxBi[ind_cedo,t] .- 5
        auxAi[ind_tarde,t] .= auxAi[ind_tarde,t] .+ 5
    end
    ############################

    auxARi=1000*ones(Restaurantes*maximum(T),maximum(T))
    for t in 1:maximum(T)
        auxARi[Restaurantes*(t-1)+1:Restaurantes*t,t] .= 11 + 24*(t-1)
    end
    auxBRi=1000*ones(Restaurantes*maximum(T),maximum(T))
    for t in 1:maximum(T)
        auxBRi[Restaurantes*(t-1)+1:Restaurantes*t,t] .= 15 + 24*(t-1)
    end

    auxAHi=1000*ones(Hoteis*maximum(T)-Hoteis,maximum(T))
    for t in 1:maximum(T)-1
        auxAHi[Hoteis*(t-1)+1:Hoteis*t,t] .= 19 + 24*(t-1)
    end
    auxBHi=1000*ones(Hoteis*maximum(T)-Hoteis,maximum(T))
    for t in 1:maximum(T)-1
        auxBHi[Hoteis*(t-1)+1:Hoteis*t,t] .= 23 + 24*(t-1)
    end

    a=[
    zeros(1,maximum(T)) ;
    auxAi;
    auxARi;
    auxAHi
    ]

    b=[
    maximum(T)*24*ones(1,maximum(T)) ;
    auxBi;
    auxBRi;
    auxBHi
    ]

    R=[]
    for t in 1:maximum(T)
        Rt=atracoes+2 .+ (t-1)*Restaurantes:atracoes+1+Restaurantes .+ (t-1)*Restaurantes
        R = [R;[Rt]]
    end

    H=[]
    for t in 1:maximum(T)-1
        Ht=atracoes+2+maximum(T)*Restaurantes .+ (t-1)*Hoteis:atracoes+1+maximum(T)*Restaurantes+Hoteis .+ (t-1)*Hoteis
        H = [H;[Ht]]
    end

    M=maximum(T)*24
    w0 = 7

    n_min=Int.(3*ones(maximum(T),1))
    n_max=Int.(5*ones(maximum(T),1))

    return (n,L,a,b,c,d,e,P,Q,R,H,M,w0,n_min,n_max)
end

function ExatoFinal(TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit", TempoE)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    #
    @variable(modelo, x[i in 1:n,j in 1:n], Bin) # Se o arco (i,j) é percorrido
    @variable(modelo, y[k in 1:n,i in 1:n], Bin) # Se o nó i é designado para ser percorrido na ordem k
    @variable(modelo, Φ[i in 1:n, t in T], Bin)  # Se o nó i é percorrido no instante t
    @variable(modelo, w[i in 1:n] >=0)           # Inicio do tempo para percorrer o nó i
    @variable(modelo, z[i in 1:n, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
    #
    @objective(modelo, Max, sum( λ*(P[k,i] + Q[i])*y[k,i] for k in 1:n, i in 1:n) - (1-λ)*sum(x[i,j]*c[i,j] for i in 1:n, j in 1:n ))
    #
    @constraint(modelo,sum(x[i,1] for i in 1:n) == 1 )
    @constraint(modelo,sum(x[1,j] for j in 1:n) == 1 )
    @constraint(modelo,[i=1:n], x[i,i] == 0 )
    @constraint(modelo,[j=2:n],sum(x[i,j] for i in 1:n) == sum(y[k,j] for k in 1:n) )
    @constraint(modelo,[i=2:n],sum(x[i,j] for j in 1:n) == sum(y[k,i] for k in 1:n) )
    @constraint(modelo,[j=2:n-1], x[1,j] == y[1,j] )
    @constraint(modelo,[i=1:n,j in 1:n,k=2:n], x[i,j] >= y[k-1,i] + y[k,j] -1 )
    @constraint(modelo,[k=1:n-1],sum(y[k,i] for i in 1:n) >= sum(y[k+1,j] for j in 1:n) )
    @constraint(modelo,[i=1:n],sum(y[k,i] for k in 1:n) <= 1 )
    @constraint(modelo,[k=1:n],sum(y[k,i] for i in 1:n) <= 1 )
    @constraint(modelo,[j in 2:n], w[j] >= w0 + e[1,j]*x[1,j] - M*(1-x[1,j]) )
    @constraint(modelo,[i in 2:n, j in 2:n], w[j] >= w[i] + (e[i,j] + d[i])*x[i,j] - M*(1-x[i,j]) )
    #
    @constraint(modelo,[i in 2:n], sum(Φ[i,t] for t in T) == sum(y[k,i] for k in 1:n) )
    @constraint(modelo,[i in 1:n,t in T],z[i,t] <= M*Φ[i,t] )
    @constraint(modelo,[i in 1:n],w[i] == sum(z[i,t] for t in T) )
    @constraint(modelo,[i in 1:n,t in T], a[i,t]*Φ[i,t] <= z[i,t] )
    @constraint(modelo,[i in 1:n,t in T], b[i,t]*Φ[i,t] >= z[i,t] )
    @constraint(modelo,[t in T], sum(Φ[i,t] for i in R[t]) == 1 )
    @constraint(modelo,[t in setdiff(T,maximum(T))], sum(Φ[i,t] for i in H[t]) == 1 )
    @constraint(modelo,[t in T], sum(Φ[i,t] for i in 2:atracoes) <=  n_max[t])
    @constraint(modelo,[t in T], sum(Φ[i,t] for i in 2:atracoes) >=  n_min[t])
    #
    optimize!(modelo)
    Cputime = MOI.get(modelo, MOI.SolveTime())
    if has_values(modelo) == true
        Gap = MOI.get(modelo, MOI.RelativeGap())
        #
        Path = zeros(n,n)
        Ordem=zeros(n,n)
        Visita=zeros(n,maximum(T))
        Horas = zeros(n)
        for i in 1:n, j in 1:n
            Path[i,j] = value(x[i,j])
        end
        for k in 1:n, i in 1:n
            Ordem[k,i] = value(y[k,i])
        end
        for i in 1:n, t in T
            Visita[i,t] = value(Φ[i,t])
        end
        for i in 1:n
            Horas[i] = value(w[i])
        end

        Lower_bound = objective_value(modelo)
        Upper_bound = dual_objective_value(modelo)
    else
        Path=-1
        Ordem=-1
        Visita=-1
        Horas=-1
        Lower_bound=-1
        Upper_bound=-1
        Gap=-1
    end
    #
    return (Path,Ordem,Visita,Horas,Lower_bound,Upper_bound,Cputime,Gap)
end

function ExatoFinal_compromisso(TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,z10,z20,z11,z21)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit", TempoE)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    #
    @variable(modelo, x[i in 1:n,j in 1:n], Bin) # Se o arco (i,j) é percorrido
    @variable(modelo, y[k in 1:n,i in 1:n], Bin) # Se o nó i é designado para ser percorrido na ordem k
    @variable(modelo, Φ[i in 1:n, t in T], Bin)  # Se o nó i é percorrido no instante t
    @variable(modelo, w[i in 1:n] >=0)           # Inicio do tempo para percorrer o nó i
    @variable(modelo, z[i in 1:n, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
    @variable(modelo, u >=0)
    #
    @objective(modelo, Min, u)
    #
    @constraint(modelo, λ*(sum(z11) - sum((P[k,i] + Q[i])*y[k,i] for k in 1:n, i in 1:n) )/(sum(z11) - sum(z10))  <= u)
    @constraint(modelo, (1-λ)*(sum(x[i,j]*c[i,j] for i in 1:n, j in 1:n) - sum(z20))/(sum(z21) - sum(z20))  <= u)
    #
    @constraint(modelo,sum(x[i,1] for i in 1:n) == 1 )
    @constraint(modelo,sum(x[1,j] for j in 1:n) == 1 )
    @constraint(modelo,[i=1:n], x[i,i] == 0 )
    @constraint(modelo,[j=2:n],sum(x[i,j] for i in 1:n) == sum(y[k,j] for k in 1:n) )
    @constraint(modelo,[i=2:n],sum(x[i,j] for j in 1:n) == sum(y[k,i] for k in 1:n) )
    @constraint(modelo,[j=2:n-1], x[1,j] == y[1,j] )
    @constraint(modelo,[i=1:n,j in 1:n,k=2:n], x[i,j] >= y[k-1,i] + y[k,j] -1 )
    @constraint(modelo,[k=1:n-1],sum(y[k,i] for i in 1:n) >= sum(y[k+1,j] for j in 1:n) )
    @constraint(modelo,[i=1:n],sum(y[k,i] for k in 1:n) <= 1 )
    @constraint(modelo,[k=1:n],sum(y[k,i] for i in 1:n) <= 1 )
    @constraint(modelo,[j in 2:n], w[j] >= w0 + e[1,j]*x[1,j] - M*(1-x[1,j]) )
    @constraint(modelo,[i in 2:n, j in 2:n], w[j] >= w[i] + (e[i,j] + d[i])*x[i,j] - M*(1-x[i,j]) )
    #
    @constraint(modelo,[i in 2:n], sum(Φ[i,t] for t in T) == sum(y[k,i] for k in 1:n) )
    @constraint(modelo,[i in 1:n,t in T],z[i,t] <= M*Φ[i,t] )
    @constraint(modelo,[i in 1:n],w[i] == sum(z[i,t] for t in T) )
    @constraint(modelo,[i in 1:n,t in T], a[i,t]*Φ[i,t] <= z[i,t] )
    @constraint(modelo,[i in 1:n,t in T], b[i,t]*Φ[i,t] >= z[i,t] )
    @constraint(modelo,[t in T], sum(Φ[i,t] for i in R[t]) == 1 )
    @constraint(modelo,[t in setdiff(T,maximum(T))], sum(Φ[i,t] for i in H[t]) == 1 )
    @constraint(modelo,[t in T], sum(Φ[i,t] for i in 2:atracoes) <=  n_max[t])
    @constraint(modelo,[t in T], sum(Φ[i,t] for i in 2:atracoes) >=  n_min[t])
    #
    optimize!(modelo)
    Cputime = MOI.get(modelo, MOI.SolveTime())
    if has_values(modelo) == true
        Gap = MOI.get(modelo, MOI.RelativeGap())
        #
        Path = zeros(n,n)
        Ordem=zeros(n,n)
        Visita=zeros(n,maximum(T))
        Horas = zeros(n)
        for i in 1:n, j in 1:n
            Path[i,j] = value(x[i,j])
        end
        for k in 1:n, i in 1:n
            Ordem[k,i] = value(y[k,i])
        end
        for i in 1:n, t in T
            Visita[i,t] = value(Φ[i,t])
        end
        for i in 1:n
            Horas[i] = value(w[i])
        end

        Lower_bound = objective_value(modelo)
        Upper_bound = dual_objective_value(modelo)
    else
        Path=-1
        Ordem=-1
        Visita=-1
        Horas=-1
        Lower_bound=-1
        Upper_bound=-1
        Gap=-1
    end

    Lower_bound = objective_value(modelo)
    Upper_bound = dual_objective_value(modelo)
    #
    return (Path,Ordem,Visita,Horas,Lower_bound,Upper_bound,Cputime,Gap)
end

function RangeK(tmin,T,Ordem_por_periodo)
    tmax = maximum(T) - tmin + 1
    T_escolhidos = setdiff(T,[tmin,tmax])
    Ktodos=[]
    #
    for t in T_escolhidos
        if t == 1
            ki = 1:sum(Ordem_por_periodo[1])
        else
            ki = sum(Ordem_por_periodo[1:t-1]) + 1:sum(Ordem_por_periodo[1:t])
        end
        Ktodos=[Ktodos;collect(ki)]
    end
    return (Ktodos)
end

function RangeK_Tres(tmin,T,Ordem_por_periodo)
    tmax = tmin+1
    T_escolhidos = setdiff(T,[tmin,tmax])
    Ktodos=[]
    #
    for t in T_escolhidos
        if t == 1
            ki = 1:sum(Ordem_por_periodo[1])
        else
            ki = sum(Ordem_por_periodo[1:t-1]) + 1:sum(Ordem_por_periodo[1:t])
        end
        Ktodos=[Ktodos;collect(ki)]
    end
    return (Ktodos)
end

function Valores_obj_Exato(Path,Ordem,P,Q,c,n)
    z1 = sum(Ordem[k,i]*(P[k,i] + Q[i]) for i in 1:n, k in 1:n)
    z2 = sum(Path[i,j]*c[i,j] for i in 1:n, j in 1:n)
    return (z1,z2)
end

#--------------------------------------------------------------------------
Random.seed!(0)

T=[1,2,3]
atracoes = 15
Restaurantes=2
Hoteis=2

TempoE=3600

(n,L,a,b,c,d,e,P,Q,R,H,M,w0,n_min,n_max) = Dados(T,atracoes,Restaurantes,Hoteis)

z1C=zeros(3)
z2C=zeros(3)
GAPC=zeros(3)
CPUC=zeros(3)

#--------------------------------------------------------------------------
λ = 0.0 # =1 maximizo o benefício, 0 minimizo custo e 0.5 equilibra os dois
#
(Path,Ordem,Visita,W,LB,UB,Cputime0E,Gap0E) = ExatoFinal(2*TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H)
(z10E,z20E) = Valores_obj_Exato(Path,Ordem,P,Q,c,n)

#Pat = Int.((round.(Path,digits=0)))
#p=Plotasolucao(n,L,Pat,atracoes,R,H)
#p[end]


#--------------------------------------------------------------------------
λ = 1.0 # =1 maximizo o benefício, 0 minimizo custo e 0.5 equilibra os dois
#
(Path,Ordem,Visita,W,LB,UB,Cputime1E,Gap1E) = ExatoFinal(TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H)
(z11E,z21E) = Valores_obj_Exato(Path,Ordem,P,Q,c,n)

#Pat = Int.((round.(Path,digits=0)))
#p=Plotasolucao(n,L,Pat,atracoes,R,H)
#p[end]

#--------------------------------------------------------------------------
global k=0
for λ ∈ [0.2,0.5,0.8]
    k=k+1
    (Path,Ordem,Visita,W,LB,UB,CputimeE,GapE) = ExatoFinal_compromisso(TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,z10E,z20E,z11E,z21E)
    (z1,z2) = Valores_obj_Exato(Path,Ordem,P,Q,c,n)
    #
    z1C[k]=sum(z1)
    z2C[k]=sum(z2)
    CPUC[k]=CputimeE
    GAPC[k]=GapE
end

Z1 = [z10E;z1C;z11E]
Z2 = [z20E;z2C;z21E]
CPU = [Cputime0E;CPUC;Cputime1E]
GAP = [Gap0E;GAPC;Gap1E]

FileIO.save("Intancia4_Z1.jld2","RH",Z1)
FileIO.save("Intancia4_Z2.jld2","RH",Z2)
FileIO.save("Intancia4_CPU.jld2","RH",CPU)
FileIO.save("Intancia4_GAP.jld2","RH",GAP)

#W1 = [z10E;z1C;z11E]
#W2 = [z20E;z2C;z21E]

#plot(W1,W2,Axes(key = :off),linewidth="2", pointtype="fsquare", plotstyle="linespoints")

using Gaston, JuMP, Gurobi, Combinatorics, LightGraphs, Random, LinearAlgebra, Printf,MathOptInterface, FileIO, JLD2, Statistics

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
        auxPR=[auxPR;[10;5;2;1;0.5]]
        if t<maximum(T)
            AuxH=[AuxH;auxH]
            auxTH=[auxTH;8*ones(Restaurantes)]
            auxPH=[auxPH;[0.5;1;2;5;10]]
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
    cedo=Int(ceil(0.1*atracoes))
    tarde=Int(ceil(0.1*atracoes))

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
    n_max=Int.(6*ones(maximum(T),1))

    return (n,L,a,b,c,d,e,P,Q,R,H,M,w0,n_min,n_max)
end

function Primeira_heuristica(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H)
    Path=zeros(n,n)
    Ordem=zeros(n,n)
    Visita=zeros(n,maximum(T))
    Horarios=zeros(n)
    UV = zeros(maximum(T))
    #
    CPUtime=0
    for dia=1:maximum(T)
        if dia == 1
            ind = Candidatos(c,Q,dia,n_max,[],1,atracoes)
            modelo = Model(Gurobi.Optimizer)
            set_optimizer_attribute(modelo, "TimeLimit", TempoH)
            set_optimizer_attribute(modelo, "MIPGap", 0.01)
            #set_optimizer_attribute(modelo, "Heuristics", 0.2)
            #
            @variable(modelo, x[i in ind,j in ind], Bin) # Se o arco (i,j) é percorrido
            @variable(modelo, y[k in 1:n_max[1]+3,i in ind], Bin) # Se o nó i é designado para ser percorrido na ordem k
            @variable(modelo, Φ[i in ind, t in T], Bin)  # Se o nó i é percorrido no instante t
            @variable(modelo, w[i in ind] >=0)           # Inicio do tempo para percorrer o nó i
            @variable(modelo, z[i in ind, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
            #
            @objective(modelo, Max, λ*sum((P[k,i] + Q[i])*y[k,i] for k in 1:n_max[1]+3, i in setdiff(ind,1)) - (1-λ)*sum(x[i,j]*c[i,j] for i in ind, j in setdiff(ind,1) ))
            #
            @constraint(modelo,sum(x[i,1] for i in ind) == 1 )
            @constraint(modelo,sum(x[1,j] for j in ind) == 1 )
            @constraint(modelo,[i=ind], x[i,i] == 0 )
            @constraint(modelo,[j=ind],sum(x[i,j] for i in ind) == sum(y[k,j] for k in 1:n_max[1]+3) )
            @constraint(modelo,[i=ind],sum(x[i,j] for j in ind) == sum(y[k,i] for k in 1:n_max[1]+3) )
            @constraint(modelo,[j=setdiff(ind,1)], x[1,j] == y[1,j] )
            @constraint(modelo,[i=ind,j in ind,k=2:n_max[1]+3], x[i,j] >= y[k-1,i] + y[k,j] -1 )
            @constraint(modelo,[k=1:n_max[1]+2],sum(y[k,i] for i in ind) >= sum(y[k+1,j] for j in ind) )
            @constraint(modelo,[i=ind],sum(y[k,i] for k in 1:n_max[1]+3) <= 1 )
            @constraint(modelo,[k=1:n_max[1]+3],sum(y[k,i] for i in ind) <= 1 )
            @constraint(modelo,[j in setdiff(ind,1)], w[j] >= w0 + e[1,j]*x[1,j] - M*(1-x[1,j]) )
            @constraint(modelo,[i in setdiff(ind,1), j in setdiff(ind,1)], w[j] >= w[i] + (e[i,j] + d[i])*x[i,j] - M*(1-x[i,j]) )
            #
            @constraint(modelo,[i in ind],sum(Φ[i,t] for t in 1:dia) == sum(y[k,i] for k in 1:n_max[1]+3) )
            @constraint(modelo,[i in ind,t in 1:dia],z[i,t] <= M*Φ[i,t] )
            @constraint(modelo,[i in ind],w[i] == sum(z[i,t] for t in 1:dia) )
            @constraint(modelo,[i in ind,t in 1:dia], a[i,t]*Φ[i,t] <= z[i,t] )
            @constraint(modelo,[i in ind,t in 1:dia], b[i,t]*Φ[i,t] >= z[i,t] )
            @constraint(modelo, sum(y[k,i] for k in 1:n_max[1]+3 for i in R[1]) == 1 )
            @constraint(modelo, sum(y[k,i] for k in 1:+n_max[1]+3 for i in H[1]) == 1 )
            #
            @constraint(modelo,[t in 1:dia], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) <=  n_max[t])
            @constraint(modelo,[t in 1:dia], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) >=  n_min[t])
            #
            @constraint(modelo,[k in 1:n_min[1]+2], sum(y[k,i] for i in ind) == 1 )
            #
            optimize!(modelo)
            X = value.(x)
            Y = value.(y)
            Phi = value.(Φ)
            W = value.(w)
            Z = value.(z)
            CPUtime = CPUtime + MOI.get(modelo, MOI.SolveTime())

            for i in ind, j in setdiff(ind,1)
                Path[i,j] = X[i,j]
            end
            for k in 1:n_max[1]+2, i in setdiff(ind,1)
                Ordem[k,i] = Y[k,i]
            end
            for i in setdiff(ind), t in [dia]
                Visita[i,t] = Phi[i,t]
            end
            for i in ind
                Horarios[i] = W[i]
            end
            #
            global visitados, ultima_visita, Kmax, W0
            Kmax = maxima_ordem(Y,ind,1,n_max[1]+3)
            visitados = Visitados(n,ind,Phi,dia)
            ultima_visita = Int(sum(i*value(x[i,1]) for i in ind))
            UV[dia] = ultima_visita
            W0 = W[ultima_visita]

        else
            if dia < maximum(T)
                ind = Candidatos(c,Q,dia,n_max,visitados,ultima_visita,atracoes)
            else
                ind = union(1,Candidatos(c,Q,dia,n_max,visitados,ultima_visita,atracoes))
            end

            modelo = Model(Gurobi.Optimizer)
            set_optimizer_attribute(modelo, "TimeLimit", TempoH)
            set_optimizer_attribute(modelo, "MIPGap", 0.01)
            set_optimizer_attribute(modelo, "Heuristics", 0.2)
            #
            @variable(modelo, x[i in ind,j in ind], Bin) # Se o arco (i,j) é percorrido
            @variable(modelo, y[k in Kmax+1:Kmax+n_max[dia]+3,i in ind], Bin) # Se o nó i é designado para ser percorrido na ordem k
            @variable(modelo, Φ[i in ind, t in T], Bin)  # Se o nó i é percorrido no instante t
            @variable(modelo, w[i in ind] >=0)           # Inicio do tempo para percorrer o nó i
            @variable(modelo, z[i in ind, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
            #
            if dia < maximum(T)
                @objective(modelo, Max, λ*sum((P[k,i] + Q[i])*y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3, i in ind) -
                (1-λ)*sum(x[i,j]*c[i,j] for i in ind, j in setdiff(ind,ultima_visita)) )
                @constraint(modelo,sum(x[i,ultima_visita] for i in ind) == 1 )
                @constraint(modelo,sum(x[ultima_visita,j] for j in ind) == 1 )
            else
                @objective(modelo, Max, λ*sum( (P[k,i] + Q[i])*y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3, i in ind)
                - (1-λ)*sum(x[i,j]*c[i,j] for i in setdiff(ind,1), j in setdiff(ind,ultima_visita) ))
                @constraint(modelo,sum(x[i,ultima_visita] for i in ind) == 1 )
                @constraint(modelo,x[1,ultima_visita] == 1 )
                @constraint(modelo,sum(x[ultima_visita,j] for j in ind) == 1 )
            end
            #
            @constraint(modelo,[j=ind],sum(x[i,j] for i in ind) == sum(y[k,j] for k in Kmax+1:Kmax+n_max[dia]+3) )
            @constraint(modelo,[i=ind],sum(x[i,j] for j in ind) == sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3) )
            @constraint(modelo,[j=setdiff(ind,ultima_visita)], x[ultima_visita,j] == y[Kmax+1,j] )
            @constraint(modelo,[i=ind,j in ind,k=Kmax+2:Kmax+n_max[dia]+3], x[i,j] >= y[k-1,i] + y[k,j] -1 )
            @constraint(modelo,[k=Kmax+1:Kmax+n_max[dia]+2],sum(y[k,i] for i in ind) >= sum(y[k+1,j] for j in ind) )
            @constraint(modelo,[i=ind],sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3) <= 1 )
            @constraint(modelo,[k=Kmax+1:Kmax+n_max[dia]+3],sum(y[k,i] for i in ind) <= 1 )
            @constraint(modelo,[j in setdiff(ind,ultima_visita)], w[j] >= W0 + (d[ultima_visita] +  e[ultima_visita,j])*x[ultima_visita,j] - M*(1-x[ultima_visita,j]) )
            @constraint(modelo,[i in setdiff(ind,ultima_visita), j in setdiff(ind,ultima_visita)], w[j] >= w[i] + (e[i,j] + d[i])*x[i,j] - M*(1-x[i,j]) )
            #
            @constraint(modelo,[i in ind],sum(Φ[i,t] for t in [dia]) == sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3) )
            @constraint(modelo,[i in setdiff(ind,ultima_visita),t in [dia]],w[i] <= M*Φ[i,t] )
            @constraint(modelo,[i in ind,t in [dia]], a[i,t]*Φ[i,t] <= w[i] )
            @constraint(modelo,[i in ind,t in [dia]], b[i,t]*Φ[i,t] >= w[i] )
            @constraint(modelo,[t in [dia]], sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3 for i in R[t]) == 1 )
            if dia < maximum(T)
               @constraint(modelo,[t in [dia]], sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3 for i in H[t])  == 1 )
            end
            if dia == maximum(T)
                @constraint(modelo,[i in ind], y[Kmax+n_max[dia]+3,i] == 0 )
            end
            @constraint(modelo,[t in [dia]], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) <=  n_max[t])
            @constraint(modelo,[t in [dia]], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) >=  n_min[t])
            #
            optimize!(modelo)
            CPUtime = CPUtime + MOI.get(modelo, MOI.SolveTime())
            X = value.(x)
            Y = value.(y)
            Phi = value.(Φ)
            W = value.(w)
            Z = value.(z)

            if dia < maximum(T)
                for i in ind, j in setdiff(ind,union(1,ultima_visita))
                    Path[i,j] = X[i,j]
                end
                for k in Kmax+1:Kmax+n_max[dia]+2, i in setdiff(ind,union(1,ultima_visita))
                    Ordem[k,i] = Y[k,i]
                end
                for i in setdiff(ind,union(1,ultima_visita)), t in [dia]
                    Visita[i,t] = Phi[i,t]
                end
                for i in setdiff(ind,ultima_visita)
                    Horarios[i] = W[i]
                end
            else
                for i in ind, j in setdiff(ind,ultima_visita)
                    Path[i,j] = X[i,j]
                end
                for k in Kmax+1:Kmax+n_max[dia]+2, i in setdiff(ind,union(1,ultima_visita))
                    Ordem[k,i] = Y[k,i]
                end
                for i in setdiff(ind,ultima_visita), t in [dia]
                    Visita[i,t] = Phi[i,t]
                end
                for i in setdiff(ind,ultima_visita)
                    Horarios[i] = W[i]
                end
            end
            #
            global visitados, ultima_visita, Kmax, W0
            Kmax = maxima_ordem(Y,ind,Kmax+1,Kmax+n_max[dia]+3)
            visitados = union(Visitados(n,ind,Phi,dia),visitados)
            ultima_visita = Int(sum(i*value(x[i,ultima_visita]) for i in ind))
            UV[dia] = ultima_visita
            W0 = W[ultima_visita]
        end
    end
    UV = Int.(UV)
    return (Path,Ordem,Visita,UV,CPUtime,Horarios)
end

function TerceiraHeuristica(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Path,Ordem,Visita,No_por_periodo,No_por_ordem,Ordem_por_periodo,tmin)
    #
    if tmin < maximum(T)
        tmax = tmin + 1
    else
        tmin = tmin-1
        tmax = tmin + 1
    end
    Ktodos = RangeK_Tres(tmin,T,Ordem_por_periodo) #Ordem das variáveis que vou fixar
    #
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit", TempoH)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "Heuristics", 0.2)
    #
    @variable(modelo, x[i in 1:n,j in 1:n], Bin) # Se o arco (i,j) é percorrido
    @variable(modelo, y[k in 1:n,i in 1:n], Bin) # Se o nó i é designado para ser percorrido na ordem k
    @variable(modelo, Φ[i in 1:n, t in T], Bin)  # Se o nó i é percorrido no instante t
    @variable(modelo, w[i in 1:n] >=0)           # Inicio do tempo para percorrer o nó i
    @variable(modelo, z[i in 1:n, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
    #
    @objective(modelo, Max, sum( λ*(P[k,i] + Q[i])*y[k,i] for k in 1:n, i in 1:n) - (1-λ)*sum(x[i,j]*c[i,j] for i in 1:n, j in 1:n ))
    #

    for i in 1:n, j in 1:n
        if Int(round(Path[i,j],digits=0)) > 0.5
            set_start_value(x[i,j], 1)
        end
    end
    for k in 1:n, i in 1:n
        if Int(round(Ordem[k,i],digits=0)) > 0.5
            set_start_value(y[k,i], 1)
        end
    end
    for i in 2:n, t in T
        if Int(round(Visita[i,t],digits=0)) > 0.5
            set_start_value(Φ[i,t], 1)
        end
    end

    for t in setdiff(T,[tmin,tmax]), i in No_por_periodo[t]
        if Int.(round.(Visita[i,t],digits=0)) > 0.5
            @constraint(modelo, Φ[i,t] == 1 )
        end
    end
    for k in Ktodos, i in No_por_ordem[k][1]
        if Int.(round.(Ordem[k,i],digits=0)) > 0.5
            @constraint(modelo, y[k,i] == 1 )
        end
    end

    #
    #for i in 1:n, j in 1:n
    #    @constraint(modelo,x[i,j] == Path[i,j] )
    #end
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
    #
    Path = zeros(n,n)
    Ordem=zeros(n,n)
    Visita=zeros(n,maximum(T))
    Horas=zeros(n)
    #
    Cputime = MOI.get(modelo, MOI.SolveTime())
    Gap = MOI.get(modelo, MOI.RelativeGap())
    #
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
    #
    Lower_bound = objective_value(modelo)
    Upper_bound = dual_objective_value(modelo)
    return (Path,Ordem,Visita,Horas,Lower_bound,Upper_bound,Cputime,Gap)
end

function Primeira_heuristica_compromisso(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,z10,z20,z11,z21)
    Path=zeros(n,n)
    Ordem=zeros(n,n)
    Visita=zeros(n,maximum(T))
    Horarios=zeros(n)
    UV = zeros(maximum(T))
    #
    CPUtime=0
    for dia=1:maximum(T)
        if dia == 1
            ind = Candidatos(c,Q,dia,n_max,[],1,atracoes)
            modelo = Model(Gurobi.Optimizer)
            set_optimizer_attribute(modelo, "TimeLimit", TempoH)
            set_optimizer_attribute(modelo, "MIPGap", 0.01)
            #set_optimizer_attribute(modelo, "Heuristics", 0.2)
            #
            @variable(modelo, x[i in ind,j in ind], Bin) # Se o arco (i,j) é percorrido
            @variable(modelo, y[k in 1:n_max[1]+3,i in ind], Bin) # Se o nó i é designado para ser percorrido na ordem k
            @variable(modelo, Φ[i in ind, t in T], Bin)  # Se o nó i é percorrido no instante t
            @variable(modelo, w[i in ind] >=0)           # Inicio do tempo para percorrer o nó i
            @variable(modelo, z[i in ind, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
            @variable(modelo, u >=0)
            #
            #@objective(modelo, Max, sum( λ*(P[k,i] + Q[i])*y[k,i] for k in 1:n_max[1]+3, i in setdiff(ind,1))
            #- (1-λ)*sum(x[i,j]*c[i,j] for i in ind, j in setdiff(ind,1) ))
            @objective(modelo, Min, u)
            #
            @constraint(modelo, λ*((z11[dia] - sum((P[k,i] + Q[i])*y[k,i] for k in 1:n_max[1]+3, i in setdiff(ind,1)))/(z11[dia] - z10[dia])) <= u )
            @constraint(modelo, (1-λ)*((sum(x[i,j]*c[i,j] for i in ind, j in setdiff(ind,1)) - z20[dia])/(z21[dia] - z20[dia])) <= u )
            #
            @constraint(modelo,sum(x[i,1] for i in ind) == 1 )
            @constraint(modelo,sum(x[1,j] for j in ind) == 1 )
            @constraint(modelo,[i=ind], x[i,i] == 0 )
            @constraint(modelo,[j=ind],sum(x[i,j] for i in ind) == sum(y[k,j] for k in 1:n_max[1]+3) )
            @constraint(modelo,[i=ind],sum(x[i,j] for j in ind) == sum(y[k,i] for k in 1:n_max[1]+3) )
            @constraint(modelo,[j=setdiff(ind,1)], x[1,j] == y[1,j] )
            @constraint(modelo,[i=ind,j in ind,k=2:n_max[1]+3], x[i,j] >= y[k-1,i] + y[k,j] -1 )
            @constraint(modelo,[k=1:n_max[1]+2],sum(y[k,i] for i in ind) >= sum(y[k+1,j] for j in ind) )
            @constraint(modelo,[i=ind],sum(y[k,i] for k in 1:n_max[1]+3) <= 1 )
            @constraint(modelo,[k=1:n_max[1]+3],sum(y[k,i] for i in ind) <= 1 )
            @constraint(modelo,[j in setdiff(ind,1)], w[j] >= w0 + e[1,j]*x[1,j] - M*(1-x[1,j]) )
            @constraint(modelo,[i in setdiff(ind,1), j in setdiff(ind,1)], w[j] >= w[i] + (e[i,j] + d[i])*x[i,j] - M*(1-x[i,j]) )
            #
            @constraint(modelo,[i in ind],sum(Φ[i,t] for t in 1:dia) == sum(y[k,i] for k in 1:n_max[1]+3) )
            @constraint(modelo,[i in ind,t in 1:dia],z[i,t] <= M*Φ[i,t] )
            @constraint(modelo,[i in ind],w[i] == sum(z[i,t] for t in 1:dia) )
            @constraint(modelo,[i in ind,t in 1:dia], a[i,t]*Φ[i,t] <= z[i,t] )
            @constraint(modelo,[i in ind,t in 1:dia], b[i,t]*Φ[i,t] >= z[i,t] )
            @constraint(modelo, sum(y[k,i] for k in 1:n_max[1]+3 for i in R[1]) == 1 )
            @constraint(modelo, sum(y[k,i] for k in 1:+n_max[1]+3 for i in H[1]) == 1 )
            #
            @constraint(modelo,[t in 1:dia], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) <=  n_max[t])
            @constraint(modelo,[t in 1:dia], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) >=  n_min[t])
            #
            @constraint(modelo,[k in 1:n_min[1]+2], sum(y[k,i] for i in ind) == 1 )
            #
            optimize!(modelo)
            X = value.(x)
            Y = value.(y)
            Phi = value.(Φ)
            W = value.(w)
            Z = value.(z)
            #
            CPUtime = CPUtime + MOI.get(modelo, MOI.SolveTime())

            for i in ind, j in setdiff(ind,1)
                Path[i,j] = X[i,j]
            end
            for k in 1:n_max[1]+2, i in setdiff(ind,1)
                Ordem[k,i] = Y[k,i]
            end
            for i in setdiff(ind), t in [dia]
                Visita[i,t] = Phi[i,t]
            end
            for i in ind
                Horarios[i] = W[i]
            end
            #
            global visitados, ultima_visita, Kmax, W0
            Kmax = maxima_ordem(Y,ind,1,n_max[1]+3)
            visitados = Visitados(n,ind,Phi,dia)
            ultima_visita = Int(sum(i*value(x[i,1]) for i in ind))
            UV[dia] = ultima_visita
            W0 = W[ultima_visita]
        else
            if dia < maximum(T)
                ind = Candidatos(c,Q,dia,n_max,visitados,ultima_visita,atracoes)
            else
                ind = union(1,Candidatos(c,Q,dia,n_max,visitados,ultima_visita,atracoes))
            end

            modelo = Model(Gurobi.Optimizer)
            set_optimizer_attribute(modelo, "TimeLimit",  TempoH)
            set_optimizer_attribute(modelo, "MIPGap", 0.01)
            set_optimizer_attribute(modelo, "Heuristics", 0.2)
            #
            @variable(modelo, x[i in ind,j in ind], Bin) # Se o arco (i,j) é percorrido
            @variable(modelo, y[k in Kmax+1:Kmax+n_max[dia]+3,i in ind], Bin) # Se o nó i é designado para ser percorrido na ordem k
            @variable(modelo, Φ[i in ind, t in T], Bin)  # Se o nó i é percorrido no instante t
            @variable(modelo, w[i in ind] >=0)           # Inicio do tempo para percorrer o nó i
            @variable(modelo, z[i in ind, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
            @variable(modelo, u >=0)

            @objective(modelo, Min, u)

            #
            if dia < maximum(T)
                @constraint(modelo, λ*((z11[dia] - sum((P[k,i] + Q[i])*y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3, i in ind))/(z11[dia] - z10[dia])) <= u )
                @constraint(modelo, (1-λ)*((sum(x[i,j]*c[i,j] for i in ind, j in setdiff(ind,ultima_visita)) - z20[dia])/(z21[dia] - z20[dia])) <= u )

                @constraint(modelo,sum(x[i,ultima_visita] for i in ind) == 1 )
                @constraint(modelo,sum(x[ultima_visita,j] for j in ind) == 1 )
            else
                @constraint(modelo, λ*((z11[dia] - sum((P[k,i] + Q[i])*y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3, i in ind))/(z11[dia] - z10[dia])) <= u )
                @constraint(modelo, (1-λ)*((sum(x[i,j]*c[i,j] for i in setdiff(ind,1), j in setdiff(ind,ultima_visita)) - z20[dia])/(z21[dia] - z20[dia])) <= u )

                @constraint(modelo,sum(x[i,ultima_visita] for i in ind) == 1 )
                @constraint(modelo,x[1,ultima_visita] == 1 )
                @constraint(modelo,sum(x[ultima_visita,j] for j in ind) == 1 )
            end
            #
            @constraint(modelo,[j=ind],sum(x[i,j] for i in ind) == sum(y[k,j] for k in Kmax+1:Kmax+n_max[dia]+3) )
            @constraint(modelo,[i=ind],sum(x[i,j] for j in ind) == sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3) )
            @constraint(modelo,[j=setdiff(ind,ultima_visita)], x[ultima_visita,j] == y[Kmax+1,j] )
            @constraint(modelo,[i=ind,j in ind,k=Kmax+2:Kmax+n_max[dia]+3], x[i,j] >= y[k-1,i] + y[k,j] -1 )
            @constraint(modelo,[k=Kmax+1:Kmax+n_max[dia]+2],sum(y[k,i] for i in ind) >= sum(y[k+1,j] for j in ind) )
            @constraint(modelo,[i=ind],sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3) <= 1 )
            @constraint(modelo,[k=Kmax+1:Kmax+n_max[dia]+3],sum(y[k,i] for i in ind) <= 1 )
            @constraint(modelo,[j in setdiff(ind,ultima_visita)], w[j] >= W0 + e[ultima_visita,j]*x[ultima_visita,j] - M*(1-x[ultima_visita,j]) )
            @constraint(modelo,[i in setdiff(ind,ultima_visita), j in setdiff(ind,ultima_visita)], w[j] >= w[i] + (e[i,j] + d[i])*x[i,j] - M*(1-x[i,j]) )
            #
            @constraint(modelo,[i in ind],sum(Φ[i,t] for t in [dia]) == sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3) )
            @constraint(modelo,[i in setdiff(ind,ultima_visita),t in [dia]],w[i] <= M*Φ[i,t] )
            @constraint(modelo,[i in ind,t in [dia]], a[i,t]*Φ[i,t] <= w[i] )
            @constraint(modelo,[i in ind,t in [dia]], b[i,t]*Φ[i,t] >= w[i] )
            @constraint(modelo,[t in [dia]], sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3 for i in R[t]) == 1 )
            if dia < maximum(T)
               @constraint(modelo,[t in [dia]], sum(y[k,i] for k in Kmax+1:Kmax+n_max[dia]+3 for i in H[t])  == 1 )
            end
            if dia == maximum(T)
                @constraint(modelo,[i in ind], y[Kmax+n_max[dia]+3,i] == 0 )
            end
            @constraint(modelo,[t in [dia]], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) <=  n_max[t])
            @constraint(modelo,[t in [dia]], sum(Φ[i,t] for i in intersect(2:atracoes,ind)) >=  n_min[t])
            #
            optimize!(modelo)
            CPUtime = CPUtime + MOI.get(modelo, MOI.SolveTime())
            X = value.(x)
            Y = value.(y)
            Phi = value.(Φ)
            W = value.(w)
            Z = value.(z)

            if dia < maximum(T)
                for i in ind, j in setdiff(ind,union(1,ultima_visita))
                    Path[i,j] = X[i,j]
                end
                for k in Kmax+1:Kmax+n_max[dia]+2, i in setdiff(ind,union(1,ultima_visita))
                    Ordem[k,i] = Y[k,i]
                end
                for i in setdiff(ind,union(1,ultima_visita)), t in [dia]
                    Visita[i,t] = Phi[i,t]
                end
                for i in setdiff(ind,ultima_visita)
                    Horarios[i] = W[i]
                end
            else
                for i in ind, j in setdiff(ind,ultima_visita)
                    Path[i,j] = X[i,j]
                end
                for k in Kmax+1:Kmax+n_max[dia]+2, i in setdiff(ind,union(1,ultima_visita))
                    Ordem[k,i] = Y[k,i]
                end
                for i in setdiff(ind,ultima_visita), t in [dia]
                    Visita[i,t] = Phi[i,t]
                end
                for i in setdiff(ind,ultima_visita)
                    Horarios[i] = W[i]
                end
            end
            #
            global visitados, ultima_visita, Kmax, W0
            Kmax = maxima_ordem(Y,ind,Kmax+1,Kmax+n_max[dia]+3)
            visitados = union(Visitados(n,ind,Phi,dia),visitados)
            ultima_visita = Int(round(sum(i*value(x[i,ultima_visita]) for i in ind),digits=0))
            UV[dia] = ultima_visita
            W0 = W[ultima_visita]
        end
    end
    UV = Int.(UV)
    return (Path,Ordem,Visita,UV,CPUtime)
end

function TerceiraHeuristica_compromisso(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Path,Ordem,Visita,No_por_periodo,No_por_ordem,Ordem_por_periodo,tmin,z10,z20,z11,z21)
    #
    if tmin < maximum(T)
        tmax = tmin + 1
    else
        tmin = tmin-1
        tmax = tmin + 1
    end
    Ktodos = RangeK_Tres(tmin,T,Ordem_por_periodo) #Ordem das variáveis que vou fixar
    #
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit", TempoH)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "Heuristics", 0.2)
    #
    @variable(modelo, x[i in 1:n,j in 1:n], Bin) # Se o arco (i,j) é percorrido
    @variable(modelo, y[k in 1:n,i in 1:n], Bin) # Se o nó i é designado para ser percorrido na ordem k
    @variable(modelo, Φ[i in 1:n, t in T], Bin)  # Se o nó i é percorrido no instante t
    @variable(modelo, w[i in 1:n] >=0)           # Inicio do tempo para percorrer o nó i
    @variable(modelo, z[i in 1:n, t in T] >=0)   # Instante de inicio para visitar o nó i no período t
    @variable(modelo, u >=0)
    #
    #@objective(modelo, Max, sum( λ*(P[k,i] + Q[i])*y[k,i] for k in 1:n, i in 1:n) - (1-λ)*sum(x[i,j]*c[i,j] for i in 1:n, j in 1:n ))
    #
    @objective(modelo, Min, u)
    #
    @constraint(modelo, λ*(sum(z11) - sum((P[k,i] + Q[i])*y[k,i] for k in 1:n, i in 1:n) )/(sum(z11) - sum(z10))  <= u)
    @constraint(modelo, (1-λ)*(sum(x[i,j]*c[i,j] for i in 1:n, j in 1:n) - sum(z20))/(sum(z21) - sum(z20))  <= u)

    for i in 1:n, j in 1:n
        if Int(round(Path[i,j],digits=0)) > 0.5
            set_start_value(x[i,j], 1)
        end
    end
    for k in 1:n, i in 1:n
        if Int(round(Ordem[k,i],digits=0)) > 0.5
            set_start_value(y[k,i], 1)
        end
    end
    for i in 2:n, t in T
        if Int(round(Visita[i,t],digits=0)) > 0.5
            set_start_value(Φ[i,t], 1)
        end
    end
    #
    for t in setdiff(T,[tmin,tmax]), i in No_por_periodo[t]
        if Int.(round.(Visita[i,t],digits=0)) > 0.5
            @constraint(modelo, Φ[i,t] == 1 )
        end
    end
    for k in Ktodos, i in No_por_ordem[k][1]
        if Int.(round.(Ordem[k,i],digits=0)) > 0.5
            @constraint(modelo, y[k,i] == 1 )
        end
    end
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
    #
    Path = zeros(n,n)
    Ordem=zeros(n,n)
    Visita=zeros(n,maximum(T))
    Horas=zeros(n)
    #
    Cputime = MOI.get(modelo, MOI.SolveTime())
    Gap = MOI.get(modelo, MOI.RelativeGap())
    #
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
    #
    return (Path,Ordem,Visita,Horas,Cputime,Gap)
end

function ExatoFinal(TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Path,Ordem,Visita)
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
    for i in 1:n, j in 1:n
        if Int(round(Path[i,j],digits=0)) > 0.5
            set_start_value(x[i,j], 1)
        end
    end
    for k in 1:n, i in 1:n
        if Int(round(Ordem[k,i],digits=0)) > 0.5
            set_start_value(y[k,i], 1)
        end
    end
    for i in 2:n, t in T
        if Int(round(Visita[i,t],digits=0)) > 0.5
            set_start_value(Φ[i,t], 1)
        end
    end
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
    #
    return (Path,Ordem,Visita,Horas,Lower_bound,Upper_bound,Cputime,Gap)
end

function ExatoFinal_compromisso(TempoE,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Path,Ordem,Visita,z10,z20,z11,z21)
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
    for i in 1:n, j in 1:n
        if Int(round(Path[i,j],digits=0)) > 0.5
            set_start_value(x[i,j], 1)
        end
    end
    for k in 1:n, i in 1:n
        if Int(round(Ordem[k,i],digits=0)) > 0.5
            set_start_value(y[k,i], 1)
        end
    end
    for i in 2:n, t in T
        if Int(round(Visita[i,t],digits=0)) > 0.5
            set_start_value(Φ[i,t], 1)
        end
    end
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
    Gap = MOI.get(modelo, MOI.RelativeGap())

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

function Valores_obj(T,No_por_periodo,Ordem_por_periodo,UV,Path,Ordem,P,Q,c)
    z1=zeros(maximum(T))
    z2=zeros(maximum(T))
    for t in T
        if t==1
            indt = union(No_por_periodo[t],1)
            z2[t] = sum(Path[i,j]*c[i,j] for i in indt, j in indt)
            z1[t] = sum(Ordem[k,i]*(P[k,i] + Q[i]) for i in indt, k in 1:Ordem_por_periodo[1])
        elseif t > 1 && t < maximum(T)
            indt = union(No_por_periodo[t],UV[t-1])
            z2[t] = sum(Path[i,j]*c[i,j] for i in indt, j in indt)
            z1[t] = sum(Ordem[k,i]*(P[k,i] + Q[i]) for i in indt, k in sum(Ordem_por_periodo[1:t-1])+1:sum(Ordem_por_periodo[1:t]))
        elseif t == maximum(T)
            indt = union(No_por_periodo[t],UV[t-1],1)
            z2[t] = sum(Path[i,j]*c[i,j] for i in indt, j in indt)
            z1[t] = sum(Ordem[k,i]*(P[k,i] + Q[i]) for i in indt, k in sum(Ordem_por_periodo[1:t-1])+1:sum(Ordem_por_periodo[1:t]))
        end
    end
    return (z1,z2)
end


#--------------------------------------------------------------------------
Random.seed!(0)

T=[1,2,3,4,5,6,7,8,9,10,11,12]
atracoes =100
Hoteis = 5
Restaurantes = 5

TempoH=120

GAP=zeros(5)
CPU1=zeros(5)
CPU2=zeros(5)

(n,L,a,b,c,d,e,P,Q,R,H,M,w0,n_min,n_max) = Dados(T,atracoes,Hoteis,Restaurantes)
Z10=zeros(maximum(T))
Z20=zeros(maximum(T))
Z11=zeros(maximum(T))
Z21=zeros(maximum(T))
#
for λ in [0,1]
    (Path,Ordem,Visita,UV,Cputime1,Horas) = Primeira_heuristica(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H)
    (No_por_periodo,No_por_ordem,Ordem_por_periodo) = No_por_periodo_e_ordem(Visita,Ordem,n,T)
    (z1,z2) = Valores_obj(T,No_por_periodo,Ordem_por_periodo,UV,Path,Ordem,P,Q,c)
    #
    g=[]
    tem=[]
    global No_por_periodo,No_por_ordem,Ordem_por_periodo,z10,z20, CPU, GAP
    for tmin in 1:2:maximum(T)
        (Path,Ordem,Visita,Horas,LB,UB,Cputime2,Gap) = TerceiraHeuristica(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Path,Ordem,Visita,No_por_periodo,No_por_ordem,Ordem_por_periodo,tmin)
        (No_por_periodo,No_por_ordem,Ordem_por_periodo) = No_por_periodo_e_ordem(Visita,Ordem,n,T)
        (z1,z2) = Valores_obj(T,No_por_periodo,Ordem_por_periodo,UV,Path,Ordem,P,Q,c)
        push!(g,Gap)
        push!(tem,Cputime2)
    end
    #
    if λ < 0.5
        GAP[1]=mean(g)
        CPU1[1]=Cputime1
        CPU2[1]=sum(tem)
        Z10 .= z1
        Z20 .= z2
        FileIO.save("50_Path0.jld2","RH",Path)
        FileIO.save("50_Ordem0.jld2","RH",Ordem)
        FileIO.save("50_Visita0.jld2","RH",Visita)
        FileIO.save("50_Horas0.jld2","RH",Horas)
    else
        GAP[2]=mean(g)
        CPU1[2]=Cputime1
        CPU2[2]=sum(tem)
        Z11 .= z1
        Z21 .= z2
        FileIO.save("50_Path1.jld2","RH",Path)
        FileIO.save("50_Ordem1.jld2","RH",Ordem)
        FileIO.save("50_Visita1.jld2","RH",Visita)
        FileIO.save("50_Horas1.jld2","RH",Horas)
    end
end

#
k=2
compr=0
for λ ∈ [0.2,0.5,0.8]
    k=k+1
    compr=Int(λ*10)
    (Path,Ordem,Visita,UV,Cputime1) = Primeira_heuristica_compromisso(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Z10,Z20,Z11,Z21)
    (No_por_periodo,No_por_ordem,Ordem_por_periodo) = No_por_periodo_e_ordem(Visita,Ordem,n,T)
    (z1C,z2C) = Valores_obj(T,No_por_periodo,Ordem_por_periodo,UV,Path,Ordem,P,Q,c)
    #
    g=[]
    tem=[]
    global No_por_periodoC,No_por_ordemC,Ordem_por_periodoC,z1C,z2C
    for tmin in 1:2:maximum(T)
        (Path,Ordem,Visita,Horas,Cputime2,Gap) = TerceiraHeuristica_compromisso(TempoH,λ,a,b,c,d,e,P,Q,n_min,n_max,w0,atracoes,R,H,Path,Ordem,Visita,No_por_periodo,No_por_ordem,Ordem_por_periodo,tmin,Z10,Z20,Z11,Z21)
        (No_por_periodo,No_por_ordem,Ordem_por_periodo) = No_por_periodo_e_ordem(Visita,Ordem,n,T)
        (z1C,z2C) = Valores_obj(T,No_por_periodo,Ordem_por_periodo,UV,Path,Ordem,P,Q,c)
        push!(g,Gap)
        push!(tem,Cputime2)
    end
    #
    GAP[k]=mean(g)
    CPU1[k]=Cputime1
    CPU2[k]=sum(tem)
    FileIO.save("50_Path$compr.jld2","RH",Path)
    FileIO.save("50_Ordem$compr.jld2","RH",Ordem)
    FileIO.save("50_Visita$compr.jld2","RH",Visita)
    FileIO.save("50_Horas$compr.jld2","RH",Horas)
end

Pat = Int.((round.(Path,digits=0)))
p=Plotasolucao(n,L,Pat,atracoes,R,H)
p[end]

FileIO.save("50_CPU1.jld2","RH",CPU1)
FileIO.save("50_CPU2.jld2","RH",CPU2)
FileIO.save("50_GAP.jld2","RH",GAP)


#-------------------------------------------------------
Path0 = FileIO.load("50_Path0.jld2","RH")
Path1 = FileIO.load("50_Path1.jld2","RH")
Path2 = FileIO.load("50_Path2.jld2","RH")
Path5 = FileIO.load("50_Path5.jld2","RH")
Path8 = FileIO.load("50_Path8.jld2","RH")

Ordem0 = FileIO.load("50_Ordem0.jld2","RH")
Ordem1 = FileIO.load("50_Ordem1.jld2","RH")
Ordem2 = FileIO.load("50_Ordem2.jld2","RH")
Ordem5 = FileIO.load("50_Ordem5.jld2","RH")
Ordem8 = FileIO.load("50_Ordem8.jld2","RH")

Visita0 = FileIO.load("50_Visita0.jld2","RH")
Visita1 = FileIO.load("50_Visita1.jld2","RH")
Visita2 = FileIO.load("50_Visita2.jld2","RH")
Visita5 = FileIO.load("50_Visita5.jld2","RH")
Visita8 = FileIO.load("50_Visita8.jld2","RH")

Horas0 = FileIO.load("50_Horas0.jld2","RH")
Horas1 = FileIO.load("50_Horas1.jld2","RH")
Horas2 = FileIO.load("50_Horas2.jld2","RH")
Horas5 = FileIO.load("50_Horas5.jld2","RH")
Horas8 = FileIO.load("50_Horas8.jld2","RH")

#----------------------------------------------------------
z10 = sum(Ordem0[k,i]*(P[k,i] + Q[i]) for i in 1:n, k in 1:n)
z20 = sum(Path0[i,j]*c[i,j] for i in 1:n, j in 1:n)

z11 = sum(Ordem1[k,i]*(P[k,i] + Q[i]) for i in 1:n, k in 1:n)
z21 = sum(Path1[i,j]*c[i,j] for i in 1:n, j in 1:n)

z12 = sum(Ordem2[k,i]*(P[k,i] + Q[i]) for i in 1:n, k in 1:n)
z22 = sum(Path2[i,j]*c[i,j] for i in 1:n, j in 1:n)

z15 = sum(Ordem5[k,i]*(P[k,i] + Q[i]) for i in 1:n, k in 1:n)
z25 = sum(Path5[i,j]*c[i,j] for i in 1:n, j in 1:n)

z18 = sum(Ordem8[k,i]*(P[k,i] + Q[i]) for i in 1:n, k in 1:n)
z28 = sum(Path8[i,j]*c[i,j] for i in 1:n, j in 1:n)

W1 = [z10;z12;z15;z18;z11]
W2 = [z20;z22;z25;z28;z21]


plot(W1,W2,Axes(key = :off),linewidth="2", pointtype="fsquare", plotstyle="linespoints")


@printf("---------------------------------------------\n")
caminho=[]
for k=1:n
    if sum(i*Ordem[k,i] for i in 1:n) > 0
        push!(caminho,sum(i*Ordem[k,i] for i in 1:n))
    end
end
caminho = Int.([1;caminho;1])
@printf("Caminho realizado: \n")
for i in 1:length(caminho)
    if i < length(caminho)
        @printf("%1d,",caminho[i])
    else
        @printf("%1d\n",caminho[i])
    end
end
@printf("---------------------------------------------\n")

@printf("---------------------------------------------\n")
for i in 1:n
    if Horas[i] > 1e-6
        @printf("O nó i=%1d foi visitado no tempo w[i]=%2.2f\n",i,Horas[i])
    end
end
@printf("---------------------------------------------\n")

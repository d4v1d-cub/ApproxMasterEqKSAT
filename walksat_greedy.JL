using Random
using SimpleHypergraphs
using Base.Filesystem

mutable struct Spin
    greedy::UInt16
    estado::Bool
end

function Nodos(nodos::Dict{Int,Bool},N::Int)::Int
    n=rand(1:N,1)[1]
    while haskey(nodos,n)
        n=rand(1:N,1)[1]
    end
    return n
end

function CrearHGraph(N::Int,M::Int,K::Int,s::Int)::Hypergraph{Bool,Spin,UInt8}
    Random.seed!(s)
    h=Hypergraph{Bool,Spin,UInt8}(N,M)
    for i in 1:M
        valor=rand(0:1,K)
        for j in 1:K
            h[Nodos(h.he2v[i],N),i]=valor[j]
        end
    end
    return h
end


function read_graph_from_file()::Hypergraph{Bool,Spin,UInt8}
    N, M = [parse(Int, x) for x in split(readline())]
    println("Reading graph with N=$(N), M=$(M)")
    flush(stdout)
    h=Hypergraph{Bool,Spin,UInt8}(N,M)
    for he in 1:M
        v = [parse(Int, x) for x in split(readline())]
        for n in v
            h[abs(n), he] = signbit(n)
        end
    end
    return h
end


function AsignaValores!(h::Hypergraph{Bool,Spin,UInt8},s::Int)::Dict{Int,Bool}
    Random.seed!(s)
    for k in keys(h.v_meta)
        h.v_meta[k]=Spin(0,(rand()<0.5))
    end
    CUnsats=Dict{Int,Bool}()
    for k in keys(h.he_meta)
        h.he_meta[k]=0
        for n in keys(h.he2v[k])
            h.he_meta[k]+= !xor(h.v_meta[n].estado,h.he2v[k][n])
        end
        if h.he_meta[k]==0
            CUnsats[k]=1
        else
            for n in keys(h.he2v[k])
                h.v_meta[n].greedy+=1
            end
        end
    end
    return CUnsats
end

function Invierte!(h::Hypergraph{Bool,Spin,UInt8},v::Int,CUnsats::Dict{Int,Bool})::Dict{Int,Bool}
    h.v_meta[v].estado=!h.v_meta[v].estado
    for i in keys(h.v2he[v])
        s1=1-2*xor(h.v_meta[v].estado,h.he2v[i][v])
        s2=((h.he_meta[i]==0)||((h.he_meta[i]+s1)==0))*s1
        h.he_meta[i]+=s1
        if s2!=0
            if s2==1
                delete!(CUnsats,i)
            else
                CUnsats[i]=true
            end
            for j in keys(h.he2v[i])
                h.v_meta[j].greedy+=s2
            end
        end
    end
    return CUnsats
end

function Greedy(h::Hypergraph{Bool,Spin,UInt8},C::Int)::Int
    greedy=length(h.he_meta)
    for i in keys(h.he2v[C])
        if h.v_meta[i].greedy<greedy
            greedy=h.v_meta[i].greedy
        end
    end
    return rand(filter( i -> h.v_meta[i].greedy==greedy,keys(h.he2v[C])),1)[1]
end

function Paso!(h::Hypergraph{Bool,Spin,UInt8},CUnsats::Dict{Int,Bool},q::Float64,s::Int)::Dict{Int,Bool}
    Random.seed!(s)
    C=rand(keys(CUnsats),1)[1]
    if rand()<q
        v=rand(keys(h.he2v[C]),1)[1]
        return Invierte!(h,v,CUnsats)
    end
    return Invierte!(h,Greedy(h,C),CUnsats)
end

function Walk!(h::Hypergraph{Bool,Spin,UInt8},alpha_u::Vector{Float64},q::Float64,s::Int)
    CUnsats=AsignaValores!(h,s)
    alpha_u[1]+=length(CUnsats)/nhv(h)
    for t in 2:length(alpha_u)
        if length(CUnsats)==0
            return 0
        end
        CUnsats=Paso!(h,CUnsats,q,t)
        alpha_u[t]+=length(CUnsats)/nhv(h)
    end
    return length(CUnsats)
end

function Imprime(fname_dyn,alphas,N, div)
    f=open(fname_dyn,"w")
    t=0
    while t <= min(N, length(alphas))
        println(f,Float64(t)/N,"\t",alphas[t + 1])
        t=t+div
    end
    while t <= length(alphas)
        println(f,Float64(t)/N,"\t",alphas[t + 1])
        t = max(Int(round(round(t * 1.05))), t + div)
    end
    close(f)
end


function Imprime_sat(fname_Psol,sat, cumul_ener, cumul_ener_sqr, nsamples,time)
    av_e = cumul_ener / nsamples
    av_e_sqr = cumul_ener_sqr / nsamples
    err = sqrt((av_e_sqr - av_e * av_e) / nsamples)
    f=open(fname_Psol,"w")
    println(f, "# t   nsamples   Psol   E   std(E)")
    println(f,time,"\t",nsamples,"\t",sat/nsamples, "\t", av_e, "\t", err)
    close(f)
end


function Datos(K::Int, N::Int, time::Float64, Exp::Int, Hist::Int, M::Int, q::Float64, 
               save_every::Int64, div::Int64, fname_dyn::String, fname_Psol::String)
    println("Begining runs")
    alpha_u=zeros(Float64,round(Int, N*time+1))
    counter = 0
    sat = 0
    cumul_ener = 0
    cumul_ener_sqr = 0
    for e in 1:Exp
        println("instance No. ", e)
        flush(stdout)
        h=CrearHGraph(N,M,K,e)
        solved=false
        hist = 0
        while hist < Hist && !solved
            unsat_final = Walk!(h,alpha_u,q,hist)
            cumul_ener += unsat_final
            cumul_ener_sqr += unsat_final * unsat_final
            if unsat_final == 0
                solved = true
                sat += 1
            end
            counter += 1
            hist += 1
        end
        if e % save_every == 0
            Imprime(fname_dyn,alpha_u/counter,N, div)
            Imprime_sat(fname_Psol, sat, cumul_ener, cumul_ener_sqr, e, time)
        end
    end
    Imprime(fname_dyn,alpha_u/counter,N, div)
    Imprime_sat(fname_Psol, sat, cumul_ener, cumul_ener_sqr, Exp, time)
end


function SingleGraph(h::Hypergraph{Bool,Spin,UInt8}, time::Float64, Hist::Int, q::Float64, 
                     div::Int64, fname_dyn::String, fname_Psol::String, save_every::Int64)
    println("Begining runs")
    N = nhv(h)
    alpha_u=zeros(Float64,round(Int, N*time+1))
    sat = 0
    cumul_ener = 0
    cumul_ener_sqr = 0
    hist = 0
    for hist in 1:Hist
        unsat_final = Walk!(h,alpha_u,q,hist)
        cumul_ener += unsat_final
        cumul_ener_sqr += unsat_final * unsat_final
        if unsat_final == 0
            sat += 1
        end
        if hist % save_every == 0
            Imprime(fname_dyn,alpha_u/hist,N, div)
            Imprime_sat(fname_Psol, sat, cumul_ener, cumul_ener_sqr, hist, time)
        end
    end
    Imprime(fname_dyn,alpha_u/Hist,N, div)
    Imprime_sat(fname_Psol, sat, cumul_ener, cumul_ener_sqr, Hist, time)
end

if ARGS[1] == "file"
    q = parse(Float64, ARGS[2])
    Hist=parse(Int64, ARGS[3])
    time=parse(Float64, ARGS[4])
    ruta=ARGS[5]
    save_every=parse(Int64, ARGS[6])

    h=read_graph_from_file()
    N = nhv(h)
    M = nhe(h)
    K = length(h.v2he[1])

    dtmin = 0.01
    div = round(Int64, max(dtmin * N, 1))

    fname_dyn = ruta*"WalkSAT_dynamics_K_$(K)_N_$(N)_alpha_$(M/N)_q_$(q)_t_$(time)_dtmin_$(dtmin)_nhist_$(Hist).txt"
    fname_Psol = ruta*"WalkSAT_Psol_K_$(K)_N_$(N)_alpha_$(M/N)_q_$(q)_t_$(time)_dtmin_$(dtmin)_nhist_$(Hist).txt"

    SingleGraph(h, time, Hist, q, div, fname_dyn, fname_Psol, save_every)

elseif ARGS[1] == "inside"

    N=parse(Int64, ARGS[2])
    M=parse(Int64, ARGS[3])
    q=parse(Float64, ARGS[4])
    K=parse(Int64, ARGS[5])
    Exp=parse(Int64, ARGS[6])
    Hist=parse(Int64, ARGS[7])
    time=parse(Float64, ARGS[8])
    # ruta="Walk-SAT/1_Algoritmo/Results/"
    ruta=ARGS[9]

    save_every=parse(Int64, ARGS[10])
    dtmin = 0.01
    div = round(Int64, max(dtmin * N, 1))

    fname_dyn = ruta*"WalkSAT_dynamics_K_$(K)_N_$(N)_alpha_$(M/N)_q_$(q)_t_$(time)_dtmin_$(dtmin)_ngraphs_$(Exp)_nhist_$(Hist).txt"
    fname_Psol = ruta*"WalkSAT_Psol_K_$(K)_N_$(N)_alpha_$(M/N)_q_$(q)_t_$(time)_dtmin_$(dtmin)_ngraphs_$(Exp)_nhist_$(Hist).txt"

    Datos(K,N,time,Exp,Hist,M,q, save_every, div, fname_dyn, fname_Psol)


else
    println("Error: first argument must be 'file' or 'inside'")
    exit(1)
end


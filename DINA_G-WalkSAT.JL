using Random
using DifferentialEquations

function Comb(size::Int,K::Int)::Array{Float64,2}
    C=Array{Float64,2}(undef,size+1,size+1)
    T=2^K-1
    for u in 0:size
        C[(u)+1,(0)+1]=(1/T)^u
        for k in 0:(size-1)
            C[(u)+1,(k+1)+1]=C[(u)+1,(k)+1]*(1+u/(k+1))*(1-1/T)
        end
    end
    return C
end

function Op(P::Array{Float64,2},size::Int)::Dict{String,Float64}
    D=Dict{String,Float64}()
    D["u"]=0
    D["s"]=0
    for s in 0:size
        for u in 0:size
            D["u"]+=u*P[(s)+1,(u)+1]
            D["s"]+=s*P[(s)+1,(u)+1]
        end
    end
    return D
end

function Op!(D::Dict{String,Float64},P::Array{Float64,2},size::Int,greedy::Array{Float64,1},q::Float64)::Dict{String,Float64}
    D["uu"]=0
    D["su"]=0
    for s in 0:size
        for u in 0:size
            D["uu"]+=u*u*P[(s)+1,(u)+1]*greedy[(s)+1]
            D["su"]+=s*u*P[(s)+1,(u)+1]*greedy[(s)+1]
        end
    end
    return D
end

function Pt(s::Int,u::Int,P::Array{Float64,2},size::Int)::Float64
    if s<0||u<0||s>size||u>size
        return 0.0
    end
    return P[(s)+1,(u)+1]
end

function TablaGreedy(P::Array{Float64,2},size::Int,q::Float64,D::Dict{String,Float64},K::Int)::Array{Float64,1}
    S2=0.0
    greedy=Array{Float64,1}()
    for s in 0:size
        S1=0.0
        for u in 0:size
            S1+=u*P[(s)+1,(u)+1]
        end
        S1/=(2*D["u"])
        push!(greedy,3*(1-S1-S2)^2+S1^2)
        S2+=(2*S1)
    end
    for i in 1:length(greedy)
        greedy[i]=q+(1-q)*greedy[i]
    end

    return greedy
end

function Derivadas!(dP::Array{Float64,2},P::Array{Float64,2},p::Tuple{Int,Int,Int,Int,Float64,Array{Float64,2}},t::Float64)
    K=p[1]
    M=p[2]
    N=p[3]
    size=p[4]
    q=p[5]
    C=p[6]
    D=Op(P,size)
    greedy=TablaGreedy(P,size,q,D,K)
    Op!(D,P,size,greedy,q)
    if D["u"]<=0
        dP=zeros(size+1,size+1)
    else
        for s in 0:size
            for u in 0:size
                U1=0.0
                for k in 0:s
                    if u+k<=size
                        U1+=C[(u)+1,(k)+1]*(s-k)*P[(u+k)+1,(s-k)+1]/D["u"]*greedy[(u+k)+1]
                    end
                end
                U2=D["su"]/D["s"]/D["u"]*(-s*P[(s)+1,(u)+1]+(s+1)*Pt(s+1,u-1,P,size))*(K-1)/(2^K-1)
                U3=D["uu"]/D["u"]/D["u"]*(-u*P[(s)+1,(u)+1]+(u+1)*Pt(s-1,u+1,P,size))*(K-1)
                dP[(s)+1,(u)+1]=(-u*P[(s)+1,(u)+1]/D["u"]*greedy[(s)+1]+U1+U2+U3)
            end
        end
    end
end

function P0(M::Int,N::Int,K::Int,size::Int)::Array{Float64,2}
    P=Array{Float64,2}(undef, size+1, size+1)
    F=Array{Float64,1}()
    push!(F,1)
    for c in 1:size
        push!(F,K*M/N/2^K/c*F[c])
    end
    for s in 0:size
        for u in 0:size
            P[(s)+1,(u)+1]=exp(-K*M/N)*(2^K-1)^s*F[(u)+1]*F[(s)+1]
        end
    end
    return P
end

function Resuelve(Ms::Vector{Int},q::Float64,N::Int,K::Int,t::Float64,tol::Float64,size::Int)::Dict{Float64,Dict{Float64,ODESolution}}
    D=Dict{Float64,Dict{Float64,ODESolution}}()
    for M in Ms
        D[M/N]=Dict{Float64,ODESolution}()
        prob=ODEProblem(Derivadas!,P0(M,N,K,size),(0.0,t),(K,M,N,size,q,Comb(size,K)))
        sol=solve(prob,Tsit5(),abstol=tol,reltol=tol)
        D[M/N][q]=sol
    end
    return D
end

function Alpha(sol::ODESolution,size::Int,K::Int)::Vector{Float64}
    alphas=Vector{Float64}()
    noCero=true
    for i in keys(sol.t)
        push!(alphas,0)
        if noCero
            for s in 0:size
                for u in 0:size
                    alphas[end]+=u*sol[i][(s)+1,(u)+1]
                end
            end
            if alphas[end]<=0
                alphas[end]=0
                noCero=false
            end
        end
    end
    return alphas ./ K
end


function Imprime(D::Dict{Float64,Dict{Float64,ODESolution}},N::Int,K::Int,size::Int,name::String)
    ruta="D/"
    for alpha in keys(D)
        for q in keys(D[alpha])
            f=open(ruta*name*"_N-$(N) _alpha-$(alpha) _q-$(q) _K-$(K) _size-$(size) .txt","w")
            sol=Alpha(D[alpha][q],size,K)
            for i in 1:length(sol)
                println(f,D[alpha][q].t[i],"\t",sol[i])
            end
            close(f)
        end
    end
end

Q=parse(Int,ARGS[1])
N=10000
K=3
size=25
tol=5e-8
t=10.0
name="RE_"
#Ms=collect(1500:500:1500)


if Q==1
    q=0.0
    Ms=[38800]
    D=Resuelve(Ms,q,N,K,t,tol,size)
    Imprime(D,N,K,size,name)
elseif Q==2
    q=0.0
    Ms=[39100]
    D=Resuelve(Ms,q,N,K,t,tol,size)
    Imprime(D,N,K,size,name)
end

struct QPweights
    Q::AbstractArray
    R::AbstractArray
    function QPweights(Q::Vector, Qf::Vector, R::Vector, H::Int)
        Qf = diagm(Qf)
        n = length(Q)
        Q = diagm(repeat(Q, inner=[1,1], outer=[1,H]) |> vec)
        R = diagm(repeat(R, inner=[1,1], outer=[1,H]) |> vec)
        Q[end-n+1:end,end-n+1:end] .= Qf    
        return new(Q,R)
    end 
end 

struct PredictionMatrices
    Abig::AbstractArray
    Bbig::AbstractArray
    function PredictionMatrices(Â::AbstractArray, B̂::AbstractArray, H::Int)
        n = size(Â,1)
        m = size(B̂,2)        
        Abig=zeros(H*n,n);
        Bbig=zeros(H*n,H*m);
        for i=1:H
            Abig[(1+(i-1)*n):i*n,1:n]=Â^i;
            for j=1:H            
                if(i>=j)                
                    Bbig[(1+(i-1)*n):i*n,(1+(j-1)*m):j*m] = Â^(i-j)*B̂;                
                end
            end
        end
        return new(Abig,Bbig)
    end 
end

mutable struct LTImodel
    A::AbstractArray
    B::AbstractArray
    n::Int
    m::Int 
    x0::AbstractArray
    u0::AbstractArray
    function LTImodel(A::AbstractArray,B::AbstractArray)
        n = size(A,1)
        m = size(B,2)
        return new(A,B,n,m,zeros(n),zeros(m))
    end 
end 

struct Constraints 
    umax::AbstractArray
    umin::AbstractArray
    CΔ::AbstractArray
    pMat::PredictionMatrices
    function Constraints(umax::AbstractArray,umin::AbstractArray,H::Int)
        m = length(umin)
        umax  = repeat(umax,H)
        umin  = repeat(umin,H)
        CΔ    = kron(ones(H,H)|>tril, I(m))
        return new(umax,umin,CΔ)
    end 
end
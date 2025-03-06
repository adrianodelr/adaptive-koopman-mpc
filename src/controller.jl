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

function build_predmat(Â::AbstractArray, B̂::AbstractArray, H::Int)
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
    return Abig,Bbig
end 

# struct PredictionMatrices
#     Abig::AbstractArray
#     Bbig::AbstractArray
#     function PredictionMatrices(Â::AbstractArray, B̂::AbstractArray, H::Int)
#         n = size(Â,1)
#         m = size(B̂,2)        
#         Abig=zeros(H*n,n);
#         Bbig=zeros(H*n,H*m);
#         for i=1:H
#             Abig[(1+(i-1)*n):i*n,1:n]=Â^i;
#             for j=1:H            
#                 if(i>=j)                
#                     Bbig[(1+(i-1)*n):i*n,(1+(j-1)*m):j*m] = Â^(i-j)*B̂;                
#                 end
#             end
#         end
#         return new(Abig,Bbig)
#     end 
# end

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
    function Constraints(umax::AbstractArray,umin::AbstractArray,H::Int)
        m = length(umin)
        umax  = repeat(umax,H)
        umin  = repeat(umin,H)
        CΔ    = kron(ones(H,H)|>tril, I(m))
        return new(umax,umin,CΔ)
    end 
end

function build_QP(ẑ0::Vector, rbig::Vector, w::QPweights, Abig::AbstractArray, Bbig::AbstractArray)
    P = Bbig'*w.Q*Bbig + w.R
    q = 2Bbig'*w.Q*(Abig*ẑ0 - rbig)
    return P,q    
end 

function get_dims(model::EDMDModel)
    return model.Ψ.p, model.buffer.m 
end 

function augment_model(model::EDMDModel)
    p,m = get_dims(model)
    Â = [model.A model.B; zeros(m,p) I(m)]
    B̂ = [model.B; I(m)]
    return Â,B̂
end 

mutable struct adaptiveKMPC
    model::EDMDModel
    Abig::AbstractArray
    Bbig::AbstractArray
    weights::QPweights
    constr::Constraints 
    rbig::Vector
    H::Int
    n̂::Int 
    solver 
    function adaptiveKMPC(model::EDMDModel,Q::Vector,Qf::Vector,R::Vector,r::Vector,H::Int,constr::Constraints)
        p,m = get_dims(model)
        Ψ = model.dictionary.Ψ
        n̂ = p+m
        # from eq. (7) 
        Â,B̂ = augment_model(model)
        weights = QPweights([Q;zeros(m)], [Qf;zeros(m)], R, H)
        Abig, Bbig = build_predmat(Â, B̂, H)               

        rbig = repeat([lifting(r,Ψ);zeros(m)],H)
        P,q = build_QP(zeros(n̂), rbig, weights, Abig, Bbig)
        
        ul,uu,CΔ = build_constraints(zeros(p), constr, model, H)

        solver = OSQP.Model()
        OSQP.setup!(solver, P=sparse(P), q=vec(q), A=SparseMatrixCSC(CΔ),l=ul, u=uu, verbose=false, warm_start=true)  
        return new(model, Abig, Bbig, weights, H, n̂, QPw, rbig, constr, solver)
    end 
end

function get_control(x0::Vector, ctrl::adaptiveKMPC)
        
    ẑ0 = lifting(x0, ctrl.dictionary.Ψ)
    lifting_and_regression!(ctrl.model)

    Â,B̂ = augment_model(model)
    
    Abig, Bbig = build_predmat(Â, B̂, H)

    P,q = build_QP([ẑ0;ctrl.model.u0], rbig, ctrl.weights, Abig, Bbig)

    ul,uu,CΔ = build_constraints(ẑ0, ctrl.constr, ctrl.model, ctrl.H)

    OSQP.update!(ctrl.solver, Px=triu(sparse(P)).nzval, q=vec(q), l=ul, u=uu, Ax=sparse(CΔ).nzval)

    Δu = OSQP.solve!(ctrl.solver).x
    u = ctrl.model.u0 + Δu[1:ctrl.model.m]

    ctrl.model.u0 = u 
    ctrl.model.x0 = x
    return u 
end  
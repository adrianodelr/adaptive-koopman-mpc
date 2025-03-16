struct Weights
    Q::AbstractArray
    R::AbstractArray
    function Weights(Q::Vector, Qf::Vector, R::Vector, H::Int)
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
    A_bold=zeros(H*n,n);
    B_bold=zeros(H*n,H*m);
    for i=1:H
        A_bold[(1+(i-1)*n):i*n,1:n]=Â^i;
        for j=1:H            
            if(i>=j)                
                B_bold[(1+(i-1)*n):i*n,(1+(j-1)*m):j*m] = Â^(i-j)*B̂;                
            end
        end
    end
    return A_bold,B_bold
end 

struct Constraints 
    ul_bold_const::AbstractArray
    uu_bold_const::AbstractArray
    CΔ::AbstractArray
    function Constraints(ul::AbstractArray,uu::AbstractArray,H::Int)
        m = length(umin)
        ul_bold_const  = kron(ones(H), ul)              # constant part of lower bound on u 
        uu_bold_const  = kron(ones(H), uu)              # constant part of upper bound on u 
        CΔ    = kron(ones(H,H)|>tril, I(m))
        return new(ul_bold_const,uu_bold_const,CΔ)
    end 
end

function build_qp(A::AbstractArray, B::AbstractArray, Ψ_r::AbstractArray, weights::Weights, H::Int,z0::Vector, uₖ₋₁::Vector)

    p,m = length(z0), length(uₖ₋₁)

    Â,B̂ = augment_model(A, B)
    A_bold, B_bold = build_predmat(Â, B̂, H)
    
    ẑ0 = [z0;uₖ₋₁]
    
    r_bold = zeros((p+m)*H)
    for i in 1:H
        istart = (i-1)*(p+m)+1
        r_bold[istart:istart+p-1] .= Ψ_r[:,i]
    end 
    Q_bold, R_bold = weights.Q, weights.R

    P = B_bold'*Q_bold*B_bold + R_bold
    q = 2B_bold'*Q_bold*(A_bold*ẑ0 - r_bold)

    return P, q
end  

function get_dims(param::EDMDParameters)
    return param.dict.p, param.buffer.m 
end  

function augment_model(A::AbstractArray, B::AbstractArray)
    p,m = size(B) 
    Â = [A B; zeros(m,p) I(m)]
    B̂ = [B; I(m)]
    return Â,B̂
end  


mutable struct AdaptiveKMPC
    edmd_params::EDMDParameters                    
    weights::Weights                  
    Ψ_r::AbstractArray
    H::Int
    constr::Constraints 
    uₖ₋₁::AbstractArray     
    solver 
    function AdaptiveKMPC(edmd_params::EDMDParameters,Q::Vector,Qf::Vector,R::Vector,r::AbstractArray,H::Int,constr::Constraints)
        
        # initialize QP (required by OSQP, since sparsity pattern is 'locked in' when calling OSQP setup function) 
        p,m = get_dims(edmd_params)
        _,N = size(r)
        weights = Weights([Q;zeros(m)], [Qf;zeros(m)], R, H)
        Ψ_r = zeros(p,N)
        for i in 1:N
            Ψ_r[:,i] .= lifting(r[:,i],edmd_params.dict) 
        end
        A = ones(p,p)
        B = ones(p,m)
        P,q = build_qp(A, B, Ψ_r, weights, H, ones(p), ones(m))
        ul_bold, uu_bold, CΔ = constr.umin, constr.umax, constr.CΔ

        solver = OSQP.Model()
        OSQP.setup!(solver, P=sparse(P), q=vec(q), A=SparseMatrixCSC(CΔ),l=ul_bold, u=uu_bold, verbose=false, warm_start=true)  

        return new(edmd_params, weights, Ψ_r, H, constr, zeros(m), solver)
    end 
end

function update_buffer!(x::AbstractArray,u::AbstractArray,t::Union{AbstractArray, Float64},ctrl::AdaptiveKMPC)
    update_buffer!(x,u,t,ctrl.edmd_params.buffer)
end 


function get_control(x0::Vector, k::Int, ctrl::AdaptiveKMPC)
    
    # get internal model 
    z0 = lifting(x0, ctrl.edmd_params.dict)
    A,B = EDMD(ctrl.edmd_params)
    
    # build QP     
    H = ctrl.H 
    p, m = get_dims(ctrl.edmd_params)
    Ψ_r = zeros(p,H)
    _,N = size(ctrl.Ψ_r) 
    for i in k:k+H-1
        if i <= N 
            Ψ_r[:,i-k+1] = ctrl.Ψ_r[:,i]
        else
            Ψ_r[:,i-k+1] = ctrl.Ψ_r[:,end]
        end 
    end
    P, q = build_qp(A, B, Ψ_r, ctrl.weights, H, z0, ctrl.uₖ₋₁)
    CΔ, ul_bold, uu_bold = build_constraints(ctrl)
    OSQP.update!(ctrl.solver, Px=triu(sparse(P)).nzval, q=vec(q), l=ul_bold, u=uu_bold, Ax=sparse(CΔ).nzval)
    
    # solve QP 
    δu_bold = OSQP.solve!(ctrl.solver).x
    uk = ctrl.uₖ₋₁ + δu_bold[1:m]
    ctrl.uₖ₋₁ = uk 

    return uk 
end   


""" 
    build_constraints(umax,umin,hx,m)

Build input Constraints of the form lu <= Au * u <= uu.
"""
function build_constraints(ctrl::AdaptiveKMPC)
    ul_bold = ctrl.constr.ul_bold_const - repeat(ctrl.uₖ₋₁, ctrl.H)  
    uu_bold = ctrl.constr.uu_bold_const - repeat(ctrl.uₖ₋₁, ctrl.H)      
    return ctrl.constr.CΔ, ul_bold, uu_bold
end 
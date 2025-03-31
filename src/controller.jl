## Code for adaptive and static KMPC 
"""
    Weights(Q::Vector, Qf::Vector, R::Vector, H::Int)

Holds MPC weight matrices. 

# Arguments
- `Q::Vector`: Main diagonal of tracking error-weight matrix, see eq. (6a)  
- `Qf::Vector`: Main diagonal of tracking error-weight matrix for the final state 
- `R::Vector`: Main diagonal of control effort weight matrix, see eq. (6a)
- `H::Int`: Prediction horizon 
"""
struct Weights
    Q::AbstractArray
    R::AbstractArray
    function Weights(Q::Vector, Qf::Vector, R::Vector, H::Int)
        Qf = diagm(Qf)
        n = length(Q)
        Q = diagm(repeat(Q, inner=[1,1], outer=[1,H]) |> vec)   # equivaltent to  kron(I_H, diagm(Q))
        R = diagm(repeat(R, inner=[1,1], outer=[1,H]) |> vec)   # equivaltent to  kron(I_H, diagm(R))
        Q[end-n+1:end,end-n+1:end] .= Qf    
        return new(Q,R)
    end 
end 

"""
    build_predmat(Â::AbstractArray, B̂::AbstractArray, H::Int)

Builds prediction matrices, see eq.(8). 

# Arguments
- `Â::AbstractArray`: Augmented state transition matrix, see eq. (7)  
- `B̂::AbstractArray`: Augmented control matrix, see eq. (7) 
- `H::Int`: Prediction horizon 

# Returns 
- `A_bold::AbstractArray`: state transition matrix, see eq. (8)   
- `B_bold::AbstractArray`: control matrix, see eq. (8)
"""
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

"""
    Constraints(ul::AbstractArray,uu::AbstractArray,H::Int)

Builds the constant part of the input constraints, shown at the end of section II b).

# Arguments
- `ul::AbstractArray`: lower limit on controls, see eq. (6a)  
- `uu::AbstractArray`: upper limit on controls, see eq. (6a) 
- `H::Int`: Prediction horizon 
"""
struct Constraints 
    ul_bold_const::AbstractArray
    uu_bold_const::AbstractArray
    CΔ::AbstractArray
    function Constraints(ul::AbstractArray,uu::AbstractArray,H::Int)
        m = length(ul)
        ul_bold_const  = kron(ones(H), ul)              # constant part of lower bound on u 
        uu_bold_const  = kron(ones(H), uu)              # constant part of upper bound on u 
        CΔ    = kron(ones(H,H)|>tril, I(m))
        return new(ul_bold_const,uu_bold_const,CΔ)
    end 
end

"""
    build_cost_function(A::AbstractArray, B::AbstractArray, Ψ_r::AbstractArray, weights::Weights, H::Int,z0::Vector, uₖ₋₁::Vector)

Builds the cost function, see eq (9) .

# Arguments
- `A::AbstractArray`: state transition matrix, see eq. (5)  
- `B::AbstractArray`: control matrix, see eq. (5)  
- `Ψ_r::AbstractArray`: lifted reference over next H steps, appears in the error term in eq. (6a)  
- `weights::Weights`: struct with weight matrices   
- `H::Int`: Prediction horizon 
- `z0::Vector`: lifted state at current operating point, see eq. (8)
- `uₖ₋₁::Vector`: control vector applied at previous time step 

# Returns 
- `P::AbstractArray`: term with quadratic dependence on δu  
- `q::AbstractArray`: term with linear dependence on δu  
"""
function build_cost_function(A::AbstractArray, B::AbstractArray, Ψ_r::AbstractArray, weights::Weights, H::Int, z0::Vector, uₖ₋₁::Vector)

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

"""
    get_dims(param::EDMDParameters)

get the state dimension of the lifted state space model (p -> size of state vector, m -> size of control vector).
"""
function get_dims(param::EDMDParameters)
    return param.dict.p, param.buffer.m 
end  

"""
    augment_model(A::AbstractArray, B::AbstractArray)

Augments the Koopman state space model, see eq. (7)

# Arguments
- `A::AbstractArray`: state transition matrix, see eq. (5)  
- `B::AbstractArray`: control matrix, see eq. (5)  

# Returns 
- `Â::AbstractArray`: augmented state transition matrix  
- `B̂::AbstractArray`: augmented control matrix  
"""
function augment_model(A::AbstractArray, B::AbstractArray)
    p,m = size(B) 
    Â = [A B; zeros(m,p) I(m)]
    B̂ = [B; I(m)]
    return Â,B̂
end  

"""
    AdaptiveKMPC(edmd_params::EDMDParameters,Q::Vector,Qf::Vector,R::Vector,r::AbstractArray,H::Int,constr::Constraints)

Constructor for Adaptive KMPC struct, which holds parameters for building and solving the MPC problem. 

# Arguments
- `edmd_params::EDMDParameters`: struct holding parameters needed for EDMD (dictionary..)  
- `Q::Vector`: Main diagonal of tracking error-weight matrix, see eq. (6a)  
- `Qf::Vector`: Main diagonal of tracking error-weight matrix for the final state 
- `R::Vector`: Main diagonal of control effort weight matrix, see eq. (6a)
- `r::AbstractArray`: reference trajectory 
- `H::Int`: Prediction horizon 
- `constr::Constraints`: struct holding constraint matrices  

# Fields: 
- `edmd_params::EDMDParameters`: struct holding parameters needed for EDMD (dictionary..)  
- `weights::Weights`: struct holding weight matrices  
- `Ψ_r::AbstractArray`: lifted reference, appears in the error term in eq. (6a)  
- `H::Int`: Prediction horizon
- `constr::Constraints`: struct holding constraint matrices
- `uₖ₋₁::Vector`: control vector applied at previous time step   
- `solver`: OSQP solver object 

"""
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
        P,q = build_cost_function(A, B, Ψ_r, weights, H, ones(p), ones(m))
        ul_bold, uu_bold, CΔ = constr.ul_bold_const, constr.uu_bold_const, constr.CΔ

        solver = OSQP.Model()
        OSQP.setup!(solver, P=sparse(P), q=vec(q), A=SparseMatrixCSC(CΔ),l=ul_bold, u=uu_bold, verbose=false, warm_start=true)  

        return new(edmd_params, weights, Ψ_r, H, constr, zeros(m), solver)
    end 
end

"""
    update_buffer!(x::AbstractArray,u::AbstractArray,t::Union{AbstractArray, Float64},ctrl::AdaptiveKMPC)

Updates the circular buffer with new data.  

# Arguments
- `x::AbstractArray`: state vector   
- `u::AbstractArray`: control vector   
- `t::Union{AbstractArray, Float64}`: time 
- `ctrl::AdaptiveKMPC`: Adaptive KMPC struct  
"""
function update_buffer!(x::AbstractArray,u::AbstractArray,t::Union{AbstractArray, Float64},ctrl::AdaptiveKMPC)
    update_buffer!(x,u,t,ctrl.edmd_params.buffer)
end 

"""
    get_control(x0::Vector, k::Int, ctrl::AdaptiveKMPC)

Computes next control action.  

# Arguments
- `x0::Vector`: current state   
- `k::Int`: index in reference trajectors, corresponding to current time  
- `ctrl::AdaptiveKMPC`: Adaptive KMPC struct   

# Returns 
- `uk::Vector`: next control action 
"""
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
    P, q = build_cost_function(A, B, Ψ_r, ctrl.weights, H, z0, ctrl.uₖ₋₁)
    CΔ, ul_bold, uu_bold = build_constraints(ctrl)
    OSQP.update!(ctrl.solver, Px=triu(sparse(P)).nzval, q=vec(q), l=ul_bold, u=uu_bold, Ax=sparse(CΔ).nzval)
    
    # solve QP 
    δu_bold = OSQP.solve!(ctrl.solver).x
    uk = ctrl.uₖ₋₁ + δu_bold[1:m]
    ctrl.uₖ₋₁ = uk 

    return uk 
end   

## Code for linearization MPC
"""
    build_predmat(Â::AbstractArray, B̂::AbstractArray, H::Int)

Builds prediction matrices, for linearization MPC, equivalent to eq.(8), but with a matrix which will be multiplied with a constant offset term which comes into play due to linearization,
see for example https://www.researchgate.net/publication/319122854_Successive_Linearization_Based_Model_Predictive_Control_of_Variable_Stiffness_Actuated_Robots .

# Arguments
- `Â::AbstractArray`: Augmented state transition matrix, equivalent to eq. (7)  
- `B̂::AbstractArray`: Augmented control matrix, equivalent to eq. (7)  
- `H::Int`: Prediction horizon 

# Returns 
- `A_bold::AbstractArray`: state transition matrix, equivalent to eq. (8)   
- `B_bold::AbstractArray`: control matrix, equivalent to eq. (8)
- `E_bold::AbstractArray`: matrix which will be multiplied with constant offset vector 
"""
function build_predmat_linearized(Â::AbstractArray, B̂::AbstractArray, H::Int)
    n = size(Â,1)
    m = size(B̂,2)        
    A_bold=zeros(H*n,n);
    B_bold=zeros(H*n,H*m);
    E_bold=zeros(H*n,n);
    S = zero(Â)
    for i=1:H
        S = S + Â^(i-1);
        E_bold[1+n*(i-1):n*i,:]=S 
        A_bold[(1+(i-1)*n):i*n,1:n]=Â^i;
        for j=1:H            
            if(i>=j)                
                B_bold[(1+(i-1)*n):i*n,(1+(j-1)*m):j*m] = Â^(i-j)*B̂;                
            end
        end
    end
    return A_bold,B_bold,E_bold
end

"""
    build_cost_function(A::AbstractArray, B::AbstractArray, K::AbstractArray, Ψ_r::AbstractArray, weights::Weights, H::Int,z0::Vector, uₖ₋₁::Vector)

Builds the cost function, equivalent to eq (9), but for linearized dynamics.  

# Arguments
- `A::AbstractArray`: state transition matrix  
- `B::AbstractArray`: control matrix  
- `K::AbstractArray`: constant offset vector   
- `r::AbstractArray`: reference over next H steps
- `weights::Weights`: struct with weight matrices   
- `H::Int`: Prediction horizon 
- `x0::Vector`: state at current operating point, see eq. (8)
- `uₖ₋₁::Vector`: control vector applied at previous time step 

# Returns 
- `P::AbstractArray`: term with quadratic dependence on δu  
- `q::AbstractArray`: term with linear dependence on δu  
"""
function build_cost_function(A::AbstractArray, B::AbstractArray, K::AbstractArray, r::AbstractArray, weights::Weights, H::Int, x0::Vector, uₖ₋₁::Vector)

    n,m = length(x0), length(uₖ₋₁)
    
    Â,B̂,K̂ = augment_model(A, B, K)
    A_bold, B_bold, E_bold = build_predmat_linearized(Â, B̂, H)

    x̂0 = [x0;uₖ₋₁]
    
    r_bold = zeros((n+m)*H)
    for i in 1:H
        istart = (i-1)*(n+m)+1
        r_bold[istart:istart+n-1] .= r[:,i]
    end 
    Q_bold, R_bold = weights.Q, weights.R

    P = B_bold'*Q_bold*B_bold + R_bold
    q = 2B_bold'*Q_bold*(A_bold*x̂0 - r_bold + E_bold*K̂)

    return P, q
end  


"""
    augment_model(A::AbstractArray, B::AbstractArray)

Augments the state space model of linearized dynamics, equivalent to eq. (7).

# Arguments
- `A::AbstractArray`: state transition matrix  
- `B::AbstractArray`: control matrix  
- `K::AbstractArray`: constant offset vector   

# Returns 
- `Â::AbstractArray`: augmented state transition matrix  
- `B̂::AbstractArray`: augmented control matrix
- `K̂::AbstractArray`: augmented version of the offset vector 
"""
function augment_model(A::AbstractArray, B::AbstractArray, K::AbstractArray)
    n,m = size(B) 
    Â = [A B; zeros(m,n) I(m)]
    B̂ = [B; I(m)]
    K̂ = [K; zeros(m,1)]
    return Â,B̂,K̂ 
end  


"""
    linearizationMPC(model::nPendulum,Q::Vector,Qf::Vector,R::Vector,r::AbstractArray,H::Int,constr::Constraints,h::Float64)

Constructor for linearization MPC, with parameters for building and solving the MPC problem. 

# Arguments
- `model::nPendulum`: holds parameters for nonlinear pendulum model   
- `Q::Vector`: Main diagonal of tracking error-weight matrix, equivalent to eq. (6a)  
- `Qf::Vector`: Main diagonal of tracking error-weight matrix for the final state 
- `R::Vector`: Main diagonal of control effort weight matrix, equivalent to eq. (6a)
- `r::AbstractArray`: reference trajectory 
- `H::Int`: Prediction horizon 
- `constr::Constraints`: struct holding constraint matrices  
- `h::Float64`: discretization step size    

# Fields: 
- `model::nPendulum`: holds parameters for nonlinear pendulum model   
- `weights::Weights`: holds weight matrices  
- `r::AbstractArray`: reference trajectory  
- `H::Int`: Prediction horizon
- `constr::Constraints`: struct holding constraint matrices
- `uₖ₋₁::Vector`: control vector applied at previous time step  
- `h::Float64`: discretization step size  
- `solver`: OSQP solver object 

"""
mutable struct linearizationMPC
    model::nPendulum                    
    weights::Weights                  
    r::AbstractArray
    H::Int
    constr::Constraints 
    uₖ₋₁::AbstractArray
    h::Float64     
    solver 
    function linearizationMPC(model::nPendulum,Q::Vector,Qf::Vector,R::Vector,r::AbstractArray,H::Int,constr::Constraints,h::Float64)
        
        # initialize QP (required by OSQP, since sparsity pattern is 'locked in' when calling OSQP setup function) 
        n,m = model.n, Int(model.n/2)

        weights = Weights([Q;zeros(m)], [Qf;zeros(m)], R, H)
        A = ones(n,n)
        B = ones(n,m)
        P,q = build_cost_function(A, B, ones(n), r, weights, H, ones(n), ones(m))
        ul_bold, uu_bold, CΔ = constr.ul_bold_const, constr.uu_bold_const, constr.CΔ

        solver = OSQP.Model()
        OSQP.setup!(solver, P=sparse(P), q=vec(q), A=SparseMatrixCSC(CΔ),l=ul_bold, u=uu_bold, verbose=false, warm_start=true)  

        return new(model, weights, r, H, constr, zeros(m), h, solver)
    end 
end

"""
    get_control(x0::Vector, k::Int, ctrl::linearizationMPC)

Computes next control action. 

# Arguments
- `x0::Vector`: current state   
- `k::Int`: index in reference trajectors, corresponding to current time  
- `ctrl::linearizationMPC`: linearization MPC parameters    

# Returns 
- `uk::Vector`: next control action 
"""
function get_control(x0::Vector, k::Int, ctrl::linearizationMPC)
    
    # obtain A and B matrix by linearization + constant offset term from linearization 
    A,B = linearize_discretize_dynamics(x0, ctrl.uₖ₋₁,ctrl.h,ctrl.model)
    K = simulate(x0,ctrl.uₖ₋₁,ctrl.h, ctrl.model) - A*x0 - B*ctrl.uₖ₋₁

    H = ctrl.H 
    n, m = ctrl.model.n, Int(ctrl.model.n/2)
    r = zeros(n,H)
    _,N = size(ctrl.r) 
    for i in k:k+H-1
        if i <= N 
            r[:,i-k+1] = ctrl.r[:,i]
        else
            r[:,i-k+1] = ctrl.r[:,end]
        end 
    end
    P, q = build_cost_function(A, B, K, r, ctrl.weights, H, x0, ctrl.uₖ₋₁)
    CΔ, ul_bold, uu_bold = build_constraints(ctrl)
    OSQP.update!(ctrl.solver, Px=triu(sparse(P)).nzval, q=vec(q), l=ul_bold, u=uu_bold, Ax=sparse(CΔ).nzval)
    
    # solve QP 
    δu_bold = OSQP.solve!(ctrl.solver).x
    uk = ctrl.uₖ₋₁ + δu_bold[1:m]
    ctrl.uₖ₋₁ = uk 

    return uk 
end   

"""
    build_constraints(ctrl::Union{AdaptiveKMPC,linearizationMPC})

Builds the input constraints, shown at the end of section II b). 

# Arguments
- `ctrl::Union{AdaptiveKMPC,linearizationMPC}`: Adaptive KMPC or linearization MPC  

# Returns 
- `CΔ::AbstractArray`: inequality constraint matrix for input constraints   
- `ul_bold::Vector`: lower limits on controls 
- `uu_bold::Vector`: upper limits on controls 
"""
function build_constraints(ctrl::Union{AdaptiveKMPC,linearizationMPC})
    ul_bold = ctrl.constr.ul_bold_const - repeat(ctrl.uₖ₋₁, ctrl.H)  
    uu_bold = ctrl.constr.uu_bold_const - repeat(ctrl.uₖ₋₁, ctrl.H)      
    return ctrl.constr.CΔ, ul_bold, uu_bold
end 

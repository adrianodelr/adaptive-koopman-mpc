"""
    Databuffer(length_buffer, n_joints)

Temporal (fixed size) data container that holds last last recorded joint positions, velocities and controls ('θ','v','u')
and is used for the by the data driven model for training purposes.   

# Arguments
- `length_buffer::Int`: Buffer size or memory of the system
- `n_joints::Int`: Number of joint positions/angles/actuated joint
"""
mutable struct Databuffer
    θ::Vector{CircularBuffer{Float64}}
    v::Vector{CircularBuffer{Float64}}
    u::Vector{CircularBuffer{Float64}}
    t::CircularBuffer{Float64}    
    N::Int
    m::Int
    function Databuffer(N,m)
        θ = [CircularBuffer{Float64}(N)   for _ in 1:m]
        v = [CircularBuffer{Float64}(N)   for _ in 1:m]
        u = [CircularBuffer{Float64}(N-1) for _ in 1:m]
        t = CircularBuffer{Float64}(N)
        new(θ,v,u,t,N,m)
    end 
end 

function updateBuffer!(x,u,t,buffer::Databuffer)
    N = model.buffer.N
    m = model.buffer.m
    for _ in 1:N 
        push!(buffer.t, t)        
        for j in 1:m
            push!(buffer.θ[j], x[j])
            push!(buffer.v[j], x[j+m])            
            push!(buffer.u[j], u[j])            
        end 
    end
end 

"""
    Dictionary(Ψ::Vector{Function})
    

# Arguments
- `parent::Union{DoublePendulum, SinglePendulum}`: the system which is beeing for datasampling.
- `length_buffer::Int`: bufferlength which determined how many sampling points in the past the 'memory' includes.
"""
# TODO : Write documentation 
mutable struct Dictionary
    Ψ::Vector{Function}
    p::Int 
    function Dictionary(Ψ::Vector{Function})
        new(Ψ,length(Ψ))
    end 
end  


"""
    EDMDModel(parent,bsize,p,h)

Construct a Data Driven Linear Model, inheriting properties from 'parent' which can be 
of type 'DoublePendulum' or 'SinglePendulum'.

# Arguments
- `parent::Union{DoublePendulum, SinglePendulum}`: the system which is beeing for datasampling.
- `length_buffer::Int`: bufferlength which determined how many sampling points in the past the 'memory' includes.
- `p::Int`: Lifting dimension depends on the dictionary being used for the EDMD.
- `h::Float64`: The sampling intervall length a.k.a discretization time.
"""
# TODO : update documentation
mutable struct EDMDModel
    const parent::Union{DoublePendulum, SinglePendulum}
    buffer::Databuffer
    A::Matrix{Float64}
    B::Matrix{Float64}
    Ψ::Dictionary
    function DDMLinearModel(parent,N,Ψ)
        n,m = get_dims(parent)
        # p = out of dict
        A=zeros(p,p)
        B=zeros(p,m)
        buffer = Databuffer(N,m)
        new(parent,buffer,A,B,Ψ)
    end 
end  

function get_buffer_data(buffer::Databuffer)
    X,X⁺ = zeros(2buffer.m, buffer.N-1)
    U = zeros(buffer.m, buffer.N-1)
    return X,X⁺,U
end 

"""
    regression!(X̃lift, Ỹlift, model::DDMLinearModel)

Perform a regression on the lifted snapshot matrices 'X̃lift' and 'Ỹlift' to obtain the best fit linear operators 
A,B relating both matrices in the controlled setting
"""
function lifting_and_regression!(model::EDMDModel)
    Ψ = model.Ψ
    p = model.p
    N = model.buffer.N
    m = model.buffer.m

    X,X⁺,U=get_buffer_data(buffer::Databuffer)

    Z,Z⁺ = zeros(p,N-1),zeros(p,N-1)
    U = zeros(m, N-1)

    for i in 1:m
        for j in 1:N
            Z[i,:].= vec([Ψ[i](X[i,k]) for k in eachindex(Ψ)])
            Z⁺[i,:].= vec([Ψ[i](X⁺[i,k]) for k in eachindex(Ψ)])            
        end 
    end

    for i in 1:m
        U[i,:].=model.buffer.u[i]
    end

    Ω = [Z;U]
    K = Z⁺*pinv(Ω)
    model.A = K[:,1:p]
    model.B = K[:,p+1:end]
end     
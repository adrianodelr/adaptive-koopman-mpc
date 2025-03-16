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

function update_buffer!(x::AbstractArray, u::AbstractArray, t::Union{AbstractArray, Float64}, buffer::Databuffer)
    Nx = size(x,2)
    Nu = size(u,2)
    m = buffer.m
    for i in 1:Nx 
        if isone(Nx)  
            push!(buffer.t, t)
        else 
            push!(buffer.t, t[i])
        end 
        for j in 1:m
            push!(buffer.θ[j], x[j,i])
            push!(buffer.v[j], x[j+m,i])            
        end 
    end 
    for i in 1:Nu
        for j in 1:m
            push!(buffer.u[j], u[j,i])
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

function lifting(x::Vector, dict::Dictionary)
    Ψ = dict.Ψ
    z = vec([Ψ[k](x) for k in eachindex(Ψ)])
    return z 
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
mutable struct EDMDParameters
    buffer::Databuffer
    dict::Dictionary
    function EDMDParameters(m,N,dict)
        p = dict.p
        buffer = Databuffer(N,m)
        new(buffer,dict)
    end 
end  

function get_buffer_data(buffer::Databuffer)
    X,X⁺ = zeros(2buffer.m, buffer.N-1),zeros(2buffer.m, buffer.N-1)
    U = zeros(buffer.m, buffer.N-1)
    
    for i in 1:buffer.m        
        X[i,:]   = buffer.θ[i][1:end-1]'
        X[i+m,:] = buffer.v[i][1:end-1]'

        X⁺[i,:]   = buffer.θ[i][2:end]'
        X⁺[i+m,:] = buffer.v[i][2:end]'

        U[i,:] = buffer.u[i]'
    end 

    return X,X⁺,U
end 

"""
    regression!(X̃lift, Ỹlift, model::DDMLinearModel)

Perform a regression on the lifted snapshot matrices 'X̃lift' and 'Ỹlift' to obtain the best fit linear operators 
A,B relating both matrices in the controlled setting
"""
function EDMD(param::EDMDParameters)
    p = param.dict.p
    m = param.buffer.m
    N = param.buffer.N

    X,X⁺,U=get_buffer_data(param.buffer)

    Z,Z⁺ = zeros(p,N-1),zeros(p,N-1)

    for i in 1:N-1
        Z[:,i] .= lifting(X[:,i],param.dict) 
        Z⁺[:,i].= lifting(X⁺[:,i],param.dict)             
    end

    Ω = [Z;U]
    K = Z⁺*pinv(Ω)
    A = K[:,1:p]
    B = K[:,p+1:end]
    return A, B
end     

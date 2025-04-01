"""
    Databuffer(N,m)

Fixed size data buffer that holds joint positions, velocities and controls ('θ','v','u').

# Arguments
- `N::Int`: Buffer size or memory of the system
- `m::Int`: DOF

# Field 
- `θ::Vector{CircularBuffer{Float64}}`: Holds angular positions
- `ω::Vector{CircularBuffer{Float64}}`: Holds angular velocities
- `u::Vector{CircularBuffer{Float64}}`: Holds controls 
- `t::CircularBuffer{Float64}`: Holds time stamps 
- `N::Int`: Buffer size or memory of the system
- `m::Int`: DOF
"""
mutable struct Databuffer
    θ::Vector{CircularBuffer{Float64}}
    ω::Vector{CircularBuffer{Float64}}
    u::Vector{CircularBuffer{Float64}}
    t::CircularBuffer{Float64}    
    N::Int
    m::Int
    function Databuffer(N,m)
        θ = [CircularBuffer{Float64}(N)   for _ in 1:m]
        ω = [CircularBuffer{Float64}(N)   for _ in 1:m]
        u = [CircularBuffer{Float64}(N-1) for _ in 1:m]
        t = CircularBuffer{Float64}(N)
        new(θ,ω,u,t,N,m)
    end 
end 

"""
    update_buffer!(x::AbstractArray, u::AbstractArray, t::Union{AbstractArray, Float64}, buffer::Databuffer)

Updates the circular buffer with a time series of new data 

# Arguments
- `x::AbstractArray`: state trajectory matrix with dimensions n x N (N = length of the trajectory)  
- `u::AbstractArray`: control trajectory matrix with dimensions m x N-1  
- `t::Union{AbstractArray, Float64}`: time 
- `buffer::Databuffer` : buffer to update 
"""
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
            push!(buffer.ω[j], x[j+m,i])            
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

Dictionary which holds lifting functions  

# Arguments
- `Ψ::Vector{Function}`: Vector of lifting functions   

# Field 
- `Ψ::Vector{Function}`: Vector of lifting functions
- `p::Int`: dimension of the lifted state  
"""
mutable struct Dictionary
    Ψ::Vector{Function}
    p::Int 
    function Dictionary(Ψ::Vector{Function})
        new(Ψ,length(Ψ))
    end 
end  

"""
    lifting(x::Vector, dict::Dictionary)

Applies lifting dictionary to state vector x.
"""
function lifting(x::Vector, dict::Dictionary)
    Ψ = dict.Ψ
    z = vec([Ψ[k](x) for k in eachindex(Ψ)])
    return z 
end


"""
    EDMDParameters(m,N,dict)

Holds parameters for performing extended dynamic mode decomposition, see section II a).

# Arguments
- `N::Int`: Buffer size or memory of the system
- `m::Int`: DOF
- `dict::Dictionary`: bufferlength which determined how many sampling points in the past the 'memory' includes.

# Fields 
- `buffer::Databuffer`: the system which is beeing for datasampling.
- `dict::Dictionary`: bufferlength which determined how many sampling points in the past the 'memory' includes.
"""
mutable struct EDMDParameters
    buffer::Databuffer
    dict::Dictionary
    function EDMDParameters(m,N,dict)
        p = dict.p
        buffer = Databuffer(N,m)
        new(buffer,dict)
    end 
end  


"""
    get_buffer_data(buffer::Databuffer)

Constructs snapshot matrices X, X⁺,U, see section II a), from the data buffer. 

"""
function get_buffer_data(buffer::Databuffer)
    X,X⁺ = zeros(2buffer.m, buffer.N-1),zeros(2buffer.m, buffer.N-1)
    U = zeros(buffer.m, buffer.N-1)
    
    for i in 1:buffer.m        
        X[i,:]   = buffer.θ[i][1:end-1]'
        X[i+buffer.m,:] = buffer.ω[i][1:end-1]'

        X⁺[i,:]   = buffer.θ[i][2:end]'
        X⁺[i+buffer.m,:] = buffer.ω[i][2:end]'

        U[i,:] = buffer.u[i]'
    end 

    return X,X⁺,U
end 

"""
    EDMD(param::EDMDParameters)

Performs extended dynamic mode decomposition and returns a lifted state space model, see eq. (5).

# Arguments
- `param::EDMDParameters`: holds parameters for EDMD, including data. 

# Returns 
- `A::AbstractArray`: state transition matrix  
- `B::AbstractArray`: control matrix  
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

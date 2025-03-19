# TODO : write documentation
mutable struct nPendulum
    const m::Vector
    const l::Vector
    const lcom::Vector
    const I::Vector
    const n::Int 
    function nPendulum(m,l,lcom,I)
        @assert length(m)==length(l)==length(lcom)==length(I) "Model parameter vectors need to be of same dimension!"
        @assert length(m)!= 1 || length(m)!= 2 "Current implementation only single and double pendulum dynamics."
        if length(m)== 1 println("building simple pendulum") end 
        if length(m)== 2 println("building double pendulum") end 
        new(m,l,lcom,I,2length(m))
    end 
end 

# TODO : write documentation
function forward_dynamics_single_pendulum(x::Vector, u::Vector, model::nPendulum; g::Float64=9.81)
    lcom1,m1,I1=model.lcom[1],model.m[1],model.I[1]
    θ1d = x[2]
    ω1d = 1/I1 * (u[1] - m1*g*lcom1*sin(x[1])) 
    return [θ1d;ω1d]    
end 

# TODO : write documentation
function forward_dynamics_double_pendulum(x::Vector, u::Vector, model::nPendulum; g::Float64=9.81)
    l1,lcom1,m1,I1=model.l[1],model.lcom[1],model.m[1],model.I[1]
    l2,lcom2,m2,I2=model.l[2],model.lcom[2],model.m[2],model.I[2]
    q=x[1:2]
    v=x[3:4]
    m11 = I1 + I2 + l1^2*m2 + 2*l1*lcom2*m2*cos(q[2])  
    m22 = I2 
    m12 = I2 + l1*lcom2*m2*cos(q[2])
    m21 = I2 + l1*lcom2*m2*cos(q[2])     
    M = [m11 m12; m21 m22]      # intertia matrix 
    c1 = -2*l1*lcom2*m2*sin(q[2])*v[1]*v[2] - l1*lcom2*m2*sin(q[2])*v[2]^2
    c2 = l1*lcom2*m2*sin(q[2])*v[1]^2
    C = [c1; c2]                # coriolis matrix + centrigual forces 
    g1 = g*lcom1*m1*sin(q[1]) + g*m2*(l1*sin(q[1]) + lcom2*sin(sum(q)))    
    g2 = g*m2*lcom2*sin(sum(q))
    G = [g1;g2]                 # gravity vector 

    θd = v
    ωd = M\(u - C - G)    
    return [θd;ωd]    
end

""" 
    rk4(x, u, h, dynamics::Function)

Runge Kutta 4 discretization of the forward dynamic model of a dynamic system, given current state 'x', controls 'u', and discretization step 'h',
assuming a ZOH on controls. 
"""
function rk4(x, u, h, dynamics::Function)
    k1 = dynamics(x, u)
    k2 = dynamics(x + 0.5h*k1, u)
    k3 = dynamics(x + 0.5h*k2, u)
    k4 = dynamics(x + h*k3, u)
    return x + h/6*(k1 + 2k2 + 2k3 + k4)
end  

# TODO : write documentation
function simulate(x::Vector, u::Vector, h::Float64, model::nPendulum)
    if model.n==2
        return rk4(x, u, h, (x,u) -> forward_dynamics_single_pendulum(x, u, model))
    elseif model.n==4 
        return rk4(x, u, h, (x,u) -> forward_dynamics_double_pendulum(x, u, model))
    end 
end 

# TODO : write documentation
function linearize_discretize_dynamics(x::Vector,u::Vector, h, model::nPendulum)
    A = ForwardDiff.jacobian(x_-> simulate(x_, u, h, model),x)
    B = ForwardDiff.jacobian(u_-> simulate(x, u_, h, model),u)    
    return A,B
end  





# build model 
sp = nPendulum([1],[1],[1],[1]);
x0 = [π,0]
h = 0.01 
N = 400
m = 1 

# preceding experiment data 
xhist = zeros(sp.n, N)
uhist = zeros(Int(sp.n/2), N-1)    
xhist[:,1] = x0   
for i in 1:N-1
    uhist[:,i] .= 0.5*sin(i)
    xhist[:,i+1] .= simulate(xhist[:,i], uhist[:,i], h, sp)
end 

# Basis functions 
θ1_func(x) = x[1]
ω1_func(x) = x[2]
s1_func(x) = sin(x[1])
c1_func(x) = cos(x[1])
ω1s1_func(x) =  x[2]*sin(x[1])
ω1c1_func(x) =  x[2]*cos(x[1])

dict = Dictionary([θ1_func, ω1_func, s1_func, c1_func, ω1s1_func, ω1c1_func]) 

## build controller 

Q = [ones(sp.n);zeros(4)]                               # weights 
Qf = [ones(sp.n);zeros(4)]                             
R = ones(m)*0.2                             
r = [π, 0]                                  # reference 
H = 20                                      # prediction horizon 
umax = [6.0]                                  # constraints on controls 
umin = [-6.0]
constr = Constraints(umax, umin, H);
edmd_model = EDMDModel(m, N, dict);         # circular buffer is created within EDMD model 

ctrl = adaptiveKMPC(edmd_model, Q, Qf, R, r, H, constr);



# preceding experiment data 
xhistn = zeros(sp.n, N)
uhistn = zeros(Int(sp.n/2), N-1)    
xhistn[:,1] = xhist[:,end]   
for i in 1:N-1
    uhistn[:,i] = get_control(xhistn[:,i], ctrl)
    xhist[:,i+1] .= simulate(xhist[:,i], uhist[:,i], h, sp)
end 

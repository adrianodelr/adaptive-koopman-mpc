module AdaptiveKoopmanMPC

greet() = print("Have fun experimenting with adaptive Koopman MPC!")

# dependencies 
using LinearAlgebra
using OSQP 
using BlockArrays, SparseArrays
using DataStructures: CircularBuffer
using PCHIPInterpolation
using JLD2 
using LaTeXStrings  
using ForwardDiff 
using CairoMakie

# Ensure standard backend Cairo is activated  
CairoMakie.activate!()  

include("model.jl") 
export nPendulum, simulate

include("EDMD.jl")  
export update_buffer!, Dictionary, EDMDParameters 

include("controller.jl")       
export Constraints, AdaptiveKMPC, get_control, linearizationMPC 

include("plotting.jl")       
export plot_tracking_results, plot_performance_metrics

end # module AdaptiveKoopmanMPC

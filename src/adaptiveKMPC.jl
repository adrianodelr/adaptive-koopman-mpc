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
import GLMakie as GL 

# Ensure standard backend Cairo is activated  
CairoMakie.activate!()  

# TODO: add references to equations and parts of paper 
include("model.jl") 
include("EDMD.jl")  
include("controller.jl")       
include("plotting.jl")       


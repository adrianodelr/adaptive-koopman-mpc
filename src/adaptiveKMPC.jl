# dependencies 
using LinearAlgebra
using OSQP 
using BlockArrays, SparseArrays

using CairoMakie
# this is used to selectively display interactive plots with OpenGL
import GLMakie as GL
# Ensure standard backend Cairo is activated  
CairoMakie.activate!()  

include("EDMD.jl")
include("controller.jl")     
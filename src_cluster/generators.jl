## 11/15/2023
# Define the generators of SU(4) group
# Generalize the SU(3) Gell-Mann matrices

using ITensors
using LinearAlgebra

# let 
lambda₁ = [
    0  1  0  0
    1  0  0  0 
    0  0  0  0
    0  0  0  0
]  

lambda₂ = [
    0 -im 0  0
    im 0  0  0
    0  0  0  0
    0  0  0  0
]

lambda₃ = [
    1  0  0  0
    0 -1  0  0
    0  0  0  0
    0  0  0  0
]

lambda₄ = [
    0  0  1  0
    0  0  0  0
    1  0  0  0
    0  0  0  0
]

lambda₅ = [
    0  0 -im  0
    0  0  0   0
    im 0  0   0
    0  0  0   0
]

lambda₆ = [
    0  0  0  0
    0  0  1  0
    0  1  0  0
    0  0  0  0
]

lambda₇ = [
    0  0  0  0
    0  0 -im 0
    0  im 0  0
    0  0  0  0
]

lambda₈ = 1/sqrt(3)*[
    1  0  0  0
    0  1  0  0
    0  0 -2  0
    0  0  0  0
]

lambda₉ = [
    0  0  0  1
    0  0  0  0
    0  0  0  0
    1  0  0  0
]

lambda₁₀ = [
    0  0  0 -im
    0  0  0  0
    0  0  0  0
    im 0  0  0
]

lambda₁₁ = [
    0  0  0  0
    0  0  0  1
    0  0  0  0
    0  1  0  0
]

lambda₁₂ = [
    0  0  0  0
    0  0  0 -im
    0  0  0  0
    0  im 0  0
]

lambda₁₃ = [
    0  0  0  0
    0  0  0  0
    0  0  0  1
    0  0  1  0
]

lambda₁₄ = [
    0  0  0  0
    0  0  0  0 
    0  0  0 -im
    0  0 im  0
]

lambda₁₅ = 1/sqrt(6) * [
    1  0  0  0
    0  1  0  0
    0  0  1  0
    0  0  0 -3
]


#     return
# endx
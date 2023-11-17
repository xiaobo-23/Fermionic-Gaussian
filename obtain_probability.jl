## 11/15/2023
# Sample a two-site MPS to obtain the probability distribution of the local Hilbtert space
# Compute the probablility of one specific state as the cost function

using ITensors
using ITensors: orthocenter, sites, copy, complex, real 

# Sample a two-site MPS to Sz for the moment
function project_probability(m :: MPS, j :: Int, observable_type :: AbstractString)
    mpsLength = length(m)
    probablility = zeros(Float64, (2, 2))


    # Move the orthogonality center of the MPS to site j
    orthogonalize!(m, j)
    if orthocenter(m) != j
        error("sample: MPS m must have orthocenter(m) == j")
    end
    
    # Check the normalization of the MPS
    if abs(1.0 - norm(m[j])) > 1E-8
        error("sample: MPS is not normalized, norm=$(norm(m[j]))")
    end

    if observable_type == "Sz"
        tmp_projn = [[1, 0], [0, 1]]
    else
        error("sample: measurement type doesn't exist")
    end

    # Sample the target observables
    result = zeros(Int, 2)
    A = m[j]
    
    for ind in j:j+1
        tmpS = siteind(m, ind)
        d = dim(tmpS)
        pdisc = 0.0
        r = rand()

        n = 1 
        An = ITensor()
        pn = 0.0

        while n <= d
            projn = ITensor(tmpS)
            projn[tmpS => 1] = tmp_projn[n][1]
            projn[tmpS => 2] = tmp_projn[n][2]
        
            An = A * dag(projn)
            pn = real(scalar(dag(An) * An))
            pdisc += pn

            # Store the probablility of a state on each site inside a two-site unit cell
            if n == 1
                probablility[ind, 1] = pn
                probablility[ind, 2] = 1 - pn
            end

            (r < pdisc) && break
            n += 1
        end
        result[ind - j + 1] = n

        if ind < mpsLength
            A = m[ind + 1] * An
            A *= (1. / sqrt(pn))
        end

    end
    return probablility
end 
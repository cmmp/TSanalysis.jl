module TSanalysis

using Debug

export ami_calc

function range(x::Vector{Float64})
    return [min(x), max(x)]
end

function double_hist(series::Vector{Float64}, lag::Int64, partitions::Int64)
    # return a double histogram of the probabilities of observations,
    # considering the specified lag.
    hist = zeros(Float64,partitions, partitions)
    len = length(series)
    for i = 1:(len - lag)
        j = i + lag
        binx = max(ceil(series[i]*partitions),1)
        biny = max(ceil(series[j]*partitions),1)
        binx = min(binx, partitions)
        biny = min(biny, partitions)
        hist[binx,biny] += 1.
    end
    return hist / sum(hist)
end

function ami_calc(dataset::Matrix{Float64}; partitions::Int64 = 16, lagmax::Int64 = 20)
    # we will operate on columns:
    D = dataset'
    n = size(D)[2]
    A = zeros(lagmax,n)
    for i in 1:n
        series = D[:,i]
        series = (series - min(series)) / (max(series) - min(series))
        for j in 1:lagmax
            hist = double_hist(series, j, partitions)
            histx = nonzeros(vec(sum(hist,1)))
            hist = nonzeros(vec(hist))
            A[j,i] = sum(hist .* log(hist)) - 2 * sum(histx .* log(histx))
        end
    end
    return A'
end

srand(1234)

X = rand(10,100)
T = rand(10, 100)

A = ami_calc(X)

#writedlm("ami.dat", X) 

end

module TSanalysis

export ami_calc

function range(x::Vector{Float64})
    return [min(x), max(x)]
end

function double_hist(series::Vector{Float64}, lag::Int64, partitions::Int64)
    # return a double histogram of the probabilities of observations,
    # considering the specified lag.
end

function ami_calc(dataset::Matrix{Float64}; partitions::Int64 = 16, lagmax::Int64 = 20)
    # we will operate on columns:
    D = dataset'
    n = size(D)[2]
    A = zeros(lagmax,n)
    for i in 1:n
        series = D[:,i]
        series = (series - min(series)) / diff(range(series))
        corr = zeros(lagmax)
        for j in 1:lagmax
            hist = zeros(partitions, partitions)

        end
        break
    end
    #@show A
end

srand(1234)

X = rand(10,100)
T = rand(10, 100)

A = ami_calc(X)

#writedlm("ami.dat", X) 

end

module TSanalysis

#using Debug
import Stats.autocor

export ami_calc, acf_calc, adx_calc

function sma(x::Vector{Float64}; nperiods::Int64 = 14)
    ret::Vector{Float64} = zeros(length(x) - nperiods + 1)
    j = 1
    for i = 1:length(ret)
        ret[i] = sum(x[j:(j+nperiods-1)]) / nperiods
        j += 1
    end
    return ret
end

function ewma(x::Vector{Float64}; nperiods::Int64 = 14)
    y::Vector{Float64} = sma(x, nperiods = nperiods)

    mult = 2. / (nperiods + 1)

    j = nperiods + 1

    for i = 2:length(y)
        y[i] = (x[j] - y[i-1]) * mult + y[i-1]
        j += 1
    end
   
    return y
end

function range(x::Vector{Float64})
    return [min(x), max(x)]
end

function double_hist(series::Vector{Float64}, lag::Int64, partitions::Int64)
    # return a double histogram of the probabilities of observations,
    # considering the specified lag.
    hist = zeros(Float64, partitions, partitions)
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

function adx_calc(dataset::Matrix{Float64}; nperiods::Int64 = 14)
    D = dataset
    nseries, ncol = size(D)
    A = zeros(Float64, nseries, 1)

    # each row of D consists of [high | low | close]
    nlen = ncol / 3

    # reference: http://en.wikipedia.org/wiki/Average_directional_movement_index

    for i = 1:nseries
        # get the separate series:
        highs = D[i,1:nlen]
        lows = D[i,(nlen+1):(2*nlen)]
        closes = D[i,(2*nlen+1):ncol]

        upmove = float64([highs[i] - highs[i-1] for i = 2:nlen])
        downmove = float64([lows[i-1] - lows[i] for i = 2:nlen])

        # compute +DM and -DM:
        mask = (upmove .> downmove) & (upmove .> 0.)
        pDM = copy(upmove) # +DM
        pDM[!mask] = 0.

        mask = (downmove .> upmove) & (downmove .> 0.)
        mDM = copy(downmove) # -DM
        mDM[!mask] = 0.

        # compute the Average True Range (ATR): 
        # reference: http://en.wikipedia.org/wiki/Average_true_range


        # compute +DI and -DI:

    end

    return A
end

function acf_calc(dataset::Matrix{Float64}; ncoef::Int64 = 20)
    D = dataset'
    n = size(D)[2]
    A = zeros(Float64,ncoef,n)
    for i in 1:n
        A[:,i] = autocor(D[:,i], 1:ncoef)
    end
    return A'
end

function ami_calc(dataset::Matrix{Float64}; partitions::Int64 = 16, ncoef::Int64 = 20)
    # we will operate on columns:
    lagmax = ncoef - 1 # so we also compute lag = 0 
    D = dataset'
    n = size(D)[2]
    A = zeros(Float64,lagmax+1,n)
    for i in 1:n
        series = D[:,i]
        series = (series - min(series)) / (max(series) - min(series))
        for j in 1:(lagmax+1)
            hist = double_hist(series, j-1, partitions)
            histx = nonzeros(vec(sum(hist,1)))
            hist = nonzeros(vec(hist))
            A[j,i] = sum(hist .* log(hist)) - 2 * sum(histx .* log(histx))
        end
    end
    return A'
end

#srand(1234)

#X = rand(10,100)
#T = rand(10, 100)

#A = ami_calc(X)

#writedlm("ami.dat", X) 

end

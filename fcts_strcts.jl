using Distributions, Random, LinearAlgebra, Parameters, DataFrames, CSV, Documenter, StatsBase, Plots

"""
model parameters
- `𝔼P` = expected total phenotypic variance
- `h²` = expected proportion of additive genetic variance
- `L` = num host gen loci (set to zero for no host genetics)
- `SL` = num lin trans microbes (set to zero for no lin micr's)
- `SN` = num non-lin trans microbes (set to zero for no non-lin micr's)
- `SE` = num env acqr microbes (set to zero for no env micr's)
"""
@with_kw mutable struct Par

    # expected total phenotypic variance
    𝔼P::Float64 = 100.0

    # expected proportion of additive genetic variance
    h²::Float64 = 0.2

    # expected proportion of additive lineal microbial variance
    λ²::Float64 = 0.2

    # expected proportion of additive non-lineal vertical microbial variance
    ν²::Float64 = 0.2

    # expected proportion of additive external microbial variance
    ϵ²::Float64 = 0.2

    # probability of mutant allele
    p::Float64 = 0.5

    # number of genetic loci
    L::Int64 = 100

    # number of lineally transmitted taxa
    SL::Int64 = 100

    # number of non-lineally vertically transmitted taxa
    SN::Int64 = 100

    # number of uninheritable environmentally acquired taxa
    SE::Int64 = 100

    # number of hosts in each generation
    n::Int64 = 10^3

    # strength of directional selection
    s::Float64 = 1e-1

    # probability of lineal transmission per lineal microbe (given presence of that microbe)
    λ::Float64 = 0.25 / SL

    # probability of non-lineal vert transmission per non-lineal microbe (given presence of that microbe)
    ν::Float64 = 0.25 / (SN*n)

    #
    # numbers of times to repeat things
    #

    # number of times to initialize parental population
    nr::Int64 = 30
    
    # number of times to repeat selection experiments
    sr::Int64 = 5
    
    # number of times to repeat production of offspring generation
    or::Int64 = 5

    # total number of simulation runs = nr * sr * or

    # expected abundance of each microbial taxa in each host parent
    KL = 100
    KN = 100
    KE = 100

    # whether nonlineal microbes are sourced from selected parents
    Nsel = true

end

"""
parameters for generating population data (`PP` = `P`opulation `P`arameters)
- determines trait architecture shared by parents and offspring
- `γ` = additive genetic fx
- `ωL` = additive lineal microbial fx
- `ωN` = additive non-lineal microbial fx
- `ωE` = additive environmental microbial fx
"""
@with_kw mutable struct PP
    γ::Vector{Float64}
    ωL::Vector{Float64}
    ωN::Vector{Float64}
    ωE::Vector{Float64}
end

"""
population data (`PD` = `P`opulation `D`ata)
- `g` = host (diploid) genotype array
- `mL` = lineal microbe abundance matrix
- `mN` = non-lineal microbe abundance matrix
- `mE` = external microbe abundance matrix
- `z` = host trait vector
"""
@with_kw mutable struct PD
    g::Array{Bool,3}
    mL::Matrix{Int64}
    mN::Matrix{Int64}
    mE::Matrix{Int64}
    z::Vector{Float64}
end

"""
initialize parents
- returns `PP` and `PD`
"""
function init(par::Par)
    @unpack_Par par

    # expected additive genetic variance
    𝔼GA = h²*𝔼P

    # expected additive lineal microbial variance
    𝔼ML = λ²*𝔼P

    # expected additive non-lineal microbial variance
    𝔼MN = ν²*𝔼P

    # expected additive environmental microbial variance
    𝔼ME = ϵ²*𝔼P

    # expected developmental noise / environmental deviation
    𝔼E = (1 - h² - λ² - ν² - ϵ²) * 𝔼P

    # total number of microbe taxa across categories
    S = SL + SN + SE

    # unscaled variance of additive genetic fx
    vΓ = 𝔼GA / (2*p*(1-p))

    # unscaled variance of additive lineal microbial fx
    vL = 𝔼ML / (KL)

    # unscaled variance of additive non-lineal microbial fx
    vN = 𝔼MN / (KN)

    # unscaled variance of additive non-lineal microbial fx
    vE = 𝔼ME / (KE)

    #
    # parental parameters
    #

    # additive fx of genetic loci
    if L > 0
        γ = rand(Normal(0,√(vΓ/L)),L)
    else
        γ = 0
    end

    # additive fx of lineal microbes
    if SL > 0
        ωL = rand(Normal(0,√(vL/SL)),SL)
    else
        ωL = 0
    end

    # additive fx of non-lineal heritable microbes
    if SN > 0
        ωN = rand(Normal(0,1√(vN/SN)),SN)
    else
        ωN = 0
    end

    # additive fx of environmental microbes
    if SE > 0
        ωE = rand(Normal(0,√(vE/SE)),SE)
    else
        ωE = 0
    end

    #
    # parental data
    #

    # genetic array
    g = rand(Bernoulli(p),n,L,2)

    # lineal microbiome matrix
    mL = rand(Poisson(KL),n,SL)

    # nonlineal microbiome matrix
    mN = rand(Poisson(KN),n,SN)

    # environmental microbiome matrix
    mE = rand(Poisson(KE),n,SE)

    #
    # collecting parental parameters and data
    #

    # collect parental parameters
    pp = PP(
        γ = γ,
        ωL = ωL,
        ωN = ωN,
        ωE = ωE
    )

    # collect parental data
    pd = PD(
        g = g,
        mL = mL,
        mN = mN,
        mE = mE,
        z = zeros(n)
    )

    #
    # traits
    #

    # genetic-microbic values
    zgm = gmvals(pp,pd,"γλνε")

    # (generative) developmental noise / environmental deviation
    e = rand(Normal(0,√𝔼E),n)

    # compute parental trait values from (generative) model
    pd.z = zgm .+ e

    return pp, pd
end

"""
compute genetic-microbic value of host traits
- used during formation of host parents and host offspring data
- `e` will need to be added later if noise is desired
"""
function gmvals(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    # allele count matrix
    c = g[:,:,1] .+ g[:,:,2]

    # (generative) additive genetic value
    gA = c*γ

    # (generative) lineal microbial value
    ml = mL*ωL

    # (generative) nonlineal microbial value
    mn = mN*ωN

    # (generative) environmental microbial value
    me = mE*ωE

    # (generative) additive microbial value
    mA = ml .+ mn .+ me

    # compute parental trait values from (generative) model
    zgm = gA .+ mA

    return zgm
end

"""
compute genetic-microbic value of host traits
- used for the analysis of variance
- just uses model of `zgm` instead of actually averaging
- `fs` stands for *factors*
    - `fs = "γλνε"` for estimating zgm from all microbes
    - `fs = "γλν"` for estimating zgm excluding env microbes
    - `fs = "γλν"` for estimating zgm including only lineal microbes
    - `fs = "γ"` for estimating zgm from only genetic data
"""
function gmvals(pp::PP,pd::PD,fs::String)
    @unpack_PP pp
    @unpack_PD pd

    # allele count matrix
    c = g[:,:,1] .+ g[:,:,2]

    # (mechanistic) additive genetic value
    gA = c*γ

    if fs == "γλνε"

        # (mechanistic) additive lineal microbial value
        ml = mL*ωL
        
        # (mechanistic) additive nonlineal microbial value
        mn = mN*ωN

        # (mechanistic) additive environmental microbial value
        me = mE*ωE

    elseif fs == "γλν"
    
        # (mechanistic) additive lineal microbial value
        ml = mL*ωL
    
        # (mechanistic) additive nonlineal microbial value
        mn = mN*ωN

        # (mechanistic) additive environmental microbial value
        me = 0

    elseif fs == "γλ"
        
        # (mechanistic) additive lineal microbial value
        ml = mL*ωL
    
        # (mechanistic) additive nonlineal microbial value
        mn = 0

        # (mechanistic) additive environmental microbial value
        me = 0

    elseif fs == "γ"
            
        # (mechanistic) additive lineal microbial value
        ml = 0

        # (mechanistic) additive nonlineal microbial value
        mn = 0

        # (mechanistic) additive environmental microbial value
        me = 0

    else

        print("invalid factor string")

        return 

    end

    # (mechanistic) additive microbial value
    mA = ml .+ mn .+ me

    # compute parental trait values from (generative) model
    zgm = gA .+ mA
    
    return zgm
end

"""
All Additive FX
- includes genetic, lineal, non-lineal, and ext
"""
function AAFX(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    # compute genetic-microbic values
    zgm = gmvals(pp,pd,"γλνε")
 
    # allele count matrix
    c = g[:,:,1] .+ g[:,:,2]

    # microbiome matrix
    m = [mL mN mE]

    # all the thingis
    a = [c m]

    # gene-microbe covariance matrix
    Σ = cov(a)

    # additive fx
    α = inv(Σ)*cov(a,zgm)

    return α[:,1]

end

"""
All Additive Variance
"""
function AA(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    α = AAFX(pp,pd)

    c = g[:,:,1] .+ g[:,:,2]

    a = [c mL mN mE]

    Aₐ = var(a*α)

    return Aₐ
end

"""
Genetic and Lineal Additive FX
"""
function LAFX(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    # compute genetic-microbic values
    zgm = gmvals(pp,pd,"γλ")
 
    # allele count matrix
    c = g[:,:,1] .+ g[:,:,2]

    # genes and lineal microbes
    a = [c mL]

    # gene-microbe covariance matrix
    Σ = cov(a)

    # additive fx
    α = inv(Σ)*cov(a,zgm)

    return α[:,1]

end

"""
Genetic and Lineal Additive Variance
"""
function LA(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    α = LAFX(pp,pd)
    
    c = g[:,:,1] .+ g[:,:,2]

    a = [c mL]

    Lₐ = var(a*α)

    return Lₐ
end

"""
Additive FX of Transmitted Materials
"""
function TAFX(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    # compute genetic-microbic values
    zgm = gmvals(pp,pd,"γλν")

    # allele count matrix
    c = g[:,:,1] .+ g[:,:,2]

    # genes and lineal microbes
    a = [c mL mN]

    # gene-microbe covariance matrix
    Σ = cov(a)

    # additive fx
    α = inv(Σ)*cov(a,zgm)

    return α[:,1]

end

"""
Additive Variance of Transmitted Materials
"""
function TA(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    α = TAFX(pp,pd)

    c = g[:,:,1] .+ g[:,:,2]

    a = [c mL mN]

    Tₐ = var(a*α)

    return Tₐ
end

"""
Additive Genetic FX
"""
function AGFX(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    # compute genetic-microbic values
    zgm = gmvals(pp,pd,"γ")

    # allele count matrix
    c = g[:,:,1] .+ g[:,:,2]

    # genetic covariance matrix
    Σ = cov(c)

    # additive fx for each locus
    γ = inv(Σ)*cov(c,zgm)

    return γ[:,1]

end

"""
Additive Genetic Variance
"""
function GA(pp::PP,pd::PD)
    @unpack_PP pp
    @unpack_PD pd

    γ = AGFX(pp,pd)

    c = g[:,:,1] .+ g[:,:,2]

    Gₐ = var(c*γ)

    return Gₐ
end

"""
selection experiment
- returns indices of selected pairs and num unique pairs
"""
function selection(par::Par, pd::PD)
    @unpack_Par par
    @unpack_PD pd

    # probability of being a member of a parental pair chosen by an offspring
    # proportional to fitness (both types of relative fitness, and absolute fitness)
    w = normalize(exp.(s*z),1)

    # select mating pairs
    wv = Weights(w)
    pinds = collect(1:n)
    p1 = sample(pinds,wv,n)
    p2 = sample(pinds,wv,n)
    pairwents = [p1 p2]
    
    # count number of times each individual in parent population occurs in a mating pair
    W = zeros(Int64, n)
    for i in 1:n
        W[i] += length(findall(x->x==i,pairwents[:,1]))
        W[i] += length(findall(x->x==i,pairwents[:,2]))
    end

    # realized absolute fitness
    W /= 2

    # expected absolute fitness
    𝔼W = n*w

    return pairwents, W, 𝔼W
end

"""
generate offspring traits
- includes transmission of genes and microbes
"""
function offspring(par::Par, pp::PP, pd::PD, pairs)
    @unpack_Par par
    @unpack_PP pp
    @unpack_PD pd

    # whether nonlineal microbes are sourced from selected parents
    if Nsel
        # indices of selected parents        
        vp = vec(pairs)
    else
        # indices of all parents
        vp = 1:n
    end

    # pick nonlineal donors
    nld = sample(vp,n)

    # non-lineal vert trans
    mNₚ = rand.(Poisson.(mN[nld,:]))

    # transmission of genetic material and lineal microbes
    gₚ = Array{Bool}(undef,n,L,2)
    mLₚ = Array{Int64}(undef,n,SL)
    for i in 1:n
        
        p1 = pairs[i,1]
        p2 = pairs[i,2]

        #
        # inheritance of genetic material
        #

        # decide which loci recombine
        free_recomb = rand(Bernoulli(),L)

        # pick one chromosome from each parent
        p1ch = rand(Bernoulli()) .+ 1
        p2ch = rand(Bernoulli()) .+ 1

        # form offspring chromosomes
        ch1 = Array{Bool}(undef,L)
        ch2 = Array{Bool}(undef,L)
        for l in 1:L
            if free_recomb[l]
                ch1[l] = g[p1,l,p1ch]
                ch2[l] = g[p2,l,p2ch]
            else
                ch1[l] = g[p1,l,p2ch]
                ch2[l] = g[p2,l,p1ch]
            end
        end

        gₚ[i,:,1] = ch1
        gₚ[i,:,2] = ch2

        # poisson distr abundances centered on mid-parent for each taxa
        mLₘ = 0.5 .* (mL[p1,:] .+ mL[p2,:])
        mLₚ[i,:] = rand.(Poisson.(mLₘ))

    end

    # offspring allele count matrix
    cₚ = gₚ[:,:,1] .+ gₚ[:,:,2]

    # offspring environmental microbiome matrix
    mEₚ = rand(Poisson(KE),n,SE)

    #
    # offspring trait components
    #

    # additive genetic value
    gAₚ = cₚ*γ

    # lineal microbial value
    mlₚ = mLₚ*ωL

    # nonlineal microbial value
    mnₚ = mNₚ*ωN

    # environmental microbial value
    meₚ = mEₚ*ωE

    # additive microbial value
    mAₚ = mlₚ .+ mnₚ .+ meₚ

    # expected developmental noise / environmental deviation
    # defined for parental pop, but stays same here for consistency
    𝔼E = (1 - h² - λ² - ν² - ϵ²)*𝔼P

    # offspring dev noise / env dev
    eₚ = rand(Normal(0,√𝔼E),n)

    # compute offspring trait values
    zₚ = gAₚ .+ mAₚ .+ eₚ

    # offspring data
    pdₚ = PD(
        g = gₚ,
        mL = mLₚ,
        mN = mNₚ,
        mE = mEₚ,
        z = zₚ
    )

    # offspring mean trait
    z̄ₚ = mean(zₚ)

    return z̄ₚ, pdₚ

end

function timeseries(par::Par,pp::PP,pd₀::PD,T::Int64)
    @unpack_PP pp
    @unpack_PD pd₀

    pdₜ = Vector{PD}(undef, T)

    pdₜ[1] = pd₀

    for t in 1:(T-1)

        pairs, W, 𝔼W = selection(par, pdₜ[t])

        z̄ₚ, pdₜ[t+1] = offspring(par, pp, pdₜ[t], pairs)

    end

    return pdₜ

end

function runs(par::Par)
    @unpack_Par par

    pp = Vector{PP}(undef,nr)
    pd = Vector{PD}(undef,nr)

    μΔz̄ = zeros(nr)
    μAᵧβ = zeros(nr)
    μAₗβ = zeros(nr)
    μAₙβ = zeros(nr)
    μAᵩβ = zeros(nr)

    Aᵧ = zeros(nr)
    Aₗ = zeros(nr)
    Aₙ = zeros(nr)
    Aᵩ = zeros(nr)

    Δz̄ = zeros(sr)
    Aᵧβ = zeros(sr)
    Aₗβ = zeros(sr)
    Aₙβ = zeros(sr)
    Aᵩβ = zeros(sr)

    for i in 1:nr

        pp[i], pd[i] = init(par)

        Aᵧ[i] = GA(pp[i],pd[i])

        Aₗ[i] = LA(pp[i],pd[i])

        Aₙ[i] = TA(pp[i],pd[i])

        Aᵩ[i] = AA(pp[i],pd[i])

        z̄ = mean(pd[i].z)

        P = var(pd[i].z)
        
        for j in 1:sr

            pairs, W, 𝔼W = selection(par, pd[i])

            z̄ₚ = zeros(or)
            for k in 1:or
                z̄ₚ[k], pdₚ = offspring(par, pp[i], pd[i], pairs)
            end

            Δz̄[j] = mean(z̄ₚ) - z̄

            β = cov(W,pd[i].z) / P

            Aᵧβ[j] = Aᵧ[i]*β

            Aₗβ[j] = Aₗ[i]*β

            Aₙβ[j] = Aₙ[i]*β

            Aᵩβ[j] = Aᵩ[i]*β

        end

        μΔz̄[i] = mean(Δz̄)
        μAᵧβ[i] = mean(Aᵧβ)
        μAₗβ[i] = mean(Aₗβ)
        μAₙβ[i] = mean(Aₙβ)
        μAᵩβ[i] = mean(Aᵩβ)

    end

    return μΔz̄, μAᵧβ, μAₗβ, μAₙβ, μAᵩβ

end

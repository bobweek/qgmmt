include("fcts_strcts.jl")

theme(:vibrant)

# number host gens
T = 10

# number of time series to repeat
ts_rep = 20

nr = 50 # number times to draw background parameters
sr = 10 # number times to draw selected parents
or = 10 # number times to form offspring


function ts_ensemble(par::Par,T::Int64,ts_rep::Int64)

    z̄ = Matrix{Float64}(undef,T,ts_rep)
    P = Matrix{Float64}(undef,T,ts_rep)

    for i in 1:ts_rep

        pp, pd = init(par)

        ts = timeseries(par,pp,pd,T)
        
        for t in 1:T
            z̄[t,i] = mean(ts[t].z)
            P[t,i] = var(ts[t].z)
        end
        
        z̄[:,i] = z̄[:,i] .- z̄[:,i][1]

    end

    return z̄, P

end

#
# nonlineals contribute to selection response
#

# all

parγλνε∅ = Par(nr = nr, 
    𝔼P = 100.0,
    h² = 0.25, λ² = 0.25, ν² = 0.25, ϵ² = 0.25, 
    L = 100, SL = 100, SN = 100, SE = 100, 
    sr = sr, or = or)

z̄γλνε∅, Pγλνε∅ = ts_ensemble(parγλνε∅,T,ts_rep)
μz̄γλνε∅ = vec(mean(z̄γλνε∅,dims=2))
σz̄γλνε∅ =  vec(std(z̄γλνε∅,dims=2))

ymx = maximum(μz̄γλνε∅)

# no novel

parγλν∅ = Par(nr = nr, 
    𝔼P = 75.0,
    h² = 0.33, λ² = 0.33, ν² = 0.33, ϵ² = 0.001, 
    L = 100, SL = 100, SN = 100, SE = 1, 
    sr = sr, or = or)

z̄γλν∅, Pγλν∅ = ts_ensemble(parγλν∅,T,ts_rep)
μz̄γλν∅ = vec(mean(z̄γλν∅,dims=2))
σz̄γλν∅ =  vec(std(z̄γλν∅,dims=2))

ymx = max(ymx,maximum(μz̄γλν∅))

# just genes + lineal

parγλ∅ = Par(nr = nr, 
    𝔼P = 50.0,
    h² = 0.49, λ² = 0.49, ν² = 0.001, ϵ² = 0.001, 
    L = 100, SL = 100, SN = 1, SE = 1, 
    sr = sr, or = or)

z̄γλ∅, Pγλ∅ = ts_ensemble(parγλ∅,T,ts_rep)
μz̄γλ∅ = vec(mean(z̄γλ∅,dims=2))
σz̄γλ∅ =  vec(std(z̄γλ∅,dims=2))

# just genes

parγ∅ = Par(nr = nr, 
    𝔼P = 25.0,
    h² = 0.99, λ² = 0.001, ν² = 0.001, ϵ² = 0.001, 
    L = 100, SL = 1, SN = 1, SE = 1, 
    sr = sr, or = or)

z̄γ∅, Pγ∅ = ts_ensemble(parγ∅,T,ts_rep)
μz̄γ∅ = vec(mean(z̄γ∅,dims=2))
σz̄γ∅ =  vec(std(z̄γ∅,dims=2))

#
# nonlineals dont contribute to selection response
#


# all

parγλνεψ = Par(nr = nr, 
    𝔼P = 100.0,
    h² = 0.25, λ² = 0.25, ν² = 0.25, ϵ² = 0.25, 
    L = 100, SL = 100, SN = 100, SE = 100, 
    sr = sr, or = or,
    Nsel=false)

z̄γλνεψ, Pγλνεψ = ts_ensemble(parγλνεψ,T,ts_rep)
μz̄γλνεψ = vec(mean(z̄γλνεψ,dims=2))
σz̄γλνεψ =  vec(std(z̄γλνεψ,dims=2))

# no ext

parγλνψ = Par(nr = nr, 
    𝔼P = 75.0,
    h² = 0.33, λ² = 0.33, ν² = 0.33, ϵ² = 0.001, 
    L = 100, SL = 100, SN = 100, SE = 1, 
    sr = sr, or = or,
    Nsel=false)

z̄γλνψ, Pγλνψ = ts_ensemble(parγλνψ,T,ts_rep)
μz̄γλνψ = vec(mean(z̄γλνψ,dims=2))
σz̄γλνψ =  vec(std(z̄γλνψ,dims=2))

# just genes + lineal

parγλψ = Par(nr = nr, 
    𝔼P = 50.0,
    h² = 0.49, λ² = 0.49, ν² = 0.001, ϵ² = 0.001, 
    L = 100, SL = 100, SN = 1, SE = 1, 
    sr = sr, or = or,
    Nsel=false)

z̄γλψ, Pγλψ = ts_ensemble(parγλψ,T,ts_rep)
μz̄γλψ = vec(mean(z̄γλψ,dims=2))
σz̄γλψ =  vec(std(z̄γλψ,dims=2))

# just genes

parγψ = Par(nr = nr, 
    𝔼P = 25.0,
    h² = 0.99, λ² = 0.001, ν² = 0.001, ϵ² = 0.001, 
    L = 100, SL = 1, SN = 1, SE = 1, 
    sr = sr, or = or,
    Nsel=false)

z̄γψ, Pγψ = ts_ensemble(parγψ,T,ts_rep)
μz̄γψ = vec(mean(z̄γψ,dims=2))
σz̄γψ =  vec(std(z̄γψ,dims=2))

#
# plots plots plots ...
#

tspl∅ = plot(1:T,μz̄γλνε∅,ribbon=σz̄γλνε∅,fillalpha=0.3,
                legend=false, 
                title="Post-Selected\n\n",
                ylim=(0,ymx),
                color=colorant"#ffb86c");
plot!(1:T,μz̄γλν∅,ribbon=σz̄γλν∅,fillalpha=0.3,color=colorant"#56a0d3");
plot!(1:T,μz̄γλ∅,ribbon=σz̄γλ∅,fillalpha=0.3,color=colorant"#50fa7b");
plot!(1:T,μz̄γ∅,ribbon=σz̄γ∅,fillalpha=0.3,color=colorant"#ff79c6");
xlabel!("Host Generations");
ylabel!("\nMean Trait Value");

tsplψ = plot(1:T,μz̄γλνεψ, ribbon=σz̄γλνεψ,fillalpha=0.3,
                label="GLNV", 
                title="Pre-Selected\n\n",
                ylim=(0,ymx),
                color=colorant"#ffb86c");
plot!(1:T,μz̄γλνψ,ribbon=σz̄γλνψ,fillalpha=0.3, label="GLN",color=colorant"#56a0d3");
plot!(1:T,μz̄γλψ,ribbon=σz̄γλψ,fillalpha=0.3, label="GL",color=colorant"#50fa7b");
plot!(1:T,μz̄γψ,ribbon=σz̄γψ,fillalpha=0.3, label="G",color=colorant"#ff79c6");
xlabel!("Host Generations\n");

tspl = plot(tspl∅,tsplψ, size=(650,250), dpi=400, linewidth=2)

savefig(tspl,"fig3.png")


#
# analyze gene-microbe covariance dynamics ...
#

function InfNaN2Zero(ρ) # ρ is a matrix
    naninds = findall(x->isnan(x),ρ)            
    for i in naninds
        ρ[i[1],i[2]] = 0
    end
    infinds = findall(x->isinf(x),ρ)
    for i in infinds
        ρ[i[1],i[2]] = 0
    end
    return ρ
end

function InfNan2Bye(ρ)
    ρ = ρ[.!isnan.(ρ)]
    ρ = ρ[.!isinf.(ρ)]
    return ρ
end

function ts_ensemble_detailed(par::Par,T::Int64,ts_rep::Int64)

    mρᵧₑ = Matrix{Float64}(undef,T,ts_rep)
    sρᵧₑ = Matrix{Float64}(undef,T,ts_rep)

    mρᵧₗ = Matrix{Float64}(undef,T,ts_rep)
    sρᵧₗ = Matrix{Float64}(undef,T,ts_rep)

    mρᵧₙ = Matrix{Float64}(undef,T,ts_rep)
    sρᵧₙ = Matrix{Float64}(undef,T,ts_rep)

    c̄  = Array{Float64}(undef, ts_rep, par.L,  T)    
    m̄ₗ = Array{Float64}(undef, ts_rep, par.SL, T)
    m̄ₙ = Array{Float64}(undef, ts_rep, par.SN, T)
    m̄ₑ = Array{Float64}(undef, ts_rep, par.SE, T)

    μRᵧₗ = Vector{Float64}(undef,T)
    μRᵧₙ = Vector{Float64}(undef,T)
    μRᵧₑ = Vector{Float64}(undef,T)

    σRᵧₗ = Vector{Float64}(undef,T)
    σRᵧₙ = Vector{Float64}(undef,T)
    σRᵧₑ = Vector{Float64}(undef,T)

    pp, pd = init(par)

    for i in 1:ts_rep        

        ts = timeseries(par,pp,pd,T)
        
        for t in 1:T

            # allele counts
            c = ts[t].g[:,:,1] .+ ts[t].g[:,:,2]

            # mean allele count at each locus across hosts (at time t for rep i)
            c̄[i,:,t] = vec(mean(c,dims=1))

            # lineal microbe abundances
            mₗ = ts[t].mL

            # mean lineal abundance for at each taxa across hosts (at time t for rep i)
            m̄ₗ[i,:,t] = vec(mean(mₗ,dims=1))

            # non-lineal microbe abundances
            mₙ = ts[t].mN

            # mean lineal abundance for at each taxa across hosts (at time t for rep i)
            m̄ₙ[i,:,t] = vec(mean(mₙ,dims=1))

            # novel microbe abundances
            mₑ = ts[t].mE

            # mean lineal abundance for at each taxa across hosts (at time t for rep i)
            m̄ₑ[i,:,t] = vec(mean(mₑ,dims=1))

            # correlation between allele counts and lineals
            ρᵧₗ = abs.(vec(cor(c,mₗ)))
            ρᵧₗ = InfNan2Bye(ρᵧₗ)
            mρᵧₗ[t,i] = mean(ρᵧₗ)
            sρᵧₗ[t,i] = std(ρᵧₗ)
            
            # correlation between allele counts and non-lineals
            ρᵧₙ = abs.(vec(cor(c,mₙ)))
            ρᵧₙ = InfNan2Bye(ρᵧₙ)
            mρᵧₙ[t,i] = mean(ρᵧₙ)
            sρᵧₙ[t,i] = std(ρᵧₙ)

            # correlation between allele counts and novels
            ρᵧₑ = abs.(vec(cor(c,mₑ)))
            ρᵧₑ = InfNan2Bye(ρᵧₑ)
            mρᵧₑ[t,i] = mean(ρᵧₑ)
            sρᵧₑ[t,i] = std(ρᵧₑ)

        end        

    end

    # gene-microbe covs among replicates
    for t in 1:T

        Rᵧₗ = abs.(vec(cor(c̄[:,:,t], m̄ₗ[:,:,t])))
        Rᵧₙ = abs.(vec(cor(c̄[:,:,t], m̄ₙ[:,:,t])))
        Rᵧₑ = abs.(vec(cor(c̄[:,:,t], m̄ₑ[:,:,t])))

        Rᵧₗ = InfNan2Bye(Rᵧₗ)
        Rᵧₙ = InfNan2Bye(Rᵧₙ)
        Rᵧₑ = InfNan2Bye(Rᵧₑ)

        μRᵧₗ[t] = mean(Rᵧₗ)
        μRᵧₙ[t] = mean(Rᵧₙ)
        μRᵧₑ[t] = mean(Rᵧₑ)

        σRᵧₗ[t] = std(Rᵧₗ)
        σRᵧₙ[t] = std(Rᵧₙ)
        σRᵧₑ[t] = std(Rᵧₑ)
    
    end

    μρᵧₗ = vec(mean(mρᵧₗ,dims=2)) 
    μρᵧₙ = vec(mean(mρᵧₙ,dims=2))
    μρᵧₑ = vec(mean(mρᵧₑ,dims=2)) 
    σρᵧₗ = vec(std(mρᵧₗ,dims=2))  
    σρᵧₙ = vec(std(mρᵧₙ,dims=2))
    σρᵧₑ = vec(std(mρᵧₑ,dims=2))  
    vρᵧₗ = vec(mean(sρᵧₗ,dims=2))  
    vρᵧₙ = vec(mean(sρᵧₙ,dims=2))
    vρᵧₑ = vec(mean(sρᵧₑ,dims=2))

    return  μρᵧₗ, μρᵧₙ, μρᵧₑ, 
            σρᵧₗ, σρᵧₙ, σρᵧₑ, 
            vρᵧₗ, vρᵧₙ, vρᵧₑ,
            
            μRᵧₗ, μRᵧₙ, μRᵧₑ, 
            σRᵧₗ, σRᵧₙ, σRᵧₑ

end

T = 10

ts_rep = 100

# these are the same as for parγλνε∅,
#      but sr and or not applicable
par∅ = Par(nr = nr, 
    𝔼P = 100.0,
    h² = 0.25, λ² = 0.25, ν² = 0.25, ϵ² = 0.25, 
    L = 100, SL = 100, SN = 100, SE = 100, 
    Nsel=false)

μρᵧₗ∅, μρᵧₙ∅, μρᵧₑ∅, 
σρᵧₗ∅, σρᵧₙ∅, σρᵧₑ∅,
vρᵧₗ∅, vρᵧₙ∅, vρᵧₑ∅,
μRᵧₗ∅, μRᵧₙ∅, μRᵧₑ∅, 
σRᵧₗ∅, σRᵧₙ∅, σRᵧₑ∅ = ts_ensemble_detailed(par∅,T,ts_rep)

# these are the same as for parγλνεψ
parψ = Par(nr = nr, 
    𝔼P = 100.0,
    h² = 0.25, λ² = 0.25, ν² = 0.25, ϵ² = 0.25, 
    L = 100, SL = 100, SN = 100, SE = 100, 
    sr = sr, or = or)

μρᵧₗψ, μρᵧₙψ, μρᵧₑψ, 
σρᵧₗψ, σρᵧₙψ, σρᵧₑψ,
vρᵧₗψ, vρᵧₙψ, vρᵧₑψ,
μRᵧₗψ, μRᵧₙψ, μRᵧₑψ, 
σRᵧₗψ, σRᵧₙψ, σRᵧₑψ = ts_ensemble_detailed(parψ,T,ts_rep)

plt_vals = [
    μρᵧₑ∅ .+ 0.5 .* σρᵧₑ∅; 
    μρᵧₗ∅ .+ 0.5 .* σρᵧₗ∅; 
    μρᵧₙ∅ .+ 0.5 .* σρᵧₙ∅;
    μρᵧₑψ .+ 0.5 .* σρᵧₑψ; 
    μρᵧₗψ .+ 0.5 .* σρᵧₗψ; 
    μρᵧₙψ .+ 0.5 .* σρᵧₙψ];

ρmx = maximum(InfNan2Bye(plt_vals))

ργl∅ₚ = plot(1:T,μρᵧₗ∅,ribbon=0.5 .* σρᵧₗ∅,fillalpha=0.3,
    legend=false,
    title="Gene - Lineal",
    ylim=(0,ρmx),
    color=colorant"#ffb86c",
    showaxis=:y);
ylabel!("\n|Correlation|");

ργn∅ₚ = plot(1:T,μρᵧₙ∅,ribbon=0.5 .* σρᵧₙ∅,fillalpha=0.3,
    legend=false, 
    title="Gene - Non-Lineal",
    ylim=(0,ρmx),
    color=colorant"#ffb86c",
    showaxis=false);

ργe∅ₚ = plot(1:T,μρᵧₑ∅,ribbon=0.5 .* σρᵧₑ∅,fillalpha=0.3,
    legend=false, 
    title="Gene - Novel",
    ylim=(0,ρmx),
    color=colorant"#ffb86c",
    guide_position=:right, guidefonthalign=:left, 
    ylabel="Pre-Selection", showaxis=false);

ργlψₚ = plot(1:T,μρᵧₗψ,ribbon=0.5 .* σρᵧₗψ,fillalpha=0.3,
    legend=false,
    ylim=(0,ρmx),
    color=colorant"#ffb86c");
xlabel!("Host Generations\n");
ylabel!("\n|Correlation|");

ργnψₚ = plot(1:T,μρᵧₙψ,ribbon=0.5 .* σρᵧₙψ,fillalpha=0.3,
    legend=false, 
    ylim=(0,ρmx),
    color=colorant"#ffb86c",
    showaxis=:x);
xlabel!("Host Generations\n");

ργeψₚ = plot(1:T,μρᵧₑψ,ribbon=0.5 .* σρᵧₑψ,fillalpha=0.3,
    legend=false, 
    ylim=(0,ρmx),
    color=colorant"#ffb86c",
    showaxis=:x,
    guide_position=:right, guidefonthalign=:left, ylabel="Post-Selection");
xlabel!("Host Generations\n");


tspl = plot(ργl∅ₚ,ργn∅ₚ,ργe∅ₚ, ργlψₚ,ργnψₚ,ργeψₚ, 
                layout=(2,3), size=(650,300), dpi=400, linewidth=2)

savefig(tspl,"fig4.png")

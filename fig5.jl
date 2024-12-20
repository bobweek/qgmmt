include("fcts_strcts.jl")

theme(:vibrant)

using LaTeXStrings

# number of repetitions at each nested step

nr = 100 # number times to draw background parameters
sr = 20 # number times to draw selected parents
or = 20 # number times to form offspring

#
# simulating data using pre-selected nonlineal source
#

prpar = Par(nr = nr, 
    h² = 0.25, λ² = 0.25, ν² = 0.25, ϵ² = 0.25, 
    L = 100, SL = 100, SN = 100, SE = 100, 
    sr = sr, or = or,
    Nsel=false)

prμΔz̄, prμAᵧβ, prμAₗβ, prμAₙβ, prμAᵩβ = runs(prpar)

#
# simulating data using post-selected nonlineal source
#

popar = Par(nr = nr, 
    h² = 0.25, λ² = 0.25, ν² = 0.25, ϵ² = 0.25, 
    L = 100, SL = 100, SN = 100, SE = 100, 
    sr = sr, or = or)

poμΔz̄, poμAᵧβ, poμAₗβ, poμAₙβ, poμAᵩβ = runs(popar)

mn = floor(minimum([poμΔz̄; prμΔz̄; poμAᵧβ; poμAₗβ; poμAₙβ; poμAᵩβ; prμAᵧβ; prμAₗβ; prμAₙβ; prμAᵩβ]))
mx = ceil(maximum([poμΔz̄; prμΔz̄; poμAᵧβ; poμAₗβ; poμAₙβ; poμAᵩβ; prμAᵧβ; prμAₗβ; prμAₙβ; prμAᵩβ]))

#
# plotting simulated data
#

# pre-sel plots

prpls = []

scpl = scatter(prμΔz̄,prμAᵧβ, legend=false, title="G", mc=colorant"#ff79c6", ms=2, showaxis=:y);
ylabel!("\n\nPredicted Δz̄");
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(prpls,scpl)

scpl = scatter(prμΔz̄,prμAₗβ, legend=false, title="GL", mc=colorant"#50fa7b", ms=2, showaxis=false);
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(prpls,scpl)

scpl = scatter(prμΔz̄,prμAₙβ, legend=false, title="GLN", mc=colorant"#bd93f9", ms=2, showaxis=false);
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(prpls,scpl)

scpl = scatter(prμΔz̄,prμAᵩβ, legend=false, title="GLNE", mc=colorant"#ffb86c", ms=2, showaxis=false);
ylabel!("\nPre-Selection", guide_position=:right);
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(prpls,scpl)

# post-sel plots

popls = []

scpl = scatter(poμΔz̄,poμAᵧβ, legend=false, mc=colorant"#ff79c6", ms=2, showaxis=:y);
ylabel!("\n\nPredicted Δz̄");
xlabel!("Observed Δz̄\n");
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(popls,scpl)

scpl = scatter(poμΔz̄,poμAₗβ, legend=false, mc=colorant"#50fa7b", ms=2, showaxis=false);
xlabel!("Observed Δz̄\n");
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(popls,scpl)

scpl = scatter(poμΔz̄,poμAₙβ, legend=false, mc=colorant"#bd93f9", ms=2, showaxis=false);
xlabel!("Observed Δz̄\n");
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(popls,scpl)

scpl = scatter(poμΔz̄,poμAᵩβ, legend=false, mc=colorant"#ffb86c", ms=2, guide_position=:right, guidefonthalign=:left, ylabel="\nPost-Selection", showaxis=false);
xlabel!("Observed Δz̄\n");
plot!([mn,mx],[mn,mx], color=colorant"#8be9fd");
push!(popls,scpl)

plpl = plot(prpls[1],prpls[2],prpls[3],prpls[4],
                popls[1],popls[2],popls[3],popls[4], 
                    xlims=(mn,mx), ylims=(mn,mx),
                        layout=(2,4), size=(800,300), dpi=400)

savefig(plpl,"fig5.png")
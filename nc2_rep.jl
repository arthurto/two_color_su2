begin
    using Optim
    using Plots 
    using FiniteDifferences
    include("base.jl") 
    MeV = 1e-3 # 1 MeV em unidades de GeV
    fm = 1/(197MeV)
    toMeVfm3 = fm^3/MeV
end

begin 
    plot(x,Omega.(x,0.79,μ=0.0,μi=0.2))
end

begin 
    plot(x,Omega.(x,0.657,7.23))
end

begin
    x = 0:1MeV:0.657
    plot(x,Omega.(x,0.657,7.23,Δ = 0.0,μ=0.20,μi=0.3)*toMeVfm3)
end


begin
    function MΔ(μ,μi,Λ = Λ, G = G;g1)
        optimize(x->Omega(x[1],Λ,G,Δ=x[2],μ=μ,μi=μi),g1,method=LBFGS())
    end

    nG = 7.23
    nΛ = 657MeV
    function MΔ(μ,μi)
        sol1 = MΔ(μ,μi,nΛ,nG,g1 =[Λ/2,0.0])
        sol2 = MΔ(μ,μi,nΛ,nG,g1 = [0.0,0.0])
        if Optim.minimum(sol1) < Optim.minimum(sol2)
            return Optim.minimizer(sol1)
        end 
        return Optim.minimizer(sol2)
    end
    M0 = MΔ(0.0,0.0)[1]
    Ω0 = Omega(0.0,0.0,M0,0.0)
    function thermo_stuff(μ,μi=0.0)::Vector{Number}
        # μi = 0.0
        M,Δ = MΔ(μ,μi)
        P = -Omega(μ,μi,M,Δ)+Ω0 
        n1 = sum(n(μ,μi,M,Δ))
        chimumu = Xμμ(μ,μi,M,Δ)
        cs2 = n1/(chimumu*μ)
        ε = -P+μ*n1
        return [M,Δ,P,ε,n1,cs2]
    end
end

begin
    plot(0.0:1MeV:790MeV,x->MΔ(x,0.0)[1],label="\$\\mu_I = 0\$ MeV",
    xlabel = "\$\\mu\$ [GeV]",ylabel="\$M\$ [GeV]") 
    plot!(0.0:1MeV:790MeV,x->MΔ(x,0.1)[1],label="\$\\mu_I = 100\$ MeV") 
    plot!(0.0:1MeV:790MeV,x->MΔ(x,0.2)[1],label="\$\\mu_I = 200\$ MeV")
    # plot!(0.0:1MeV:790MeV,x->MΔ(x,0.3)[1],label="\$\\mu_I = 300\$ MeV")
    # plot!(0.0:1MeV:790MeV,x->MΔ(x,0.4)[1],label="\$\\mu_I = 400\$ MeV")
    # savefig("/home/arthurp/docs/phd/NJLSU2/ShuShengXu/nc2/Mμ.png")
end

# Até aqui está tudo funcionando muito bem. 
# Agora vamos ver como eu posso fazer as contas do que falta. 
begin
    plot(0.0:1MeV:657MeV,x->thermo_stuff(x)[1],label="\$\\mu_I = 0\$ MeV",
    xlabel = "\$\\mu\$ [GeV]",ylabel="\$M\$ [GeV]")
    plot!(0.0:1MeV:657MeV,x->thermo_stuff(x,0.1)[1],label="\$\\mu_I = 0\$ MeV")
end

begin
    plot(0.0:1MeV:657MeV,x->thermo_stuff(x)[6],label="\$\\mu_I = 0\$ MeV",
    xlabel = "\$\\mu\$ [GeV]",ylabel="\$c_s^2\$")
    plot!(0.0:1MeV:657MeV,x->thermo_stuff(x,0.05)[6],label="\$\\mu_I = 50\$ MeV")
    plot!(0.0:1MeV:657MeV,x->thermo_stuff(x,0.1)[6],label="\$\\mu_I = 100\$ MeV")
    plot!(0.0:1MeV:657MeV,x->thermo_stuff(x,0.2)[6],label="\$\\mu_I = 200\$ MeV")
    savefig("cs2_μ.png")
end

# agora plotando como função da densidade
begin
    range = 0.0:1MeV:450MeV
    plotting = (i,μi=0.0) -> (x->thermo_stuff(x,μi)[i]).(range)

    # A desindade vai estar em GeV^4 
    # precisamos converter para MeV/fm^3
    plot(plotting(5)*toMeVfm3,plotting(6),label="\$\\mu_I = 0\$ MeV",
    xlabel = "\$n_q\$ [MeV/fm\$^3\$]",ylabel="\$c_s^2\$",xlims=(0,1500))
    plot!(plotting(5,0.05)*toMeVfm3,plotting(6,0.05),label="\$\\mu_I = 50\$ MeV")
    plot!(plotting(5,0.1)*toMeVfm3,plotting(6,0.1),label="\$\\mu_I = 100\$ MeV")
    plot!(plotting(5,0.2)*toMeVfm3,plotting(6,0.2),label="\$\\mu_I = 200\$ MeV")
    savefig("cs2_n.png")
end

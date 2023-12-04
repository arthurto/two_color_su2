# This program will solve the NJL ShuSheng
# minimizing the thermodynamical potential. 

begin
    include("base.jl")
    using Optim
    using Plots
    using ForwardDiff
	using NLsolve

    MeV = 1e-3 # MeV in GeV
    fm = 1/(197.3MeV)

    function solveΩ(μ,μi,M = 0.1,Δ = 0.1,T=nothing)
        initial_x = [M,Δ]
        if isnothing(T)
            res = optimize(x-> Omega(μ,μi,x[1],x[2]),initial_x,BFGS())
        else
            res = optimize(x-> OmegaT(μ,μi,x[1],x[2],T),initial_x,BFGS())
        end
        return [Optim.minimizer(res),Optim.minimum(res)]
    end

    function multisolveΩ(μ,μi,T=nothing)

        guesses = [
            [0.2,0.0],
            [0.4,0.4],
            [0.2,0.3],
            [0.25,0.0]
            ]

        dummy_sol = [[0.0,0.0],Inf]

        for guess in guesses
            if isnothing(T)
                sol = solveΩ(μ,μi,guess[1],guess[2])
            else
                sol = solveΩ(μ,μi,guess[1],guess[2],T)
            end

            if sol[2] < dummy_sol[2]
                dummy_sol = sol
            end
        end

        return dummy_sol
    end
end 

# multisolveΩ(0.1,0.1)
# solveΩ(0.1,0.1)

# Criar uma forma de resolver não linear as equações de 
# gap como derivada do omega

# solveΩ(0.0,0.1)

# begin
#     x = LinRange(0.0,0.5,50)
#     y = LinRange(0.,0.5,50)
#     # μ, μi, M, Δ
#     μ = 0.34
#     μi = 0.0*0.135
#     # plot(x,Omega.(0.0,0.2,x,0.0))
#     f = (x,y) -> Omega(μ,μi,x,y)
#     z = @. f(x',y)
#     contour(x,y,z,levels = 1000)
#     s = solveΩ(μ,μi,0.2,0.0)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩ(μ,μi,0.4,0.0)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩ(μ,μi,0.4,0.4)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩ(μ,μi,0.2,0.3)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩ(μ,μi,0.25,0.0)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
# end

# begin
#     x = LinRange(0.0,0.5,50)
#     y = LinRange(0.,0.5,50)
#     # μ, μi, M, Δ
#     μ = 0.3
#     μi = 1.3*0.135
#     # plot(x,Omega.(0.0,0.2,x,0.0))
#     f = (x,y) -> OmegaT(μ,μi,x,y)
#     z = @. f(x',y)
#     contour(x,y,z,levels = 1000)
#     s = solveΩT(μ,μi,0.2,0.0)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩT(μ,μi,0.4,0.0)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩT(μ,μi,0.4,0.4)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩT(μ,μi,0.2,0.3)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
#     s = solveΩT(μ,μi,0.25,0.0)
#     scatter!([s[1][1]],[s[1][2]],label = s[2])
# end

# begin
#     x = LinRange(0.0,2.0,100)
#     μ = 0.3
#     plot(x,(x->multisolveΩ(μ,x*0.135)[1][2]).(x),label="Δ")
#     plot!(x,(x->multisolveΩ(μ,x*0.135)[1][1]).(x),label="M")
# end

# begin
#     x = LinRange(0.0,0.5,100)
#     # x = LinRange(0.37,0.4,100)
#     μ = 1.5
#     T = 1e-3
#     # scatter(x,(x->multisolveΩ(x,μ*0.135,T)[1][2]).(x),label="Δ")
#     # scatter!(x,(x->multisolveΩ(x,μ*0.135,T)[1][1]).(x),label="M")
#     scatter(x,(x->multisolveΩ(x,μ*0.135)[1][2]).(x),label="Δ",
#     title = "\$ \\mu_I/m_\\pi = $μ\$",
#     ylabel = "[GeV]",xlabel = "\$\\mu_q\$ [GeV]")
#     scatter!(x,(x->multisolveΩ(x,μ*0.135)[1][1]).(x),label="M")
# end

# multisolveΩ(0.35,0.0)

# This part of the code is responsible for the 
# beta equilíbrium and charge neutrality 

begin
    # The Electron number density, pressure and energy density are defined
    function n_e(μ)
        return μ^3/(3*π^2)
    end

    function p_e(μ)
        return μ^4/(4*π^2)
    end

    function ε_e(μ)
        return μ^4/(12*π^2)
    end

    # Defining the sistem of equations to solve
    function eq_to_solve(μ,μi) 
        # Solving these equations as a function of μ = μ_q

        # Step 1: Use μ and μi to find 
        # M and Δ
        M,Δ = multisolveΩ(μ,μi)[1]

        # Step 2: Use μ, μi, M and Δ to find the equation for beta equilibrium.
        # Solving μi(μ)

        nd,nu = n(μ,μi/2,M,Δ)
        ne = n_e(-μi/2) # Here we impose that μ_e = -μ_i

        # The equation for eletrical neutrality:
        eqq = 2*nu/3 - nd/3 - ne 
    
        return eqq
    end

    function μi_μ(μ)
        return nlsolve(x -> eq_to_solve(μ,x[1]),[0.0]).zero[1]/2
    end

    # This functions returns the thermodynamical 
    # quantities like Pressure, Energy Density,
    # Number density
    #
    function thermo_stuff(μ)

        # This part solves for the variables
        #
        μi = μi_μ(μ)
        M,Δ = multisolveΩ(μ,μi)[1]

        # This part calculates the thermodynamics
        # 
        nd,nu = n(μ,μi/2,M,Δ)
        ne = n_e(-μi/2) # Here we impose that μ_e = -μ_i
        nb = (nd + nu)/3
        Pq = -Omega(μ,μi,M,Δ)-0.029477621219439525 
        P = Pq
        ε = -Pq +nu*(μ+μi/2) + nd*(μ-μi/2) + ε_e(-μi/2)

        return [μi,M,Δ,nu,nd,ne,nb,P,ε]
    end
end



begin 
    x = 0.0MeV:1MeV:Λ
    data = hcat(thermo_stuff.(x)...)
end

data[2,:]
begin 
    plot(data[9,:],data[8,:].-data[8,1])
end
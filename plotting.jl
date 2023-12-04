begin 
    include("/home/arthurp/docs/phd/NJLSU2/ShuShengXu/nc2/solving.jl")
    # include("../../../general_use/tov/tov_calc.jl")
    using Plots
    using DelimitedFiles
    using Interpolations
end

# cs2b(0.01,0.01)

begin
    x = LinRange(0.0,3.0,200)
    μ = 0.0
    plot(x,(x->multisolveΩ(μ,x*0.135)[1][2]).(x),
    xlims = (0.0,3.0),ylims=(0.0,0.33),
label = "\$M\$ [GeV]",xlabel = "\$μ_I /m_π\$",ylabel="[GeV]")
    plot!(x,(x->multisolveΩ(μ,x*0.135)[1][1]).(x),
    label = "\$ Δ \$[GeV]",legend=:right)
    # savefig("/home/arthurp/docs/phd/NJLSU2/ShuShengXu/nc2/MxDelta.png")
end

begin
    x1 = LinRange(0.25,0.3315,100)
    x2 = LinRange(0.3325,0.5,100)
    μi = 0.0
    # plot(x,(x->multisolveΩ(x,μi*0.135)[1][2]).(x))
    scatter(x1,(x->multisolveΩ(x,μi*0.135)[1][1]).(x1),label=false,c=1,xlims=(0.25,0.5),ylims=(0.0,0.33),
    xlabel = "\$M\$ [GeV]",
    ylabel = "\$μ_q\$ [GeV]",
    title = "\$μ_I = 0\$ ")
    scatter!(x2,(x->multisolveΩ(x,μi*0.135)[1][1]).(x2),c=1,label = "\$M\$ [GeV]")
    plot!([0.332,0.332],[0.0,0.33],c=:black,lw=2,label = false,line=:dash)
    # savefig("effective_massμ.png")
end

begin
    x = LinRange(0.0,2.0,50)
    y = LinRange(0.0,0.4,300)

    f = (x,y) -> multisolveΩ(y,-x*0.135)[1][1]
    z = @. f(x',y)
    contourf(x,y,z,levels = 100,
            xlabel = "\$μ_I/m_π\$",
            ylabel = "\$ μ_q\$[GeV]",
            linewidth=0.0,
            ms = 6,
            title = "Effective Mass - \$M\$ [GeV]",
            fmt = :png,
            c=:rainbow)
    # savefig("effective_mass.png")
end

begin
    x = LinRange(0.0,2.0,100)
    y = LinRange(0.0,0.4,100)

    f = (x,y) -> multisolveΩ(y,x*0.135)[1][2]
    z = @. f(x',y)
    contourf(x,y,z,levels = 100,xlabel = "\$μ_I\$")
end

# Now plotting the μ_I as a function of μ_q
# 
begin 
    scatter(LinRange(0.3,0.5,100),x->μi_μ(x),
    xlabel = "\$μ_q\$ [GeV]",
    ylabel = "\$μ_I\$ [GeV]",label = nothing)
end

# Agora eu vou resolver a pressão como função de μ_q 
#
begin
   x = LinRange(0.3,0.5,100)
   #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
   #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]
   scatter(x,x->thermo_stuff(x)[4],label = "\$μ_u\$",
   xlabel = "\$μ_q\$ [GeV]",ylabel = "\$ [GeV]\$") 
   scatter!(x,x->thermo_stuff(x)[5],label = "\$μ_d\$") 
end

begin 
    x = LinRange(0.3,0.4,200)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    P0 = (((1e3/197)^3)/346.255)*1e3
    P = (x->thermo_stuff(x)[8]*P0).(x)
    ε = (x->thermo_stuff(x)[9]*P0).(x)
    scatter(P,ε,xlims=(0.0,0.13),ylims=(0.0,1.3)) 
end

begin 
    x = LinRange(0.3,0.9,200)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    P0 = (((1.0e3/197.0)^3.0)/346.255)*1.0e3
    P = (x->thermo_stuff(x)[8]*P0).(x)
    ε = (x->thermo_stuff(x)[9]*P0).(x)
    scatter(P,ε,xlims=(0.0,2.0),ylims=(0.0,8.0)) 
end

begin 
    x = LinRange(0.3,0.6,100)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    P0 = (((1.0e3/197.0)^3.0)/346.255)*1.0e3
    # P = (x->thermo_stuff(x)[8]*P0).(x)
    ε = (x->thermo_stuff(x)[9]*P0).(x)
    μi = (x->thermo_stuff(x)[1]).(x)

    scatter(ε,x,xlims=(0.0,6.0),ylims=(0.0,0.6),
    label = "\$μ_q\$",xlabel = "\$ε/P_0\$",ylabel = "[GeV]")
    scatter!(ε,μi,label = "\$μ_I\$")
end

begin
    x = LinRange(0.0,2.0,25)
    y = LinRange(0.0,0.6,150)

    f = (x,y) -> multisolveΩ(y,-x*0.135)[1][1]
    z = @. f(x',y)
    contourf(x,y,z,levels = 100,
            xlabel = "\$μ_I/m_π\$",
            ylabel = "\$ μ_q\$[GeV]",
            linewidth=0.0,
            ms = 6,
            title = "Effective Mass - \$M\$ [GeV]",
            fmt = :png,
            c=:rainbow)

    # Plotting the path for 
    # the equilibrium solution
    x = LinRange(0.3,0.6,100)
            # thermo_stuff
            #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
            #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]
        
    μi = -(x->thermo_stuff(x)[1]).(x)
    scatter!(μi/0.135,x,label = "EoS path")
    # savefig("effective_mass.png")
end


# EOS 
begin 
    x = LinRange(0.3,0.66,300)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    P0 = (1.0e12)/(197.0^3.0) # convertendo GeV^4 -> MeV/fm^3
    P = (x->thermo_stuff(x)[8]*P0).(x)
    ε = (x->thermo_stuff(x)[9]*P0).(x)
    # plot(P,ε,xlims=(-0.1,1.0),ylims=(0.0,3.0))
    
    # Saving the eos
    open("eos.dat", "w") do io
        writedlm(io, [P ε])
    end

    plot(ε,P)
end

begin 
    x = LinRange(0.3,0.6,100)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    P0 = (1.0e12)/(197.0^3.0) # convertendo GeV^4 -> MeV/fm^3
    P = (x->thermo_stuff(x)[8]*P0).(x)
    ε = (x->thermo_stuff(x)[9]*P0).(x)
    # plot(P,ε,xlims=(-0.1,1.0),ylims=(0.0,3.0))
    plot(ε,P)
end

begin 
    x = LinRange(0.3,0.4,100)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    # P0 = (1.0e12)/(197.0^3.0) # convertendo GeV^4 -> MeV/fm^3
    # P = (x->thermo_stuff(x)[8]*P0).(x)
    nb = (x->thermo_stuff(x)[7]*(1e9/(197.0^3))).(x)
    EA = (x->thermo_stuff(x)[9]*1e3/thermo_stuff(x)[7]).(x)
    # plot(P,ε,xlims=(-0.1,1.0),ylims=(0.0,3.0))
    scatter(nb/0.15,EA,ylims = (900,1000))
end

begin 
    N = 200
    x = LinRange(0.3,0.66,N)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    matrix = Matrix(undef,N+1,10)

    matrix[1,:] = ["μ","μi","M","Δ","nu","nd","ne","nb","P","ε"]
    for i ∈ 1:N
        matrix[i+1,:] = append!([x[i]],thermo_stuff(x[i]))
    end
    # plot(P,ε,xlims=(-0.1,1.0),ylims=(0.0,3.0))
    
    # Saving the eos
    # open("research/hotqcdmui/Ours/ShuShengXu/eos_full.csv", "w") do io
    open("eos_full.csv", "w") do io
        writedlm(io, matrix,',')
    end
end

begin 
    x = LinRange(0.3,0.4,100)
    # thermo_stuff
    #[μi,M,Δ,nu,nd,ne,nb,P,ε ]
    #[1 ,2,3,4 ,5 ,6 ,7 ,8,9]

    # P0 = (1.0e12)/(197.0^3.0) # convertendo GeV^4 -> MeV/fm^3
    # P = (x->thermo_stuff(x)[8]*P0).(x)
    nb = (x->thermo_stuff(x)[7]*(1e9/(197.0^3))).(x)
    nu = (x->thermo_stuff(x)[4]/thermo_stuff(x)[7]).(x)
    nd = (x->thermo_stuff(x)[5]/thermo_stuff(x)[7]).(x)
    # plot(P,ε,xlims=(-0.1,1.0),ylims=(0.0,3.0))
    plot(nb/0.15,nu/3,ylims=(0,1))
    plot!(nb/0.15,nd/3)
end

begin
    x = LinRange(0.32,0.34,200)
    μ = 0.0
    plot(x,(x->multisolveΩ(x,μ*0.135)[1][1]).(x),
label = "\$μ_I = 0\$",xlabel = "\$μ [GeV]\$",ylabel="M [GeV]")
    plot!(x,(x->thermo_stuff(x)[2]).(x),
    label = "\$μ_I \\neq 0\$",legend=:right)
end

# This block will solve the TOV equations
# for a single star with center density
begin 
    # Creating the EOS 

    # Range of μ_q 
    x = LinRange(0.3,0.6,100)

    # Converting factors
    conv = (1e12/(197^3))# GeV^4 -> MeV/fm³

    # Bag constant 
    B = 1 # MeV/fm³

    P = ((x->thermo_stuff(x)[8]*conv - B).(x))*nutogu
    ε = (x->thermo_stuff(x)[9]*conv + B).(x)*nutogu
    interp = linear_interpolation(P,ε)
    function eos(x)
        if x < P[1] return 0.0 end
        return interp(x)
    end

    M,R = MR_diag(P[15:5:end],EOS=eos)
    # M,R = MR_diag()
    plot(R/1000,M,xlims=(0,10),ylims=(0,2))
end
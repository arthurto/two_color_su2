# This file will keep the basic functions
# to build the thermodynamical potential
# of the NJL SU(2) with pion condensation 
# in T = 0. 
begin
    using QuadGK
    using NLsolve
    using ForwardDiff
    # Under development...
    # Missing Features:


    # Constantes
    # Shu-Sheng Xu - Nuclear Physics B 971 (2021) 115540
    m_c = 0.000
    Λ = 0.657
    G = 7.23
    Nc = 2.0

    # Constantes (nossas)
    # m_c = 4.93561e-3
    # Λ = 659.325e-3
    # G = 2.07691/(Λ^2)
    # Nc = 3.0

    # Parâmetros para a resolução numérica
    Integration_max_evals = 16
end 


begin
    function Ek(k,M)
        return sqrt(k^2 + M^2)
    end

    function Ekp(k = 0.0,
        M = 0.0,
        μi = 0.0,
        Δ = 0.0)
        # return sqrt((Ek(k,M) + μi)^2+(2*G*Δ)^2)
        return sqrt((Ek(k,M) + μi)^2+Δ^2)
    end
end

begin 
    function IOmega_TRS(M = 0.0,
        μi = 0.0,
        Δ = 0.0)
        return quadgk(k->k*k*(Ekp(k,M,μi,Δ)+Ekp(k,M,-μi,Δ)),0.0,Λ,maxevals = Integration_max_evals)[1]/(2.0*π*π)
    end

    function IOmega_TRS(M = 0.0,
        μi = 0.0,
        Δ = 0.0,Λ=Λ)
        return quadgk(k->k*k*(Ekp(k,M,μi,Δ)+Ekp(k,M,-μi,Δ)),0.0,Λ,maxevals = Integration_max_evals)[1]/(2.0*π*π)
    end
end 

begin
    function fermi_momentum(μ,M)
        if abs2(M) < abs2(μ)
            return sqrt(μ^2 - M^2)
        end
        return 0.0 
    end

    function pff(M,μ,μi,Δ)
        # if μ < Ekp(k,M,μi,Δ) return 0.0 end # dumbass
        Ek0 =  μi + fermi_momentum(μ,Δ)
        return fermi_momentum(Ek0,M)
    end
end 

begin
    function tl(k,M,μ,μi,Δ)
        a = -(Ekp(k,M,μi,Δ)-μ)
        a = a*(sign(a)+1)/2
        return a
    end

    function kf_solve(M,μ,μi,Δ)
        # Encontra o momento de fermi numericamente 
        # por motivos de preguiça
        if μ < Ekp(0.0,M,μi,Δ)
            return 0.0
        end

        sol = nlsolve(k->μ - Ekp(k[1],M,μi,Δ),[μ])
        return sol.zero[1]
    end

    # kf_solve(0.1,0.3,0.1,0.0)

    function IOmega_mui(M,μ,μi,Δ)

        # The fermi momentum from the up quark
        pf = pff(M,μ,-μi,Δ)
        ip = quadgk(k->k^2*(tl(k,M,μ,μi,Δ)),0.0,pf,maxevals= Integration_max_evals)[1]
        
        # The fermi momentum from the down quark
        pf = pff(M,μ,μi,Δ)
        ip += quadgk(k->k^2*(tl(k,M,μ,-μi,Δ)),0.0,pf,maxevals= Integration_max_evals)[1]

        return ip/(2*(π^2))
    end

    function Omega(μ,μi,M,Δ)

        σ = M - m_c

        Ω = (σ^2 + Δ^2)/(4*G) - 2*Nc*(IOmega_TRS(M,μi/2,Δ)+IOmega_mui(M,μ,μi/2,Δ))

        return Ω
    end

    function Omega(M,Λ,G=G;μ=0.0,μi=0.0,Δ = 0.0)
        σ = M - m_c

        Ω = (σ^2 + Δ^2)/(4*G) - 2*Nc*(IOmega_TRS(M,μi/2,Δ,Λ)+IOmega_mui(M,μ,μi/2,Δ))

        return Ω
    end
    function Omega(M,Λ,G)
        μ = 0.0
        μi = 0.0 
        Δ = 0.0
        σ = M - m_c

        Ω = (σ^2 + Δ^2)/(4*G) - 2*Nc*(IOmega_TRS(M,μi/2,Δ,Λ)+IOmega_mui(M,μ,μi/2,Δ))

        return Ω
    end

    function OmegaT(μ,μi,M,Δ,T=1e-3)

        σ = M - m_c

        Ω = (σ^2 + Δ^2)/(4*G) - 2*Nc*(IOmega_TRS(M,μi/2,Δ)+IOmega_muiT(M,μ,μi/2,Δ,T))

        return Ω
    end

    function n(μ,μi,M,Δ) # The densities are defined is therms of this
    
        pf = pff(M,μ,μi,Δ)
        ip = pf^3/3
        pf = pff(M,μ,-μi,Δ)
        im = pf^3/3

        return [ 
            Nc*(im)/(π^2), # n_u 
            Nc*(ip)/(π^2)  # n_d 
            ] 
    end

    function subXμμ(μ,μi,M,Δ)
        twoGΔ = Δ
        sq1 = fermi_momentum(μ,twoGΔ)
        # Return zero if fermimomentum is zero 
        Ek0 = μi + sq1
        pf = fermi_momentum(Ek0,M) # Remember to always return to the real values
        if sq1 == 0.0 return Nc*pf*μ/(π^2) end 
        return Nc*pf*μ*(1+μi/sq1)/(π^2)
    end

    function Xμμ(μ,μi,M,Δ)
        # δ = 1e-5
        # np = sum(n(μ+δ,μi,M,Δ))
        # nm = sum(n(μ-δ,μi,M,Δ))
        # return (np-nm)/(2*δ)
        return subXμμ(μ,μi,M,Δ)+subXμμ(μ,-μi,M,Δ)
    end

    function epsilon(μ,μi,M,Δ)
        Pq = -Omega(μ,μi,M,Δ)-0.029477621219439525
        nu,nd = n(μ,μi/2,M,Δ)
        ε = -Pq +nu*(μ+μi/2) + nd*(μ-μi/2)
        return ε
    end
end
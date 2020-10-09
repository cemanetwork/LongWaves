# NONLINEAR 1-D LONG WAVE MODEL FOR PROPAGATION OVER CONSTANT DEPTH
#CEMAN: COASTAL ECOSYSTEM MANAGEMENT NETWORK
# BY: GERMAN RIVILLAS-OSPINA
# DATE: JUNIO  DE 2020
#LUGAR: BARRANQUILLA-COLOMBIA

using Plots

# CONSTANT VALUES
const g=9.81
const T=100.0
const dt=0.5
const dx=100.0
const Hₒ=50.0

function LW()

    # INITIAL CONDITIONS
    C=(g*Hₒ)^0.5
    L=C*T

    # DEFINING THE DOMAIN
    xmax=0:20.0:L
    xin=length(xmax)
    tin=round.(Int,L)

    #DEFINING VALUES TO PHYSICAL VARIABLES
    η=zeros(1,xin)
    HₒA=fill(Hₒ, (1,xin))
    H=fill(Hₒ, (1,xin))
    Hₙ=zeros(1,xin)
    H₁=zeros(1,xin)
    H₂=zeros(1,xin)
    Hᵣ=zeros(1,xin)
    U=zeros(1,xin)
    U₁=zeros(1,xin)
    U₂=zeros(1,xin)
    Uₙ=zeros(1,xin)

    η₀=1.5
    dtmax=round.(Int,tin/2)
    N=0
    η₁=0.0
    η₂=0.0
    pi=4*atan(1.)
println(xmax)
    while N < dtmax

        #TEMPORAL INCREMENTS
        N=N+1

        #COURANT NUMBER
        TERM1=(N-1)*dt
        TERM2=dx/C

        # FREE SURFACE ELEVATION CALCULATION
        if TERM1<TERM2

            D=sin(2*pi*(N)*(dt/T))
            η[1]=η₀*D+η₁

        else

            A=sin(2*pi*(N-1)*dt/T)
            B=sin((2*pi*(N-1)*dt)/T-(dx/L))
            η₁=η[1]-η₀*A
            η₂=η[2]-η₀*B

            η₁=η₁+(dt/dx)*C*(η₂-η₁)
            D=sin(2*pi*(N)*(dt/T))
            η[1]=η₀*D+η₁

        end

        for i = 2:xin-1
                # DETERMINING WATER DEPTH VARIATION
                Hₙ[1]=HₒA[1]+η[1]

                H₁[i]=(H[i+1]+H[i])*U[i+1]

                H₂[i]=(H[i]+H[i-1])*U[i]

                Hᵣ[i]=H₁[i]-H₂[i]

                Hₙ[i]=H[i]-(dt/(2*dx))*Hᵣ[i]

                # WATER VELOCITIES FIELD CALCULATION CENTRAL SCHEME
                η[i]=Hₙ[i]-HₒA[i]

                U₁[i]=(U[i+1]+U[i])^2-(U[i]+U[i-1])^2
                U₂[i]=g*(η[i]-η[i-1])

                Uₙ[i]=U[i]-(dt/(8dx))U₁[i]-(dt/dx)U₂[i]
                #println(Uₙ[i],"   ", η[i])
        end

            Uₙ[1]=Uₙ[2]
            Uₙ[xin]=Uₙ[xin-1]*(g/Hₙ[xin-1])^0.5

        for i = 1:xin
            H[i]=Hₙ[i]
            U[i]=Uₙ[i]

        end
    end
    #CONVERT RANGE TO VECTORS
    # RETURN U AND η AS VECTORS, NOT 2D ARRAYS
    return collect(xmax), U[1,:], η[1,:]

end
# DEFINING THE RESULTS IN VECTORS
x, u, η = LW()

# PLOTTING RESULTS

    plot(
        x,
        [η u],
         title = ["Surface Elevation" "Velocity field"],
         layout = (2, 1),
         color = ["blue" "green"],
         legend = false
         )
    savefig("G:\\My Drive\\Uninorte\\Modelacion\\Julia\\LongWaves\\result.png")



##_____________________

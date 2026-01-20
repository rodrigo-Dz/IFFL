# --- Main script for IFFL analysis ---
# --- Rodrigo Aguilar


# --- Required libraries ---
using DifferentialEquations, Plots, Printf, OrderedCollections, LinearAlgebra, LaTeXStrings, Statistics

# --- Model ---
M = "Phosporilation"  # Options: "Sequestration", "Enzimatic", "Phosporilation"
mm = include(string("Library/$M.jl"))



#######   Frecuency    #########################################

n = 10 
f_range = 10.0 .^ range(-3, 2.0, length=n) 

m1_focus = 1  # Elegir un m1 

fig3d = plot(
    xlabel=L"Time",
    ylabel=L"log_{10}(Frecuency)",
    zlabel=L"[U_1]",
    legend=false,
    size=(1400, 900),
    camera=(40, 35),
    grid=true,
    foreground_color_legend=nothing
)

p_current = deepcopy(mm.p)

# seq
#p_current[:m1] = m1_focus
# enz
# phos
p_current[:k1] = m1_focus

color_grad = cgrad(:viridis, n, categorical=true)
all_times_min = Float64[]
all_times_max = Float64[]
all_U1_max = Float64[]


f_range_sorted = sort(f_range, rev=true)

for (f_idx, f_val) in enumerate(f_range_sorted)    
    p_current[:f] = f_val

    # seq
    p_values = [
        p_current[:m1], p_current[:m2], p_current[:mx],
        p_current[:g], p_current[:d], p_current[:g0],
        p_current[:f], p_current[:A]
    ]

    # enz
    #=p_values = [
        p_current[:m1], p_current[:m2], p_current[:mx],
        p_current[:g], p_current[:d], p_current[:r],
        p_current[:f], p_current[:A]
    ]=#

    # phos
    #=p_values = [
        p_current[:k1], p_current[:k2], p_current[:mU],
        p_current[:mX], p_current[:d], p_current[:f],
        p_current[:A]
    ]=#
    
    prob = ODEProblem(mm.IFFL, mm.u0, tspan, p_values)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6, saveat=0.05)
    
    times = sol.t
    U1_concentration = sol[1, :]
    
    mask = times .> 50.0
    times_steady = times[mask]
    U1_steady = U1_concentration[mask]
    
    time_window = 500

    mask2 = times_steady .<= (minimum(times_steady) + time_window)
    times_plot = times_steady[mask2]
    U1_plot = U1_steady[mask2]
    freq_y = fill(log10(f_val), length(times_plot))
    
    push!(all_times_min, minimum(times_plot))
    push!(all_times_max, maximum(times_plot))
    push!(all_U1_max, maximum(U1_plot))
    
    alpha_val = 0.5 + 0.5 * (f_idx / n)  # Más opaco al frente
    
    plot!(
        fig3d,
        times_plot,
        freq_y,
        U1_plot,
        color=color_grad[f_idx],
        linewidth=2.0,
        alpha=alpha_val
    )
end


p_current[:f] = 0
# seq
p_values = [
    p_current[:m1], p_current[:m2], p_current[:mx],
    p_current[:g], p_current[:d], p_current[:g0],
    p_current[:f], p_current[:A]
]

# enz
#=p_values = [
    p_current[:m1], p_current[:m2], p_current[:mx],
    p_current[:g], p_current[:d], p_current[:r],
    p_current[:f], p_current[:A]
]=#

# phos
#=p_values = [
    p_current[:k1], p_current[:k2], p_current[:mU],
    p_current[:mX], p_current[:d], p_current[:f],
    p_current[:A]
]=#

prob = ODEProblem(mm.IFFL, mm.u0, tspan, p_values)
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6, saveat=0.05)

times = sol.t
U1_concentration = sol[1, :]

# Estado estacionario (después de 50s)
mask = times .> 50.0
times_steady = times[mask]
U1_steady = U1_concentration[mask]


time_window = 500

mask2 = times_steady .<= (minimum(times_steady) + time_window)
times_plot = times_steady[mask2]
U1_plot = U1_steady[mask2]
freq_y = fill(log10(1000), length(times_plot))

push!(all_times_min, minimum(times_plot))
push!(all_times_max, maximum(times_plot))
push!(all_U1_max, maximum(U1_plot))

# Graficar con transparencia ajustada
alpha_val = 0.5 

plot!(
    fig3d,
    times_plot,
    freq_y,
    U1_plot,
    color=color_grad[1],
    linewidth=2.0,
    alpha=alpha_val
)


# Planos de referencia
t_min, t_max = minimum(all_times_min), maximum(all_times_max)
f_min, f_max = extrema(log10.(f_range))
u_min, u_max = 0.0, maximum(all_U1_max) * 1.1

plot!(fig3d,
    [t_min, t_max, t_max, t_min, t_min],
    [f_min, f_min, f_max, f_max, f_min],
    [u_min, u_min, u_min, u_min, u_min],
    color=:gray, alpha=0.05, linewidth=0.5)

save_path_3d = "Sequestration_FrecuencyRespose_M$(m1_focus)_$M.png"
savefig(fig3d, save_path_3d)



##########    PARTE 2: Análisis de Respuesta en distintas M1

# Crear figura
fig_bode_magnitude = plot(
    xlabel=L"log_{10}(Frecuencia [Hz])",
    ylabel=L"[dB]",
    legend=:topright,
    size=(1200, 700),
    grid=true,
    xlim=(log10(minimum(f_range))-0.5, log10(maximum(f_range))+0.5)
)

m1_range = [ 0.5, 1, 3]
colors_m1 = [:blue, :red, :green]

for (m1_idx, m1_val) in enumerate(m1_range)
    p_current = deepcopy(mm.p)
    p_current[:m1] = m1_val
    
    println("  Analizando m1 = $m1_val...")
    
    amplitudes = Float64[]
    dB_values = Float64[]
    
    for (f_idx, f_val) in enumerate(f_range)
        p_current[:f] = f_val
        # seq
        p_values = [
            p_current[:m1], p_current[:m2], p_current[:mx],
            p_current[:g], p_current[:d], p_current[:g0],
            p_current[:f], p_current[:A]
        ]

        # enz
        #=p_values = [
            p_current[:m1], p_current[:m2], p_current[:mx],
            p_current[:g], p_current[:d], p_current[:r],
            p_current[:f], p_current[:A]
        ]=#

        # phos
        #=p_values = [
            p_current[:k1], p_current[:k2], p_current[:mU],
            p_current[:mX], p_current[:d], p_current[:f],
            p_current[:A]
        ]=#
        
        prob = ODEProblem(mm.IFFL, mm.u0, tspan, p_values)
        sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6, saveat=0.005)
        
        times = sol.t

        mask = times .> 50.0
        A_signal = (maximum(sol[1, mask]) - minimum(sol[1, mask]))/2
        A_ref = p_current[:A]
        dB = 20*log10(A_signal/A_ref)

        push!(dB_values, dB)
        
        println(" dB=$(round(dB, digits=2))")
    end
    
    plot!(
        fig_bode_magnitude,
        log10.(f_range),
        dB_values,
        label="m1 = $m1_val",
        color=colors_m1[m1_idx],
        linewidth=2.5,
        alpha=0.8
    )
    
end  


save_path_bode_mag = "Sequestration_FrequencyRespose_$M.png"
savefig(fig_bode_magnitude, save_path_bode_mag)






############ Dose Response 3D

n = 25
A_range = 10.0 .^ range(-3, 2.0, length=n)  

m1_focus = 1 
f_focus = 0.001 # Elegir una frecuencia fija para dosis-respuesta

fig_dose_3D = plot(
    xlabel=L"Time",
    ylabel=L"log_{10}(Amplitude)",
    zlabel=L"[U_1]",
    legend=false,
    size=(1400, 900),
    camera=(40, 35),
    grid=true,
    foreground_color_legend=nothing
)

p_current = deepcopy(mm.p)
p_current[:m1] = m1_focus
p_current[:f] = f_focus

color_grad = cgrad(:viridis, n, categorical=true)
all_times_min = Float64[]
all_times_max = Float64[]
all_U1_max = Float64[]

# Ordenar frecuencias de mayor a menor para dibujar correctamente
A_range_sorted = sort(A_range, rev=true)

for (A_idx, A_val) in enumerate(A_range_sorted)
    
    p_current[:A] = A_val
    # seq
    p_values = [
        p_current[:m1], p_current[:m2], p_current[:mx],
        p_current[:g], p_current[:d], p_current[:g0],
        p_current[:f], p_current[:A]
    ]

    # enz
    #=p_values = [
        p_current[:m1], p_current[:m2], p_current[:mx],
        p_current[:g], p_current[:d], p_current[:r],
        p_current[:f], p_current[:A]
    ]=#

    # phos
    #=p_values = [
        p_current[:k1], p_current[:k2], p_current[:mU],
        p_current[:mX], p_current[:d], p_current[:f],
        p_current[:A]
    ]=#
    
    prob = ODEProblem(mm.IFFL, mm.u0, tspan, p_values)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6, saveat=0.05)
    
    times = sol.t
    U1_concentration = sol[1, :]
    
    # Estado estacionario (después de 50s)
    mask = times .> 0.0
    times_steady = times[mask]
    U1_steady = U1_concentration[mask]
    
    # Segmento de tiempo para visualización
    if A_val > 0
        period = 1.0 / A_val
        time_window = min(30.0, 8 * period)  # Mostrar hasta 8 periodos
    else
        time_window = 30.0
    end
    
    time_window = 500

    mask2 = times_steady .<= (minimum(times_steady) + time_window)
    times_plot = times_steady[mask2]
    U1_plot = U1_steady[mask2]
    A_y = fill(log10(A_val), length(times_plot))
    
    push!(all_times_min, minimum(times_plot))
    push!(all_times_max, maximum(times_plot))
    push!(all_U1_max, maximum(U1_plot))
    
    # Graficar con transparencia ajustada
    alpha_val = 0.5 + 0.5 * (A_idx / n) 
    
    plot!(
        fig_dose_3D,
        times_plot,
        A_y,
        log10.(U1_plot),
        color=color_grad[A_idx],
        linewidth=2.0,
        alpha=alpha_val
    )
end


# Planos de referencia
t_min, t_max = minimum(all_times_min), maximum(all_times_max)
f_min, f_max = extrema(log10.(A_range))
u_min, u_max = 0.0, maximum(all_U1_max) * 1.1

plot!(fig_dose_3D,
    [t_min, t_max, t_max, t_min, t_min],
    [f_min, f_min, f_max, f_max, f_min],
    [u_min, u_min, u_min, u_min, u_min],
    color=:gray, alpha=0.05, linewidth=0.5)

save_path_3d = "Sequestration_DoseRespose_M$(m1_focus)_F$(f_focus)_$M.png"
savefig(fig_dose_3D, save_path_3d)
println("  ✓ Gráfica 3D guardada: $save_path_3d")



######### Análisis de Dosis-Respuesta

fig_dose = plot(
    xlabel=L"log_{10}(Amplitud)",
    ylabel=L"[dB]",
    legend=:topright,
    size=(1200, 700),
    grid=true,
    xlim=(log10(minimum(A_range))-0.5, log10(maximum(A_range))+0.5)
)


colors_m1 = [:blue, :red, :green]

for (m1_idx, m1_val) in enumerate(m1_range)
    p_current = deepcopy(mm.p)
    p_current[:m1] = m1_val
    p_current[:f] = f_focus

    println("  Analizando m1 = $m1_val...")
    
    amplitudes = Float64[]
    dB_values = Float64[]
    
    for (A_idx, A_val) in enumerate(A_range)
        p_current[:A] = A_val
        # seq
        p_values = [
            p_current[:m1], p_current[:m2], p_current[:mx],
            p_current[:g], p_current[:d], p_current[:g0],
            p_current[:f], p_current[:A]
        ]

        # enz
        #=p_values = [
            p_current[:m1], p_current[:m2], p_current[:mx],
            p_current[:g], p_current[:d], p_current[:r],
            p_current[:f], p_current[:A]
        ]=#

        # phos
        #=p_values = [
            p_current[:k1], p_current[:k2], p_current[:mU],
            p_current[:mX], p_current[:d], p_current[:f],
            p_current[:A]
        ]=#
        
        prob = ODEProblem(mm.IFFL, mm.u0, tspan, p_values)
        sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6, saveat=0.005)
        
        times = sol.t

        mask = times .> 50.0
        A_signal = (maximum(sol[1, mask]) - minimum(sol[1, mask]))/2
        A_ref = p_current[:A]
        dB = 20*log10(A_signal/A_ref)

        push!(dB_values, dB)
        
        println(" dB=$(round(dB, digits=2))")
    end
    
    plot!(
        fig_dose,
        log10.(A_range),
        dB_values,
        label="m1 = $m1_val",
        color=colors_m1[m1_idx],
        linewidth=2.5,
        alpha=0.8
    )
    
end  

save_dose = "Sequestration_DoseRespose_F$(f_focus)_$M.png"
savefig(fig_dose, save_dose)






# --- Main script for IFFL analysis ---
# --- Rodrigo Aguilar


# --- Required libraries ---
using DifferentialEquations, Plots, Printf, OrderedCollections, LinearAlgebra, LaTeXStrings, Statistics

# --- Model ---
M = "Phosporilation"  # Options: "Sequestration", "Enzimatic", "Phosporilation"
mm = include(string("Library/$M.jl"))
  


# ===========================================================================

n_points = 10  
m1_range_log = range(-1.0, 2.0, length=n_points)
m2_range_log= range(-1.0, 2.0, length=n_points)  
m1_range = 10.0 .^ m1_range_log
m2_range = 10.0 .^ m2_range_log

tspan = (0.0, 20000.0) 

 # Frecuencia de referencia (Hz)
f_ref = 0.001 

# Frecuencias a explorar (Hz)
n_frequencies = 10
f_search_range = 10.0 .^ range(-3.0, 2.0, length=10)  

# almacenar resultados
DB_ref_grid = zeros(n_points, n_points)      
res_freq_grid = zeros(n_points, n_points)    # Frecuencia de resonancia (log10)
DB_res_grid = zeros(n_points, n_points)      # DB en frecuencia de resonancia
DB_diff_grid = zeros(n_points, n_points)     # Diferencia DB_res - DB_ref
resonance_found_grid = zeros(Bool, n_points, n_points)  # Si se encontró resonancia


# CALCULAR DB
function calculate_db_response(f_val, m1_val, m2_val, p_base, tspan, u0, IFFL)
    p_current = deepcopy(p_base)
    p_current[:f] = f_val

    # Seq
    #=p_current[:m1] = m1_val
    p_current[:m2] = m2_val
    p_values = [
        p_current[:m1], p_current[:m2], p_current[:mx],
        p_current[:g], p_current[:d], p_current[:g0],
        p_current[:f], p_current[:A]
    ]=#

    # Enz 
    #=p_current[:m1] = m1_val
    p_current[:m2] = m2_val
    p_values = [
        p_current[:m1], p_current[:m2], p_current[:mx],
        p_current[:g], p_current[:d], p_current[:r],
        p_current[:f], p_current[:A]
    ]=#

    # Phos
    p_current[:k1] = m1_val
    p_current[:k2] = m2_val
    p_values = [
        p_current[:k1], p_current[:k2], p_current[:mU],
        p_current[:mX], p_current[:d], p_current[:f],
        p_current[:A]
    ]

    prob = ODEProblem(IFFL, u0, tspan, p_values)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6, saveat=0.01)
    
    times = sol.t
    U1_concentration = sol[1, :]
    
    mask = times .> 50.0
    
    if sum(mask) > 0
        A_signal = (maximum(sol[2, mask]) - minimum(sol[2, mask]))/2
        A_ref = p_current[:A]
        dB = 20*log10(A_signal/A_ref)
        return dB, A_signal, A_ref
    else 
        return NaN, NaN, NaN

    end

end

# ===========================================================================
#  RESONANCIA PARA CADA COMBINACIÓN
# ===========================================================================

total_combinations = n_points^2

for (i, m1_val) in enumerate(m1_range)
    for (j, m2_val) in enumerate(m2_range)
        println("   Analizando combinación $i, $j    ")
        try
            # 1. Calcular DB a frecuencia de referencia
            db_ref, _, _ = calculate_db_response(f_ref, m1_val, m2_val, mm.p, tspan, mm.u0, mm.IFFL)
            println("     • DB a f_ref=$(f_ref) Hz: $(round(db_ref, digits=2)) dB")
            
            # 2. Buscar frecuencia de resonancia (máxima amplitud)
            db_values = Float64[]
            f_values = Float64[]
            
            for f_test in f_search_range
                println("       - Probando f=$(round(f_test, digits=4)) Hz")
                db_test, amplitude, mean_val = calculate_db_response(f_test, m1_val, m2_val, mm.p, tspan, mm.u0, mm.IFFL)
                push!(db_values, db_test)
                push!(f_values, f_test)
            end
            
            println("     • DB values: ", db_values)
            println("     • F values: ", f_values)
            # Extraer primer valor de DB para referencia
            db_ref_search = db_values[1]
            DB_ref_grid[i, j] = db_values[1]
            # Encontrar máxima respuesta (resonancia)
            if !isempty(db_values)
                max_db_idx = argmax(db_values)
                max_db = db_values[max_db_idx]
                res_freq = f_values[max_db_idx]
                
                # Verificar si es una resonancia real (umbral de 3 dB)
                if max_db > 3.0 && max_db_idx > 1 && max_db_idx < length(db_values)
                    resonance_found_grid[i, j] = true
                    res_freq_grid[i, j] = log10(res_freq)
                    DB_res_grid[i, j] = max_db
                    DB_diff_grid[i, j] = max_db - db_ref_search
                else
                    resonance_found_grid[i, j] = false
                    res_freq_grid[i, j] = -Inf  # No resonancia
                    DB_res_grid[i, j] = db_ref_search
                    DB_diff_grid[i, j] = 0.0
                end
            else
                resonance_found_grid[i, j] = false
                res_freq_grid[i, j] = -Inf
                DB_res_grid[i, j] = db_ref_search
                DB_diff_grid[i, j] = 0.0
            end
        catch e
            println("   Error en (m1=$(round(m1_val,3)), m2=$(round(m2_val,3))): $e")
            DB_ref_grid[i, j] = 0.0
            resonance_found_grid[i, j] = false
            res_freq_grid[i, j] = -Inf
            DB_res_grid[i, j] = 0.0
            DB_diff_grid[i, j] = 0.0
        end
    end
end

# ===========================================================================
# GRÁFICA 1: DB A FRECUENCIA DE REFERENCIA (0.001 Hz)
# ===========================================================================

fig1 = plot(
    xlabel=L"log_{10}(m_1)",
    ylabel=L"log_{10}(m_2)",
    title="Respuesta DB a f = $(f_ref) Hz",
    size=(900, 800),
    aspect_ratio=1,
    colorbar_title="DB (dB)"
)

DB_masked = copy(DB_ref_grid)
DB_masked[DB_ref_grid .== 0.0] .= NaN

heatmap!(
    fig1,
    m1_range_log,
    m2_range_log,
    DB_masked,
    clims=(-40, 20),
    c=:balance,
    colorbar=true,
    nan_color=:white
)


# Guardar gráfica
save_path1 = "./MultipleFunctions_IFFL/DB_at_ref_frequency_$M.png"
savefig(fig1, save_path1)

# =====××××====10==××××========10===××××===========××××==========××××=========××××======××××===========××××===========××××==========××××===============×=======×==========
# GRÁFICA 2: FRECUENCIA DE RESONANCIA (log10)
# ====×======××××======×===========××××=========××××==========××××============××××============××××=====

# Crea×r máscara para regiones con resonancia
res_freq_masked = copy(res_freq_grid)
res_freq_masked[.!resonance_found_grid] .= NaN

fig2 = plot(
    xlabel=L"log_{10}(m_1)",
    ylabel=L"log_{10}(m_2)",
    title="Frecuencia de Resonancia (log₁₀ Hz)",
    size=(900, 800),
    aspect_ratio=1,
    colorbar_title="log₁₀(f_res)"
)

# Heatmap de frecuencia de resonancia
heatmap!(fig2, 
         m1_range_log, 
         m2_range_log, 
         res_freq_masked,
         clims=(-3, 2),
         c=:viridis,
         colorbar=true)

# Contornos
contour!(fig2, 
         m1_range_log, 
         m2_range_log, 
         res_freq_grid,
         levels=10,
         color=:white,
         linewidth=1,
         alpha=0.7)

# Línea diagonal
plot!(fig2, 
      [-1, 2], [-1, 2], 
      linewidth=2, 
      linestyle=:dash, 
      color=:white)

save_path2 = "./MultipleFunctions_IFFL/Resonance_frequency_$M.png"
savefig(fig2, save_path2)

# ===========================================================================
# GRÁFICA 3: DIFERENCIA DB (Resonancia - Referencia)
# ===========================================================================

fig3 = plot(
    xlabel=L"log_{10}(m_1)",
    ylabel=L"log_{10}(m_2)",
    title="Ganancia de Resonancia: DB_res - DB_ref",
    size=(900, 800),
    aspect_ratio=1,
    colorbar_title="ΔDB (dB)"
)

# Enmascarar regiones sin resonancia
DB_diff_masked = copy(DB_diff_grid)
DB_diff_masked[.!resonance_found_grid] .= 0.0

heatmap!(fig3, 
         m1_range_log, 
         m2_range_log, 
         DB_diff_masked,
         clims=(0, 30),
         c=:hot,
         colorbar=true)

contour!(fig3, 
         m1_range_log, 
         m2_range_log, 
         DB_diff_masked,
         levels=10,
         color=:black,
         linewidth=1,
         alpha=0.7)

# Línea diagonal
plot!(fig3, 
      [-1, 2], [-1, 2], 
      linewidth=2, 
      linestyle=:dash, 
      color=:white)

save_path3 = "./MultipleFunctions_IFFL/Resonance_gain_$M.png"
savefig(fig3, save_path3)

# Barrido del laboratorio de la sesión 7: razón de participación inversa Y₂
# del modelo de Bouchaud-Mézard para varias intensidades de intercambio J y
# tamaños N, en campo medio y en una red aleatoria con 4 socios en promedio.
# Guarda Data/condensacion_bm.csv (modelo, J, N, replica, Y2).
#
# Uso (desde la raíz del repositorio; tarda unos minutos):
#     julia --project=Code -t auto Code/simular_condensacion.jl
#
# Las simulaciones grandes son demasiado lentas para el cuaderno, así que se
# corren aquí una vez. `red_aleatoria` y `simular_bm` son las mismas funciones
# del apéndice de Code/s07_condensacion.jl.

using CSV, DataFrames, Random, Statistics

"Red aleatoria con grado promedio c: se emparejan al azar c 'medias aristas' por agente."
function red_aleatoria(N, c; rng = Xoshiro(7))
    vecinos = [Int[] for _ in 1:N]
    extremos = shuffle(rng, repeat(1:N, c))
    for k in 1:2:length(extremos)-1
        a, b = extremos[k], extremos[k+1]
        a != b && (push!(vecinos[a], b); push!(vecinos[b], a))
    end
    for i in 1:N
        isempty(vecinos[i]) && push!(vecinos[i], mod1(i + 1, N))   # nadie queda aislado
    end
    vecinos
end

"Modelo de Bouchaud-Mézard con σ = 1; ver el apéndice del cuaderno de la sesión 7."
function simular_bm(N, J; vecinos = nothing, φ = 0.0, dt = 0.01, T = 200.0,
                    rng = Xoshiro(1), Y2 = false)
    w = ones(N)
    nuevo = similar(w)
    pasos = round(Int, T / dt)
    suma_Y2, cuenta = 0.0, 0
    for k in 1:pasos
        w .*= exp.(sqrt(dt) .* randn(rng, N) .- dt / 2)
        if vecinos === nothing
            w .+= J * dt .* (mean(w) .- w)
        else
            for i in 1:N
                s = 0.0
                for j in vecinos[i]
                    s += w[j] - w[i]
                end
                nuevo[i] = w[i] + J * dt * s / length(vecinos[i])
            end
            w .= nuevo
        end
        φ > 0 && (w .+= φ * dt .* (mean(w) .- w))
        w ./= mean(w)
        if Y2 && k > pasos ÷ 2 && k % 10 == 0
            suma_Y2 += sum(w .^ 2) / N^2
            cuenta += 1
        end
    end
    Y2 ? (w, suma_Y2 / cuenta) : w
end

const JS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8]
const NS = [1000, 5000, 20000]
const REPLICAS = 4
const T = 400.0

tareas = vec(collect(Iterators.product(("campo medio", "red"), JS, NS, 1:REPLICAS)))
Y2 = zeros(length(tareas))
println("Corriendo ", length(tareas), " simulaciones con ", Threads.nthreads(), " hilos...")
Threads.@threads for k in eachindex(tareas)
    modelo, J, N, r = tareas[k]
    vecinos = modelo == "red" ? red_aleatoria(N, 4; rng = Xoshiro(100r + N)) : nothing
    Y2[k] = simular_bm(N, J; vecinos, T, rng = Xoshiro(r), Y2 = true)[2]
end

resultado = DataFrame(modelo = first.(tareas), J = getindex.(tareas, 2),
                      N = getindex.(tareas, 3), replica = last.(tareas), Y2 = Y2)
sort!(resultado, [:modelo, :J, :N, :replica])
CSV.write(joinpath(@__DIR__, "..", "Data", "condensacion_bm.csv"), resultado)
println(combine(groupby(resultado, [:modelo, :J, :N]), :Y2 => median => :Y2_mediana))

# Genera las figuras de las notas de repaso (Class_notes/repaso/figuras/).
# Uso, desde la raíz del repositorio:  julia --project=Code Code/figuras_repaso.jl

ENV["GKSwstype"] = "100"          # GR sin ventana
using Plots, Distributions, Random, Statistics, StatsBase, LinearAlgebra
gr()

const SALIDA = joinpath(@__DIR__, "..", "Class_notes", "repaso", "figuras")
mkpath(SALIDA)

# Paleta del curso (la misma del cuaderno de la sesión 1)
const AZUL    = "#2a78d6"
const NARANJA = "#eb6834"
const AQUA    = "#1baf7a"
const GRIS    = "#52514e"
const TINTA   = "#0b0b0b"
const FONDO   = "#fcfcfb"
const RAMPA   = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial, claro → oscuro

default(size = (720, 320), dpi = 150, legend = false, grid = :y, gridalpha = 0.15,
        framestyle = :axes, linewidth = 2, titlefontsize = 11, guidefontsize = 9,
        tickfontsize = 8, legendfontsize = 8, background_color = FONDO,
        foreground_color_legend = nothing, background_color_legend = nothing,
        foreground_color_axis = GRIS, foreground_color_border = GRIS,
        foreground_color_guide = TINTA, foreground_color_text = GRIS,
        titlefontcolor = TINTA, margin = 4Plots.mm)

guardar(p, nombre) = (savefig(p, joinpath(SALIDA, nombre)); println("  ", nombre))
etiqueta(p, x, y, s, color = GRIS; pos = :left, tam = 8) =
    annotate!(p, x, y, text(s, tam, pos, color))

println("Generando figuras en ", normpath(SALIDA))

# ── 1. Escalas: exponencial vs. ley de potencias ──────────────────────────────
let x = range(1, 20, length = 400)
    expo = exp.(-(x .- 1))           # CCDF exponencial (normalizada en x = 1)
    pot  = x .^ -2.0                 # CCDF de Pareto con α = 2
    paneles = []
    for (titulo, xs, ys) in [("Escala lineal", :identity, :identity),
                             ("Semilog: eje y logarítmico", :identity, :log10),
                             ("Log-log: ambos ejes logarítmicos", :log10, :log10)]
        p = plot(x, expo; color = AZUL, xscale = xs, yscale = ys, title = titulo,
                 xlabel = "x", ylabel = "P(X > x)", ylims = ys == :log10 ? (1e-9, 1.5) : (0, 1.05),
                 label = "exponencial")
        plot!(p, x, pot; color = NARANJA, label = "ley de potencias")
        ys == :log10 && plot!(p; yticks = 10.0 .^ (0:-2:-8))
        xs == :log10 && plot!(p; xticks = ([1, 2, 5, 10, 20], ["1", "2", "5", "10", "20"]))
        push!(paneles, p)
    end
    plot!(paneles[1]; legend = :topright)
    etiqueta(paneles[2], 6, 1e-1, "ley de potencias", NARANJA)
    etiqueta(paneles[2], 2.0, 1e-7, "exponencial: recta", AZUL)
    etiqueta(paneles[3], 1.1, 3e-3, "ley de potencias: recta", NARANJA)
    etiqueta(paneles[3], 3.2, 1e-6, "exponencial", AZUL)
    guardar(plot(paneles...; layout = (1, 3), size = (1080, 330)), "escalas.png")
end

# ── 2. Teorema del límite central: sumas de exponenciales ─────────────────────
let z = range(-3.5, 4.5, length = 500)
    p = plot(z, pdf.(Normal(), z); color = GRIS, ls = :dash, lw = 1.5,
             title = "Suma estandarizada de n exponenciales", xlabel = "z", ylabel = "densidad",
             ylims = (0, 1.05))
    for (i, n) in enumerate((1, 2, 5, 30))
        # S_n ~ Gamma(n, 1); estandarizada: z = (S_n − n)/√n
        dens = [pdf(Gamma(n, 1), n + sqrt(n) * zi) * sqrt(n) for zi in z]
        plot!(p, z, dens; color = RAMPA[i])
    end
    etiqueta(p, -0.55, 0.98, "n = 1", RAMPA[1])
    etiqueta(p, 0.75, 0.52, "n = 2", RAMPA[2])
    etiqueta(p, 1.75, 0.22, "n = 5", RAMPA[3])
    etiqueta(p, -2.6, 0.17, "n = 30", RAMPA[4])
    etiqueta(p, 2.3, 0.07, "normal", GRIS)
    guardar(p, "tlc.png")
end

# ── 3. Caminata aleatoria y la regla √n ───────────────────────────────────────
let rng = Xoshiro(3), n = 1_000
    p = plot(title = "30 caminatas aleatorias simples", xlabel = "paso n", ylabel = "posición Sₙ")
    for _ in 1:30
        plot!(p, 0:n, [0; cumsum(rand(rng, (-1, 1), n))]; color = AZUL, alpha = 0.25, lw = 1)
    end
    t = 0:n
    plot!(p, t, 2 .* sqrt.(t); color = NARANJA); plot!(p, t, -2 .* sqrt.(t); color = NARANJA)
    etiqueta(p, 1010, 2sqrt(n), "+2√n", NARANJA); etiqueta(p, 1010, -2sqrt(n), "−2√n", NARANJA)
    guardar(plot(p; xlims = (0, 1120)), "caminata.png")
end

# ── 4. Gibrat: promedio de ensamble vs. trayectoria típica ────────────────────
let rng = Xoshiro(4), T = 100, M = 10_000
    X = ones(M, T + 1)
    for t in 1:T
        X[:, t+1] = X[:, t] .* rand(rng, (1.5, 0.6), M)
    end
    p = plot(title = "Crecimiento multiplicativo: ×1.5 o ×0.6 con igual probabilidad",
             xlabel = "periodo t", ylabel = "riqueza (escala log)", yscale = :log10,
             ylims = (1e-8, 1e5), yticks = 10.0 .^ (-8:2:4))
    for i in 1:25
        plot!(p, 0:T, X[i, :]; color = GRIS, alpha = 0.25, lw = 1)
    end
    plot!(p, 0:T, 1.05 .^ (0:T); color = AZUL)
    plot!(p, 0:T, vec(median(X; dims = 1)); color = NARANJA)
    etiqueta(p, 58, 3e3, "promedio: crece 5% por periodo", AZUL)
    etiqueta(p, 50, 1e-4, "mediana: decrece", NARANJA)
    etiqueta(p, 2, 3e-7, "líneas grises: 25 trayectorias individuales", GRIS)
    guardar(p, "gibrat.png")
end

# ── 5. Campo medio: solución gráfica y bifurcación ────────────────────────────
let m = range(-1.2, 1.2, length = 400)
    p1 = plot(m, m; color = GRIS, ls = :dash, lw = 1.5, title = "Solución gráfica de m = tanh(βJ m)",
              xlabel = "m", ylabel = "", aspect_ratio = :equal, xlims = (-1.2, 1.2), ylims = (-1.2, 1.2),
              grid = :both, legend = :topleft, label = "y = m")
    for (i, K) in enumerate((0.5, 1.0, 1.5))
        plot!(p1, m, tanh.(K .* m); color = RAMPA[i+1], label = "βJ = $K")
    end
    mstar(K) = (x = 0.9; for _ in 1:2_000; x = tanh(K * x); end; x)   # iteración de punto fijo
    ms = mstar(1.5)
    scatter!(p1, [-ms, 0, ms], [-ms, 0, ms]; color = RAMPA[4], ms = 5, msw = 0, label = "")
    etiqueta(p1, -0.35, -1.08, "puntos: soluciones con βJ = 1.5", RAMPA[4])

    K = range(0, 3, length = 300)
    rama = [k <= 1 ? 0.0 : mstar(k) for k in K]
    p2 = plot(K, rama; color = AZUL, title = "Magnetización de equilibrio", xlabel = "βJ",
              ylabel = "m*", ylims = (-1.1, 1.1))
    plot!(p2, K, -rama; color = AZUL)
    plot!(p2, K[K .> 1], zeros(count(K .> 1)); color = AZUL, ls = :dot, lw = 1.5)
    vline!(p2, [1]; color = GRIS, ls = :dash, lw = 1)
    etiqueta(p2, 0.93, -0.5, "transición: βJ = 1", GRIS; pos = :right)
    etiqueta(p2, 1.9, 0.15, "m = 0 inestable", AZUL)
    etiqueta(p2, 0.1, 0.15, "desorden", GRIS); etiqueta(p2, 2.0, 0.75, "consenso", GRIS)
    guardar(plot(p1, p2; layout = (1, 2), size = (900, 380)), "campo_medio.png")
end

# ── 6. Histéresis en campo medio ──────────────────────────────────────────────
let K = 1.5
    m = range(-0.999, 0.999, length = 1000)
    h = atanh.(m) .- K .* m                 # inversa de m = tanh(K m + h), con β = 1
    mc = sqrt(1 - 1 / K)                     # bordes de las ramas estables
    hc = atanh(mc) - K * mc                  # (negativo); el salto ocurre en ±|hc|
    est = abs.(m) .> mc
    p = plot(title = "Histéresis en campo medio (βJz = 1.5)", xlabel = "incentivo común h",
             ylabel = "adopción neta m", xlims = (-0.8, 0.8), ylims = (-1.1, 1.1))
    plot!(p, h[m .< -mc], m[m .< -mc]; color = AZUL)
    plot!(p, h[m .> mc], m[m .> mc]; color = AZUL)
    plot!(p, h[.!est], m[.!est]; color = AZUL, ls = :dot, lw = 1.5)
    plot!(p, [-hc, -hc], [-mc, 0.93]; color = NARANJA, arrow = true, lw = 1.5)
    plot!(p, [hc, hc], [mc, -0.93]; color = NARANJA, arrow = true, lw = 1.5)
    etiqueta(p, -hc + 0.03, 0.0, "salto al subir h", NARANJA)
    etiqueta(p, hc - 0.03, 0.0, "salto al bajar h", NARANJA; pos = :right)
    etiqueta(p, -0.78, 0.55, "punteado: equilibrio inestable", AZUL)
    guardar(p, "histeresis.png")
end

# ── 7. Retrato de fase en una dimensión ───────────────────────────────────────
let f(x) = x * (1 - x) * (x - 0.3), x = range(-0.1, 1.1, length = 400)
    p = plot(x, f.(x); color = AZUL, title = "Retrato de fase de ẋ = x(1 − x)(x − 0.3)",
             xlabel = "x", ylabel = "ẋ = f(x)", ylims = (-0.09, 0.09))
    hline!(p, [0]; color = GRIS, lw = 1)
    scatter!(p, [0, 1], [0, 0]; color = AZUL, ms = 7, msw = 0)                       # estables
    scatter!(p, [0.3], [0]; color = FONDO, ms = 7, msc = AZUL, msw = 2)              # inestable
    for (xa, dir) in [(0.13, -1), (0.6, 1), (1.06, -1), (-0.06, 1)]
        plot!(p, [xa - 0.035dir, xa + 0.035dir], [0, 0]; color = NARANJA, arrow = true, lw = 2)
    end
    etiqueta(p, 0.0, -0.02, "estable", AZUL; pos = :center)
    etiqueta(p, 1.0, -0.02, "estable", AZUL; pos = :center)
    etiqueta(p, 0.3, 0.02, "inestable", AZUL; pos = :center)
    etiqueta(p, 0.33, -0.065, "flechas: hacia dónde se mueve x", NARANJA)
    guardar(p, "retrato_fase.png")
end

# ── 8. Diagrama de bifurcación del mapa logístico ─────────────────────────────
let rs = 2.5:0.0015:4.0
    xs = Float64[]; ys = Float64[]
    for r in rs
        x = 0.5
        for _ in 1:600; x = r * x * (1 - x); end
        for _ in 1:150; x = r * x * (1 - x); push!(xs, r); push!(ys, x); end
    end
    p = scatter(xs, ys; ms = 0.6, msw = 0, color = AZUL, alpha = 0.5,
                title = "Mapa logístico xₜ₊₁ = r xₜ (1 − xₜ)", xlabel = "r",
                ylabel = "valores visitados", ylims = (0, 1.05), grid = false)
    vline!(p, [3.0, 3.5699]; color = GRIS, ls = :dash, lw = 1)
    etiqueta(p, 2.55, 0.95, "punto fijo", GRIS)
    etiqueta(p, 3.02, 0.12, "periodo 2", GRIS)
    etiqueta(p, 3.59, 0.04, "caos (r ≈ 3.57)", GRIS)
    guardar(p, "bifurcacion.png")
end

# ── 9. Eigenvectores: el círculo se vuelve elipse ─────────────────────────────
let A = [2 1; 1 2], θ = range(0, 2π, length = 300)
    c = [cos.(θ) sin.(θ)]'; e = A * c
    p = plot(c[1, :], c[2, :]; color = GRIS, lw = 1.5, aspect_ratio = :equal,
             title = "A = [2 1; 1 2] transforma el círculo unitario", grid = :both,
             xlims = (-3.6, 3.6), ylims = (-3.3, 3.3), xlabel = "", ylabel = "")
    plot!(p, e[1, :], e[2, :]; color = AZUL)
    v1, v2 = [1, 1] / sqrt(2), [1, -1] / sqrt(2)
    for (v, λ) in [(v1, 3), (v2, 1)]
        plot!(p, [0, v[1]], [0, v[2]]; color = NARANJA, arrow = true, lw = 2)
        plot!(p, [0, λ * v[1]], [0, λ * v[2]]; color = NARANJA, ls = :dash, lw = 1.5)
    end
    x = [1.0, 0.0]; Ax = A * x
    plot!(p, [0, x[1]], [0, x[2]]; color = GRIS, arrow = true, lw = 1.5)
    plot!(p, [0, Ax[1]], [0, Ax[2]]; color = GRIS, arrow = true, lw = 1.5, ls = :dash)
    etiqueta(p, 2.2, 2.45, "v₁: estirado ×3", NARANJA)
    etiqueta(p, 0.8, -1.0, "v₂: ×1", NARANJA)
    etiqueta(p, -3.5, 2.9, "gris: x → Ax cambia de dirección", GRIS)
    etiqueta(p, -3.5, -2.9, "círculo gris → elipse azul", GRIS)
    guardar(plot(p; size = (560, 520)), "eigen.png")
end

# ── 10. Marchenko-Pastur: ruido vs. datos con estructura ──────────────────────
let rng = Xoshiro(10), N = 200, T = 800
    Q = T / N; λm, λp = (1 - 1 / sqrt(Q))^2, (1 + 1 / sqrt(Q))^2
    mp(λ) = λm <= λ <= λp ? Q / (2π) * sqrt((λp - λ) * (λ - λm)) / λ : 0.0
    ruido = eigvals(Symmetric(cor(randn(rng, T, N))))
    # Modelo con un factor de mercado y 4 sectores de 50 activos
    f = randn(rng, T); s = randn(rng, T, 4)
    R = 0.35 .* f .+ 0.45 .* s[:, repeat(1:4, inner = 50)] .+ randn(rng, T, N)
    estr = eigvals(Symmetric(cor(R)))
    λs = range(0, 3, length = 400)
    p1 = histogram(ruido; normalize = :pdf, bins = range(0, 3, length = 46), color = GRIS,
                   alpha = 0.45, lw = 0, title = "Ruido puro (N = 200, T = 800)",
                   xlabel = "eigenvalor λ", ylabel = "densidad", xlims = (0, 3.2))
    plot!(p1, λs, mp.(λs); color = NARANJA)
    etiqueta(p1, 1.6, 0.85, "Marchenko-Pastur", NARANJA)
    p2 = histogram(estr; normalize = :pdf, bins = range(0, 12, length = 97), color = GRIS,
                   alpha = 0.45, lw = 0, title = "Mercado + 4 sectores", xlabel = "eigenvalor λ",
                   ylabel = "", xlims = (0, 12.5), ylims = (0, 1.6))
    plot!(p2, λs, mp.(λs); color = NARANJA)
    fuera = filter(>(λp), estr)
    scatter!(p2, min.(fuera, 12.2), fill(0.05, length(fuera)); color = AZUL, ms = 6, msw = 0)
    etiqueta(p2, 3.5, 0.45, "$(length(fuera)) eigenvalores sobre λ₊: señal", AZUL)
    maximum(fuera) > 12.2 && etiqueta(p2, 12.3, 0.16, "λ₁ ≈ $(round(Int, maximum(fuera))) →", AZUL; pos = :right)
    guardar(plot(p1, p2; layout = (1, 2), size = (980, 340)), "marchenko_pastur.png")
end

# ── 11. Ising 2D: simulación vs. solución exacta ──────────────────────────────
function ising_m(L, T; barridos = 3_000, rng = Xoshiro(11))
    s = ones(Int, L, L); β = 1 / T; ms = Float64[]
    for b in 1:barridos
        for _ in 1:L^2
            i, j = rand(rng, 1:L), rand(rng, 1:L)
            v = s[mod1(i+1, L), j] + s[mod1(i-1, L), j] + s[i, mod1(j+1, L)] + s[i, mod1(j-1, L)]
            ΔE = 2 * s[i, j] * v
            (ΔE <= 0 || rand(rng) < exp(-β * ΔE)) && (s[i, j] = -s[i, j])
        end
        b > barridos ÷ 3 && push!(ms, abs(mean(s)))
    end
    mean(ms)
end
let Ts = 1.0:0.1:3.6, Tc = 2 / log(1 + sqrt(2))
    sim = [ising_m(32, T) for T in Ts]
    Tf = range(1.0, 3.6, length = 400)
    exacta = [T < Tc ? (1 - sinh(2 / T)^-4)^(1 / 8) : 0.0 for T in Tf]
    p = plot(Tf, exacta; color = GRIS, title = "Modelo de Ising en la red cuadrada",
             xlabel = "temperatura T (en unidades de J)", ylabel = "|m|", ylims = (0, 1.05))
    scatter!(p, Ts, sim; color = AZUL, ms = 5, msw = 0)
    vline!(p, [Tc]; color = GRIS, ls = :dash, lw = 1)
    etiqueta(p, Tc + 0.05, 0.9, "Tc ≈ 2.27", GRIS)
    etiqueta(p, 1.05, 0.25, "línea: solución exacta (Onsager)", GRIS)
    etiqueta(p, 1.05, 0.15, "puntos: Metropolis, red 32 × 32", AZUL)
    etiqueta(p, 2.6, 0.3, "efecto de tamaño finito", AZUL)
    guardar(p, "ising.png")
end

# ── 12. Intervalos de confianza: qué significa "95%" ──────────────────────────
let rng = Xoshiro(12), K = 50, n = 30
    p = plot(title = "50 intervalos de 95% para la media (μ = 0, muestras de 30 datos)",
             xlabel = "intervalo", ylabel = "estimación ± 1.96 SE", ylims = (-1.1, 1.1))
    hline!(p, [0]; color = GRIS, lw = 1)
    fallas = 0
    for k in 1:K
        x = randn(rng, n); μ̂ = mean(x); se = std(x) / sqrt(n)
        falla = abs(μ̂) > 1.96se; fallas += falla
        plot!(p, [k, k], [μ̂ - 1.96se, μ̂ + 1.96se]; color = falla ? NARANJA : AZUL, lw = falla ? 3 : 2)
        scatter!(p, [k], [μ̂]; color = falla ? NARANJA : AZUL, ms = 3, msw = 0)
    end
    etiqueta(p, 1, -0.95, "$fallas de $K no contienen μ = 0 (se esperan ≈ 2.5)", NARANJA)
    guardar(p, "intervalos.png")
end

# ── 13. Intercambio de dinero: exponencial y gamma ────────────────────────────
function intercambio(N, pasos, λ; rng = Xoshiro(13))
    m = ones(N)
    for _ in 1:pasos
        i, j = rand(rng, 1:N), rand(rng, 1:N); i == j && continue
        ε = rand(rng); bolsa = (1 - λ) * (m[i] + m[j])
        m[i], m[j] = λ * m[i] + ε * bolsa, λ * m[j] + (1 - ε) * bolsa
    end
    m
end
let N = 10_000, pasos = 4_000_000
    m0 = intercambio(N, pasos, 0.0); m5 = intercambio(N, pasos, 0.5)
    bordes = range(0, 6, length = 49); c = (bordes[1:end-1] .+ bordes[2:end]) ./ 2
    dens(m) = fit(Histogram, m, bordes).weights ./ (length(m) * step(bordes))
    x = range(0.01, 6, length = 300)
    k = 1 + 3 * 0.5 / (1 - 0.5)            # forma aproximada de la gamma con ahorro
    p = plot(x, exp.(-x); color = GRIS, ls = :dash, lw = 1.5, yscale = :log10,
             ylims = (1e-3, 1.5), title = "Intercambio aleatorio de dinero (N = 10 000)",
             xlabel = "dinero m (promedio = 1)", ylabel = "densidad (escala log)")
    plot!(p, x, pdf.(Gamma(k, 1 / k), x); color = GRIS, ls = :dash, lw = 1.5)
    scatter!(p, c, max.(dens(m0), 1e-4); color = AZUL, ms = 4, msw = 0)
    scatter!(p, c, max.(dens(m5), 1e-4); color = NARANJA, ms = 4, msw = 0)
    etiqueta(p, 3.3, 0.06, "sin ahorro: exponencial (recta)", AZUL)
    etiqueta(p, 0.75, 0.02, "ahorro 50%: gamma", NARANJA)
    etiqueta(p, 0.15, 1.6e-3, "líneas punteadas: teoría", GRIS)
    guardar(p, "intercambio.png")
end

println("Listo.")

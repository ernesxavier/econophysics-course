#%% hide
md"""
# Sesión 3: Procesos estocásticos

**Objetivos de la sesión**
1. Simular los procesos aleatorios básicos: caminata aleatoria, cadena de
   Markov, ecuación de Langevin y movimiento browniano geométrico.
2. Entender la simulación de Euler-Maruyama como una regla de actualización:
   "tendencia + ruido".
3. Construir un proceso con volatilidad que cambia en el tiempo (GARCH) y
   ajustarlo a datos reales.
4. Responder la pregunta del laboratorio: **¿reproduce GARCH los hechos
   estilizados de la sesión 1?**

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s03_procesos_estocasticos.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B1–E.B3, E.I1–E.I5
(procesos estocásticos), P.B2 (varianza de una suma), C.B3 y C.I4
(derivadas y ecuaciones diferenciales). Ruta avanzada: E.A1–E.A4.
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: ¿qué proceso genera estas series?

En las sesiones 1 y 2 medimos tres propiedades de los retornos: casi no
tienen autocorrelación, su volatilidad se agrupa en el tiempo y sus colas
siguen una ley de potencias con ``\alpha \approx 3``. Hoy cambiamos la
pregunta: **¿qué mecanismo aleatorio produce series así?**

Compara los retornos reales con retornos normales independientes de la misma
volatilidad. Los dos oscilan alrededor de cero, pero sólo uno tiene periodos
de calma y periodos de turbulencia.

**Activo:** $(@bind activo Select(unique(datos.activo)))
"""

#%%
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame);

#%%
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    r = diff(log.(serie.cierre))           # retornos logarítmicos diarios
    σ_r = std(r)
end;

#%% hide
let
    n = min(length(r), 2000)
    real = 100 .* r[end-n+1:end]
    iid = 100σ_r .* randn(Xoshiro(3), n)
    lim = maximum(abs, real) * 1.05
    p1 = plot(real; color = AZUL, lw = 0.6, ylims = (-lim, lim), ylabel = "r (%)",
              title = "$activo: últimos $n días")
    p2 = plot(iid; color = GRIS, lw = 0.6, ylims = (-lim, lim), ylabel = "r (%)",
              xlabel = "día", title = "Normales independientes con la misma volatilidad")
    plot(p1, p2; layout = (2, 1), size = (720, 460))
end

#%% hide
md"""
La estrategia de hoy es la de un físico: proponer el modelo más simple,
simularlo, compararlo con los datos y, cuando falla, agregar sólo lo
necesario. Empezamos por el más simple de todos.

## 1. Repaso: sumas y tasas de cambio

Dos herramientas sostienen toda la sesión.

**La varianza de una suma.** Si ``\xi_1, \dots, \xi_n`` son independientes,

```math
\mathbb{E}\Big[\sum_i \xi_i\Big] = \sum_i \mathbb{E}[\xi_i], \qquad
\mathrm{Var}\Big[\sum_i \xi_i\Big] = \sum_i \mathrm{Var}[\xi_i].
```

**La derivada como regla de actualización.** Si ``x`` cambia a una tasa
``f(x)``, en un intervalo pequeño ``\Delta t`` cambia aproximadamente
``f(x)\,\Delta t``:

```math
x_{t + \Delta t} \approx x_t + f(x_t)\,\Delta t.
```

Así se simula una ecuación diferencial: aplicando la regla paso a paso. En la
sección 4 le agregaremos ruido.
"""

#%% hide
md"""
## 2. La caminata aleatoria

En cada paso, un caminante sube 1 con probabilidad ``p`` o baja 1 con
probabilidad ``1 - p``. Su posición después de ``n`` pasos es
``S_n = \xi_1 + \dots + \xi_n``. Cada paso tiene media ``2p - 1`` y varianza
``4p(1-p)``, así que, por la sección 1,

```math
\mathbb{E}[S_n] = n(2p - 1), \qquad \mathrm{Var}[S_n] = 4np(1-p).
```

La media crece como ``n``, pero la dispersión crece como ``\sqrt{n}``: es la
firma de la **difusión**. Con ``p = 1/2`` es el modelo de Bachelier (1900)
para los precios.

Pasos: $(@bind n_pasos Slider(100:100:5000, default = 1000, show_value = true))
  Probabilidad de subir ``p``: $(@bind p_subir Slider(0.40:0.01:0.60, default = 0.50, show_value = true))
"""

#%% hide
let
    rng = Xoshiro(4)
    S = [cumsum(rand(rng, n_pasos) .< p_subir) .* 2 .- (1:n_pasos) for _ in 1:200]
    t = 0:n_pasos
    μ = t .* (2p_subir - 1); σ = sqrt.(4 .* t .* p_subir .* (1 - p_subir))
    p = plot(xlabel = "paso n", ylabel = "posición Sₙ", legend = :topleft,
             title = "200 caminatas aleatorias con p = $(p_subir)")
    for (i, s) in enumerate(S)
        plot!(p, t, [0; s]; color = AZUL, alpha = 0.12, lw = 0.8, label = i == 1 ? "Trayectorias" : "")
    end
    plot!(p, t, μ; color = NARANJA, lw = 2, label = "Media teórica")
    plot!(p, t, μ .+ 2σ; color = NARANJA, lw = 1.5, ls = :dash, label = "Media ± 2 desviaciones")
    plot!(p, t, μ .- 2σ; color = NARANJA, lw = 1.5, ls = :dash, label = "")
end

#%% hide
md"""
En el paso final, la varianza simulada de las 200 caminatas es
$(round(var_caminata.simulada, digits = 0)) y la teórica es
$(round(var_caminata.teorica, digits = 0)). Duplica el número de pasos: la
dispersión no se duplica, crece en un factor ``\sqrt{2} \approx 1.41``.
"""

#%% hide
var_caminata = let rng = Xoshiro(4)
    finales = [sum(2 .* (rand(rng, n_pasos) .< p_subir) .- 1) for _ in 1:200]
    (simulada = var(finales), teorica = 4n_pasos * p_subir * (1 - p_subir))
end;

#%% hide
md"""
## 3. Cadenas de Markov: un mercado con dos estados

Una **cadena de Markov** salta entre estados con probabilidades que sólo
dependen del estado actual. Imagina un mercado que puede estar en **calma**
(volatilidad baja) o en **turbulencia** (volatilidad tres veces mayor). Cada
día, desde la calma se pasa a turbulencia con probabilidad ``p``, y desde la
turbulencia se vuelve a la calma con probabilidad ``q``:

```math
P = \begin{pmatrix} 1-p & p \\ q & 1-q \end{pmatrix}.
```

A la larga, la fracción de días en turbulencia es ``\pi_T = p/(p+q)`` y un
episodio turbulento dura en promedio ``1/q`` días.

Probabilidad de entrar en turbulencia ``p``: $(@bind p_markov Slider(0.002:0.002:0.05, default = 0.01, show_value = true))
  de volver a la calma ``q``: $(@bind q_markov Slider(0.01:0.01:0.2, default = 0.05, show_value = true))
"""

#%% hide
mercado = simular_regimenes(p_markov, q_markov; n = 3000, σ_calma = 0.008, factor = 3);

#%% hide
let
    m = mercado
    p1 = plot(100 .* m.r; color = AZUL, lw = 0.6, ylabel = "r (%)",
              title = "Retornos simulados con dos regímenes")
    turb = findall(m.estado .== 2)
    isempty(turb) || vspan!(p1, rangos(turb); color = NARANJA, alpha = 0.15, lw = 0)
    p2 = plot(1:60, autocor(abs.(m.r), 1:60); color = AZUL, lw = 2, xlabel = "rezago τ (días)",
              ylabel = "ρ(τ)", title = "Autocorrelación de |r|", ylims = (-0.05, 0.5))
    plot!(p2, 1:60, autocor(m.r, 1:60); color = GRIS, lw = 1.5)
    hline!(p2, [0]; color = :black, lw = 0.5)
    annotate!(p2, 30, 0.42, text("azul: |r| · gris: r", 8, GRIS))
    plot(p1, p2; layout = (2, 1), size = (720, 480))
end

#%% hide
md"""
Fracción de días en turbulencia: $(round(100 * mean(mercado.estado .== 2), digits = 1))%
en la simulación contra $(round(100p_markov / (p_markov + q_markov), digits = 1))% de la
teoría (zonas sombreadas).

Este modelo tan simple ya reproduce dos hechos estilizados: los retornos no
tienen autocorrelación, pero ``|r|`` sí, porque un día turbulento suele seguir
a otro. Y la mezcla de dos normales con varianzas distintas tiene colas más
pesadas que una normal (curtosis en exceso de
$(round(kurtosis(mercado.r), digits = 1))). Pero la autocorrelación de ``|r|``
decae como ``(1 - p - q)^\tau``: exponencialmente, más rápido que en los datos.
"""

#%% hide
md"""
## 4. Tendencia + ruido: la ecuación de Langevin

Una **ecuación diferencial estocástica** escribe el cambio de ``X`` como una
tendencia más un ruido:

```math
dX_t = A(X_t)\,dt + B(X_t)\,dW_t .
```

Se simula con el **método de Euler-Maruyama**, la regla de actualización de
la sección 1 con un término aleatorio:

```math
X_{k+1} = X_k + A(X_k)\,\Delta t + B(X_k)\,\sqrt{\Delta t}\,Z_k, \qquad Z_k \sim \mathcal{N}(0, 1).
```

El ruido escala con ``\sqrt{\Delta t}``, no con ``\Delta t``: es la difusión de
la sección 2.

**Ejemplo: Ornstein-Uhlenbeck.** Con ``A(x) = -\theta x`` y ``B = \sigma``, el
proceso es atraído hacia 0 con fuerza ``\theta``, mientras el ruido lo aleja.
Es el modelo de Vasicek para tasas de interés y sirve para cualquier variable
que "regresa a su nivel normal", como un diferencial de precios.

Fuerza de regreso ``\theta``: $(@bind θ_ou Slider(0.1:0.1:3.0, default = 1.0, show_value = true))
  Ruido ``\sigma``: $(@bind σ_ou Slider(0.1:0.1:2.0, default = 0.5, show_value = true))
"""

#%% hide
let
    T, Δt, x0 = 6.0, 0.01, 3.0
    t = 0:Δt:T
    rng = Xoshiro(5)
    p = plot(xlabel = "tiempo t", ylabel = "X", legend = :topright,
             title = "Ornstein-Uhlenbeck desde X₀ = 3")
    for i in 1:25
        x = euler_maruyama_ref(x -> -θ_ou * x, x -> σ_ou, x0, Δt, length(t) - 1, rng)
        plot!(p, t, x; color = AZUL, alpha = 0.25, lw = 0.8, label = i == 1 ? "Trayectorias" : "")
    end
    m = x0 .* exp.(-θ_ou .* t)
    s = sqrt.(σ_ou^2 .* (1 .- exp.(-2θ_ou .* t)) ./ (2θ_ou))
    plot!(p, t, m; color = NARANJA, lw = 2, label = "Media exacta")
    plot!(p, t, m .+ 2s; color = NARANJA, lw = 1.5, ls = :dash, label = "Media ± 2 desviaciones")
    plot!(p, t, m .- 2s; color = NARANJA, lw = 1.5, ls = :dash, label = "")
end

#%% hide
md"""
El proceso olvida su punto de partida en un tiempo del orden de ``1/\theta``
y después fluctúa con varianza ``\sigma^2/(2\theta)`` =
$(round(σ_ou^2 / (2θ_ou), digits = 3)). Más fuerza de regreso significa
fluctuaciones más pequeñas; más ruido, más grandes.
"""

#%% hide
md"""
## 5. El movimiento browniano geométrico

El modelo estándar para el precio de una acción (y la base de Black-Scholes)
es

```math
dP_t = \mu P_t\,dt + \sigma P_t\,dW_t .
```

El cambio es proporcional al precio: lo que es aleatorio es el **retorno**.
Por eso el precio nunca es negativo y sus retornos logarítmicos son normales
e independientes.

Tendencia anual ``\mu``: $(@bind μ_gbm Slider(-0.1:0.02:0.3, default = 0.08, show_value = true))
  Volatilidad anual ``\sigma``: $(@bind σ_gbm Slider(0.05:0.05:0.8, default = 0.3, show_value = true))
"""

#%% hide
gbm_trayectorias = simular_gbm(μ_gbm, σ_gbm; P0 = 100.0, años = 10, M = 500);

#%% hide
let
    g = gbm_trayectorias
    p = plot(yscale = :log10, xlabel = "años", ylabel = "precio (escala log)", legend = :topleft,
             title = "Movimiento browniano geométrico (μ = $(μ_gbm), σ = $(σ_gbm))")
    for i in 1:40
        plot!(p, g.t, g.P[:, i]; color = GRIS, alpha = 0.25, lw = 0.7, label = i == 1 ? "Trayectorias" : "")
    end
    plot!(p, g.t, 100 .* exp.(μ_gbm .* g.t); color = AZUL, lw = 2, label = "Promedio: crece a tasa μ")
    plot!(p, g.t, vec(median(g.P; dims = 2)); color = NARANJA, lw = 2,
          label = "Mediana simulada: crece a tasa μ − σ²/2")
end

#%% hide
md"""
El precio **promedio** crece a la tasa ``\mu``, pero el precio **típico** (la
mediana) crece a ``\mu - \sigma^2/2`` = $(round(μ_gbm - σ_gbm^2 / 2, digits = 3)) por año.
Con volatilidad alta, la mayoría de las trayectorias termina por debajo del
promedio, que lo sostienen unas pocas con mucha suerte. (Es el mismo efecto
del crecimiento multiplicativo del repaso E.I6, y lo explica el lema de Itô,
repaso E.A2.)

¿Reproduce este modelo los hechos estilizados? No: sus retornos logarítmicos
son normales e independientes, así que no tienen colas pesadas ni
agrupamiento de volatilidad. Hace falta que la volatilidad **cambie en el
tiempo**.
"""

#%% hide
md"""
## 6. GARCH: volatilidad con memoria

### 6.1 El modelo

El modelo GARCH(1,1) (Bollerslev, 1986) hace que la varianza de hoy dependa
del retorno de ayer y de la varianza de ayer:

```math
r_t = \sigma_t\,\varepsilon_t, \qquad
\sigma_t^2 = \omega + a\,r_{t-1}^2 + b\,\sigma_{t-1}^2, \qquad \varepsilon_t \sim \mathcal{N}(0, 1).
```

Un día con un retorno grande sube la volatilidad del día siguiente, y ``b``
controla cuánto tarda en olvidarse. El proceso es estacionario si
``a + b < 1``; la persistencia ``a + b`` define la **vida media** de un choque
de volatilidad: ``\ln(1/2)/\ln(a + b)`` días.

Reacción a las noticias ``a``: $(@bind a_garch Slider(0.0:0.01:0.3, default = 0.08, show_value = true))
  Memoria ``b``: $(@bind b_garch Slider(0.0:0.01:0.99, default = 0.90, show_value = true))
"""

#%% hide
if a_garch + b_garch >= 1
    md"""
    !!! warning "Proceso no estacionario"
        Con ``a + b \ge 1`` la varianza no tiene un nivel de largo plazo. Baja ``a`` o ``b``.
    """
else
    let
        g = simular_garch(3000, a_garch, b_garch, 0.01^2; rng = Xoshiro(6))
        vida = log(0.5) / log(a_garch + b_garch)
        p1 = plot(100 .* g.r; color = AZUL, lw = 0.6, ylabel = "r (%)",
                  title = "GARCH(1,1): a + b = $(a_garch + b_garch), vida media ≈ $(round(vida, digits = 1)) días")
        p2 = plot(100 .* sqrt.(g.σ2); color = NARANJA, lw = 1.2, ylabel = "σₜ (%)", xlabel = "día",
                  title = "Volatilidad condicional")
        plot(p1, p2; layout = (2, 1), size = (720, 480))
    end
end

#%% hide
md"""
### 6.2 Laboratorio: ¿reproduce GARCH los hechos estilizados?

Ajustamos el GARCH(1,1) a los retornos de $activo por máxima verosimilitud
(búsqueda en una cuadrícula de valores de ``a`` y ``b``, con ``\omega`` elegido
para que la varianza de largo plazo coincida con la de los datos). Después
simulamos 20 series del modelo ajustado, del mismo largo que los datos, y las
comparamos con la serie real.
"""

#%% hide
ajuste_garch = ajustar_garch(r);

#%% hide
md"""
Parámetros estimados: ``a = `` $(ajuste_garch.a), ``b = `` $(ajuste_garch.b),
persistencia ``a + b = `` $(round(ajuste_garch.a + ajuste_garch.b, digits = 4)),
vida media de un choque ≈ $(round(Int, log(0.5) / log(ajuste_garch.a + ajuste_garch.b))) días.
"""

#%% hide
simulaciones_garch = [simular_garch(length(r), ajuste_garch.a, ajuste_garch.b, ajuste_garch.s2;
                                    rng = Xoshiro(100 + i)).r for i in 1:20];

#%% hide
let
    L = 250
    acf_real = autocor(abs.(r), 1:L)
    acfs = reduce(hcat, [autocor(abs.(s), 1:L) for s in simulaciones_garch])
    p1 = plot(1:L, maximum(acfs; dims = 2); fillrange = minimum(acfs; dims = 2), color = NARANJA,
              fillalpha = 0.2, lw = 0, label = "GARCH: rango de 20 simulaciones",
              xlabel = "rezago τ (días)", ylabel = "ρ(τ) de |r|", legend = :topright,
              title = "Agrupamiento de volatilidad", ylims = (-0.05, 0.4))
    plot!(p1, 1:L, vec(mean(acfs; dims = 2)); color = NARANJA, lw = 2, label = "GARCH: promedio")
    plot!(p1, 1:L, acf_real; color = AZUL, lw = 2, label = "$activo")
    hline!(p1, [0]; color = :black, lw = 0.5, label = "")
    p2 = plot(xscale = :log10, yscale = :log10, xlims = (0.5, 30), ylims = (1e-4, 1.2),
              xticks = ([0.5, 1, 2, 5, 10, 20], ["0.5", "1", "2", "5", "10", "20"]),
              xlabel = "x (desviaciones estándar)", ylabel = "P(|z| ≥ x)", legend = :bottomleft,
              title = "Colas")
    zs = reduce(vcat, [abs.((s .- mean(s)) ./ std(s)) for s in simulaciones_garch])
    s, P = ccdf_empirica(zs)
    plot!(p2, s, P; color = NARANJA, lw = 2, label = "GARCH (20 simulaciones)")
    s, P = ccdf_empirica(abs.((r .- mean(r)) ./ std(r)))
    scatter!(p2, s, P; color = AZUL, ms = 2.5, msw = 0, label = "$activo")
    xr = exp10.(range(log10(0.5), log10(30), length = 200))
    plot!(p2, xr, 2 .* ccdf.(Normal(), xr); color = GRIS, ls = :dash, lw = 1.5, label = "Normal")
    plot(p1, p2; layout = (2, 1), size = (720, 620))
end

#%% hide
md"""
**Lo que GARCH sí reproduce.** Los retornos no tienen autocorrelación (por
construcción) y la autocorrelación de ``|r|`` en los primeros días es muy
parecida a la de los datos: el agrupamiento de volatilidad de corto plazo.

**Lo que no reproduce.** Compara las curvas lejos del origen: en el rezago
100, la autocorrelación de ``|r|`` es
$(round(comparacion_garch.acf100_real, digits = 3)) en los datos y
$(round(comparacion_garch.acf100_garch, digits = 3)) en GARCH. La del modelo
decae **exponencialmente**, como ``(a+b)^\tau``; la de los datos decae mucho
más despacio (memoria larga, sesión 1). Y en la gráfica de colas, el modelo
suele quedar por debajo de los datos en los eventos extremos: sus colas son
pesadas, pero en general menos que las reales (cambia de activo para
comprobarlo).

GARCH es un buen modelo de **corto plazo** y por eso se usa tanto para
pronosticar volatilidad y medir riesgo. Pero no es el mecanismo completo.
Modelos con memoria larga, innovaciones con colas pesadas o, en la semana 14,
agentes que interactúan, intentan cerrar esa brecha.
"""

#%% hide
comparacion_garch = (acf100_real = autocor(abs.(r), [100])[1],
                     acf100_garch = mean(autocor(abs.(s), [100])[1] for s in simulaciones_garch));

#%% hide
md"""
## 7. [Avanzado] De las trayectorias a las distribuciones: Fokker-Planck

En lugar de seguir trayectorias, podemos preguntar cómo evoluciona la
**densidad** ``p(x, t)`` del proceso. Para ``dX = A\,dt + B\,dW`` la respuesta
es la **ecuación de Fokker-Planck** (repaso E.A4):

```math
\frac{\partial p}{\partial t} = -\frac{\partial}{\partial x}\big[A(x)\,p\big]
+ \frac{1}{2}\frac{\partial^2}{\partial x^2}\big[B(x)^2\,p\big].
```

Para Ornstein-Uhlenbeck desde ``X_0 = x_0``, su solución es una normal con
media ``x_0 e^{-\theta t}`` y varianza ``\frac{\sigma^2}{2\theta}(1 - e^{-2\theta t})``.
La gráfica compara el histograma de 5000 trayectorias de Euler-Maruyama con
esa solución (usa ``\theta`` y ``\sigma`` de la sección 4).

Tiempo ``t``: $(@bind t_fp Slider(0.05:0.05:4.0, default = 0.5, show_value = true))
"""

#%% hide
let
    x0, Δt = 3.0, 0.01
    n = round(Int, t_fp / Δt)
    rng = Xoshiro(7)
    finales = [euler_maruyama_ref(x -> -θ_ou * x, x -> σ_ou, x0, Δt, n, rng)[end] for _ in 1:5000]
    m = x0 * exp(-θ_ou * t_fp)
    v = σ_ou^2 * (1 - exp(-2θ_ou * t_fp)) / (2θ_ou)
    xs = range(-3, 4, length = 400)
    p = histogram(finales; normalize = :pdf, bins = range(-3, 4, length = 90), color = AZUL,
                  alpha = 0.5, lw = 0, label = "5000 trayectorias (Euler-Maruyama)",
                  xlabel = "x", ylabel = "densidad", legend = :topleft,
                  title = "Distribución de X en t = $(t_fp)")
    plot!(p, xs, pdf.(Normal(m, sqrt(v)), xs); color = NARANJA, lw = 2,
          label = "Solución de Fokker-Planck")
end

#%% hide
md"""
**Black-Scholes es una ecuación de difusión.** El precio ``V(S, t)`` de una
opción sobre una acción que sigue un movimiento browniano geométrico cumple

```math
\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}
+ r_f S \frac{\partial V}{\partial S} - r_f V = 0,
```

con ``r_f`` la tasa libre de riesgo. Con el cambio de variables
``x = \ln S`` y un reescalamiento del tiempo y de ``V``, se convierte en la
ecuación del calor, ``\partial_\tau u = \partial_x^2 u``: la misma que describe
la difusión de la caminata aleatoria. Las notas formales de la sesión
desarrollan el argumento.
"""

#%% hide
md"""
## Resumen

| Proceso | Regla | Qué reproduce | Qué no |
|---|---|---|---|
| Caminata aleatoria | ``S_{n+1} = S_n \pm 1`` | Difusión: dispersión ``\propto \sqrt{n}`` | Colas, agrupamiento |
| Cadena de Markov (dos regímenes) | Volatilidad salta entre dos niveles | Agrupamiento de corto plazo, colas algo pesadas | Memoria larga |
| Ornstein-Uhlenbeck | ``dX = -\theta X\,dt + \sigma\,dW`` | Regreso a la media | (No es un modelo de retornos) |
| Browniano geométrico | ``dP = \mu P\,dt + \sigma P\,dW`` | Precios positivos; media ≠ mediana | Colas, agrupamiento |
| GARCH(1,1) | ``\sigma_t^2 = \omega + a r_{t-1}^2 + b \sigma_{t-1}^2`` | Agrupamiento de corto plazo | Memoria larga, colas tan pesadas como las reales |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `euler_maruyama(A, B, x0, Δt, n, rng)`, que devuelva un vector
con ``X_0, X_1, \dots, X_n`` simulados con la regla de la sección 4. `A` y `B`
son funciones. Usa `randn(rng)` una vez por paso, en orden.
"""

#%%
function euler_maruyama(A, B, x0, Δt, n, rng)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> euler_maruyama(x -> -2x, x -> 0.5, 1.0, 0.01, 100, Xoshiro(1)),
          euler_maruyama_ref(x -> -2x, x -> 0.5, 1.0, 0.01, 100, Xoshiro(1));
          pista = md"Crea `X = zeros(n + 1)` con `X[1] = x0`. En un ciclo `for k in 1:n`, calcula `X[k+1] = X[k] + A(X[k]) * Δt + B(X[k]) * sqrt(Δt) * randn(rng)`.")

#%% hide
md"""
**2.** Escribe `garch_varianza(r, ω, a, b)`, que reciba una serie de retornos
y devuelva el vector de varianzas condicionales ``\sigma_1^2, \dots, \sigma_n^2``
del GARCH(1,1). Empieza con la varianza de largo plazo,
``\sigma_1^2 = \omega/(1 - a - b)``.
"""

#%%
function garch_varianza(r, ω, a, b)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> garch_varianza(RETORNOS_PRUEBA, 1e-6, 0.1, 0.85),
          garch_varianza_ref(RETORNOS_PRUEBA, 1e-6, 0.1, 0.85);
          pista = md"Crea `σ2 = zeros(length(r))` con `σ2[1] = ω / (1 - a - b)`. Para `t` desde 2, usa `σ2[t] = ω + a * r[t-1]^2 + b * σ2[t-1]`.")

#%% hide
md"""
### Ruta básica: interpretación

3. En la caminata aleatoria con ``p = 1/2``, ¿cuánto hay que aumentar el
   número de pasos para que la dispersión típica se triplique?
4. En el modelo de dos regímenes con ``p = 0.02`` y ``q = 0.1``, ¿qué fracción
   del tiempo pasa el mercado en turbulencia y cuánto dura un episodio
   turbulento? ¿Por qué este modelo produce agrupamiento de volatilidad?
5. En el browniano geométrico, ¿por qué la mediana crece más despacio que el
   promedio? ¿Qué le pasa a un inversionista "típico" cuando la volatilidad
   aumenta y la tendencia ``\mu`` no cambia?
6. Cambia de activo en la sección 0 y compara la persistencia ``a + b`` y la
   vida media de los choques. ¿Qué mercado olvida más rápido sus crisis?
7. Resume en dos oraciones qué hechos estilizados reproduce GARCH y cuáles
   no, citando las gráficas de la sección 6.2.

### Ruta avanzada

8. Usa la ecuación de Fokker-Planck para demostrar que la varianza
   estacionaria de Ornstein-Uhlenbeck es ``\sigma^2/(2\theta)``, y verifícalo
   con la simulación de la sección 7 para ``t`` grande.
9. Para GARCH(1,1) con innovaciones normales, la curtosis es
   ``3(1 - (a+b)^2)/(1 - (a+b)^2 - 2a^2)`` cuando el denominador es positivo.
   Compárala con la curtosis de series simuladas para varios ``(a, b)``. ¿Qué
   pasa cerca de ``(a+b)^2 + 2a^2 = 1``?
10. Repite el laboratorio con innovaciones ``t`` de Student: ajusta también
    los grados de libertad ``\nu`` por máxima verosimilitud. ¿Mejoran las
    colas? ¿Y la memoria larga?
11. Verifica numéricamente que en el browniano geométrico ``\ln P_t`` es normal
    con media ``\ln P_0 + (\mu - \sigma^2/2)t`` y varianza ``\sigma^2 t``.

## Para la próxima sesión

- Lectura: Slanina (2014), secciones 2.1 y 2.2; Requião da Cunha (2021),
  secciones 2.3.5 y 3.1.
- Ruta avanzada: Requião da Cunha (2021), secciones 3.3 y 4.3; Mantegna &
  Stanley, capítulos 10 y 14.
- En la semana 4 estudiaremos el **crecimiento multiplicativo**. Una pista:
  la ecuación de GARCH, ``\sigma_t^2 = \omega + (a\varepsilon_{t-1}^2 + b)\,\sigma_{t-1}^2``,
  es un proceso de Kesten, y eso explica por qué genera colas de ley de
  potencias.
"""

#%% hide
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

#%%
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Statistics, StatsBase, Distributions,
          Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

#%%
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie o modelo
    const GRIS    = "#52514e"   # referencias
end;

#%%
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

#%%
"Simulación de Euler-Maruyama de dX = A(X) dt + B(X) dW; devuelve X₀, …, Xₙ."
function euler_maruyama_ref(A, B, x0, Δt, n, rng)
    X = zeros(n + 1)
    X[1] = x0
    for k in 1:n
        X[k+1] = X[k] + A(X[k]) * Δt + B(X[k]) * sqrt(Δt) * randn(rng)
    end
    X
end

#%%
"Mercado con dos regímenes de volatilidad que siguen una cadena de Markov."
function simular_regimenes(p, q; n = 3000, σ_calma = 0.008, factor = 3, rng = Xoshiro(8))
    estado = ones(Int, n)
    for t in 2:n
        u = rand(rng)
        estado[t] = estado[t-1] == 1 ? (u < p ? 2 : 1) : (u < q ? 1 : 2)
    end
    σ = ifelse.(estado .== 1, σ_calma, factor * σ_calma)
    (estado = estado, r = σ .* randn(rng, n))
end

#%%
"Agrupa índices consecutivos en pares (inicio, fin) para sombrear intervalos."
function rangos(idx)
    out = Float64[]
    inicio = idx[1]
    for i in 2:length(idx)
        if idx[i] != idx[i-1] + 1
            append!(out, [inicio - 0.5, idx[i-1] + 0.5]); inicio = idx[i]
        end
    end
    append!(out, [inicio - 0.5, idx[end] + 0.5])
end

#%%
"Trayectorias diarias de un browniano geométrico (solución exacta paso a paso)."
function simular_gbm(μ, σ; P0 = 100.0, años = 10, M = 2000, rng = Xoshiro(9))
    Δt = 1 / 252
    n = round(Int, años / Δt)
    incrementos = (μ - σ^2 / 2) * Δt .+ σ * sqrt(Δt) .* randn(rng, n, M)
    P = P0 .* exp.(vcat(zeros(1, M), cumsum(incrementos; dims = 1)))
    (t = (0:n) .* Δt, P = P)
end

#%%
"""
Simula un GARCH(1,1) con varianza de largo plazo `s2`, es decir,
ω = s2 (1 − a − b). Descarta 500 pasos iniciales.
"""
function simular_garch(n, a, b, s2; rng = Xoshiro(10))
    ω = s2 * (1 - a - b)
    σ2 = s2
    for _ in 1:500
        σ2 = ω + (a * randn(rng)^2 + b) * σ2
    end
    r, v = zeros(n), zeros(n)
    for t in 1:n
        v[t] = σ2
        r[t] = sqrt(σ2) * randn(rng)
        σ2 = ω + a * r[t]^2 + b * σ2
    end
    (r = r, σ2 = v)
end

#%%
"Varianzas condicionales de un GARCH(1,1) dada una serie de retornos."
function garch_varianza_ref(r, ω, a, b)
    σ2 = zeros(length(r))
    σ2[1] = ω / (1 - a - b)
    for t in 2:length(r)
        σ2[t] = ω + a * r[t-1]^2 + b * σ2[t-1]
    end
    σ2
end

#%%
"Log-verosimilitud normal del GARCH(1,1) con varianza de largo plazo s2."
function garch_loglik(x, a, b, s2)
    ω = s2 * (1 - a - b)
    σ2 = s2
    ℓ = 0.0
    for t in eachindex(x)
        ℓ -= 0.5 * (log(σ2) + x[t]^2 / σ2)
        σ2 = ω + a * x[t]^2 + b * σ2
    end
    ℓ
end

#%%
"""
Ajuste del GARCH(1,1) por máxima verosimilitud en una cuadrícula de (a, b),
con la varianza de largo plazo fijada en la varianza muestral.
"""
function ajustar_garch(r; as = 0.005:0.005:0.4, bs = 0.5:0.0025:0.9975)
    x = r .- mean(r)
    s2 = var(x)
    mejor = (ℓ = -Inf, a = NaN, b = NaN, s2 = s2)
    for a in as, b in bs
        a + b < 0.9995 || continue
        ℓ = garch_loglik(x, a, b, s2)
        ℓ > mejor.ℓ && (mejor = (ℓ = ℓ, a = a, b = b, s2 = s2))
    end
    mejor
end

#%%
"Muestra si la respuesta de un ejercicio es correcta."
function verificar(f, esperado; rtol = 1e-8, pista = md"")
    caja(tipo, titulo, texto) = Markdown.MD(Markdown.Admonition(tipo, titulo, [texto]))
    resultado = try
        f()
    catch e
        return caja("danger", "Tu código produce un error", md"$(sprint(showerror, e))")
    end
    if resultado === missing || (resultado isa AbstractArray && any(ismissing, resultado))
        caja("warning", "Pendiente", md"Reemplaza `missing` por tu código.")
    elseif resultado isa AbstractArray && esperado isa AbstractArray && size(resultado) == size(esperado) &&
           all(isapprox.(resultado, esperado; rtol, atol = 1e-12))
        caja("correct", "¡Correcto!", md"Bien hecho.")
    elseif resultado isa typeof(esperado) && all(isapprox.(resultado, esperado; rtol, atol = 1e-12))
        caja("correct", "¡Correcto!", md"Bien hecho.")
    else
        caja("danger", "Todavía no", pista)
    end
end

#%%
const RETORNOS_PRUEBA = 0.01 .* randn(Xoshiro(2026), 500);

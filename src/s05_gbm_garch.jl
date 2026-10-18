#%% hide
md"""
# Sesión 5: Langevin, GBM y GARCH

**Objetivos de la sesión**
1. Entender la simulación de **Euler-Maruyama** como una regla de
   actualización: "tendencia + ruido".
2. Simular la ecuación de **Langevin** (Ornstein-Uhlenbeck) y el **movimiento
   browniano geométrico**.
3. Construir un proceso con volatilidad que cambia en el tiempo (**GARCH**) y
   ajustarlo a datos reales.
4. Responder la pregunta del laboratorio: **¿reproduce GARCH los hechos
   estilizados de las sesiones 1 y 2?**

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s05_gbm_garch.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B3 (la derivada como
tasa de cambio). Para las notas formales: C.I4 (ecuaciones diferenciales),
E.I1 (autocorrelación), E.I4–I6 (GARCH, Langevin y crecimiento
multiplicativo). Ruta avanzada: E.A1–A4 (browniano, Itô, ecuación maestra y
Fokker-Planck), S.A1 (máxima verosimilitud).

**Conceptos transversales.** Se retoman los **retornos logarítmicos** (sesión
1) y la **caminata aleatoria** (sesión 4).
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: ¿qué proceso genera estas series?

En la sesión 4 vimos que un mercado que salta entre calma y turbulencia ya
produce agrupamiento de volatilidad, pero con memoria demasiado corta. Hoy
probamos otra idea: que la volatilidad **cambie suavemente** y dependa de lo
que pasó ayer. Compara otra vez los retornos reales con retornos normales
independientes de la misma volatilidad:

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
## 1. Repaso: la derivada como regla de actualización

Si ``x`` cambia a una tasa ``f(x)``, en un intervalo pequeño ``\Delta t``
cambia aproximadamente ``f(x)\,\Delta t``:

```math
x_{t + \Delta t} \approx x_t + f(x_t)\,\Delta t.
```

Así se simula una ecuación diferencial: aplicando la regla paso a paso. En la
sección 2 le agregaremos ruido.
"""

#%% hide
md"""
## 2. Tendencia + ruido: la ecuación de Langevin

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
la caminata aleatoria (sesión 4).

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
## 3. El movimiento browniano geométrico

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
## 4. GARCH: volatilidad con memoria

### 4.1 El modelo

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
### 4.2 Laboratorio: ¿reproduce GARCH los hechos estilizados?

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
más despacio (memoria larga, sesión 2). Y en la gráfica de colas, el modelo
suele quedar por debajo de los datos en los eventos extremos: sus colas son
pesadas, pero en general menos que las reales (cambia de activo para
comprobarlo).

GARCH es un buen modelo de **corto plazo** y por eso se usa tanto para
pronosticar volatilidad y medir riesgo. Pero no es el mecanismo completo.
Modelos con memoria larga, innovaciones con colas pesadas o, en la sesión 20,
agentes que interactúan, intentan cerrar esa brecha.
"""

#%% hide
comparacion_garch = (acf100_real = autocor(abs.(r), [100])[1],
                     acf100_garch = mean(autocor(abs.(s), [100])[1] for s in simulaciones_garch));

#%% hide
md"""
## 5. [Avanzado] De las trayectorias a las distribuciones: Fokker-Planck

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
esa solución (usa ``\theta`` y ``\sigma`` de la sección 2).

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
| Ornstein-Uhlenbeck | ``dX = -\theta X\,dt + \sigma\,dW`` | Regreso a la media | (No es un modelo de retornos) |
| Browniano geométrico | ``dP = \mu P\,dt + \sigma P\,dW`` | Precios positivos; media ≠ mediana | Colas, agrupamiento |
| GARCH(1,1) | ``\sigma_t^2 = \omega + a r_{t-1}^2 + b \sigma_{t-1}^2`` | Agrupamiento de corto plazo | Memoria larga, colas tan pesadas como las reales |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `euler_maruyama(A, B, x0, Δt, n, rng)`, que devuelva un vector
con ``X_0, X_1, \dots, X_n`` simulados con la regla de la sección 2. `A` y `B`
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

3. En el browniano geométrico, ¿por qué la mediana crece más despacio que el
   promedio? ¿Qué le pasa a un inversionista "típico" cuando la volatilidad
   aumenta y la tendencia ``\mu`` no cambia?
4. Cambia de activo en la sección 0 y compara la persistencia ``a + b`` y la
   vida media de los choques. ¿Qué mercado olvida más rápido sus crisis?
5. Resume en dos oraciones qué hechos estilizados reproduce GARCH y cuáles
   no, citando las gráficas de la sección 4.2.
6. En Ornstein-Uhlenbeck, ¿qué pasa con las fluctuaciones si duplicas la
   fuerza de regreso ``\theta``? ¿Y si duplicas el ruido ``\sigma``?


### Conceptos que vuelven

Preguntas cortas para repasar ideas de sesiones anteriores; se responden con
lo que ya sabes y con el laboratorio de hoy.

- **R1.** En el browniano geométrico los retornos logarítmicos son normales.
  ¿Por qué el precio nunca se vuelve negativo, aunque un retorno logarítmico
  sí pueda serlo? (Sesión 1.)
- **R2.** El ruido de Euler-Maruyama crece con ``\sqrt{\Delta t}``. ¿Qué tiene
  que ver con la dispersión ``\sqrt{n}`` de la caminata aleatoria? (Sesión 4.)

### Ruta avanzada

7. Usa la ecuación de Fokker-Planck para demostrar que la varianza
   estacionaria de Ornstein-Uhlenbeck es ``\sigma^2/(2\theta)``, y verifícalo
   con la simulación de la sección 5 para ``t`` grande.
8. Para GARCH(1,1) con innovaciones normales, la curtosis es
   ``3(1 - (a+b)^2)/(1 - (a+b)^2 - 2a^2)`` cuando el denominador es positivo.
   Compárala con la curtosis de series simuladas para varios ``(a, b)``. ¿Qué
   pasa cerca de ``(a+b)^2 + 2a^2 = 1``?
9. Repite el laboratorio con innovaciones ``t`` de Student: ajusta también
   los grados de libertad ``\nu`` por máxima verosimilitud. ¿Mejoran las
   colas? ¿Y la memoria larga?
10. Verifica numéricamente que en el browniano geométrico ``\ln P_t`` es normal
    con media ``\ln P_0 + (\mu - \sigma^2/2)t`` y varianza ``\sigma^2 t``.

## Para la próxima sesión

- Lectura: Requião da Cunha (2021), sección 2.3.5. Complementaria: Slanina
  (2014), sección 2.2. Para la ruta avanzada, Mantegna y Stanley (1999),
  capítulo 10; Requião da Cunha (2021), secciones 3.3 y 4.3; de consulta,
  Gardiner (2021).
- En la sesión 6 estudiaremos el **crecimiento multiplicativo**. Una pista:
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

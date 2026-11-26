### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 72045a1e-ca99-11f1-90c6-d143fcae313a
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Dates, Statistics, StatsBase, Distributions,
          Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 72045938-ca99-11f1-997a-8789a555af3c
md"""
# Sesión 20: Mercados artificiales y expectativas heterogéneas

**Objetivos de la sesión**
1. Comprender el **modelo de Brock y Hommes (1998)**: cómo la competencia
   evolutiva entre **fundamentalistas** y **seguidores de tendencia** genera
   burbujas y cracs endógenos sin necesidad de choques externos.
2. Explorar los modelos estocásticos de **fundamentalistas y chartistas con
   contagio** (Lux y Marchesi, 1999; Giardina y Bouchaud, 2003; Franke y
   Westerhoff, 2012) y observar cómo los cambios de régimen en la población
   crean ráfagas de volatilidad.
3. **Cerrar el ciclo del curso** aplicando las herramientas de diagnóstico de las
   **sesiones 1 y 2** (CCDF de colas, curtosis, autocorrelación de ``r_t`` y de
   ``|r_t|``) para comparar datos reales del S&P 500, caminatas aleatorias y
   mercados artificiales.
4. **[Avanzado]** Analizar la pérdida de estabilidad del equilibrio fundamental
   en ``\beta_c = \frac{1}{C}\ln\!\left(\frac{R}{g - R}\right)``, construir el
   **diagrama de bifurcación** y medir el **exponente de Lyapunov**.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s20_mercados_artificiales.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B3 (estacionariedad y
promedios temporales), S.B5 (pruebas de permutación), E.I1 (autocorrelación),
E.I4 (ARCH y GARCH), C.I5 (puntos fijos y estabilidad en 1D), S.I2 (Ljung-Box),
S.I4 (datos dependientes). Ruta avanzada: L.A3 (jacobiano y estabilidad),
C.A5 (mapas discretos, bifurcaciones y caos), S.A3 (calibración por momentos
simulados), S.A4 (memoria larga).

**Conceptos transversales.** Se retoman los **retornos logarítmicos** (sesión
1), el **logit** (sesión 13) y los **puntos fijos y su estabilidad** (sesiones
12 y 14).
"""

# ╔═╡ 7204595e-ca99-11f1-8c3f-e7271e4f8224
TableOfContents(title = "Contenido")

# ╔═╡ 7204596a-ca99-11f1-a174-7107e6dc34fd
md"""
## 0. El fenómeno: ¿por qué los mercados alternan entre calma y tormenta?

En las **sesiones 1 y 2** abrimos el curso midiendo tres hechos estilizados universales
en los retornos financieros ``r_t = \ln(P_t/P_{t-1})``:
1. **Sin autocorrelación lineal:** ``\hat\rho_r(\tau) \approx 0``.
2. **Colas pesadas:** curtosis en exceso ``\kappa \gg 0`` y eventos extremos de
   varias desviaciones estándar.
3. **Agrupamiento de volatilidad:** ``\hat\rho_{|r|}(\tau) > 0`` durante semanas
   o meses.

Una caminata aleatoria gaussiana (izquierda) tiene la misma varianza todos los
días. En cambio, cuando en un mercado conviven **fundamentalistas** (que anclan
el precio a su valor fundamental) y **chartistas** (que siguen la tendencia y se
contagian entre sí), el propio cambio en la fracción de especuladores produce
etapas tranquilas seguidas de **burbujas y tormentas de volatilidad** (derecha):
"""

# ╔═╡ 72045974-ca99-11f1-b0b3-7b599481ad7c
md"""
## 1. El modelo de Brock y Hommes (1998): selección evolutiva de predictores

### Las reglas en palabras

Sea ``x_t = p_t - p^*`` la desviación del precio del activo respecto a su **valor
fundamental** ``p^*``, y ``R = 1 + r > 1`` el rendimiento bruto libre de riesgo.
El precio de hoy vacía el mercado promediando los pronósticos de los agentes:

```math
R\,x_t = \sum_{h=1}^H n_{h,t}\,f_{h,t} ,
```

donde ``n_{h,t}`` es la fracción de agentes que usa la regla predictiva ``h`` y
``f_{h,t} = \mathbb{E}_{h,t}[x_{t+1}]`` es su pronóstico para mañana.

En el caso clásico con ``H = 2`` reglas:
1. **Fundamentalistas (``h = 1``):** pronostican que la desviación volverá a
   cero, ``f_{1,t} = 0``, pero pagan un **costo de información** ``C \ge 0`` por
   investigar los fundamentales.
2. **Chartistas / seguidores de tendencia (``h = 2``):** extrapolan la última
   desviación con fuerza ``g > 0``, ``f_{2,t} = g\,x_{t-1}``, sin costo de
   información (``C_2 = 0``).

Al cierre de cada periodo, cada regla acumula una **ganancia realizada**
``\pi_{h,t-1} = (x_{t-1} - R\,x_{t-2})(f_{h,t-2} - R\,x_{t-2}) - C_h``, y los
agentes eligen qué regla usar mediante la probabilidad **logit** (sesión 13) con
**intensidad de elección** ``\beta \ge 0``:

```math
n_{2,t} = \frac{1}{1 + \exp\!\big(\beta\,(\pi_{1,t-1} - \pi_{2,t-1})\big)} .
```

### Laboratorio interactivo: calma, ciclos y burbujas en Brock-Hommes

Mueve la **intensidad de elección** ``\beta`` y la **extrapolación chartista**
``g`` en el régimen canónico ``R < g < R^2`` (con ``R = 1.10``, ``R^2 = 1.21`` y
costo fundamentalista ``C = 1.0``):

- **Intensidad de elección ``\beta``:** $(@bind beta_bh Slider(1.0:0.2:4.6, default = 3.8, show_value = true))
- **Fuerza de tendencia ``g``:** $(@bind g_bh Slider(1.12:0.02:1.20, default = 1.20, show_value = true))
"""

# ╔═╡ 72045988-ca99-11f1-b3d7-7134e184652c
let
    R = 1.10
    C = 1.0
    n2_eq = 1.0 / (1.0 + exp(-beta_bh * C))
    mult = g_bh * n2_eq / R
    bc = g_bh > R ? log(R / (g_bh - R)) / C : Inf
    m = round(mult, digits = 2)
    regimen = if mult < 1.0
        "**Equilibrio fundamental estable** (`g·n₂*/R = $m < 1`). Cualquier desviación se amortigua hacia `x* = 0`."
    else
        "**Equilibrio fundamental inestable** (`g·n₂*/R = $m > 1`, porque `β > β_c`). En la calma cerca de `x = 0` los agentes ahorran el costo `C` pasándose al chartismo; eso infla una burbuja hasta que `x_t` crece y hace rentable volver a ser fundamentalista."
    end
    Markdown.parse("""
    > **Diagnóstico con β = $(beta_bh) y g = $(g_bh).** En el punto fijo fundamental
    > `x* = 0`, la fracción de chartistas es `n₂* = $(round(n2_eq, digits = 3))` y el
    > umbral crítico es `β_c ≈ $(round(bc, digits = 2))`.
    >
    > $regimen
    """)
end

# ╔═╡ 72045992-ca99-11f1-a3c8-9d1d92558d24
md"""
## 2. Fundamentalistas, chartistas y contagio estocástico

En los mercados reales, además de la selección por ganancias hay **contagio de
opinión** (como en el modelo de Kirman de la sesión 15) y órdenes ruidosas. Lux
y Marchesi (1999) y Giardina y Bouchaud (2003) mostraron con modelos de agentes
que el ir y venir entre fundamentalismo y chartismo basta para generar colas
pesadas y agrupamiento de volatilidad. Aquí simulamos la versión compacta de
**Franke y Westerhoff (2012)** (modelo DCA-HPM), con los parámetros que ellos
estimaron para el S&P 500:

- El log-precio ``p_t`` responde al exceso de demanda de una fracción
  ``n^f_{t-1}`` de **fundamentalistas** y una fracción
  ``z_{t-1} = n^c_{t-1} = 1 - n^f_{t-1}`` de **chartistas**:

```math
p_t = p_{t-1} + \mu\bigl(n^f_{t-1} d^f_{t-1} + n^c_{t-1} d^c_{t-1}\bigr), \qquad
d^f_t = \phi(p^* - p_t) + \sigma_f\varepsilon^f_t, \qquad
d^c_t = \chi(p_t - p_{t-1}) + \sigma_c\varepsilon^c_t .
```

  Los fundamentalistas empujan el precio hacia ``p^*`` y los chartistas siguen
  la tendencia, con una demanda más ruidosa (``\sigma_c > \sigma_f``).
- La fracción de fundamentalistas es logit en el **atractivo relativo** del
  fundamentalismo:

```math
n^f_t = \frac{1}{1 + e^{-\beta a_{t-1}}}, \qquad
a_t = \alpha_0 + \alpha_n\,(n^f_t - n^c_t) + \alpha_p\,(p_t - p^*)^2 .
```

  ``\alpha_0`` es la predisposición, ``\alpha_n > 0`` el **contagio** (la
  estrategia mayoritaria atrae a más agentes) y ``\alpha_p > 0`` el **temor al
  desalineamiento**: cuanto más lejos está el precio del fundamental, más
  atractivo resulta apostar a que regrese.

Como ``\mathrm{Var}(r_t \mid t-1) = \mu^2\bigl[(n^f_{t-1})^2\sigma_f^2 + (n^c_{t-1})^2\sigma_c^2\bigr]``,
la volatilidad sube cuando dominan los chartistas. Explora cómo cambia la
dinámica con la **intensidad de elección** ``\beta`` (estimada: ``\beta = 1``)
y el **ruido chartista** ``\sigma_c`` (estimado: ``\sigma_c = 2.087``, con
``\sigma_f = 0.758``):

- **Intensidad de elección ``\beta``:** $(@bind beta_lm Slider(0.5:0.25:2.0, default = 1.0, show_value = true))
- **Ruido chartista ``\sigma_c``:** $(@bind sigma_c_lm Slider([1.0, 1.5, 2.087, 2.5, 3.0], default = 2.087, show_value = true))
"""

# ╔═╡ 720459a6-ca99-11f1-beb4-555e7bdce368
md"""
Cada vez que la curva naranja (fracción de chartistas ``z_t``) sube, el mercado
pierde su ancla fundamentalista y aparece un **racimo de alta volatilidad**; el
contagio ``\alpha_n`` mantiene a la mayoría en su estrategia durante un tiempo.
Cuando el desalineamiento acumulado hace atractivo el fundamentalismo
(``z_t`` baja), regresa la calma.
"""

# ╔═╡ 720459a6-ca99-11f1-b086-ef4ff44161c1
md"""
## 3. Cierre del ciclo: el análisis de las sesiones 1 y 2 en series simuladas y reales

Esta sección es el punto de partida del **taller para casa del Módulo III**
(`Code/taller3_cierre_ciclo.jl`).
En clase basta con leer la tabla; en el taller se repite el análisis con otros
modelos y parámetros.

Pongamos a prueba lo aprendido en las 20 sesiones aplicando exactamente las
gráficas de diagnóstico de las **sesiones 1 y 2** a tres series de retornos
estandarizados ``z_t = (r_t - \hat\mu)/\hat\sigma``:

1. **Datos reales (S&P 500)** de `Data/precios.csv`.
2. **Caminata aleatoria gaussiana** (movimiento browniano geométrico, sesión 5).
3. **Mercado artificial heterogéneo** (Franke y Westerhoff, sección 2), con el
   mismo número de días que el S&P 500.
"""

# ╔═╡ 720459ba-ca99-11f1-8138-ff78ddc3e4d7
md"""
**Lectura de la tabla.** En una normal sólo el 0.27% de los días tiene
``\lvert z \rvert > 3``. El mercado artificial reproduce los dos
hechos estilizados que la caminata aleatoria no puede: colas más pesadas que la
normal y agrupamiento de volatilidad que persiste decenas de días, sin memoria
lineal apreciable en ``r_t``. Pero su curtosis queda muy por debajo de la del
S&P 500: un modelo tan sencillo explica el mecanismo, no todos los números.
Calibrar un modelo de agentes para que reproduzca estos momentos es justo lo
que hacen Franke y Westerhoff (2012) con el método de momentos simulados
(repaso S.A3).
"""

# ╔═╡ 720459c4-ca99-11f1-9ac5-f3d362c08707
md"""
## 4. [Avanzado] Bifurcaciones, estabilidad y caos en Brock-Hommes

¿Por qué el equilibrio fundamental ``x^* = 0`` pierde estabilidad exactamente al
cruzar un umbral crítico ``\beta_c`` (repaso L.A3 y C.A5)?

En el punto fijo fundamental ``x_t = 0``, las ganancias estacionarias son
``\pi_1^* = -C`` (los fundamentalistas pagan el costo de información) y
``\pi_2^* = 0`` (los chartistas no pagan costo). Por tanto, la fracción
estacionaria de chartistas en el origen es:

```math
n_2^*(\beta) = \frac{1}{1 + e^{-\beta C}} .
```

Linealizando ``R\,x_t = n_{2,t}\,g\,x_{t-1}`` alrededor de ``x^* = 0``, el
multiplicador local es ``\lambda_1 = g\,n_2^*(\beta) / R``. Cuando
``R < g < 2R`` y ``C > 0``, el punto fijo fundamental se vuelve inestable
(``\lambda_1 > 1``) cuando ``\beta`` supera el **valor crítico de bifurcación**:

```math
\beta_c = \frac{1}{C}\,\ln\!\left(\frac{R}{g - R}\right) .
```

La figura izquierda muestra el **diagrama de bifurcación** de ``x_t`` al barrer
``\beta`` (descartando el transitorio; como empezamos en ``x_0 > 0`` sólo se
ve la rama positiva, y la negativa es su reflejo), y la figura derecha muestra
el **exponente de Lyapunov** ``\lambda(\beta)`` (cuando ``\lambda > 0``, hay
caos determinista). Tras la bifurcación en ``\beta_c``, los dos equilibrios no
fundamentales ``\pm x^\ast`` pierden a su vez estabilidad (bifurcación de
Neimark-Sacker, ``\lambda \approx 0``) y la dinámica se vuelve cuasiperiódica y
luego caótica:
"""

# ╔═╡ 720459ce-ca99-11f1-9916-f33b28a01215
md"""
## Resumen

| Modelo / Concepto | Ecuación o mecanismo | Qué explica |
|---|---|---|
| **Vaciado de mercado (Brock-Hommes)** | ``R\,x_t = \sum_{h=1}^H n_{h,t}\,f_{h,t}`` | El precio refleja el pronóstico medio ponderado por las cuotas de cada estrategia |
| **Elección evolutiva logit** | ``n_{h,t} = e^{\beta \pi_{h,t-1}} / \sum_k e^{\beta \pi_{k,t-1}}`` | Los agentes migran hacia las reglas que ganaron más recientemente |
| **Umbral crítico ``\beta_c``** | ``\beta_c = \frac{1}{C}\ln\!\left(\frac{R}{g - R}\right)`` para ``R < g < 2R`` | La calma hace que los agentes ahorren el costo ``C`` pasándose al chartismo y desestabilizando ``x^* = 0`` |
| **Contagio y cambio de régimen (Franke-Westerhoff)** | ``\sigma_t^2 = \mu^2\bigl[(1-z_{t-1})^2\sigma_f^2 + z_{t-1}^2\sigma_c^2\bigr]`` con ``\sigma_c > \sigma_f`` | Las olas de chartismo ``z_t`` generan colas pesadas y agrupamiento de volatilidad sin memoria lineal en ``r_t`` |
"""

# ╔═╡ 720459d8-ca99-11f1-aa7e-23c97971346b
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `fraccion_logit_dos(pi_1, pi_2, beta)`, que devuelva la fracción
``n_2`` de agentes que eligen la estrategia 2 cuando las ganancias recientes son
`pi_1` y `pi_2` y la intensidad de elección es `beta`:

```math
n_2 = \frac{1}{1 + \exp\!\big(\beta\,(\pi_1 - \pi_2)\big)} .
```
"""

# ╔═╡ 720459d8-ca99-11f1-aa52-b36c16d3aba7
function fraccion_logit_dos(pi_1, pi_2, beta)
    missing  # reemplaza esta línea
end

# ╔═╡ 720459e2-ca99-11f1-b4d4-7dcd4aeface8
md"""
**2.** *Si alcanza (margen de la clase).* Escribe `paso_precio_bh(x_ant, n2, g, R)`, que calcule la nueva
desviación del precio ``x_t`` en el modelo de Brock-Hommes con fundamentalistas
(``f_{1,t} = 0``) y seguidores de tendencia (``f_{2,t} = g\,x_{t-1}``) a partir
de la ecuación de vaciado de mercado ``R\,x_t = n_2\,g\,x_{t-1}``.
"""

# ╔═╡ 720459ec-ca99-11f1-bf99-27d5104f9942
function paso_precio_bh(x_ant, n2, g, R)
    missing  # reemplaza esta línea
end

# ╔═╡ 720459f6-ca99-11f1-b8d3-b56deb6e8ab9
md"""
### Ruta básica: interpretación

3. En el laboratorio de la sección 1, ¿por qué un costo de información ``C > 0``
   para ser fundamentalista es indispensable para que el aumento de ``\beta``
   desestabilice el equilibrio fundamental cuando ``R < g < 2R``? ¿Qué ocurre si
   ``C = 0``?
4. En el laboratorio de la sección 2, sube la intensidad de elección ``\beta``
   y después el ruido chartista ``\sigma_c``. ¿Qué pasa con las olas de
   chartismo ``z_t`` y con los racimos de volatilidad? ¿Por qué un mercado
   dominado por fundamentalistas es más tranquilo?
5. En la tabla de la sección 3, compara la caminata aleatoria gaussiana y el
   mercado artificial frente al S&P 500. ¿Qué dos hechos estilizados captura el
   mercado artificial que la caminata aleatoria no puede reproducir? ¿Cuál no
   alcanza a reproducir bien?
"""

# ╔═╡ 720459f6-ca99-11f1-92f7-07bd0433376a
md"""

### Conceptos que vuelven

Preguntas cortas para repasar ideas de sesiones anteriores; se responden con
lo que ya sabes y con el laboratorio de hoy.

- **R1.** En la sección 3 comparamos retornos logarítmicos estandarizados de
  las tres series. ¿Por qué comparamos retornos y no precios? (Sesión 1.)
- **R2.** La fracción de chartistas de Brock-Hommes es una probabilidad logit.
  ¿Qué representa ``\beta`` en el mercado y qué pasa si es muy grande? (Sesión
  13.)
- **R3.** En el diagnóstico de la sección 1, el equilibrio fundamental pierde
  estabilidad cuando el multiplicador ``g\,n_2^\ast/R`` supera 1. ¿Qué tiene
  en común con la estabilidad de los puntos fijos que vimos antes? (Sesiones
  12 y 14.)

### Ruta avanzada: completa el código

Los ejercicios 6 a 8 son opcionales.

**6.** Escribe `beta_critico_bh(g, R, C)`, que devuelva la intensidad
de elección crítica ``\beta_c = \frac{1}{C}\ln\!\left(\frac{R}{g - R}\right)``
en la que el equilibrio fundamental ``x^* = 0`` pierde estabilidad cuando
``R < g < 2R`` y ``C > 0``.
"""

# ╔═╡ 72045a02-ca99-11f1-be6b-bf846fdc646b
function beta_critico_bh(g, R, C)
    missing  # reemplaza esta línea
end

# ╔═╡ 72045a0a-ca99-11f1-84bc-a56a46d3317c
md"""
**7.** Escribe `multiplicador_origen_bh(beta, g, R, C)`, que calcule
el eigenvalor linealizado ``\lambda_1 = \frac{g\,n_2^*(\beta)}{R}`` en el punto
fijo fundamental ``x^* = 0``, donde ``n_2^*(\beta) = \frac{1}{1 + e^{-\beta C}}``.
"""

# ╔═╡ 72045a0a-ca99-11f1-a31a-a525bc67a70f
function multiplicador_origen_bh(beta, g, R, C)
    missing  # reemplaza esta línea
end

# ╔═╡ 72045a14-ca99-11f1-a16b-334b6c296970
md"""
**8.** En el diagrama de bifurcación de la sección 4, ¿qué significa
   económicamente que el exponente de Lyapunov ``\lambda(\beta)`` sea positivo
   para valores altos de ``\beta`` aunque el modelo no tenga ningún término de
   ruido aleatorio?

## Para la próxima sesión (proyecto final)

- Lectura: LeBaron (2006). Complementaria: Hommes (2013), capítulo 1; Lux y
  Marchesi (1999). Para la ruta avanzada, Hommes (2013), capítulos 2 y 6; Slanina (2014),
  sección 3.2; Giardina y Bouchaud (2003); Brock y Hommes (1998).
- **Entrega**: taller para casa del Módulo III (cierre del ciclo con los modelos
  de agentes frente a los datos), en `Code/taller3_cierre_ciclo.jl`.
- Presentación de los proyectos finales: réplica y extensión de un artículo de
  econofísica usando las herramientas de simulación y diagnóstico estadístico
  del curso (etapas y rúbrica en `programa/proyecto_final.md`).
- Repaso recomendado para el proyecto: S.B4 (transformación logarítmica), S.I1
  (bootstrap y bootstrap por bloques), S.I5 (supuestos de regresión) y S.A3
  (calibración de modelos de agentes por momentos simulados).
"""

# ╔═╡ 72045a1e-ca99-11f1-9475-0bf39586ed6a
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 72045a28-ca99-11f1-adfc-bdcade772936
begin
    const AZUL       = "#2a78d6"
    const NARANJA    = "#eb6834"
    const GRIS       = "#52514e"
    const GRIS_CLARO = "#c9c7c0"
    const BLANCO     = "#ffffff"
end;

# ╔═╡ 72045a28-ca99-11f1-bb79-b5f72d9c6406
"""
Simula el modelo determinista de Brock y Hommes (1998) con fundamentalistas
(`f_1 = 0`, costo `C`) y seguidores de tendencia (`f_2 = g * x_{t-1}`, costo `0`)
en el régimen canónico `R < g < R^2`.
"""
function simular_brock_hommes(beta, g; R = 1.10, C = 1.0, pasos = 200, calent = 500, x0 = 0.25)
    T = calent + pasos
    x = zeros(Float64, T + 3)
    n2 = fill(0.5, T + 3)
    x[1:3] .= (x0, x0 * 1.02, x0 * 0.98)

    for t in 4:(T + 3)
        ret_real = x[t - 1] - R * x[t - 2]
        pi1 = -R * x[t - 2] * ret_real - C
        pi2 = (g * x[t - 3] - R * x[t - 2]) * ret_real
        diff_pi = clamp(beta * (pi1 - pi2), -40.0, 40.0)
        n2[t] = 1.0 / (1.0 + exp(diff_pi))
        x[t] = (n2[t] * g * x[t - 1]) / R
    end

    idx = (calent + 4):(T + 3)
    (x = x[idx], n2 = n2[idx])
end

# ╔═╡ 7204597e-ca99-11f1-8553-21fb958df804
let
    R = 1.10
    C = 1.0
    sim = simular_brock_hommes(beta_bh, g_bh; R = R, C = C, pasos = 180, calent = 500)
    ts = 1:length(sim.x)

    p1 = plot(ts, sim.x; color = AZUL, lw = 1.8, label = "Desviación x_t = p_t − p*",
              xlabel = "periodo t", ylabel = "desviación x_t",
              legend = :outerbottom, legend_columns = 2,
              title = "Trayectoria de precios y cuota chartista")
    hline!(p1, [0.0]; color = GRIS_CLARO, ls = :dash, lw = 1.2, label = "Fundamental (x = 0)")
    plot!(p1, ts, sim.n2; color = NARANJA, lw = 1.4, ls = :dot,
          label = "Fracción chartista n_{2,t}")

    # rango mínimo de los ejes: si el precio converge a x* = 0 la escala no colapsa
    lo, hi = min(0.0, minimum(sim.x)), max(0.3, maximum(sim.x))
    lims = (lo - 0.05(hi - lo), hi + 0.05(hi - lo))
    p2 = scatter(sim.x[1:end-1], sim.x[2:end]; color = AZUL, alpha = 0.65, ms = 3.5, msw = 0,
                 xlabel = "desviación pasada x_{t-1}", ylabel = "desviación actual x_t",
                 xlims = lims, ylims = lims, title = "Retrato en el plano (x_{t-1}, x_t)")
    plot!(p2, sim.x[1:end-1], sim.x[2:end]; color = AZUL, alpha = 0.25, lw = 1.0)
    scatter!(p2, [0.0], [0.0]; color = NARANJA, marker = :xcross, ms = 6, msw = 2)

    plot(p1, p2; layout = (1, 2), size = (880, 370), bottom_margin = 4Plots.mm)
end

# ╔═╡ 72045a34-ca99-11f1-a139-45cd87f281da
"""
Simula el modelo DCA-HPM de Franke y Westerhoff (2012): fundamentalistas y
chartistas con contagio (herding), predisposición y temor al desalineamiento.
Los valores por omisión son los que los autores estimaron para el S&P 500.
Devuelve el log-precio `p`, los retornos `r` y la fracción de chartistas `z`.
"""
function simular_franke_westerhoff(; pasos = 2500, calent = 1000, beta = 1.0, mu = 0.01,
                                   phi = 0.12, chi = 1.5, sigma_f = 0.758, sigma_c = 2.087,
                                   alpha_0 = -0.327, alpha_n = 1.79, alpha_p = 18.43,
                                   p_fund = 0.0, rng = Xoshiro(42))
    T = calent + pasos + 1
    p = fill(p_fund, T)
    nf = fill(0.5, T)
    for t in 3:T
        d_f = phi * (p_fund - p[t - 1]) + sigma_f * randn(rng)
        d_c = chi * (p[t - 1] - p[t - 2]) + sigma_c * randn(rng)
        p[t] = p[t - 1] + mu * (nf[t - 1] * d_f + (1 - nf[t - 1]) * d_c)
        # Atractivo relativo del fundamentalismo en t - 1: contagio + predisposición + desalineamiento
        a = alpha_n * (2 * nf[t - 1] - 1) + alpha_0 + alpha_p * (p[t - 1] - p_fund)^2
        nf[t] = 1 / (1 + exp(-beta * a))
    end
    idx = (calent + 2):T
    (p = p[idx], r = p[idx] .- p[idx .- 1], z = 1 .- nf[idx])
end

# ╔═╡ 7204599c-ca99-11f1-acc9-b3c53724b071
let
    sim = simular_franke_westerhoff(; pasos = 1500, beta = beta_lm, sigma_c = sigma_c_lm, rng = Xoshiro(42))
    t = 1:length(sim.r)

    p1 = plot(t, sim.p; color = AZUL, lw = 1.4, label = "Log-precio p_t",
              xlabel = "día t", ylabel = "log-precio p_t",
              legend = :outerbottom, legend_columns = 2,
              title = "Precio simulado y valor fundamental")
    hline!(p1, [0.0]; color = GRIS, ls = :dash, lw = 1.3, label = "Fundamental p* = 0")

    escala = round(maximum(abs, sim.r), digits = 2)
    p2 = plot(t, sim.r; color = AZUL, lw = 0.7, alpha = 0.85, label = "Retorno r_t",
              xlabel = "día t", ylabel = "retorno / fracción z_t",
              legend = :outerbottom, legend_columns = 2,
              title = "Retornos r_t y olas de chartismo z_t")
    plot!(p2, t, escala .* sim.z; color = NARANJA, lw = 1.9,
          label = "z_t (×$escala)")

    plot(p1, p2; layout = (1, 2), size = (880, 360), bottom_margin = 4Plots.mm)
end

# ╔═╡ 72045a3c-ca99-11f1-9577-13a48baf26b4
"Ejemplo inicial de la sección 0: caminata aleatoria frente a mercado artificial heterogéneo."
ejemplo_calma_tormenta = let
    rng = Xoshiro(14)
    sim = simular_franke_westerhoff(; pasos = 900, rng = rng)
    z_abm = (sim.r .- mean(sim.r)) ./ std(sim.r)
    z_gbm = randn(Xoshiro(99), length(z_abm))
    (r_gbm = z_gbm, r_abm = z_abm)
end;

# ╔═╡ 72045974-ca99-11f1-a107-a11749a7d068
let
    s = ejemplo_calma_tormenta
    t = 1:length(s.r_gbm)
    p1 = plot(t, s.r_gbm; color = GRIS_CLARO, lw = 0.8,
              xlabel = "día t", ylabel = "retorno estandarizado z_t", ylims = (-6.5, 6.5),
              title = "Caminata aleatoria: volatilidad constante")
    hline!(p1, [-3, 3]; color = GRIS, ls = :dot, lw = 1)

    p2 = plot(t, s.r_abm; color = AZUL, lw = 0.8,
              xlabel = "día t", ylabel = "retorno estandarizado z_t", ylims = (-6.5, 6.5),
              title = "Mercado heterogéneo: calma y tormentas")
    hline!(p2, [-3, 3]; color = NARANJA, ls = :dot, lw = 1.2)

    plot(p1, p2; layout = (1, 2), size = (880, 320), bottom_margin = 4Plots.mm)
end

# ╔═╡ 72045a3c-ca99-11f1-8aa3-7dd8c0801684
"Calcula la CCDF empírica de |z|, la curtosis en exceso y las autocorrelaciones de z y |z|."
function resumen_hechos_estilizados(r; max_lag = 25)
    z = (r .- mean(r)) ./ std(r)
    abs_z = abs.(z)
    xs = sort(abs_z)
    n = length(xs)
    indices = unique(round.(Int, range(1, n - 1, length = min(250, n - 1))))
    x_ccdf = xs[indices]
    y_ccdf = [(n - i) / n for i in indices]
    (z = z,
     x_ccdf = x_ccdf,
     y_ccdf = y_ccdf,
     curtosis = kurtosis(z),
     extremos3 = count(>(3.0), abs_z) / n,
     acf_r = autocor(z, 1:max_lag),
     acf_abs = autocor(abs_z, 1:max_lag))
end

# ╔═╡ 72045a46-ca99-11f1-9d5f-c51ef7d94e39
"Comparación de la sección 3 entre el S&P 500 real, una caminata aleatoria y el mercado artificial."
comparacion_semana1 = let
    ruta = joinpath(@__DIR__, "..", "Data", "precios.csv")
    df = CSV.read(ruta, DataFrame)
    sp = filter(row -> row.activo == "S&P 500", df)
    sort!(sp, :fecha)
    r_sp = diff(log.(sp.cierre))

    N = length(r_sp)
    r_gbm = randn(Xoshiro(101), N)
    r_abm = simular_franke_westerhoff(; pasos = N, rng = Xoshiro(2026)).r

    (sp = resumen_hechos_estilizados(r_sp),
     gbm = resumen_hechos_estilizados(r_gbm),
     abm = resumen_hechos_estilizados(r_abm))
end;

# ╔═╡ 720459b0-ca99-11f1-b03b-9384321db7ce
let
    c = comparacion_semana1
    x_grid = range(0.2, 5.5, length = 150)
    ccdf_norm = 2.0 .* ccdf.(Normal(0, 1), x_grid)

    p1 = plot(x_grid, ccdf_norm; color = GRIS_CLARO, lw = 2.0, ls = :dash,
              yscale = :log10, ylims = (1e-4, 1.0), xlims = (0.0, 5.5),
              xlabel = "umbral x (desviaciones estándar)", ylabel = "P(|z| > x)",
              legend = :outerbottom, legend_columns = 2,
              label = "Normal estándar",
              title = "1. Colas: CCDF de |z_t|")
    plot!(p1, c.sp.x_ccdf, c.sp.y_ccdf; color = GRIS, lw = 2.2, label = "S&P 500 real")
    plot!(p1, c.gbm.x_ccdf, c.gbm.y_ccdf; color = NARANJA, lw = 1.8, ls = :dot, label = "Caminata aleatoria")
    plot!(p1, c.abm.x_ccdf, c.abm.y_ccdf; color = AZUL, lw = 2.2, label = "Mercado artificial")

    lags = 1:25
    banda = 1.96 / sqrt(length(c.sp.z))
    p2 = plot(lags, c.sp.acf_r; color = GRIS, lw = 1.8, marker = :circle, ms = 3, msw = 0,
              ylims = (-0.15, 0.25), xlabel = "rezago τ (días)", ylabel = "ACF de r_t",
              legend = :outerbottom, legend_columns = 2,
              label = "S&P 500 real",
              title = "2. Memoria lineal: ACF de r_t")
    plot!(p2, lags, c.abm.acf_r; color = AZUL, lw = 1.8, marker = :square, ms = 3, msw = 0,
          label = "Mercado artificial")
    hline!(p2, [-banda, banda]; color = NARANJA, ls = :dash, lw = 1.2, label = "Banda ±1.96/√T")
    hline!(p2, [0.0]; color = GRIS_CLARO, lw = 1, label = "")

    p3 = plot(lags, c.sp.acf_abs; color = GRIS, lw = 2.2, marker = :circle, ms = 3.5, msw = 0,
              ylims = (-0.05, 0.45), xlabel = "rezago τ (días)", ylabel = "ACF de |r_t|",
              legend = :outerbottom, legend_columns = 2,
              label = "S&P 500 real",
              title = "3. Volatilidad: ACF de |r_t|")
    plot!(p3, lags, c.abm.acf_abs; color = AZUL, lw = 2.2, marker = :square, ms = 3.5, msw = 0,
          label = "Mercado artificial")
    plot!(p3, lags, c.gbm.acf_abs; color = NARANJA, lw = 1.8, ls = :dot,
          label = "Caminata aleatoria")
    hline!(p3, [-banda, banda]; color = GRIS_CLARO, ls = :dash, lw = 1.2, label = "")

    plot(p1, p2, p3; layout = (1, 3), size = (960, 360), bottom_margin = 5Plots.mm)
end

# ╔═╡ 720459ba-ca99-11f1-aa61-9d62df6b9b0b
let
    c = comparacion_semana1
    f(x) = round(x, digits = 3)
    fila(nombre, d) = "| $nombre | $(f(d.curtosis)) | $(f(100 * d.extremos3))% | $(f(d.acf_r[1])) | $(f(d.acf_abs[5])) | $(f(d.acf_abs[25])) |"
    Markdown.parse(join([
        "| Serie | Curtosis en exceso | Días con ∣z∣ > 3 | ACF de r, rezago 1 | ACF de ∣r∣, rezago 5 | ACF de ∣r∣, rezago 25 |",
        "|---|---:|---:|---:|---:|---:|",
        fila("**S&P 500 real** (`Data/precios.csv`)", c.sp),
        fila("**Caminata aleatoria gaussiana**", c.gbm),
        fila("**Mercado artificial (Franke-Westerhoff)**", c.abm)], "\n"))
end

# ╔═╡ 72045a46-ca99-11f1-96f2-4d211ea71d63
"Diagrama de bifurcación y exponente de Lyapunov del modelo de Brock-Hommes frente a β."
diagrama_bifurcacion_bh = let
    R = 1.10
    g = 1.20
    C = 1.0
    beta_c = log(R / (g - R)) / C

    betas_grid = range(1.5, 4.6, length = 150)
    betas_pts = Float64[]
    x_pts = Float64[]

    for b in betas_grid
        sim = simular_brock_hommes(b, g; R = R, C = C, pasos = 80, calent = 800, x0 = 0.25)
        for val in sim.x
            push!(betas_pts, b)
            push!(x_pts, val)
        end
    end

    betas_lyap = range(1.5, 4.6, length = 60)
    lyap = Float64[]
    for b in betas_lyap
        T_cal = 600
        T_med = 800
        x = [0.25, 0.25 * 1.02, 0.25 * 0.98]
        dx = [1.0, 0.0, 0.0]
        eps = 1e-7
        suma_log = 0.0

        for t in 1:(T_cal + T_med)
            paso_bh(v) = begin
                ret = v[1] - R * v[2]
                pi1 = -R * v[2] * ret - C
                pi2 = (g * v[3] - R * v[2]) * ret
                n2 = 1.0 / (1.0 + exp(clamp(b * (pi1 - pi2), -40.0, 40.0)))
                [(n2 * g * v[1]) / R, v[1], v[2]]
            end
            x_sig = paso_bh(x)
            if t > T_cal
                x_pert = paso_bh(x .+ eps .* dx)
                diff_v = (x_pert .- x_sig) ./ eps
                norma = sqrt(sum(abs2, diff_v))
                if norma > 0.0
                    suma_log += log(norma)
                    dx = diff_v ./ norma
                end
            end
            x = x_sig
        end
        push!(lyap, suma_log / T_med)
    end

    (betas_pts = betas_pts, x_pts = x_pts, betas_lyap = collect(betas_lyap),
     lyap = lyap, beta_c = beta_c)
end;

# ╔═╡ 720459ce-ca99-11f1-8e0e-917599e49a41
let
    b = diagrama_bifurcacion_bh
    p1 = scatter(b.betas_pts, b.x_pts; color = AZUL, alpha = 0.22, ms = 1.5, msw = 0,
                 xlabel = "intensidad de elección β", ylabel = "estados visitados x_t",
                 title = "Diagrama de bifurcación (Brock-Hommes)")
    vline!(p1, [b.beta_c]; color = NARANJA, ls = :dash, lw = 1.8)

    p2 = plot(b.betas_lyap, b.lyap; color = AZUL, lw = 2.0, marker = :circle, ms = 3, msw = 0,
              xlabel = "intensidad de elección β", ylabel = "exponente de Lyapunov λ",
              legend = :outerbottom, legend_columns = 2,
              label = "Exponente λ(β)",
              title = "Estabilidad (λ < 0) y caos (λ > 0)")
    hline!(p2, [0.0]; color = GRIS, ls = :dot, lw = 1.4, label = "Umbral de caos (λ = 0)")
    vline!(p2, [b.beta_c]; color = NARANJA, ls = :dash, lw = 1.8,
           label = "Bifurcación β_c ≈ $(round(b.beta_c, digits=2))")

    plot(p1, p2; layout = (1, 2), size = (880, 370), bottom_margin = 4Plots.mm)
end

# ╔═╡ 72045a50-ca99-11f1-b3cb-0d25921fb76f
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

# ╔═╡ 720459e2-ca99-11f1-be3b-6173fbf2414c
verificar(() -> [fraccion_logit_dos(-1.0, 0.0, 0.0),
                 fraccion_logit_dos(-1.0, 0.0, 2.0),
                 fraccion_logit_dos(1.5, 0.5, 3.0)],
          [0.5, 1.0 / (1.0 + exp(-2.0)), 1.0 / (1.0 + exp(3.0))];
          pista = md"Usa `1.0 / (1.0 + exp(beta * (pi_1 - pi_2)))`.")

# ╔═╡ 720459ec-ca99-11f1-8c7b-29cd8b131922
verificar(() -> [paso_precio_bh(2.0, 0.85, 1.2, 1.01),
                 paso_precio_bh(-1.5, 0.5, 1.1, 1.02),
                 paso_precio_bh(0.0, 0.9, 1.3, 1.01)],
          [0.85 * 1.2 * 2.0 / 1.01, 0.5 * 1.1 * (-1.5) / 1.02, 0.0];
          pista = md"Despeja `x_t` de `R * x_t = n2 * g * x_ant`: `(n2 * g * x_ant) / R`.")

# ╔═╡ 72045a02-ca99-11f1-9973-a9d2afb79aa6
verificar(() -> [beta_critico_bh(1.20, 1.01, 1.0),
                 beta_critico_bh(1.35, 1.01, 0.5)],
          [log(1.01 / (1.20 - 1.01)), log(1.01 / (1.35 - 1.01)) / 0.5];
          pista = md"Usa `log(R / (g - R)) / C`.")

# ╔═╡ 72045a14-ca99-11f1-9079-f10fa772c78d
verificar(() -> [multiplicador_origen_bh(0.0, 1.20, 1.01, 1.0),
                 multiplicador_origen_bh(3.0, 1.20, 1.01, 1.0)],
          [1.20 * 0.5 / 1.01, (1.20 / 1.01) / (1.0 + exp(-3.0))];
          pista = md"Calcula `n2 = 1.0 / (1.0 + exp(-beta * C))` y devuelve `(g * n2) / R`.")

# ╔═╡ Cell order:
# ╟─72045938-ca99-11f1-997a-8789a555af3c
# ╟─7204595e-ca99-11f1-8c3f-e7271e4f8224
# ╟─7204596a-ca99-11f1-a174-7107e6dc34fd
# ╟─72045974-ca99-11f1-a107-a11749a7d068
# ╟─72045974-ca99-11f1-b0b3-7b599481ad7c
# ╟─7204597e-ca99-11f1-8553-21fb958df804
# ╟─72045988-ca99-11f1-b3d7-7134e184652c
# ╟─72045992-ca99-11f1-a3c8-9d1d92558d24
# ╟─7204599c-ca99-11f1-acc9-b3c53724b071
# ╟─720459a6-ca99-11f1-beb4-555e7bdce368
# ╟─720459a6-ca99-11f1-b086-ef4ff44161c1
# ╟─720459b0-ca99-11f1-b03b-9384321db7ce
# ╟─720459ba-ca99-11f1-aa61-9d62df6b9b0b
# ╟─720459ba-ca99-11f1-8138-ff78ddc3e4d7
# ╟─720459c4-ca99-11f1-9ac5-f3d362c08707
# ╟─720459ce-ca99-11f1-8e0e-917599e49a41
# ╟─720459ce-ca99-11f1-9916-f33b28a01215
# ╟─720459d8-ca99-11f1-aa7e-23c97971346b
# ╠═720459d8-ca99-11f1-aa52-b36c16d3aba7
# ╟─720459e2-ca99-11f1-be3b-6173fbf2414c
# ╟─720459e2-ca99-11f1-b4d4-7dcd4aeface8
# ╠═720459ec-ca99-11f1-bf99-27d5104f9942
# ╟─720459ec-ca99-11f1-8c7b-29cd8b131922
# ╟─720459f6-ca99-11f1-b8d3-b56deb6e8ab9
# ╟─720459f6-ca99-11f1-92f7-07bd0433376a
# ╠═72045a02-ca99-11f1-be6b-bf846fdc646b
# ╟─72045a02-ca99-11f1-9973-a9d2afb79aa6
# ╟─72045a0a-ca99-11f1-84bc-a56a46d3317c
# ╠═72045a0a-ca99-11f1-a31a-a525bc67a70f
# ╟─72045a14-ca99-11f1-9079-f10fa772c78d
# ╟─72045a14-ca99-11f1-a16b-334b6c296970
# ╟─72045a1e-ca99-11f1-9475-0bf39586ed6a
# ╠═72045a1e-ca99-11f1-90c6-d143fcae313a
# ╠═72045a28-ca99-11f1-adfc-bdcade772936
# ╠═72045a28-ca99-11f1-bb79-b5f72d9c6406
# ╠═72045a34-ca99-11f1-a139-45cd87f281da
# ╠═72045a3c-ca99-11f1-9577-13a48baf26b4
# ╠═72045a3c-ca99-11f1-8aa3-7dd8c0801684
# ╠═72045a46-ca99-11f1-9d5f-c51ef7d94e39
# ╠═72045a46-ca99-11f1-96f2-4d211ea71d63
# ╠═72045a50-ca99-11f1-b3cb-0d25921fb76f

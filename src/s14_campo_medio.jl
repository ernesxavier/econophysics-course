#%% hide
md"""
# Sesión 14: Campo medio, histéresis y crisis

**Objetivos de la sesión**
1. Resolver el **campo medio** gráficamente y compararlo con la simulación.
2. Entender cuándo aparecen **equilibrios múltiples** y cómo el
   **multiplicador social** amplifica los incentivos (Brock y Durlauf).
3. Observar la **histéresis**: cambios graduales de incentivos que producen
   saltos de comportamiento que no se revierten al deshacer la causa.
4. Explorar cómo la **heterogeneidad** entre agentes controla el tamaño de
   las **avalanchas** y las crisis (ruta avanzada).

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s14_campo_medio.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B4 (distribución
normal), C.I1 (Taylor), C.I5 (puntos fijos), F.I4 (campo medio). Ruta
avanzada: C.A4 (Stirling), F.A1 (Landau), F.A2 (exponentes críticos), F.A3
(campo aleatorio y avalanchas).

**Conceptos transversales.** Se completan los **puntos fijos y su
estabilidad** (solución gráfica del campo medio), que vuelven en las sesiones
17 y 20.
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: cambios abruptos que no se revierten al deshacer la causa

Un gobierno quiere impulsar una tecnología nueva y ofrece un subsidio ``h``
que sube poco a poco, desde ``-0.6`` hasta ``+0.6``, y después vuelve a bajar.
¿Qué pasa cuando retira el subsidio?

Depende de cuánto se imiten los ciudadanos entre sí. Abajo comparamos dos
poblaciones de 2,000 agentes: en la de la izquierda la imitación es débil
(``\beta J = 0.7``); en la de la derecha es fuerte (``\beta J = 1.5``).
"""

#%% hide
let
    ps = map(((0.7, "Imitación débil (βJ = 0.7): reversible"),
              (1.5, "Imitación fuerte (βJ = 1.5): irreversible"))) do (J, titulo)
        c = histeresis(2_000, J; hs = -0.6:0.03:0.6)
        p = plot(c.h_sube, c.m_sube; color = AZUL, lw = 2, marker = :circle, ms = 2.5, msw = 0,
                 label = "h sube", xlabel = "incentivo común h", ylabel = "adopción neta m",
                 ylims = (-1.05, 1.05), legend = :topleft, title = titulo)
        plot!(p, c.h_baja, c.m_baja; color = NARANJA, lw = 2, marker = :square, ms = 2.5, msw = 0,
              label = "h baja")
        vline!(p, [0.0]; color = GRIS_CLARO, ls = :dot, label = "")
    end
    plot(ps...; layout = (1, 2), size = (880, 360), bottom_margin = 4Plots.mm)
end

#%% hide
md"""
Con imitación débil, la adopción sube suavemente y regresa por el mismo
camino cuando se retira el subsidio. Con imitación fuerte, al principio casi
nadie adopta aunque ``h`` ya sea positivo; de pronto la población **salta** a
la adopción masiva, y cuando el subsidio vuelve a ``h = 0`` la adopción **se
queda arriba**. Hoy veremos de dónde salen ese salto y esa memoria.
"""

#%% hide
md"""
## 1. Campo medio: la solución gráfica

En la cuadrícula de la sesión 13 cada agente sólo veía a 4 vecinos. Brock y
Durlauf (2001) estudiaron el caso en que a cada agente le importa el
**promedio** de toda la población: su utilidad social es ``J\, s_i\, m``. En
ese caso, el campo que siente el agente es ``h + Jm``, el mismo para todos, y
la regla logit da (repaso F.I4)

```math
\mathbb{P}(s_i = +1) = \sigma\big(2\beta(h + Jm)\big)
\quad\Longrightarrow\quad
\mathbb{E}[s_i] = \tanh\big(\beta(h + Jm)\big) .
```

Si la población es grande, el promedio de las elecciones es ``m`` mismo:

```math
m = \tanh\big(\beta(h + Jm)\big) .
```

Es una ecuación de **punto fijo** (repaso C.I5): ``m`` aparece en los dos
lados. Se resuelve buscando dónde se cruzan la recta ``y = m`` y la curva
``y = \tanh(\beta(h + Jm))``.

Fijamos ``\beta = 1`` (así ``J`` y ``h`` se miden en unidades del ruido).

Intensidad de la interacción ``J``: $(@bind J_cm Slider([0.5, 0.8, 1.0, 1.2, 1.5, 2.0], default = 1.5, show_value = true))

Incentivo común ``h``: $(@bind h_cm Slider([-0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4], default = 0.0, show_value = true))
"""

#%% hide
let
    x = range(-1.2, 1.2, length = 400)
    sol = soluciones_cm(J_cm, h_cm)
    p = plot(x, x; color = GRIS, lw = 1.5, label = "y = m", xlims = (-1.2, 1.2), ylims = (-1.2, 1.2),
             aspect_ratio = :equal, xlabel = "adopción neta m", ylabel = "y", legend = :topleft,
             title = "J = $J_cm, h = $h_cm: $(length(sol) == 1 ? "1 solución" : "$(length(sol)) soluciones")",
             size = (520, 500))
    plot!(p, x, tanh.(h_cm .+ J_cm .* x); color = AZUL, lw = 2.5, label = "y = tanh(h + J m)")
    for m in sol
        estable = J_cm * (1 - tanh(h_cm + J_cm * m)^2) < 1
        scatter!(p, [m], [m]; ms = 7, color = estable ? NARANJA : BLANCO, msc = NARANJA, msw = 2,
                 label = "")
    end
    scatter!(p, [NaN], [NaN]; ms = 7, color = NARANJA, msw = 0, label = "Equilibrio estable")
    scatter!(p, [NaN], [NaN]; ms = 7, color = BLANCO, msc = NARANJA, msw = 2, label = "Equilibrio inestable")
end

#%% hide
md"""
- Con ``J < 1`` y ``h = 0`` hay una sola solución, ``m = 0``: la población se
  divide a la mitad. La conformidad es demasiado débil frente al ruido.
- Con ``J > 1`` y ``h = 0``, la pendiente de ``\tanh`` en el origen es mayor que
  1 y aparecen **dos soluciones más**, ``\pm m^\ast``: dos consensos posibles
  (cerca de ``J = 1``, ``m^\ast \approx \sqrt{3(J - 1)}`` por desarrollo de
  Taylor; repaso C.I1, F.A1 y F.A2). ``m = 0`` se vuelve **inestable**:
  cualquier desviación se amplifica.
- La transición está en ``\beta J = 1``. En lenguaje de Brock y Durlauf: hay
  **equilibrios múltiples** cuando la interacción social es más fuerte que el
  ruido.
- Con ``h \ne 0`` y ``J`` grande puede haber tres soluciones, y una de ellas va
  **contra** el incentivo: la población puede estar atrapada en el consenso
  "equivocado" (una **trampa de coordinación**).

**¿Funciona?** Simulamos 2,000 agentes que interactúan con el promedio de
todos (con ``h = 0``) y los comparamos con la solución gráfica.
"""

#%%
barrido_cm = let
    Js = 0.4:0.1:2.0
    (J = Js, m = [mean(abs.(campo_medio(2_000, J, 0.0; barridos = 400, rng = Xoshiro(round(Int, 10J))).m[201:end]))
                  for J in Js])
end;

#%% hide
let
    x = range(0.4, 2.0, length = 300)
    p = plot(x, [maximum(abs.(soluciones_cm(J, 0.0))) for J in x]; color = NARANJA, lw = 2.5,
             label = "Campo medio: mayor solución de m = tanh(J m)", xlabel = "interacción βJ",
             ylabel = "|m|", ylims = (0, 1.05), legend = :topleft,
             title = "2,000 agentes que miran el promedio de todos")
    scatter!(p, barrido_cm.J, barrido_cm.m; color = AZUL, ms = 5, msw = 0, label = "Simulación (baño térmico)")
    vline!(p, [1]; color = GRIS, ls = :dash, label = "βJ = 1")
end

#%% hide
md"""
Aquí la simulación sigue a la teoría en todo el rango: cuando todos
interactúan con todos, el campo medio es **exacto** para ``N`` grande. En la
red cuadrada de la sesión 13 no lo es: con ``z = 4`` vecinos el campo medio
predice ``T_c = Jz = 4``, mientras que la solución exacta es 2.269. El campo
medio ignora las fluctuaciones locales y **sobreestima** el orden (repaso
F.I4).
"""

#%% hide
md"""
## 2. Laboratorio: histéresis

Ahora fijamos la interacción y movemos **lentamente** el incentivo ``h``:
primero de ``-0.6`` a ``0.6`` (por ejemplo, un subsidio que crece) y luego de
regreso. En cada valor de ``h`` dejamos que la población se ajuste unos
barridos y medimos ``m``.

Interacción ``\beta J``: $(@bind J_hist Slider([0.8, 1.2, 1.5, 2.0], default = 1.5, show_value = true))
"""

#%%
ciclo = histeresis(2_000, J_hist; hs = -0.6:0.02:0.6);

#%% hide
let
    p = plot(ciclo.h_sube, ciclo.m_sube; color = AZUL, lw = 2, marker = :circle, ms = 2.5, msw = 0,
             label = "Simulación: h sube", xlabel = "incentivo común h", ylabel = "adopción neta m",
             ylims = (-1.05, 1.05), legend = :topleft, title = "Ciclo de histéresis con βJ = $J_hist")
    plot!(p, ciclo.h_baja, ciclo.m_baja; color = NARANJA, lw = 2, marker = :square, ms = 2.5, msw = 0,
          label = "Simulación: h baja")
    for (ramas, c) in ((rama_cm(J_hist, -0.6:0.005:0.6, -1.0), AZUL), (rama_cm(J_hist, 0.6:-0.005:-0.6, 1.0), NARANJA))
        plot!(p, ramas.h, ramas.m; color = c, lw = 1, ls = :dash, label = "")
    end
    plot!(p, [NaN], [NaN]; color = GRIS, ls = :dash, label = "Campo medio")
    if J_hist > 1
        hc = salto_cm(J_hist)
        vline!(p, [-hc, hc]; color = GRIS_CLARO, ls = :dot, label = "Saltos teóricos: h = ±$(round(hc, digits = 2))")
    end
    p
end

#%% hide
md"""
Con ``\beta J > 1``:

- Al subir ``h``, la población se queda en la rama de **no adopción** aunque
  el incentivo ya sea positivo. Sigue ahí hasta que esa rama desaparece, y
  entonces **salta** de golpe a casi la adopción total.
- Al bajar ``h``, se queda en la rama de **adopción** hasta un valor negativo
  de ``h``, y entonces salta hacia abajo.
- Entre los dos saltos, el estado de la población depende de **su historia**,
  no sólo del incentivo actual.

**Lecturas económicas.** Un cambio gradual en los incentivos (precios,
regulación, información) puede producir un cambio de comportamiento
**repentino**: adopción de tecnologías, cambios de normas sociales, corridas
bancarias. Y revertir el incentivo no basta para revertir el
comportamiento: hay que ir mucho más allá. Con ``\beta J < 1`` (mueve el
deslizador a 0.8) la respuesta es suave y reversible, aunque la interacción
social amplifica la pendiente en ``h = 0`` por un factor ``1/(1 - \beta J)``:
el **multiplicador social** de las notas formales.

Los saltos de la simulación ocurren un poco **después** que los teóricos.
Cerca del punto donde desaparece una rama, la población se ajusta muy
despacio (en física se llama *frenado crítico*), y 20 barridos por valor de
``h`` no bastan para completar el salto. Con un cambio más lento del
incentivo, los saltos se acercan a los teóricos.
"""

#%% hide
md"""
## 3. [Avanzado] Campo aleatorio, avalanchas y crisis

En la sección 2 todos los agentes eran iguales salvo por el ruido. En la
realidad, cada agente tiene su propio incentivo privado ``h_i``: unos están
casi convencidos de adoptar, otros muy lejos. Este es el **modelo de Ising
con campo aleatorio** (repaso F.A3). A temperatura cero (sin ruido), el
agente ``i`` adopta si

```math
h + h_i + J m > 0 ,
```

con ``h_i`` normal con media 0 y desviación estándar ``R`` (repaso P.B4): la
**heterogeneidad**. Cuando ``h`` sube un poco, adopta el agente más cercano
a su umbral; eso sube ``m``, lo que puede convencer al siguiente, y así: una
**avalancha**.

Heterogeneidad ``R``: $(@bind R_rfim Slider([0.4, 0.6, 0.8, 1.0, 1.5], default = 0.6, show_value = true))
"""

#%%
avalanchas = campo_aleatorio(10_000, 1.0, R_rfim);

#%% hide
let
    d = avalanchas
    p1 = plot(d.h, d.m; color = AZUL, lw = 2, seriestype = :steppost, ylims = (-1.05, 1.05),
              xlabel = "incentivo común h", ylabel = "adopción neta m", legend = false, xlims = (-2, 2),
              title = "J = 1, R = $R_rfim")
    s = sort(unique(d.tamaños))
    P = [mean(d.tamaños .>= x) for x in s]
    p2 = scatter(s, P; xscale = :log10, yscale = :log10, color = AZUL, ms = 3, msw = 0, legend = false,
                 xlabel = "tamaño de la avalancha (agentes)", ylabel = "P(tamaño ≥ s)",
                 title = "Avalancha más grande: $(miles(maximum(d.tamaños))) de 10,000")
    plot(p1, p2; layout = (1, 2), size = (920, 360), bottom_margin = 5Plots.mm)
end

#%% hide
md"""
- Con **poca heterogeneidad** (``R = 0.4`` o ``0.6``), hay una avalancha
  gigante: una fracción grande de la población cambia de golpe. Es una
  **crisis**.
- Con **mucha heterogeneidad** (``R = 1.5``), todas las avalanchas son
  pequeñas y la adopción crece suavemente.
- En medio, cerca de ``R_c = J\sqrt{2/\pi} \approx 0.8``, las avalanchas tienen
  **todos los tamaños**: su distribución es una ley de potencias.

Bouchaud (2013) usa este modelo para explicar cambios abruptos de opinión,
pánicos y caídas del consumo: la diversidad de opiniones protege contra los
saltos colectivos, y cuando la diversidad baja, el sistema se vuelve frágil.
La derivación de ``R_c`` está en las notas formales.
"""

#%% hide
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Campo medio | ``m = \tanh(\beta(h + Jm))``; equilibrios múltiples si ``\beta J > 1`` |
| Multiplicador social | La conformidad amplifica el efecto de un incentivo por ``1/(1 - \beta J(1 - m^2))`` |
| Histéresis | Saltos al mover ``h``; el estado depende de la historia |
| Campo aleatorio [Avanzado] | Poca heterogeneidad (``R < J\sqrt{2/\pi}``) → crisis; mucha → cambio suave |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `punto_fijo_cm(J, h; m0 = 1.0, pasos = 200)`, que resuelva la
ecuación de campo medio ``m = \tanh(h + J m)`` iterando `m = tanh(h + J * m)`
durante `pasos` iteraciones a partir de `m0`.
"""

#%%
function punto_fijo_cm(J, h; m0 = 1.0, pasos = 200)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [punto_fijo_cm(0.5, 0.0), punto_fijo_cm(1.5, 0.0; m0 = 1.0), punto_fijo_cm(1.5, 0.0; m0 = -1.0)],
          [0.0, rama_cm(1.5, [0.0], 1.0).m[1], -rama_cm(1.5, [0.0], 1.0).m[1]];
          rtol = 1e-5,
          pista = md"Empieza con `m = m0`, repite `m = tanh(h + J * m)` en un ciclo `for _ in 1:pasos` y devuelve `m`.")

#%% hide
md"""
**2.** Escribe `multiplicador_social(J, m)`, el factor ``1 / (1 - J(1 - m^2))``
por el que la interacción social ``J`` (con ``\beta = 1``) amplifica el efecto
de un cambio pequeño en el incentivo ``h`` alrededor de un equilibrio `m`.
"""

#%%
function multiplicador_social(J, m)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [multiplicador_social(0.0, 0.0), multiplicador_social(0.8, 0.0), multiplicador_social(1.5, 0.9)],
          [1.0, 5.0, 1 / (1 - 1.5 * (1 - 0.9^2))];
          pista = md"La fórmula es `1 / (1 - J * (1 - m^2))`.")

#%% hide
md"""
### Ruta básica: interpretación

3. Un gobierno subsidia una tecnología y, cuando la adopción ya es alta,
   retira el subsidio. Según la sección 2, ¿qué puede pasar? ¿De qué
   depende?
4. ¿Por qué el campo medio funciona perfecto en la sección 1 y no en la red
   cuadrada de la sesión 13?
5. ¿Qué significa que exista un equilibrio estable con ``m < 0`` cuando
   ``h > 0``? Da un ejemplo económico o social de esta trampa de
   coordinación.

### Ruta avanzada

6. Usa el desarrollo de Taylor de ``\tanh`` para mostrar que, con ``h = 0`` y
   ``\beta J`` un poco mayor que 1, ``m^\ast \approx \sqrt{3(\beta J - 1)}``.
   Compara con la gráfica de la sección 1.
7. En la sección 3, mide el tamaño de la avalancha más grande para varios
   valores de ``R`` y para ``N = 1{,}000`` y ``10{,}000``. ¿Cómo ubicarías
   ``R_c`` con estos datos?

## Para la próxima sesión

- Lectura: Requião da Cunha (2021), capítulo 6. Para la ruta avanzada,
  Brock y Durlauf (2001) y Bouchaud, Marsili y Nadal (2023).
- **Entrega**: taller para casa del Módulo II (desigualdad con la ENIGH y
  contagio en redes), en `Code/taller2_desigualdad_contagio.jl`.
- En la sesión 15 empieza el Módulo III (Agentes, estrategias y mercados
  artificiales): modelos de opinión (votante, Kirman) y herding en mercados.
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
    using PlutoUI, Statistics, Random, Distributions, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

#%%
begin
    const AZUL       = "#2a78d6"   # serie principal, opción +1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, opción -1
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
    const RAMPA      = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

#%%
"Función logística."
σ(x) = 1 / (1 + exp(-x))

#%%
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

#%%
"""
Modelo de campo medio con N agentes y β = 1: cada revisión elige un agente
al azar y le asigna +1 con probabilidad logística σ(2(h + J m)), donde m es
la adopción neta de los demás (baño térmico). Devuelve la última
configuración y m después de cada barrido.
"""
function campo_medio(N, J, h; barridos = 400, rng = Xoshiro(1), s = rand(rng, [-1, 1], N))
    s = copy(s)
    total = sum(s)
    m = Float64[]
    for _ in 1:barridos
        for _ in 1:N
            i = rand(rng, 1:N)
            m_otros = (total - s[i]) / (N - 1)
            nuevo = rand(rng) < σ(2 * (h + J * m_otros)) ? 1 : -1
            total += nuevo - s[i]
            s[i] = nuevo
        end
        push!(m, total / N)
    end
    (s = s, m = m)
end

#%%
"Soluciones de m = tanh(h + J m) en [-1, 1], encontradas por cambio de signo y bisección."
function soluciones_cm(J, h)
    f(m) = tanh(h + J * m) - m
    xs = range(-1.0, 1.0, length = 2001)
    raices = Float64[]
    for k in 1:length(xs)-1
        a, b = xs[k], xs[k+1]
        f(a) == 0 && push!(raices, a)
        if f(a) * f(b) < 0
            for _ in 1:60
                c = (a + b) / 2
                f(a) * f(c) <= 0 ? (b = c) : (a = c)
            end
            push!(raices, (a + b) / 2)
        end
    end
    raices
end

#%%
"Rama de campo medio que se sigue al mover h lentamente: se itera m = tanh(h + J m) desde la solución anterior."
function rama_cm(J, hs, m0)
    m = m0
    ms = map(hs) do h
        for _ in 1:5_000
            m = tanh(h + J * m)
        end
        m
    end
    (h = collect(hs), m = ms)
end

#%%
"Valor de h donde desaparece una rama en el campo medio (J > 1)."
function salto_cm(J)
    m = sqrt(1 - 1 / J)       # donde la pendiente de tanh(h + J m) es 1
    J * m - atanh(m)
end

#%%
"Ciclo de histéresis simulado: h sube y luego baja, con 20 barridos en cada valor."
function histeresis(N, J; hs = -0.6:0.02:0.6, barridos = 20, rng = Xoshiro(4))
    s = fill(-1, N)
    m_sube, m_baja = Float64[], Float64[]
    for h in hs
        r = campo_medio(N, J, h; barridos, rng, s)
        s = r.s
        push!(m_sube, mean(r.m[end-4:end]))
    end
    for h in reverse(hs)
        r = campo_medio(N, J, h; barridos, rng, s)
        s = r.s
        push!(m_baja, mean(r.m[end-4:end]))
    end
    (h_sube = collect(hs), m_sube = m_sube, h_baja = collect(reverse(hs)), m_baja = m_baja)
end

#%%
"""
Ising con campo aleatorio a temperatura cero en campo medio. Los incentivos
privados h_i son normales con desviación R. Se sube h desde -∞: los agentes
adoptan en orden de h_i decreciente, y cada adopción sube el campo de todos
en 2J/N. Devuelve m después de cada avalancha, el h en que ocurrió y el
tamaño de cada avalancha.
"""
function campo_aleatorio(N, J, R; rng = Xoshiro(5))
    hi = sort(R .* randn(rng, N); rev = true)
    arriba = 0                       # agentes que ya adoptaron
    hs, ms, tamaños = Float64[], Float64[], Int[]
    while arriba < N
        m = 2arriba / N - 1
        h = -hi[arriba+1] - J * m    # h justo para que adopte el siguiente
        inicio = arriba
        arriba += 1
        while arriba < N && h + hi[arriba+1] + J * (2arriba / N - 1) > 0
            arriba += 1
        end
        push!(hs, h)
        push!(ms, 2arriba / N - 1)
        push!(tamaños, arriba - inicio)
    end
    (h = [-5.0; hs], m = [-1.0; ms], tamaños = tamaños)
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

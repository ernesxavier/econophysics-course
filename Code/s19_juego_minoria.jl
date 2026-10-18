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

# ╔═╡ 256527d8-ca8a-11f1-8f74-d9bef392ffdd
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, Statistics, Random, Distributions, Plots, Markdown, LinearAlgebra
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 256526f2-ca8a-11f1-884c-2ff3f7ee1b1d
md"""
# Sesión 19: Juego de la minoría

**Objetivos de la sesión**
1. Entender el **problema del bar El Farol** (Arthur, 1994) y por qué en los
   juegos de congestión o especulación no puede existir una regla predictiva
   única compartida por todos.
2. Implementar las reglas del **juego de la minoría** (Challet y Zhang, 1997):
   historias de longitud ``M``, tablas de estrategias y aprendizaje por puntajes
   virtuales.
3. Interpretar el juego de la minoría como un **modelo de mercado financiero**
   donde ``A(t)`` es el exceso de demanda y ``\sigma^2/N`` mide la volatilidad y
   la ineficiencia en el uso del recurso.
4. Descubrir el colapso de escala con el parámetro de control
   ``\alpha = 2^M / N`` y distinguir la **fase hacinada (ineficiente)** de la
   **fase diluida** alrededor del mínimo de volatilidad en ``\alpha_c \approx 0.34``.
5. **[Avanzado]** Medir la **predictibilidad** ``H/N`` y la fracción de
   **agentes congelados** ``\phi(\alpha)`` para explicar la transición de fase
   en ``\alpha_c = 1 - \phi(\alpha_c)``.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s19_juego_minoria.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B2 (media y varianza de una
suma), E.B3 (estacionariedad y promedios temporales), F.B4 (transiciones de
fase), S.I4 (datos dependientes y promedios sobre realizaciones). Ruta
avanzada: F.A2 (exponentes críticos), F.A4 (vidrios de espín y frustración).
"""

# ╔═╡ 25652724-ca8a-11f1-a7b1-89bb41b9e06d
TableOfContents(title = "Contenido")

# ╔═╡ 25652724-ca8a-11f1-8229-df8198032751
md"""
## 0. El fenómeno: el bar El Farol y la trampa de pensar igual

En el modelo de Ising (sesión 13) y en los juegos de coordinación (sesión 17),
a cada agente le convenía hacer lo mismo que la mayoría. Pero en muchos
problemas económicos ocurre exactamente lo contrario: **gana quien está en la
minoría**.

W. Brian Arthur (1994) lo ilustró con el **problema del bar El Farol** en Santa
Fe:
- Cada jueves por la noche, 100 personas deciden si van al bar o se quedan en
  casa.
- El bar tiene cupo cómodo para **60 personas**. Si asisten 60 o menos, quienes
  fueron pasan una gran noche; si asisten más de 60, el bar se abarrota y habría
  sido mejor quedarse en casa.

¿Qué pasa si todos usan la **misma regla razonable** ("si la semana pasada
estuvo cómodo, hoy voy; si estuvo lleno, hoy no voy")? Todos van al mismo
tiempo, llenan el bar, a la semana siguiente casi nadie va y el bar queda casi
vacío. En cambio, cuando cada persona usa su propio par de reglas simples sobre
las últimas semanas y aprende cuál funciona mejor, la asistencia se estabiliza
sola alrededor del cupo de 60 (Challet, Marsili y Ottino, 2004):
"""

# ╔═╡ 25652738-ca8a-11f1-bd98-a516e4f10e8a
md"""
Sin que nadie coordine los turnos, un conjunto de agentes con **racionalidad
acotada e inductiva** resuelve un problema de coordinación que una regla
deductiva única e idéntica para todos vuelve imposible.
"""

# ╔═╡ 25652740-ca8a-11f1-9cfc-cfaff8882fe3
md"""
## 1. Repaso: el referente aleatorio (``\sigma^2 / N = 1``)

Para estudiar matemáticamente el problema de El Farol, Challet y Zhang (1997)
propusieron una versión simétrica (con umbral en el 50%): el **juego de la
minoría**.

En cada turno ``t``, ``N`` agentes (con ``N`` impar) eligen entre dos acciones:
``a_i(t) \in \{-1, +1\}``. La **acción agregada** es

```math
A(t) = \sum_{i=1}^N a_i(t) ,
```

y ganan quienes quedaron del lado minoritario, es decir, quienes eligieron el
signo opuesto a ``A(t)``: ``a_{\text{min}}(t) = -\mathrm{sgn}(A(t))``.

Si ``k`` agentes eligen ``+1`` y ``N - k`` eligen ``-1``, entonces
``A(t) = 2k - N`` y el número total de ganadores en ese turno es

```math
N_{\text{gan}}(t) = \frac{N - |A(t)|}{2} .
```

Cuanto más cerca de ``0`` se mantenga ``A(t)``, más agentes ganan (hasta un
máximo de ``(N-1)/2`` cuando ``|A(t)| = 1``). Por eso medimos las fluctuaciones
con la **volatilidad reescalada** ``\sigma^2 / N = \langle A^2\rangle / N``.

**¿Cuánto vale ``\sigma^2/N`` si todos eligen al azar?** Si cada agente lanza
una moneda justa independiente (``a_i = \pm 1`` con probabilidad ``1/2``), cada
término tiene media ``0`` y varianza ``1``. La varianza de la suma de ``N``
variables independientes es ``\sigma_{\text{azar}}^2 = N`` (repaso P.B2), así
que el **referente aleatorio** es exactamente:

```math
\frac{\sigma_{\text{azar}}^2}{N} = 1 .
```
"""

# ╔═╡ 2565274c-ca8a-11f1-ab95-41d48413308c
# Número de ganadores en un turno con N agentes y acción agregada A
ganadores_minoria(A, N) = (N - abs(A)) ÷ 2

# ╔═╡ 25652756-ca8a-11f1-8981-b7f20b911d8e
# Con N = 101: si A = 1 ganan 50 agentes (el máximo); si todos hacen lo mismo (A = 101) no gana nadie
ganadores_minoria(1, 101), ganadores_minoria(21, 101), ganadores_minoria(101, 101)

# ╔═╡ 25652756-ca8a-11f1-a181-a580ba7c492b
md"""
## 2. Reglas del juego de la minoría y lectura como mercado

¿Cómo deciden los agentes si no lanzan monedas?

1. **Memoria e historias (``P = 2^M``):** Todos observan la lista de los
   últimos ``M`` bandos ganadores: una cadena binaria de longitud ``M``, como
   `(-1, +1, +1)`. Como hay ``M`` pasos pasados y 2 opciones en cada uno,
   existen ``P = 2^M`` **historias posibles**, numeradas ``\mu \in \{1, \dots, P\}``.
2. **Tablas de estrategias:** Una estrategia ``s`` es una tabla completa que
   dice qué hacer (``-1`` o ``+1``) ante cada una de las ``P = 2^M`` historias.
   Al inicio del juego, cada agente recibe ``S = 2`` estrategias fijas sorteadas
   al azar entre todas las ``2^P`` tablas posibles.
3. **Puntajes virtuales y elección:** Cada agente lleva un puntaje acumulado
   ``U_{i,s}(t)`` para cada una de sus ``S = 2`` estrategias. En el turno ``t``:
   - Mira la historia actual ``\mu(t)`` y usa la estrategia que tenga **mayor
     puntaje** hasta ese momento.
   - Una vez que se conoce la acción agregada ``A(t)``, actualiza los puntajes
     de **sus dos estrategias** según habrían acertado o fallado:
     ```math
     U_{i,s}(t + 1) = U_{i,s}(t) - s\bigl(\mu(t)\bigr)\,A(t) .
     ```
     Si la estrategia recomendaba el lado minoritario (signo opuesto a
     ``A(t)``), ``-s(\mu)A(t) > 0`` y su puntaje sube; si recomendaba ir con la
     multitud, su puntaje baja.

**El juego de la minoría como mercado financiero.** Si interpretamos
``a_i(t) = +1`` como una orden de **compra** y ``a_i(t) = -1`` como una orden
de **venta**:
- ``A(t) = N_{\text{compra}}(t) - N_{\text{venta}}(t)`` es el **exceso de
  demanda**, que mueve el precio ``p(t)`` (en minúscula, para no confundirlo
  con el número de historias ``P``): ``r(t) = \ln p(t+1) - \ln p(t) \propto A(t)``.
- Cuando casi todos compran (``A(t) > 0``), empujan el precio hacia arriba y
  compran caro, mientras que la minoría que vende obtiene un precio alto.
- ``\sigma^2 / N = \langle A^2\rangle / N`` es la **volatilidad del mercado**.

### Laboratorio: del efecto manada a la coordinación espontánea

Fijamos ``N = 101`` agentes con ``S = 2`` estrategias cada uno. Mueve la
longitud de la memoria ``M`` (que cambia el número de historias ``P = 2^M`` y
el cociente ``\alpha = 2^M / N``):

Memoria ``M`` (con ``N = 101``): $(@bind M_lab Slider([2, 3, 4, 5, 6, 7, 8, 10], default = 5, show_value = true))
"""

# ╔═╡ 2565277e-ca8a-11f1-830a-edf9ff6800aa
md"""
## 3. Laboratorio: la curva universal ``\sigma^2 / N`` contra ``\alpha = 2^M / N``

¿Qué pasa si cambiamos el número de agentes ``N``? Savit, Manuca y Riolo (1999)
encontraron un hecho notable: si graficamos la volatilidad ``\sigma^2/N`` contra
la memoria ``M`` para ``N = 51``, ``101`` y ``201``, las curvas quedan
desfasadas; pero si las graficamos contra el **parámetro de control**

```math
\alpha = \frac{P}{N} = \frac{2^M}{N} ,
```

las tres curvas **colapsan sobre una sola curva universal** con un mínimo cerca
del punto crítico ``\alpha_c \approx 0.3374``. Cada punto promedia 20 sorteos
de estrategias, y en cada uno dejamos correr el juego ``50P`` turnos antes de
medir: con muchas historias, el aprendizaje tarda en llegar al estado
estacionario.
"""

# ╔═╡ 25652792-ca8a-11f1-97c9-053ccff4d8d6
md"""
La gráfica derecha resume toda la física del juego de la minoría:

1. **Izquierda del mínimo (``\alpha < \alpha_c \approx 0.34``, fase hacinada):**
   En escala log-log, ``\sigma^2 / N`` cae como una recta de pendiente ``-1``,
   es decir, ``\sigma^2 / N \propto \alpha^{-1} = N / 2^M``. Cuando muchos
   inversionistas usan el mismo conjunto pequeño de señales técnicas, sus
   intentos por ganarle al mercado crean **volatilidad endógena masiva**.
2. **En el mínimo (``\alpha \approx \alpha_c \approx 0.34``, fase eficiente):**
   La volatilidad es mínima (``\sigma^2 / N \approx 0.2``): la diversidad de
   estrategias permite que los agentes se turnen implícitamente y aprovechen el
   recurso mucho mejor que tirando monedas.
3. **Derecha del mínimo (``\alpha > \alpha_c``, fase diluida):** Con demasiadas
   historias frente al número de agentes, el aprendizaje pierde coordinación y
   ``\sigma^2 / N`` sube desde abajo hacia ``1`` (el límite aleatorio).
"""

# ╔═╡ 2565279c-ca8a-11f1-b39c-e50ac452a555
md"""
## 4. [Avanzado] Teoría del punto crítico ``\alpha_c``: predictibilidad y congelamiento

¿Por qué el mínimo de volatilidad en ``\alpha_c \approx 0.3374`` no es un simple
cruce suave, sino una verdadera **transición de fase** (repaso F.B4 y F.A4)?

Challet y Marsili (1999) demostraron que el parámetro de orden que cambia en
``\alpha_c`` es la **predictibilidad** ``H / N``: el promedio sobre las ``P``
historias del cuadrado del exceso de demanda medio condicional
``\langle A \mid \mu\rangle``:

```math
\frac{H}{N} = \frac{1}{P N}\sum_{\mu=1}^P \langle A \mid \mu\rangle^2 .
```

- **En la fase simétrica (``\alpha < \alpha_c``):** ``\langle A \mid \mu\rangle = 0``
  para toda historia ``\mu``, de modo que ``H / N = 0``. Cualquier asimetría
  predecible es inmediatamente explotada y borrada por la multitud de agentes:
  el mercado es **impredecible (sin oportunidades de arbitraje)**, pero muy
  volátil.
- **En la fase asimétrica (``\alpha > \alpha_c``):** ``H / N > 0``. Hay más
  historias ``P`` que agentes libres capaces de arbitrarlas, por lo que quedan
  sesgos permanentes ``\langle A \mid \mu\rangle \ne 0`` que una estrategia
  adicional podría explotar.

**El mecanismo microscópico: agentes congelados.** Cuando ``H > 0``, para una
fracción ``\phi(\alpha)`` de los agentes una de sus dos estrategias resulta
sistemáticamente mejor que la otra; su diferencia de puntajes crece sin límite y
esos agentes quedan **congelados** (usan siempre la misma estrategia). Los
únicos que siguen ajustando sus decisiones son los ``(1 - \phi)N`` **agentes
libres**.

Hacer ``H = 0`` requiere cancelar ``P`` promedios ``\langle A \mid \mu\rangle``
usando ``(1 - \phi)N`` agentes libres. La transición ocurre exactamente cuando
el número de variables libres iguala al número de restricciones ``P``:

```math
(1 - \phi(\alpha_c))\,N = P \quad\Longrightarrow\quad 1 - \phi(\alpha_c) = \alpha_c \approx 0.3374 .
```

**Cómo se calculan las curvas.** Challet, Marsili y Zecchina (2000) mostraron
que el estado estacionario del juego minimiza ``H`` en función de las
preferencias medias ``m_i \in [-1, 1]`` de los agentes entre sus dos
estrategias. Por eso no simulamos el juego turno a turno: minimizamos ``H``
con descenso de gradiente proyectado desde ``m_i = 0`` (puntajes iniciales
iguales, como en la simulación), con ``N = 201`` y 12 sorteos de estrategias
por punto. Un agente está congelado si termina en ``m_i = \pm 1``. Con ``N``
finito la caída de ``1 - \phi`` se suaviza y cruza la recta ``\alpha`` un poco a
la derecha de ``\alpha_c``; al crecer ``N`` el cruce se acerca a ``0.3374``.
"""

# ╔═╡ 256527a6-ca8a-11f1-bf5f-dd5016d4f421
md"""
## Resumen

| Concepto / Régimen | Regla o condición | Comportamiento del mercado |
|---|---|---|
| **Juego de la minoría** | ``N`` agentes con ``S=2`` tablas sobre ``P = 2^M`` historias; gana ``-\mathrm{sgn}(A(t))`` | Aprendizaje inductivo sin equilibrio deductivo común |
| **Referente aleatorio** | Elección independiente ``a_i = \pm 1`` con probabilidad ``1/2`` | ``\sigma_{\text{azar}}^2 / N = 1`` |
| **Fase hacinada (``\alpha = 2^M/N < \alpha_c``)** | Pocas historias frente al número de agentes (``P \ll N``) | Efecto manada: ``\sigma^2/N \sim \alpha^{-1} > 1``, pero mercado impredecible (``H/N = 0``) |
| **Punto crítico (``\alpha_c \approx 0.3374``)** | Grados de libertad libres igualan restricciones: ``1 - \phi(\alpha_c) = \alpha_c`` | Mínimo global de volatilidad (``\sigma^2/N \approx 0.25 < 1``) y máxima eficiencia |
| **Fase diluida (``\alpha > \alpha_c``)** | Muchas historias frente al número de agentes (``P \gtrsim N``) | ``\sigma^2/N < 1`` tiende a ``1``; quedan patrones predecibles (``H/N > 0``) y agentes congelados |
"""

# ╔═╡ 256527b0-ca8a-11f1-add5-554f84931308
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `parametro_alfa(M, N)`, que calcule el parámetro de control
``\alpha = 2^M / N`` del juego de la minoría para una memoria `M` y `N` agentes.
"""

# ╔═╡ 256527b0-ca8a-11f1-b208-0b8f68edb9be
function parametro_alfa(M, N)
    missing  # reemplaza esta línea
end

# ╔═╡ 256527ba-ca8a-11f1-b5d1-cb8ad5590aac
md"""
**2.** Escribe `actualizar_puntaje(U, s_mu, A)`, que devuelva el nuevo puntaje
virtual de una estrategia cuyo puntaje actual es `U`, que ante la historia
actual recomendó la acción `s_mu` (``\pm 1``), cuando la acción agregada fue `A`
(usando la regla ``U_{\text{nuevo}} = U - s(\mu)\,A``).
"""

# ╔═╡ 256527c4-ca8a-11f1-80e6-bba6e90ce24f
function actualizar_puntaje(U, s_mu, A)
    missing  # reemplaza esta línea
end

# ╔═╡ 256527ce-ca8a-11f1-bf29-cbd46c091ade
md"""
### Ruta básica: interpretación

3. En el problema del bar El Farol, explica por qué una predicción pública y
   confiable del periódico local ("este jueves irán 45 personas al bar") se
   destruye a sí misma en cuanto los clientes la creen.
4. En el laboratorio de la sección 2 con ``N = 101``, compara lo que ocurre
   cuando ``M = 2`` y cuando ``M = 5``. ¿Por qué tener **menos** memoria
   (``M = 2``) produce fluctuaciones mucho más violentas que tirar monedas al
   azar?
5. En el juego de la minoría las tablas de estrategias de cada agente se sortean
   al azar en ``t = 0`` y quedan fijas (**desorden congelado**). ¿Por qué para
   construir la curva de la sección 3 es indispensable promediar sobre varias
   realizaciones del sorteo de estrategias y no basta con una sola corrida larga
   (repaso S.I4)?

### Ruta avanzada

6. En la fase hacinada (``\alpha < \alpha_c``), la volatilidad es muy alta
   (``\sigma^2/N \gg 1``) pero la predictibilidad es nula (``H/N = 0``). Explica
   por qué estas dos propiedades no se contradicen: ¿por qué un mercado sin
   oportunidades lineales de arbitraje según la historia ``\mu`` puede aun así
   sufrir enormes oscilaciones?
7. Explica qué es un **agente congelado** en el juego de la minoría y por qué
   la fracción de agentes congelados ``\phi(\alpha)`` es casi nula cuando
   ``\alpha < \alpha_c`` pero positiva cuando ``\alpha > \alpha_c``.
8. Usando el argumento de conteo de ecuaciones y variables libres de Challet y
   Marsili, explica por qué la transición entre ``H = 0`` y ``H > 0`` ocurre
   cuando ``1 - \phi(\alpha_c) = \alpha_c``.

## Para la próxima sesión

- Lectura: Slanina (2014), sección 5.1. Complementaria: Challet y Zhang
  (1997); Requião da Cunha (2021), sección 7.1.1; Slanina (2014), sección 5.3.
  Para la ruta avanzada, Slanina (2014), sección 5.2; de consulta: Challet,
  Marsili y Zhang (2005).
- En la **sesión 20** cerramos el curso uniendo todas las piezas: construiremos
  **mercados financieros artificiales con expectativas heterogéneas**
  (**Brock-Hommes** y modelos de **fundamentalistas vs. chartistas**) y les
  aplicaremos el análisis estadístico de las **sesiones 1 y 2** para ver cómo emergen
  las colas pesadas y el agrupamiento de volatilidad.
"""

# ╔═╡ 256527ce-ca8a-11f1-89b2-15555f2a9315
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 256527e4-ca8a-11f1-997c-6bb2a8a3aeb1
begin
    const AZUL       = "#2a78d6"   # serie principal
    const NARANJA    = "#eb6834"   # segunda serie o referencia teórica
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

# ╔═╡ 256527ec-ca8a-11f1-9328-5558d0b1ce09
"""
Simula el bar El Farol con `N` agentes, cupo `L` y memoria `M` (Challet, Marsili
y Ottino, 2004). Cada agente tiene dos tablas que dicen si ir (+1) o quedarse (-1)
según si el bar estuvo cómodo en las últimas `M` semanas; usa la de mayor puntaje,
y una tabla gana un punto si acertó (ir con asistencia ≤ `L`, quedarse con > `L`).
Devuelve la asistencia de las últimas `semanas` semanas.
"""
function simular_el_farol(N, L, M; semanas = 40, calent = 500, rng = Xoshiro(19))
    P = 1 << M
    mascara = P - 1
    est1 = rand(rng, (-1, 1), P, N)
    est2 = rand(rng, (-1, 1), P, N)
    ΔU = 1e-9 .* randn(rng, N)
    μ = rand(rng, 0:mascara)
    asistencia = Int[]
    for t in 1:(calent + semanas)
        idx = μ + 1
        n = 0
        @inbounds for i in 1:N
            n += (ΔU[i] >= 0.0 ? est1[idx, i] : est2[idx, i]) == 1
        end
        g = n <= L ? 1 : -1   # +1 si ir era la buena decisión
        @inbounds for i in 1:N
            ΔU[i] += (est1[idx, i] - est2[idx, i]) * g
        end
        μ = ((μ << 1) & mascara) | (n <= L ? 1 : 0)
        t > calent && push!(asistencia, n)
    end
    asistencia
end

# ╔═╡ 256527ec-ca8a-11f1-a774-4b27fac101a5
"""
Ejemplo motivador de la sección 0 (El Farol con 100 agentes y cupo 60):
compara una población que comparte la misma regla de tendencia con una población
que aprende entre reglas heterogéneas.
"""
ejemplo_el_farol = let
    rng = Xoshiro(7)
    T = 40
    # 1) Todos usan la misma regla con pequeño ruido idiosincrático: ir si la semana pasada asistieron <= 60
    hom = zeros(Int, T)
    asist = 48
    for t in 1:T
        hom[t] = asist
        van = asist <= 60 ? 88 + rand(rng, -4:4) : 14 + rand(rng, -4:4)
        asist = clamp(van, 0, 100)
    end
    # 2) Estrategias heterogéneas inductivas: dos tablas por agente, memoria de 5 semanas
    ind = simular_el_farol(100, 60, 5; semanas = T, calent = 500, rng = Xoshiro(19))
    (homogenea = hom, inductiva = ind)
end;

# ╔═╡ 2565272e-ca8a-11f1-9108-91f549a8b38a
let
    ef = ejemplo_el_farol
    sem = 1:length(ef.homogenea)
    p1 = plot(sem, ef.homogenea; color = NARANJA, lw = 2.0, marker = :circle, ms = 3, msw = 0,
              xlabel = "semana", ylabel = "asistentes (de 100)", ylims = (0, 102),
              title = "Regla idéntica: oscila entre lleno y casi vacío")
    hline!(p1, [60]; color = GRIS, ls = :dash, lw = 1.5)

    p2 = plot(sem, ef.inductiva; color = AZUL, lw = 2.0, marker = :circle, ms = 3, msw = 0,
              xlabel = "semana", ylabel = "asistentes (de 100)", ylims = (0, 102),
              title = "Reglas diversas: se estabiliza en torno a 60")
    hline!(p2, [60]; color = GRIS, ls = :dash, lw = 1.5)
    plot(p1, p2; layout = (1, 2), size = (880, 340), bottom_margin = 4Plots.mm)
end

# ╔═╡ 256527f6-ca8a-11f1-9dcf-3f17a5d8a1ef
"""
Simula una realización del juego de la minoría con `N` agentes, memoria `M` y
`S = 2` estrategias por agente. Devuelve la serie estacionaria de `A(t)`, la
volatilidad reescalada `σ²/N`, el promedio de ganadores, la predictibilidad
`H/N` y la fracción de agentes congelados `ϕ`.
"""
function simular_mg_serie(N, M; pasos = 800, calent = 300, rng = Xoshiro(13))
    P = 1 << M
    mascara = P - 1
    est1 = rand(rng, (-1, 1), P, N)
    est2 = rand(rng, (-1, 1), P, N)
    # Semidiferencia de puntajes ΔU_i = U_{i,1} - U_{i,2} con desempate inicial infinitesimal
    ΔU = 1e-9 .* randn(rng, N)
    suma_σ = zeros(Float64, N)

    μ = rand(rng, 0:mascara)
    As = Vector{Int}(undef, pasos)

    for t in 1:(calent + pasos)
        idx = μ + 1
        A = 0
        @inbounds for i in 1:N
            if ΔU[i] >= 0.0
                A += est1[idx, i]
                t > calent && (suma_σ[i] += 1.0)
            else
                A += est2[idx, i]
                t > calent && (suma_σ[i] -= 1.0)
            end
        end
        @inbounds for i in 1:N
            ΔU[i] -= (est1[idx, i] - est2[idx, i]) * A
        end
        μ = ((μ << 1) & mascara) | (A < 0 ? 1 : 0)
        if t > calent
            As[t - calent] = A
        end
    end

    # Preferencia media m_i ∈ [-1, 1] de cada agente y predictibilidad exacta sobre las P historias
    m = suma_σ ./ pasos
    suma_H = 0.0
    @inbounds for idx in 1:P
        A_med_mu = 0.0
        for i in 1:N
            s_bar = 0.5 * (est1[idx, i] + est2[idx, i])
            xi    = 0.5 * (est1[idx, i] - est2[idx, i])
            A_med_mu += s_bar + xi * m[i]
        end
        suma_H += A_med_mu * A_med_mu
    end

    vol = mean(abs2, As) / N
    gan = mean(ganadores_minoria.(As, N))
    H_norm = suma_H / (P * N)
    congelados = count(x -> abs(x) > 0.90, m) / N
    (A = As, vol = vol, ganadores = gan, H = H_norm, congelados = congelados)
end

# ╔═╡ 2565276a-ca8a-11f1-9934-ff6c0a460d5d
sim_m = simular_mg_serie(101, M_lab; pasos = 600, calent = max(300, 50 * 2^M_lab), rng = Xoshiro(13));

# ╔═╡ 25652772-ca8a-11f1-8aed-a19a432b27ad
let
    s = sim_m
    N = 101
    α = round(2.0^M_lab / N, digits = 2)
    vol = round(s.vol, digits = 2)
    p1 = plot(1:200, s.A[1:200]; color = AZUL, lw = 1.3,
              xlabel = "turno t", ylabel = "exceso de demanda A(t)",
              ylims = (-85, 85), legend = :outerbottom, legend_columns = 2,
              label = "Juego de la minoría (σ²/N = $vol)",
              title = "M = $M_lab (P = $(2^M_lab) historias, α = 2^M/N = $α)")
    hline!(p1, [sqrt(N), -sqrt(N)]; color = NARANJA, ls = :dash, lw = 1.6,
            label = "Banda aleatoria ±√N (σ²/N = 1)")
    hline!(p1, [0.0]; color = GRIS_CLARO, lw = 1, label = "")

    xs = -80:2:80
    p2 = histogram(s.A; bins = -81:4:81, normalize = :pdf, color = AZUL, alpha = 0.55, lw = 0,
                   label = "Distribución empírica de A(t)",
                   xlabel = "exceso de demanda A", ylabel = "densidad",
                   xlims = (-85, 85), legend = :outerbottom,
                   title = α < 0.2 ? "α ≪ α_c: efecto manada (colas anchas)" :
                           α < 1.0 ? "α cerca de α_c: menor volatilidad que el azar" :
                                     "α ≫ α_c: se acerca al azar gaussiano")
    plot!(p2, xs, pdf.(Normal(0, sqrt(N)), xs); color = NARANJA, lw = 2.4,
          label = "Azar independiente Normal(0, √N)")

    plot(p1, p2; layout = (1, 2), size = (900, 390), bottom_margin = 4Plots.mm)
end

# ╔═╡ 2565277e-ca8a-11f1-86bf-91a81a07527b
let
    N = 101
    α = round(2.0^M_lab / N, digits = 2)
    vol = round(sim_m.vol, digits = 2)
    gan = round(sim_m.ganadores, digits = 1)
    rng = Xoshiro(1)
    gan_azar = round(mean(ganadores_minoria(sum(rand(rng, (-1, 1), N)), N) for _ in 1:5000), digits = 1)
    regimen = α < 0.2 ? "**Fase hacinada (`α < α_c`):** con tan pocas historias (`P = $(2^M_lab)`), muchos agentes tienen exactamente las mismas tablas y cambian de estrategia en bloque. Esas estampidas disparan la volatilidad a `σ²/N = $vol` (peor que tirar monedas, que da 1)." :
              α < 1.0 ? "**Región eficiente (`α` cerca de `α_c`):** el número de historias (`P = $(2^M_lab)`) está equilibrado con el número de agentes (`N = 101`). Los agentes se reparten entre ``+1`` y ``-1`` con volatilidad `σ²/N = $vol < 1`, logrando que ganen en promedio **$gan agentes por turno** (frente a $gan_azar al azar)." :
                           "**Fase diluida (`α > α_c`):** hay tantas historias posibles (`P = $(2^M_lab)`) que las estrategias de distintos agentes casi no se parecen y el mercado se comporta casi como si eligieran al azar (`σ²/N = $vol ≈ 1`)."
    Markdown.parse(regimen)
end

# ╔═╡ 25652800-ca8a-11f1-bec9-1ba9627d13f0
"Barrido de la volatilidad σ²/N frente a M y α = 2^M / N para varios tamaños de población N."
function barrido_escala_mg(Ns = [51, 101, 201]; Ms = 2:10, R = 20)
    series = map(Ns) do N
        vols, ses = Float64[], Float64[]
        for M in Ms
            P = 2^M
            vs = [simular_mg_serie(N, M; pasos = max(600, 50P), calent = max(250, 50P),
                                   rng = Xoshiro(1000 * N + 50 * M + r)).vol for r in 1:R]
            push!(vols, mean(vs))
            push!(ses, std(vs) / sqrt(R))
        end
        (N = N, Ms = collect(Ms), alfas = (2.0 .^ collect(Ms)) ./ N, vol = vols, se = ses)
    end
    (series = series,)
end

# ╔═╡ 25652788-ca8a-11f1-9bf2-236bf4e4a65c
barrido_mg = barrido_escala_mg([51, 101, 201]; Ms = 2:10, R = 20);

# ╔═╡ 25652788-ca8a-11f1-ac7a-a737754832a1
let
    b = barrido_mg
    colores = (AZUL, NARANJA, GRIS)
    marcadores = (:circle, :square, :utriangle)

    p1 = plot(; xlabel = "memoria M", ylabel = "volatilidad reescalada σ² / N",
              yscale = :log10, ylims = (0.18, 25), legend = :outerbottom, legend_columns = 3,
              title = "Sin reescalar: el mínimo se desplaza con N")
    hline!(p1, [1.0]; color = GRIS_CLARO, ls = :dash, lw = 1.5, label = "")

    p2 = plot(; xlabel = "parámetro de control α = 2^M / N", ylabel = "volatilidad reescalada σ² / N",
              xscale = :log10, yscale = :log10, xlims = (0.015, 25), ylims = (0.18, 25),
              legend = :outerbottom, legend_columns = 2,
              title = "Colapso de escala contra α = 2^M / N")
    hline!(p2, [1.0]; color = GRIS, ls = :dot, lw = 1.5, label = "Azar (σ²/N = 1)")
    vline!(p2, [0.3374]; color = NARANJA, alpha = 0.45, ls = :dash, lw = 1.5,
           label = "Punto crítico α_c ≈ 0.34")

    for (k, d) in enumerate(b.series)
        plot!(p1, d.Ms, d.vol; color = colores[k], lw = 2.0, marker = marcadores[k],
              ms = 4.5, msw = 0, label = "N = $(d.N)")
        plot!(p2, d.alfas, d.vol; yerror = d.se, color = colores[k], lw = 2.0, msc = colores[k],
              marker = marcadores[k], ms = 4.5, label = "N = $(d.N)")
    end

    plot(p1, p2; layout = (1, 2), size = (900, 410), bottom_margin = 4Plots.mm)
end

# ╔═╡ 25652800-ca8a-11f1-b6e2-c3c53805fb20
"""
Curvas estacionarias de predictibilidad H/N y fracción de agentes congelados ϕ(α)
frente a α = P/N. Minimiza H(m) = (1/P) Σ_μ (Ω(μ) + Σ_i ξ_i(μ) m_i)² sobre m ∈ [-1, 1]^N
por descenso de gradiente proyectado desde m = 0 (Challet, Marsili y Zecchina, 2000).
"""
curvas_transicion_mg = let
    N = 201
    Ps = [10, 20, 34, 48, 58, 64, 70, 78, 90, 120, 180, 300, 500, 800]
    R = 12
    alfas = Ps ./ N
    Hs, cong = Float64[], Float64[]
    for P in Ps
        h_acum = 0.0
        c_acum = 0.0
        for r in 1:R
            rng = Xoshiro(7000 + 31 * P + r)
            est1 = rand(rng, (-1.0, 1.0), P, N)
            est2 = rand(rng, (-1.0, 1.0), P, N)
            Ω = vec(0.5 .* sum(est1 .+ est2; dims = 2))
            ξ = 0.5 .* (est1 .- est2)
            η = 0.9 * P / opnorm(ξ)^2   # paso estable: el gradiente es 2ξ'(Ω + ξm)/P
            m = zeros(Float64, N)
            for _ in 1:4000
                m .= clamp.(m .- η .* (ξ' * (Ω .+ ξ * m)) ./ P, -1.0, 1.0)
            end
            A_med = Ω .+ ξ * m
            h_acum += sum(abs2, A_med) / (P * N)
            c_acum += count(x -> abs(x) > 1 - 1e-6, m) / N
        end
        push!(Hs, h_acum / R)
        push!(cong, c_acum / R)
    end
    (alfas = alfas, H = Hs, congelados = cong, libres = 1.0 .- cong)
end;

# ╔═╡ 2565279c-ca8a-11f1-9f49-b9778dc06f62
let
    t = curvas_transicion_mg
    p1 = plot(t.alfas, t.H; color = AZUL, lw = 2.4, marker = :circle, ms = 4.5, msw = 0,
              xscale = :log10, xlabel = "parámetro de control α = P / N", ylims = (-0.01, 0.42),
              ylabel = "predictibilidad H / N", legend = :outerbottom,
              label = "Predictibilidad estacionaria H / N",
              title = "Parámetro de orden: H / N = 0 si α < α_c")
    vline!(p1, [0.3374]; color = NARANJA, ls = :dash, lw = 1.6, label = "Umbral crítico α_c ≈ 0.34")
    hline!(p1, [0.0]; color = GRIS_CLARO, lw = 1, label = "")

    α_lin = range(0.05, 0.95, length = 100)
    p2 = plot(t.alfas, t.libres; color = AZUL, lw = 2.4, marker = :square, ms = 4.5, msw = 0,
              xscale = :log10, ylims = (0.0, 1.05),
              xlabel = "parámetro de control α = P / N",
              ylabel = "fracción", legend = :outerbottom, legend_columns = 2,
              label = "Agentes libres 1 − ϕ(α)",
              title = "Transición donde 1 − ϕ(α) = α")
    plot!(p2, t.alfas, t.congelados; color = GRIS, lw = 1.8, ls = :dot, marker = :circle, ms = 3.5, msw = 0,
          label = "Agentes congelados ϕ(α)")
    plot!(p2, α_lin, α_lin; color = NARANJA, lw = 2.0, ls = :dash, label = "Recta α = P / N")
    vline!(p2, [0.3374]; color = NARANJA, alpha = 0.4, ls = :dash, label = "")

    plot(p1, p2; layout = (1, 2), size = (900, 400), bottom_margin = 4Plots.mm)
end

# ╔═╡ 2565280a-ca8a-11f1-941d-f172b7d79f65
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

# ╔═╡ 256527ba-ca8a-11f1-9c1c-198e895188cc
verificar(() -> [parametro_alfa(2, 100), parametro_alfa(5, 100), parametro_alfa(8, 200)],
          [4 / 100, 32 / 100, 256 / 200];
          pista = md"Usa `(2.0 ^ M) / N`.")

# ╔═╡ 256527c4-ca8a-11f1-bda1-6de7e1bc1836
verificar(() -> [actualizar_puntaje(10.0, 1, 15), actualizar_puntaje(10.0, -1, 15), actualizar_puntaje(0.0, 1, -7)],
          [-5.0, 25.0, 7.0];
          pista = md"Si `s_mu` y `A` tienen signos opuestos el puntaje sube: `U - s_mu * A`.")

# ╔═╡ Cell order:
# ╟─256526f2-ca8a-11f1-884c-2ff3f7ee1b1d
# ╟─25652724-ca8a-11f1-a7b1-89bb41b9e06d
# ╟─25652724-ca8a-11f1-8229-df8198032751
# ╟─2565272e-ca8a-11f1-9108-91f549a8b38a
# ╟─25652738-ca8a-11f1-bd98-a516e4f10e8a
# ╟─25652740-ca8a-11f1-9cfc-cfaff8882fe3
# ╠═2565274c-ca8a-11f1-ab95-41d48413308c
# ╠═25652756-ca8a-11f1-8981-b7f20b911d8e
# ╟─25652756-ca8a-11f1-a181-a580ba7c492b
# ╠═2565276a-ca8a-11f1-9934-ff6c0a460d5d
# ╟─25652772-ca8a-11f1-8aed-a19a432b27ad
# ╟─2565277e-ca8a-11f1-86bf-91a81a07527b
# ╟─2565277e-ca8a-11f1-830a-edf9ff6800aa
# ╠═25652788-ca8a-11f1-9bf2-236bf4e4a65c
# ╟─25652788-ca8a-11f1-ac7a-a737754832a1
# ╟─25652792-ca8a-11f1-97c9-053ccff4d8d6
# ╟─2565279c-ca8a-11f1-b39c-e50ac452a555
# ╟─2565279c-ca8a-11f1-9f49-b9778dc06f62
# ╟─256527a6-ca8a-11f1-bf5f-dd5016d4f421
# ╟─256527b0-ca8a-11f1-add5-554f84931308
# ╠═256527b0-ca8a-11f1-b208-0b8f68edb9be
# ╟─256527ba-ca8a-11f1-9c1c-198e895188cc
# ╟─256527ba-ca8a-11f1-b5d1-cb8ad5590aac
# ╠═256527c4-ca8a-11f1-80e6-bba6e90ce24f
# ╟─256527c4-ca8a-11f1-bda1-6de7e1bc1836
# ╟─256527ce-ca8a-11f1-bf29-cbd46c091ade
# ╟─256527ce-ca8a-11f1-89b2-15555f2a9315
# ╠═256527d8-ca8a-11f1-8f74-d9bef392ffdd
# ╠═256527e4-ca8a-11f1-997c-6bb2a8a3aeb1
# ╠═256527ec-ca8a-11f1-9328-5558d0b1ce09
# ╠═256527ec-ca8a-11f1-a774-4b27fac101a5
# ╠═256527f6-ca8a-11f1-9dcf-3f17a5d8a1ef
# ╠═25652800-ca8a-11f1-bec9-1ba9627d13f0
# ╠═25652800-ca8a-11f1-b6e2-c3c53805fb20
# ╠═2565280a-ca8a-11f1-941d-f172b7d79f65

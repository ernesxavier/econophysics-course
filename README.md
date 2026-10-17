# Econofísica

Material del curso de Econofísica (licenciatura en Economía, Física o
Matemáticas), con simulaciones en Julia.

## Estructura

```
├── programa/
│   ├── programa.md            → Programa rediseñado (15 semanas)
│   └── programa-anterior.md   → Programa del semestre anterior
├── Class_notes/
│   ├── s01_hechos_estilizados.md → Notas formales, sesión 1
│   ├── s02_colas_pesadas.md   → Notas formales, sesión 2
│   ├── s03_procesos_estocasticos.md → Notas formales, sesión 3
│   ├── s04_crecimiento_multiplicativo.md → Notas formales, sesión 4
│   ├── s05_correlaciones.md   → Notas formales, sesión 5
│   ├── s06_intercambio.md     → Notas formales, sesión 6
│   ├── s07_condensacion.md    → Notas formales, sesión 7
│   ├── s08_redes.md           → Notas formales, sesión 8
│   ├── s09_contagio.md        → Notas formales, sesión 9
│   ├── s10_ising.md           → Notas formales, sesión 10
│   └── repaso/                → Repaso de matemáticas por niveles (índice por semana)
│       └── figuras/           → Figuras de las notas (generadas con Code/figuras_repaso.jl)
├── docs/
│   ├── diagnostico.md         → Revisión del curso anterior
│   └── benchmark.md           → Econofísica en otras universidades
├── Data/
│   ├── precios.csv            → IPC, USD/MXN, S&P 500 y Bitcoin (2000–sep 2026)
│   ├── sp100_precios.csv      → 95 acciones del S&P 100 (2015–sep 2026)
│   ├── sp100_sectores.csv     → Nombre y sector GICS de esas acciones
│   ├── enigh2024_hogares.csv  → Ingreso mensual por hogar, ENIGH 2024 (INEGI)
│   ├── condensacion_bm.csv    → Barrido simulado del modelo de Bouchaud-Mézard (sesión 7)
│   ├── localidades_2020.csv   → Población por localidad, Censo 2020 (INEGI), ≥ 100 hab.
│   └── municipios_2010_2020.csv → Población por municipio, Censos 2010 y 2020 (INEGI)
├── Code/                      → Todo el código Julia, un solo entorno
│   ├── s01_hechos_estilizados.jl → Cuaderno de Pluto, sesión 1 (programa nuevo)
│   ├── s02_colas_pesadas.jl   → Cuaderno de Pluto, sesión 2
│   ├── s03_procesos_estocasticos.jl → Cuaderno de Pluto, sesión 3
│   ├── s04_crecimiento_multiplicativo.jl → Cuaderno de Pluto, sesión 4
│   ├── s05_correlaciones.jl   → Cuaderno de Pluto, sesión 5
│   ├── s06_intercambio.jl     → Cuaderno de Pluto, sesión 6
│   ├── s07_condensacion.jl    → Cuaderno de Pluto, sesión 7
│   ├── s08_redes.jl           → Cuaderno de Pluto, sesión 8
│   ├── s09_contagio.jl        → Cuaderno de Pluto, sesión 9
│   ├── s10_ising.jl           → Cuaderno de Pluto, sesión 10
│   ├── descargar_datos.jl     → Regenera Data/precios.csv desde Yahoo Finance
│   ├── descargar_sp100.jl     → Regenera los datos del S&P 100 (Yahoo Finance)
│   ├── descargar_enigh.jl     → Regenera los datos de la ENIGH 2024 (INEGI)
│   ├── simular_condensacion.jl → Regenera Data/condensacion_bm.csv (tarda minutos)
│   ├── descargar_inegi.jl     → Regenera los datos de INEGI (localidades y municipios)
│   ├── figuras_repaso.jl      → Regenera las figuras de Class_notes/repaso/figuras/
│   │   Notebooks del semestre anterior:
│   ├── sesion_1.ipynb         → Caminata aleatoria, colas pesadas, TLC, datos financieros
│   ├── sesion_3.ipynb         → Redes: medidas básicas, familia Medici
│   ├── sesion_4.ipynb         → Redes aleatorias: ER, WS, BA
│   ├── sesion_5.ipynb         → Teoría de juegos clásica y evolutiva
│   ├── sesion_6.ipynb         → Juegos de coordinación en redes
│   ├── sesion_7.ipynb         → Dilema del prisionero en red
│   ├── sesion_8.ipynb         → Introducción a modelos de agentes: Schelling
│   ├── sesion_9.ipynb         → Modelo de Axelrod
│   ├── Project.toml
│   └── Manifest.toml
├── src/                       → Fuentes de los cuadernos de Pluto (formato #%%)
│   ├── construir_cuaderno.jl  → Genera y ejecuta un cuaderno a partir de su fuente
│   └── sNN_tema.jl            → Fuente de Code/sNN_tema.jl (sesiones 2 a 10)
└── README.md
```

## Requisitos

- Julia 1.13 (cuadernos de Pluto probados con Pluto 1.0.4)

## Entorno de Julia

Todo el repositorio usa **un solo entorno**: `Code/Project.toml` +
`Code/Manifest.toml`. Incluye Pluto, así que no hace falta instalar nada en el
entorno global. Los cuadernos de Pluto activan este entorno en su primera
celda de paquetes (`Pkg.activate(@__DIR__)`).

## Uso

Desde la raíz del repositorio, la primera vez:

```
julia --project=Code -e 'using Pkg; Pkg.instantiate()'
```

**Cuadernos de Pluto:**

```
julia --project=Code -e 'using Pluto; Pluto.run()'
```

y en el navegador abrir `Code/s01_hechos_estilizados.jl`.

**Notebooks de Jupyter del semestre anterior:** abrir `Code/sesion_*.ipynb` en
VS Code o Jupyter con el entorno `Code`.

**Construir un cuaderno desde su fuente.** Cada cuaderno de Pluto tiene una
fuente en `src/` con marcadores `#%%` (celda visible) y `#%% hide` (código
oculto), más fácil de editar como texto. Este comando genera el cuaderno, lo
ejecuta sin interfaz, reporta las celdas con error y guarda las gráficas en
la carpeta de salidas para revisarlas:

```
julia --project=Code src/construir_cuaderno.jl src/s02_colas_pesadas.jl Code/s02_colas_pesadas.jl /tmp/salidas
```

La fuente y el cuaderno deben coincidir. Los cambios grandes se hacen en la
fuente y se regenera el cuaderno; si se edita el cuaderno en Pluto, hay que
pasar el cambio a la fuente antes de volver a generarlo.

**Actualizar los datos:**

```
julia --project=Code Code/descargar_datos.jl
julia --project=Code Code/descargar_sp100.jl
julia --project=Code Code/descargar_enigh.jl
```

**Actualizar paquetes** (por ejemplo, al cambiar de versión de Julia):

```
julia -e 'using Pkg; Pkg.activate("Code"); Pkg.update()'
```

## Pendientes

- [x] Sesión 1: hechos estilizados (Pluto).
- [ ] Notebooks nuevos: estimación de colas, procesos
      estocásticos, crecimiento multiplicativo, matrices aleatorias, modelos
      cinéticos, Bouchaud-Mézard, contagio, Ising, votante y Kirman, minoría,
      Brock-Hommes.
- [ ] Módulo compartido con funciones comunes (`gini`, `ccdf`, `hill`, ...). Ya no puede ir en `src/`, que ahora guarda las fuentes de los cuadernos; por ejemplo, `Code/comun.jl`.
- [x] Datos en CSV para no depender de la API de Yahoo.
- [ ] Migrar el resto de las sesiones a Pluto.jl.
- [ ] Notebooks anteriores: en `sesion_1` la celda 31 grafica `Normal(...)`
      antes de cargar `StatsPlots`; en `sesion_8` la celda 16 da error
      (`agentes.feliz`), igual que en el original.
- [ ] Notebooks con huecos y ejercicios separados en ruta básica y avanzada.
- [ ] Examen diagnóstico de matemáticas (semana 1).
- [ ] Lista de artículos para presentaciones.

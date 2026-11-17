# Econofísica

Material del curso de Econofísica (licenciatura en Economía, Física o
Matemáticas), con simulaciones en Julia.

**Autor:** Ernesto X. Vera Gómez, Centro de Investigaciones Socioeconómicas (CISE),
Universidad Autónoma de Coahuila.

## Estructura

```
├── programa/
│   ├── programa.md            → Programa rediseñado (20 sesiones; semestre o trimestre)
│   └── proyecto_final.md      → Guía y rúbrica del proyecto final
├── Class_notes/
│   ├── s01_retornos.md            → Notas formales, sesión 1
│   ├── s02_volatilidad.md         → Notas formales, sesión 2
│   ├── s03_colas_pesadas.md       → Notas formales, sesión 3
│   ├── s04_caminata_markov.md     → Notas formales, sesión 4
│   ├── s05_gbm_garch.md           → Notas formales, sesión 5
│   ├── s06_crecimiento_multiplicativo.md → Notas formales, sesión 6
│   ├── s07_correlaciones.md       → Notas formales, sesión 7
│   ├── s08_intercambio.md         → Notas formales, sesión 8
│   ├── s09_condensacion.md        → Notas formales, sesión 9
│   ├── s10_redes_modelos.md       → Notas formales, sesión 10
│   ├── s11_redes_mercado.md       → Notas formales, sesión 11
│   ├── s12_contagio.md            → Notas formales, sesión 12
│   ├── s13_logit_ising.md         → Notas formales, sesión 13
│   ├── s14_campo_medio.md         → Notas formales, sesión 14
│   ├── s15_votante_kirman.md      → Notas formales, sesión 15
│   ├── s16_schelling.md           → Notas formales, sesión 16
│   ├── s17_replicador_goodwin.md  → Notas formales, sesión 17
│   ├── s18_cooperacion_redes.md   → Notas formales, sesión 18
│   ├── s19_juego_minoria.md       → Notas formales, sesión 19
│   ├── s20_mercados_artificiales.md → Notas formales, sesión 20
│   └── repaso/                → Repaso de matemáticas por niveles (índice por sesión)
│       ├── diagnostico_inicial.md → Diagnóstico de matemáticas de la sesión 1
│       └── figuras/           → Figuras de las notas (generadas con Code/figuras_repaso.jl)
├── docs/
│   └── benchmark.md           → Econofísica en otras universidades
├── Data/
│   ├── precios.csv            → IPC, USD/MXN, S&P 500 y Bitcoin (2000–sep 2026); no incluido, ver Uso
│   ├── sp100_precios.csv      → 95 acciones del S&P 100 (2015–sep 2026); no incluido, ver Uso
│   ├── sp100_sectores.csv     → Nombre y sector GICS de esas acciones
│   ├── enigh2024_hogares.csv  → Ingreso mensual por hogar, ENIGH 2024 (INEGI)
│   ├── condensacion_bm.csv    → Barrido simulado del modelo de Bouchaud-Mézard (sesión 9)
│   ├── condensacion_impuesto.csv → Barrido con impuesto redistribuido (sesión 9)
│   ├── localidades_2020.csv   → Población por localidad, Censo 2020 (INEGI), ≥ 100 hab.
│   └── municipios_2010_2020.csv → Población por municipio, Censos 2010 y 2020 (INEGI)
├── Code/                      → Todo el código Julia, un solo entorno
│   ├── s01_retornos.jl            → Cuaderno de Pluto, sesión 1
│   ├── s02_volatilidad.jl         → Cuaderno de Pluto, sesión 2
│   ├── s03_colas_pesadas.jl       → Cuaderno de Pluto, sesión 3
│   ├── s04_caminata_markov.jl     → Cuaderno de Pluto, sesión 4
│   ├── s05_gbm_garch.jl           → Cuaderno de Pluto, sesión 5
│   ├── s06_crecimiento_multiplicativo.jl → Cuaderno de Pluto, sesión 6
│   ├── s07_correlaciones.jl       → Cuaderno de Pluto, sesión 7
│   ├── s08_intercambio.jl         → Cuaderno de Pluto, sesión 8
│   ├── s09_condensacion.jl        → Cuaderno de Pluto, sesión 9
│   ├── s10_redes_modelos.jl       → Cuaderno de Pluto, sesión 10
│   ├── s11_redes_mercado.jl       → Cuaderno de Pluto, sesión 11
│   ├── s12_contagio.jl            → Cuaderno de Pluto, sesión 12
│   ├── s13_logit_ising.jl         → Cuaderno de Pluto, sesión 13
│   ├── s14_campo_medio.jl         → Cuaderno de Pluto, sesión 14
│   ├── s15_votante_kirman.jl      → Cuaderno de Pluto, sesión 15
│   ├── s16_schelling.jl           → Cuaderno de Pluto, sesión 16
│   ├── s17_replicador_goodwin.jl  → Cuaderno de Pluto, sesión 17
│   ├── s18_cooperacion_redes.jl   → Cuaderno de Pluto, sesión 18
│   ├── s19_juego_minoria.jl       → Cuaderno de Pluto, sesión 19
│   ├── s20_mercados_artificiales.jl → Cuaderno de Pluto, sesión 20
│   ├── taller1_hechos_estilizados.jl → Taller para casa del Módulo I
│   ├── taller2_desigualdad_contagio.jl → Taller para casa del Módulo II
│   ├── taller3_cierre_ciclo.jl    → Taller para casa del Módulo III
│   ├── descargar_datos.jl     → Regenera Data/precios.csv desde Yahoo Finance
│   ├── descargar_sp100.jl     → Regenera los datos del S&P 100 (Yahoo Finance)
│   ├── descargar_enigh.jl     → Regenera los datos de la ENIGH 2024 (INEGI)
│   ├── simular_condensacion.jl → Regenera los datos de condensación de la sesión 9 (tarda minutos)
│   ├── descargar_inegi.jl     → Regenera los datos de INEGI (localidades y municipios)
│   ├── figuras_repaso.jl      → Regenera las figuras de Class_notes/repaso/figuras/
│   ├── Project.toml
│   └── Manifest.toml
├── src/                       → Fuentes de los cuadernos de Pluto (formato #%%)
│   ├── construir_cuaderno.jl  → Genera y ejecuta un cuaderno a partir de su fuente
│   ├── revisar_material.jl    → Revisa claves de repaso, enlaces y referencias a sesiones
│   ├── revisar_formulas.jl    → Revisa que GitHub muestre bien las fórmulas de los .md
│   ├── sNN_tema.jl            → Fuente de Code/sNN_tema.jl (20 fuentes, sesiones 1 a 20)
│   └── tallerN_tema.jl        → Fuente de Code/tallerN_tema.jl (talleres de los tres módulos)
├── CITATION.cff           → Cómo citar el curso (botón "Cite this repository" de GitHub)
├── LICENSE                → Licencia MIT (código)
├── LICENSE-CC-BY-4.0      → Licencia CC BY 4.0 (notas, programa y documentos)
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

Desde la raíz del repositorio, la primera vez, instalar los paquetes y
descargar los precios de Yahoo Finance (tarda unos minutos):

```
julia --project=Code -e 'using Pkg; Pkg.instantiate()'
julia --project=Code Code/descargar_datos.jl
julia --project=Code Code/descargar_sp100.jl
```

Los precios no vienen en el repositorio porque los términos de Yahoo Finance
no permiten redistribuirlos. Las fechas de descarga están fijas, así que todos
obtienen los mismos rendimientos. Los datos de INEGI y las simulaciones sí
vienen incluidos en `Data/`.

**Cuadernos de Pluto:**

```
julia --project=Code -e 'using Pluto; Pluto.run()'
```

y en el navegador abrir `Code/s01_retornos.jl`.

**Construir un cuaderno desde su fuente.** Cada cuaderno de Pluto tiene una
fuente en `src/` con marcadores `#%%` (celda visible) y `#%% hide` (código
oculto), más fácil de editar como texto. Este comando genera el cuaderno, lo
ejecuta sin interfaz, reporta las celdas con error y guarda las gráficas en
la carpeta de salidas para revisarlas:

```
julia --project=Code src/construir_cuaderno.jl src/s03_colas_pesadas.jl Code/s03_colas_pesadas.jl /tmp/salidas
```

La fuente y el cuaderno deben coincidir. Los cambios grandes se hacen en la
fuente y se regenera el cuaderno; si se edita el cuaderno en Pluto, hay que
pasar el cambio a la fuente antes de volver a generarlo.

**Regenerar los datos de INEGI** (ya vienen incluidos):

```
julia --project=Code Code/descargar_enigh.jl
julia --project=Code Code/descargar_inegi.jl
```

**Revisar el material** antes de cada commit: la consistencia (claves del
repaso, enlaces y números de sesión) y las fórmulas tal como las muestra
GitHub. La segunda revisión usa `gh` con sesión iniciada y tarda unos
20 segundos; también acepta archivos sueltos.

```
julia src/revisar_material.jl
julia src/revisar_formulas.jl
julia src/revisar_formulas.jl Class_notes/s04_caminata_markov.md
```

**Actualizar paquetes** (por ejemplo, al cambiar de versión de Julia):

```
julia -e 'using Pkg; Pkg.activate("Code"); Pkg.update()'
```

## Uso de inteligencia artificial

Se utilizó Claude (Anthropic) como apoyo en la revisión del código y del
material del curso. El contenido, las decisiones y cualquier error son
responsabilidad del autor.

## Licencia

© 2026 Ernesto X. Vera Gómez, Centro de Investigaciones Socioeconómicas (CISE),
Universidad Autónoma de Coahuila.

- **Notas, programa y documentos** (`Class_notes/`, `programa/`, `docs/`,
  este README y las figuras): [CC BY 4.0](LICENSE-CC-BY-4.0). Puedes
  copiarlos, adaptarlos y redistribuirlos, incluso con fines comerciales,
  siempre que cites la fuente e indiques si hiciste cambios.
- **Código** (`Code/` y `src/`, incluidos los cuadernos de Pluto y Jupyter
  con su texto): [MIT](LICENSE).
- **Datos** (`Data/`): los descargados conservan los términos de su fuente
  (INEGI, Yahoo Finance); los generados por simulación (`condensacion_*.csv`)
  están bajo CC BY 4.0.

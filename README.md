# Econofísica

Material del curso de Econofísica (licenciatura en Economía, Física o
Matemáticas), con simulaciones en Julia.

## Estructura

```
├── programa/
│   ├── programa.md            → Programa rediseñado (15 semanas)
│   └── programa-anterior.md   → Programa del semestre anterior
├── docs/
│   ├── diagnostico.md         → Revisión del curso anterior
│   └── benchmark.md           → Econofísica en otras universidades
├── Data/
│   └── precios.csv            → IPC, USD/MXN, S&P 500 y Bitcoin (2000–sep 2026)
├── Code/                      → Todo el código Julia, un solo entorno
│   ├── s01_hechos_estilizados.jl → Cuaderno de Pluto, sesión 1 (programa nuevo)
│   ├── descargar_datos.jl     → Regenera Data/precios.csv desde Yahoo Finance
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

**Actualizar los datos:**

```
julia --project=Code Code/descargar_datos.jl
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
- [ ] Módulo compartido `src/` con funciones comunes (`gini`, `ccdf`, `hill`, ...).
- [x] Datos en CSV para no depender de la API de Yahoo.
- [ ] Migrar el resto de las sesiones a Pluto.jl.
- [ ] Notebooks anteriores: en `sesion_1` la celda 31 grafica `Normal(...)`
      antes de cargar `StatsPlots`; en `sesion_8` la celda 16 da error
      (`agentes.feliz`), igual que en el original.
- [ ] Notebooks con huecos y ejercicios separados en ruta básica y avanzada.
- [ ] Examen diagnóstico de matemáticas (semana 1).
- [ ] Lista de artículos para presentaciones.

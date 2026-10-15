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
├── Code/                      → Notebooks del semestre anterior
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

- Julia ≥ 1.10 (los notebooks se ejecutaron con 1.12)
- Paquetes en `Code/Project.toml`

## Uso

```julia
using Pkg; Pkg.activate("Code"); Pkg.instantiate()
```

Después, abrir cualquier `Code/sesion_*.ipynb` en VS Code o Jupyter.

## Pendientes

- [ ] Notebooks nuevos: hechos estilizados con estimación de colas, procesos
      estocásticos, crecimiento multiplicativo, matrices aleatorias, modelos
      cinéticos, Bouchaud-Mézard, contagio, Ising, votante y Kirman, minoría,
      Brock-Hommes.
- [ ] Módulo compartido `src/` con funciones comunes (`gini`, `ccdf`, `hill`, ...).
- [ ] Datos en CSV para no depender de la API de Yahoo.
- [ ] Migrar a Pluto.jl (deslizadores para las transiciones de fase).
- [ ] Notebooks con huecos y ejercicios separados en ruta básica y avanzada.
- [ ] Examen diagnóstico de matemáticas (semana 1).
- [ ] Lista de artículos para presentaciones.

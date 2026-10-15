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
├── sesiones/                  → Cuadernos de Pluto del programa nuevo
│   └── s01_hechos_estilizados.jl
├── Data/
│   └── precios.csv            → IPC, USD/MXN, S&P 500 y Bitcoin (2000–sep 2026)
├── Code/                      → Código Julia con entorno compartido
│   ├── descargar_datos.jl     → Regenera Data/precios.csv desde Yahoo Finance
│   ├── sesion_*.ipynb         → Notebooks del semestre anterior:
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

## Entornos de Julia

- **`Code/`**: un solo `Project.toml` + `Manifest.toml` para los scripts y los
  notebooks de Jupyter.
- **Cuadernos de Pluto**: cada `sesiones/*.jl` guarda sus propios paquetes al
  final del archivo. Pluto los instala solo; no requieren `activate` ni
  `instantiate`.

## Uso

### Cuadernos de Pluto (programa nuevo)

Sólo hace falta instalar Pluto una vez:

```julia
using Pkg; Pkg.add("Pluto")
using Pluto; Pluto.run()
```

En el navegador, abrir `sesiones/s01_hechos_estilizados.jl`. Cada cuaderno
trae sus paquetes con versiones fijas; Pluto los instala la primera vez que se
abre (tarda unos minutos).

Para actualizar los datos:

```
julia --project=Code Code/descargar_datos.jl
```

### Notebooks del semestre anterior

```julia
using Pkg; Pkg.activate("Code"); Pkg.instantiate()
```

Después, abrir cualquier `Code/sesion_*.ipynb` en VS Code o Jupyter.

## Pendientes

- [x] Sesión 1: hechos estilizados (Pluto).
- [ ] Notebooks nuevos: estimación de colas, procesos
      estocásticos, crecimiento multiplicativo, matrices aleatorias, modelos
      cinéticos, Bouchaud-Mézard, contagio, Ising, votante y Kirman, minoría,
      Brock-Hommes.
- [ ] Módulo compartido `src/` con funciones comunes (`gini`, `ccdf`, `hill`, ...).
- [x] Datos en CSV para no depender de la API de Yahoo.
- [ ] Migrar el resto de las sesiones a Pluto.jl.
- [ ] Notebooks con huecos y ejercicios separados en ruta básica y avanzada.
- [ ] Examen diagnóstico de matemáticas (semana 1).
- [ ] Lista de artículos para presentaciones.

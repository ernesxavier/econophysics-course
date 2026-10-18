# Proyecto final: réplica y extensión de un artículo

El proyecto vale el **35%** de la calificación. En equipos de dos o tres
personas, de preferencia mixtos (Economía con Física o Matemáticas), eligen un
artículo de econofísica, **reproducen** su resultado central con las
herramientas del curso y lo **extienden** en una dirección propia.

El objetivo no es descubrir algo nuevo, sino hacer bien lo que hace un
investigador antes de confiar en un resultado: volver a obtenerlo, entender
de qué depende y ver hasta dónde llega.

## Qué significa replicar y extender

- **Replicar** es reproducir una figura o una tabla central del artículo, con
  su propio código, y explicar las diferencias que queden. Una réplica
  honesta que no sale igual, con una buena explicación de por qué, vale tanto
  como una que sale idéntica.
- **Extender** es hacer una sola pregunta nueva y contestarla bien. Por
  ejemplo:
  - aplicar el método a **datos mexicanos** (los de `Data/` u otros);
  - **cambiar un supuesto** del modelo y ver qué hechos sobreviven;
  - **barrer un parámetro** que el artículo deja fijo, promediando sobre
    varias semillas;
  - comparar el modelo con un **modelo nulo** (datos barajados, caminata
    gaussiana, red aleatoria con los mismos grados).

## Etapas y fechas

Las fechas se cuentan por sesiones, para que valgan en los dos calendarios
del programa.

| Etapa | Cuándo | Qué se entrega | Peso dentro del proyecto |
|---|---|---|---|
| 1. Equipo y artículo | Al terminar la sesión 10 | Nombres y artículo elegido (de la lista o propio, con visto bueno del profesor) | — |
| 2. Propuesta | Al terminar la sesión 12 | Una página: qué figura se replica, con qué datos o modelo, y la pregunta de la extensión | 5% |
| 3. Avance | Al terminar la sesión 17 | El cuaderno con la réplica funcionando, aunque no esté pulida | 10% |
| 4. Asesorías | Semestre: las dos sesiones de asesoría del calendario. Trimestre: 15 minutos por equipo en horario de consulta, después del avance | Dudas y un plan para la extensión | — |
| 5. Entrega final | El día de la presentación | Escrito, cuaderno y diapositivas | 85% |

## Qué se entrega

**1. El escrito** (8 a 12 páginas, sin contar referencias ni apéndices):

1. *Introducción:* la pregunta del artículo y por qué importa, en una página.
2. *El modelo o el método:* lo necesario para entender la réplica, con las
   ecuaciones centrales y sin derivaciones largas.
3. *Réplica:* la figura o tabla original junto a la de ustedes, y las
   diferencias explicadas.
4. *Extensión:* la pregunta nueva, el diseño del experimento y el resultado.
5. *Discusión:* qué aprendieron, qué limita las conclusiones y qué harían
   después.
6. *Contribución de cada integrante,* en un párrafo.

**2. El cuaderno de Pluto,** que debe correr de principio a fin en el entorno
del curso (`Code/Project.toml`) sin intervención:

- semillas fijas en todas las simulaciones;
- los datos en `Data/` o un script que los descargue;
- menos de 10 minutos de ejecución (si algo tarda más, guarden el resultado en
  un CSV y expliquen cómo se generó, como en la sesión 9).

**3. La presentación:** 15 minutos y 5 de preguntas. Todos los integrantes
presentan y todos deben poder contestar preguntas sobre cualquier parte del
proyecto.

## Rúbrica de la entrega final

| Criterio | Puntos | Qué se espera para la nota completa |
|---|---:|---|
| Réplica | 30 | La figura o tabla central está reproducida y las diferencias están explicadas |
| Extensión | 25 | Una pregunta clara, un experimento que la contesta (con varias semillas o con modelo nulo) y un resultado bien interpretado |
| Interpretación y discusión | 15 | Se distingue lo que muestran los resultados de lo que no; se discuten las limitaciones |
| Reproducibilidad | 10 | El cuaderno corre sin errores y reproduce todas las figuras del escrito |
| Escrito | 10 | Claro, con figuras legibles (ejes, unidades, leyendas) y referencias completas |
| Presentación | 10 | Se entiende sin haber leído el artículo, respeta el tiempo y contesta bien las preguntas |

La ruta avanzada del curso no es obligatoria en el proyecto. Si el equipo
incluye una derivación analítica y la compara con la simulación, cuenta
dentro de la réplica o de la extensión, no como puntos extra.

## Ideas de proyecto

Cada idea dice el artículo base, las sesiones del curso que le sirven de
punto de partida y una posible extensión. También se puede proponer otro
artículo, con el visto bueno del profesor antes de la etapa 1.

| Tema | Artículo base | Sesiones | Una posible extensión |
|---|---|---|---|
| Difusión de la cultura | Axelrod (1997) | 15, 16 | ¿Cambia la polarización final en una red de mundo pequeño o con ruido? |
| Restaurante de Kolkata | Chakrabarti et al. (2009) | 19 | Comparar estrategias de aprendizaje con la del juego de la minoría |
| Libro de órdenes | Maslov (2000) | 1, 2, 5 | ¿Qué hechos estilizados reproduce? ¿Y si los agentes tienen memoria? |
| Procesos de Hawkes | Filimonov y Sornette (2012) | 2, 5 | Medir la reflexividad en los días extremos del IPC o del USD/MXN |
| DebtRank | Battiston et al. (2012) | 11, 12 | La red del mercado de la sesión 11 como red de exposiciones |
| Ingreso de la ENIGH | Drăgulescu y Yakovenko (2001) | 8, 9 | Separar ingreso laboral y transferencias; comparar entidades |
| Zipf de ciudades | Gabaix (1999) | 3, 6 | Crecimiento de los municipios de 2010 a 2020: ¿se cumple Gibrat? |
| Interacciones sociales | Brock y Durlauf (2001) | 13, 14 | Histéresis en la adopción de una tecnología en una red |

Los datos de la ENIGH, de las localidades y de los municipios ya están en
`Data/`.

## Referencias de las ideas

- Axelrod, R. (1997). The dissemination of culture: A model with local
  convergence and global polarization. *Journal of Conflict Resolution*,
  41(2), 203–226. https://doi.org/10.1177/0022002797041002001
- Battiston, S., Puliga, M., Kaushik, R., Tasca, P., & Caldarelli, G. (2012).
  DebtRank: Too central to fail? Financial networks, the FED and systemic
  risk. *Scientific Reports*, 2, 541. https://doi.org/10.1038/srep00541
- Brock, W. A., & Durlauf, S. N. (2001). Discrete choice with social
  interactions. *Review of Economic Studies*, 68(2), 235–260.
  https://doi.org/10.1111/1467-937X.00168
- Chakrabarti, A. S., Chakrabarti, B. K., Chatterjee, A., & Mitra, M. (2009).
  The Kolkata Paise Restaurant problem and resource utilization. *Physica A*,
  388(12), 2420–2426. https://doi.org/10.1016/j.physa.2009.02.039
- Drăgulescu, A., & Yakovenko, V. M. (2001). Exponential and power-law
  probability distributions of wealth and income in the United Kingdom and
  the United States. *Physica A*, 299(1–2), 213–221.
  https://doi.org/10.1016/S0378-4371(01)00298-9
- Filimonov, V., & Sornette, D. (2012). Quantifying reflexivity in financial
  markets: Toward a prediction of flash crashes. *Physical Review E*, 85(5),
  056108. https://doi.org/10.1103/PhysRevE.85.056108
- Gabaix, X. (1999). Zipf's law for cities: An explanation. *Quarterly
  Journal of Economics*, 114(3), 739–767.
  https://doi.org/10.1162/003355399556133
- Maslov, S. (2000). Simple model of a limit order-driven market. *Physica
  A*, 278(3–4), 571–578. https://doi.org/10.1016/S0378-4371(00)00067-4

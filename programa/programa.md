# Econofísica

**Nivel**: Licenciatura de Economía, Física o Matemáticas.
**Requisitos**: Álgebra lineal, probabilidad, cálculo diferencial e integral.
**Duración**: 1 semestre de 15 semanas.
**Lenguaje**: Julia.

**Objetivo**: Que las y los estudiantes comprendan cómo la interacción entre
muchos agentes heterogéneos genera regularidades estadísticas a escala
macroscópica (colas pesadas, desigualdad, crisis, consenso), y que dominen las
herramientas de la física estadística para estudiarlas: procesos estocásticos,
máxima entropía, campo medio, transiciones de fase y redes complejas. Al
finalizar, serán capaces de contrastar un modelo con datos reales y de
construir, simular y analizar modelos de agentes en Julia.

**Hilo conductor**: el curso abre midiendo los hechos estilizados en datos
reales (semana 1) y cierra construyendo modelos de agentes que los reproducen
de forma endógena (semana 14).

> **Abreviaturas de lecturas.** MS: Mantegna & Stanley · SL: Slanina ·
> RC: Requião da Cunha · SC: Sen & Chakrabarti · NW: Newman · JK: Jackson ·
> CC: Chakrabarti et al. (2013) · HM: Hommes · KI: Kirman (2011) ·
> TKH: Thurner, Klimek & Hanel. Referencias completas al final.

---

## Módulo I: Fluctuaciones y hechos empíricos

**Objetivo específico**: Caracterizar estadísticamente las fluctuaciones
económicas y financieras, y entender qué procesos estocásticos las generan.

- **Semana 1: Econofísica y hechos estilizados**
	- Qué es la econofísica: de la interacción local a los patrones globales.
	- Hipótesis de mercado eficiente. Retornos logarítmicos y escalas de tiempo.
	- Hechos estilizados: colas pesadas, ausencia de autocorrelación lineal,
	  agrupamiento de volatilidad, normalidad agregacional.
	- *Laboratorio*: IPC (BMV), USD/MXN, S&P 500 y BTC; ACF de r y |r|; cómo
	  cambia la distribución al agregar en el tiempo.
	- Lecturas: MS 1–2, 5, 7; SL 1 (1.1–1.2); RC 2.4; Cont (2001).
- **Semana 2: Colas pesadas y leyes de potencias**
	- Distribuciones estables y teorema del límite central generalizado.
	- Pareto, lognormal, ley cúbica inversa de los retornos, ley de Zipf.
	- Estimación: CCDF log-log, estimador de Hill, máxima verosimilitud,
	  prueba de Kolmogorov-Smirnov.
	- *Laboratorio*: exponente de cola de los retornos; Zipf de ciudades
	  mexicanas (INEGI).
	- Lecturas: MS 3–4, 8–9; RC 2.2.5–2.2.7; TKH 3; Clauset, Shalizi & Newman
	  (2009); Gabaix (2009).
- **Semana 3: Procesos estocásticos**
	- Caminata aleatoria, cadenas de Markov y ecuación maestra.
	- Langevin, Ornstein-Uhlenbeck, Fokker-Planck. Movimiento browniano
	  geométrico.
	- Black-Scholes como ecuación de difusión (media sesión).
	- GARCH: volatilidad con memoria.
	- *Laboratorio*: Euler-Maruyama vs. solución de Fokker-Planck; GARCH(1,1)
	  reproduce el agrupamiento de volatilidad de la semana 1.
	- Lecturas: SL 2; RC 3, 4.3, 2.3.5; MS 10, 14; Gardiner (consulta).
- **Semana 4: Crecimiento multiplicativo**
	- Ley de Gibrat y lognormal.
	- Proceso de Kesten: un mecanismo genérico para leyes de potencias.
	- Distribución de tamaños de empresas y ciudades.
	- *Laboratorio*: Gibrat puro vs. Gibrat con barrera inferior; exponente
	  resultante vs. parámetros.
	- Lecturas: SL 7 (7.1, 7.3); Gabaix (2009).
- **Semana 5: Correlaciones y matrices aleatorias**
	- Matriz de correlación de N activos. Espectro de Marchenko-Pastur: qué es
	  ruido y qué es señal.
	- Modo de mercado y sectores.
	- Distancia entre activos y árbol de expansión mínima: puente hacia redes.
	- *Laboratorio*: ~100 acciones; espectro vs. Marchenko-Pastur; árbol
	  coloreado por sector.
	- Lecturas: MS 12–13; RC 5.2; SL 1.4; Laloux et al. (1999); Mantegna (1999).
	- **Entrega**: examen para casa del Módulo I.

## Módulo II: Distribución, redes y transiciones de fase

**Objetivo específico**: Usar máxima entropía, campo medio y teoría de redes
para explicar la desigualdad, el contagio y los cambios abruptos colectivos.

- **Semana 6: Máxima entropía y modelos cinéticos de intercambio**
	- Dinero conservado → distribución de Boltzmann-Gibbs.
	- Modelo de Dragulescu-Yakovenko; ahorro uniforme (Chakraborti-Chakrabarti,
	  distribución gamma); ahorro heterogéneo (cola de Pareto).
	- Curva de Lorenz e índice de Gini.
	- *Laboratorio*: simular los tres modelos y compararlos con la ENIGH.
	- Lecturas: CC 2–5; RC 8.3; SL 7.2; Yakovenko & Rosser (2009);
	  Garibaldi & Scalas (consulta).
- **Semana 7: Condensación de riqueza y redistribución**
	- Modelo de Bouchaud-Mézard: intercambio + crecimiento multiplicativo.
	- Solución estacionaria de Fokker-Planck → Pareto. Transición de condensación.
	- Impuestos y redistribución como experimento numérico.
	- *Laboratorio*: transición de condensación vs. intensidad de intercambio;
	  efecto de un impuesto sobre el Gini.
	- Lecturas: Bouchaud & Mézard (2000); CC 4.5; Ribeiro (consulta).
- **Semana 8: Redes: estructura y modelos**
	- Medidas: grado, centralidades, agrupamiento, asortatividad, Laplaciano.
	- Erdős-Rényi, Watts-Strogatz, Barabási-Albert, modelo de configuración.
	- Percolación y robustez: fallas aleatorias vs. ataques dirigidos.
	- *Laboratorio*: comparar modelos con una red real (árbol de la semana 5 o
	  comercio internacional); ataque dirigido vs. aleatorio.
	- Lecturas: NW 6–7, 10–13, 15; JK 2–5; RC 8.1; SL 6.
- **Semana 9: Contagio y riesgo sistémico**
	- Cascadas de umbral (Watts).
	- Redes interbancarias: compensación de Eisenberg-Noe y DebtRank.
	- "Robusto pero frágil": diversificación vs. contagio.
	- *Laboratorio*: cascadas en ER vs. BA; red interbancaria sintética.
	- Lecturas: Watts (2002); Eisenberg & Noe (2001); Battiston et al. (2012);
	  Acemoglu, Ozdaglar & Tahbaz-Salehi (2015); Glasserman & Young (2016);
	  Battiston et al. (2016).
- **Semana 10: Ising y elección discreta**
	- Modelo de Ising: campo medio, transición de fase, histéresis.
	- Logit ↔ distribución de Boltzmann. Elección discreta con interacciones
	  sociales (Brock-Durlauf).
	- Campo aleatorio: avalanchas y crisis.
	- *Laboratorio*: Metropolis; campo medio vs. simulación; ciclo de histéresis.
	- Lecturas: SC 2–3; RC 6; Brock & Durlauf (2001); Bouchaud, Marsili &
	  Nadal (2023).

## Módulo III: Agentes, estrategias y mercados artificiales

**Objetivo específico**: Construir modelos de agentes con racionalidad acotada
y aprendizaje, y usarlos para reproducir patrones empíricos.

- **Semana 11: Opinión, herding y segregación**
	- Modelo de votante: tiempo de consenso, efecto de la red.
	- Modelo de hormigas de Kirman: herding y bimodalidad.
	- Schelling: preferencias leves → segregación fuerte.
	- *Laboratorio*: votante en redes; Kirman; Schelling completo con Agents.jl.
	- Lecturas: SC 3–4, 6; SL 3, 8.2; RC 8.2; KI 2–3, 6; Castellano, Fortunato
	  & Loreto (2009); Kirman (1993).
- **Semana 12: Juegos evolutivos y cooperación**
	- Dinámica del replicador como sistema dinámico.
	- Lotka-Volterra ↔ ciclo de Goodwin.
	- Dilema del prisionero en redes: reciprocidad de red, b/c > k.
	- *Laboratorio*: Halcón-Paloma; Goodwin; cooperación en retícula vs. red BA.
	- Lecturas: RC 7; SL 8.1; Nowak & May (1992); Goodwin (1967);
	  Strogatz (consulta); Schecter & Gintis (consulta).
- **Semana 13: Juego de la minoría**
	- El Farol. Reglas del juego de la minoría.
	- Volatilidad σ²/N vs. α = 2^M/N: fases eficiente e ineficiente.
	- El juego de la minoría como modelo de mercado.
	- *Laboratorio*: curva σ²/N vs. α para varios N.
	- Lecturas: SL 5; RC 7.1.1; Challet & Zhang (1997); Challet, Marsili &
	  Zhang (2005, consulta).
- **Semana 14: Mercados artificiales y expectativas heterogéneas**
	- Brock-Hommes: selección evolutiva de reglas de predicción, bifurcaciones
	  y caos.
	- Fundamentalistas vs. chartistas (Lux-Marchesi, Giardina-Bouchaud).
	- **Cierre del ciclo**: aplicar el análisis de la semana 1 a series
	  simuladas. ¿Qué hechos estilizados reproduce cada modelo?
	- *Laboratorio*: diagrama de bifurcación de Brock-Hommes; hechos
	  estilizados en series simuladas.
	- Lecturas: HM 1, 6; SL 3.2; LeBaron (2006); Lux & Marchesi (1999);
	  Giardina & Bouchaud (2003).
- **Semana 15: Proyecto final**
	- Presentaciones.

---

## Evaluación

| Componente | Peso |
|---|---|
| Laboratorios. Cada uno cierra con una pregunta cuantitativa (estimar un exponente, comparar simulación con campo medio, etc.) | 30% |
| Presentación de un artículo de investigación (15 min, desde la semana 3) | 15% |
| Examen para casa del Módulo I | 15% |
| Proyecto final: réplica y extensión de un artículo (escrito + presentación) | 40% |

**Ideas para proyecto final**: modelo de Axelrod; problema del restaurante de
Kolkata; modelos de libro de órdenes (Maslov); procesos de Hawkes; DebtRank con
datos reales; ajuste de la distribución de ingreso de la ENIGH; Zipf de
ciudades y empresas mexicanas; Brock-Durlauf con datos de adopción
tecnológica.

---

## Bibliografía

### Textos base
- Mantegna, R. N., & Stanley, H. E. (1999). *An Introduction to Econophysics:
  Correlations and Complexity in Finance*. Cambridge University Press.
  https://doi.org/10.1017/CBO9780511755767
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford University Press.
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Sen, P., & Chakrabarti, B. K. (2014). *Sociophysics: An Introduction*.
  Oxford University Press.
- Newman, M. (2018). *Networks* (2.ª ed.). Oxford University Press.
  https://doi.org/10.1093/oso/9780198805090.001.0001

### Textos complementarios
- Chakrabarti, B. K., Chakraborti, A., Chakravarty, S. R., & Chatterjee, A.
  (2013). *Econophysics of Income and Wealth Distributions*. Cambridge
  University Press.
- Hommes, C. (2013). *Behavioral Rationality and Heterogeneous Expectations in
  Complex Economic Systems*. Cambridge University Press.
- Jackson, M. O. (2008). *Social and Economic Networks*. Princeton University Press.
- Kirman, A. (2011). *Complex Economics: Individual and Collective
  Rationality*. Routledge.
- Thurner, S., Klimek, P., & Hanel, R. (2018). *Introduction to the Theory of
  Complex Systems*. Oxford University Press.
- Gardiner, C. (2021). *Elements of Stochastic Methods*. AIP Publishing.
- Garibaldi, U., & Scalas, E. (2010). *Finitary Probabilistic Methods in
  Econophysics*. Cambridge University Press.
- Ribeiro, M. B. (2020). *Income Distribution Dynamics of Economic Systems*.
  Cambridge University Press.
- Strogatz, S. (2019). *Nonlinear Dynamics and Chaos* (2.ª ed.). CRC Press.
- Schecter, S., & Gintis, H. (2016). *Game Theory in Action*. Princeton
  University Press.
- Challet, D., Marsili, M., & Zhang, Y.-C. (2005). *Minority Games*. Oxford
  University Press.

### Artículos
- Acemoglu, D., Ozdaglar, A., & Tahbaz-Salehi, A. (2015). Systemic risk and
  stability in financial networks. *American Economic Review*, 105(2), 564–608.
  https://doi.org/10.1257/aer.20130456
- Battiston, S., Puliga, M., Kaushik, R., Tasca, P., & Caldarelli, G. (2012).
  DebtRank: Too central to fail? Financial networks, the FED and systemic risk.
  *Scientific Reports*, 2, 541. https://doi.org/10.1038/srep00541
- Battiston, S., Farmer, J. D., Flache, A., Garlaschelli, D., Haldane, A. G.,
  Heesterbeek, H., Hommes, C., Jaeger, C., May, R., & Scheffer, M. (2016).
  Complexity theory and financial regulation. *Science*, 351(6275), 818–819.
  https://doi.org/10.1126/science.aad0299
- Bouchaud, J.-P., & Mézard, M. (2000). Wealth condensation in a simple model
  of economy. *Physica A*, 282(3–4), 536–545.
  https://doi.org/10.1016/S0378-4371(00)00205-3
- Bouchaud, J.-P., Marsili, M., & Nadal, J.-P. (2023). Application of spin
  glass ideas in social sciences, economics and finance. En *Spin Glass Theory
  and Far Beyond* (pp. 561–579). World Scientific.
  https://doi.org/10.1142/9789811273926_0027
- Brock, W. A., & Durlauf, S. N. (2001). Discrete choice with social
  interactions. *Review of Economic Studies*, 68(2), 235–260.
  https://doi.org/10.1111/1467-937X.00168
- Castellano, C., Fortunato, S., & Loreto, V. (2009). Statistical physics of
  social dynamics. *Reviews of Modern Physics*, 81(2), 591–646.
  https://doi.org/10.1103/RevModPhys.81.591
- Challet, D., & Zhang, Y.-C. (1997). Emergence of cooperation and organization
  in an evolutionary game. *Physica A*, 246(3–4), 407–418.
  https://doi.org/10.1016/S0378-4371(97)00419-6
- Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
  distributions in empirical data. *SIAM Review*, 51(4), 661–703.
  https://doi.org/10.1137/070710111
- Cont, R. (2001). Empirical properties of asset returns: stylized facts and
  statistical issues. *Quantitative Finance*, 1(2), 223–236.
  https://doi.org/10.1080/713665670
- Eisenberg, L., & Noe, T. H. (2001). Systemic risk in financial systems.
  *Management Science*, 47(2), 236–249.
  https://doi.org/10.1287/mnsc.47.2.236.9835
- Gabaix, X. (2009). Power laws in economics and finance. *Annual Review of
  Economics*, 1, 255–294.
  https://doi.org/10.1146/annurev.economics.050708.142940
- Giardina, I., & Bouchaud, J.-P. (2003). Bubbles, crashes and intermittency in
  agent based market models. *European Physical Journal B*, 31(3), 421–437.
  https://doi.org/10.1140/epjb/e2003-00050-6
- Glasserman, P., & Young, H. P. (2016). Contagion in financial networks.
  *Journal of Economic Literature*, 54(3), 779–831.
  https://doi.org/10.1257/jel.20151228
- Goodwin, R. M. (1967). A growth cycle. En C. H. Feinstein (Ed.),
  *Socialism, Capitalism and Economic Growth: Essays Presented to Maurice Dobb*
  (pp. 54–58). Cambridge University Press.
- Kirman, A. (1993). Ants, rationality, and recruitment. *Quarterly Journal of
  Economics*, 108(1), 137–156. https://doi.org/10.2307/2118498
- Laloux, L., Cizeau, P., Bouchaud, J.-P., & Potters, M. (1999). Noise
  dressing of financial correlation matrices. *Physical Review Letters*,
  83(7), 1467–1470. https://doi.org/10.1103/PhysRevLett.83.1467
- LeBaron, B. (2006). Agent-based financial markets: Matching stylized facts
  with style. En D. Colander (Ed.), *Post Walrasian Macroeconomics*
  (pp. 221–236). Cambridge University Press.
  https://doi.org/10.1017/CBO9780511617751.013
- Lux, T., & Marchesi, M. (1999). Scaling and criticality in a stochastic
  multi-agent model of a financial market. *Nature*, 397(6719), 498–500.
  https://doi.org/10.1038/17290
- Mantegna, R. N. (1999). Hierarchical structure in financial markets.
  *European Physical Journal B*, 11(1), 193–197.
  https://doi.org/10.1007/s100510050929
- Nowak, M. A., & May, R. M. (1992). Evolutionary games and spatial chaos.
  *Nature*, 359(6398), 826–829. https://doi.org/10.1038/359826a0
- Watts, D. J. (2002). A simple model of global cascades on random networks.
  *PNAS*, 99(9), 5766–5771. https://doi.org/10.1073/pnas.082090499
- Yakovenko, V. M., & Rosser, J. B. (2009). Colloquium: Statistical mechanics
  of money, wealth, and income. *Reviews of Modern Physics*, 81(4), 1703–1725.
  https://doi.org/10.1103/RevModPhys.81.1703

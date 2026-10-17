# Econofísica

**Nivel**: Licenciatura de Economía, Física o Matemáticas.
**Requisitos**: Probabilidad básica y cálculo diferencial. Lo demás se repasa
en clase cuando se necesita.
**Duración**: 20 sesiones de 2 horas (dos por semana). El mismo núcleo sirve
para un semestre de 15 semanas o un trimestre de 11 (ver Calendarios).
**Lenguaje**: Julia.

**Objetivo**: Que las y los estudiantes comprendan cómo la interacción entre
muchos agentes heterogéneos genera regularidades estadísticas a escala
macroscópica (colas pesadas, desigualdad, crisis, consenso), y que dominen las
herramientas de la física estadística para estudiarlas: procesos estocásticos,
máxima entropía, campo medio, transiciones de fase y redes complejas. Al
finalizar, serán capaces de contrastar un modelo con datos reales y de
construir, simular y analizar modelos de agentes en Julia.

**Hilo conductor**: el curso abre midiendo los hechos estilizados en datos
reales (sesiones 1 y 2) y cierra construyendo modelos de agentes que los
reproducen de forma endógena (sesión 20 y taller del Módulo III).

## Metodología

Cada sesión sigue el orden **fenómeno → reglas del modelo en palabras →
simular y medir → formalizar**. La matemática llega como explicación de algo
que ya se observó numéricamente, no como requisito para empezar.

**Estructura de cada sesión (2 horas)**:

| Momento | Tiempo | Qué ocurre |
|---|---|---|
| Fenómeno | 15 min | Un dato real o una paradoja que motiva la pregunta de la sesión |
| Modelo | 35 min | Reglas en palabras, la ecuación central y a lo más una derivación corta |
| Laboratorio | 50 min | En parejas, con el cuaderno de Pluto: completar código, simular, medir |
| Cierre | 20 min | Pregunta de medición o interpretación; es lo que se evalúa |

**Dos rutas.** Todo el contenido sin marca es la **ruta básica**: obligatoria,
basada en simulación, medición e interpretación económica, y es lo único que se
ve en clase. Los puntos marcados **[Avanzado]** forman la **ruta avanzada**:
derivaciones analíticas y su comparación con la simulación. Se estudia fuera
de clase en la sección "Para profundizar" de las notas, es opcional y da
puntos extra (ver Evaluación).

- **Repaso**: al inicio de la sesión que lo requiere, unos minutos de la
  herramienta matemática necesaria, con su clave del repaso por niveles
  (`Class_notes/repaso/`).
- **Notebooks con huecos**: al principio se completan 2 o 3 líneas clave; los
  huecos crecen conforme avanza el curso.
- **Simulaciones interactivas** (Pluto.jl): deslizadores para explorar
  transiciones de fase.
- **Equipos mixtos**: estudiantes de Economía junto con estudiantes de Física
  o Matemáticas.
- **Lecturas**: una lectura base por sesión, corta y en inglés; las
  complementarias son opcionales.

### Conceptos transversales

Estas ideas sostienen todo el curso. Cada una se presenta con calma una vez y
se retoma después con preguntas en los laboratorios.

| Concepto | Se presenta | Se retoma |
|---|---|---|
| Retornos logarítmicos | Sesión 1 | 2, 5, 20 |
| Distribuciones, CCDF y escalas logarítmicas | Sesiones 1 y 3 | 6, 8, 9, 10 |
| Media, varianza y correlación | Sesiones 1 y 2 | 7, 19 |
| Caminata aleatoria | Sesión 4 | 5, 15 |
| Cadenas de Markov | Sesión 4 | 13, 15 |
| Promediar simulaciones sobre varias corridas | Sesión 4 | 9, 15, 19 |
| Logit (elección discreta) | Sesión 13 | 18, 20 |
| Puntos fijos y su estabilidad | Sesiones 12 y 14 | 17, 20 |

> **Abreviaturas de lecturas.** MS: Mantegna & Stanley · SL: Slanina ·
> RC: Requião da Cunha · SC: Sen & Chakrabarti · NW: Newman · JK: Jackson ·
> CC: Chakrabarti et al. (2013) · HM: Hommes · KI: Kirman (2011) ·
> TKH: Thurner, Klimek & Hanel · WR: Wilensky & Rand.
> Las lecturas marcadas [Avanzado] corresponden a la ruta avanzada.
> Referencias completas al final.

---

## Módulo I: Fluctuaciones y hechos empíricos

**Objetivo específico**: Caracterizar estadísticamente las fluctuaciones
económicas y financieras, y entender qué procesos estocásticos las generan.

- **Sesión 1: Retornos y colas pesadas**
	- *Fenómeno*: ¿el mercado es un volado? Días de caída que "no deberían"
	  ocurrir.
	- *Diagnóstico* de matemáticas (10 preguntas, no cuenta para la calificación).
	- *Repaso*: logaritmos, media, varianza, histogramas.
	- Qué es la econofísica: de la interacción local a los patrones globales.
	- Hipótesis de mercado eficiente. Retornos logarítmicos y escalas de tiempo.
	- Hechos estilizados 1 y 2: colas pesadas y ausencia de autocorrelación
	  lineal.
	- *Laboratorio*: IPC (BMV), USD/MXN, S&P 500 y BTC; histograma de r frente
	  a la normal; ACF de r.
	- Lectura: MS 1–2. Complementaria: SL 1 (1.1–1.2).
- **Sesión 2: Volatilidad y modelos nulos**
	- *Fenómeno*: calma y tormenta; los días malos llegan juntos.
	- Hechos estilizados 3 y 4: agrupamiento de volatilidad y normalidad
	  agregacional.
	- Modelos nulos: ¿qué tan especial es lo que medimos? Barajar los datos.
	- Universalidad de los hechos estilizados.
	- *Laboratorio*: ACF de |r|; cómo cambia la distribución al agregar en el
	  tiempo; serie original frente a serie barajada.
	- Lectura: Cont (2001). Complementaria: MS 5, 7; RC 2.4.
- **Sesión 3: Colas pesadas y leyes de potencias**
	- *Fenómeno*: un mundo sin tamaño típico; las localidades de México.
	- *Repaso*: por qué una ley de potencias es una recta en escala log-log.
	- Pareto, lognormal, ley cúbica inversa de los retornos, ley de Zipf.
	- Estimación: CCDF log-log y estimador de Hill (como receta).
	- **[Avanzado]** Distribuciones estables y teorema del límite central
	  generalizado. Máxima verosimilitud y método de Clauset et al.
	- *Laboratorio*: exponente de cola de los retornos; Zipf de localidades
	  mexicanas (INEGI).
	- Lectura: RC 2.2.7. Complementaria: Gabaix (2009); TKH 3.
	  **[Avanzado]** MS 3–4; RC 2.2.5–2.2.6; Clauset, Shalizi & Newman (2009).
- **Sesión 4: Caminata aleatoria y cadenas de Markov**
	- *Fenómeno*: un mercado que alterna entre calma y crisis.
	- *Repaso*: esperanza y varianza de una suma.
	- Caminata aleatoria: por qué la dispersión crece como √t.
	- Cadenas de Markov: matriz de transición y distribución estacionaria.
	- Promediar simulaciones: una corrida no basta.
	- *Laboratorio*: caminatas simuladas frente a √t; un mercado con dos
	  estados y la duración de sus crisis.
	- Lectura: SL 2.1. Complementaria: RC 3.1.
- **Sesión 5: Langevin, GBM y GARCH**
	- *Fenómeno*: ¿qué proceso genera estas series? Precio que no puede ser
	  negativo y volatilidad con memoria.
	- *Repaso*: la derivada como tasa de cambio.
	- Ecuación de Langevin como "tendencia + ruido"; movimiento browniano
	  geométrico; simulación con Euler-Maruyama.
	- GARCH: volatilidad con memoria.
	- **[Avanzado]** Ecuación maestra y Fokker-Planck. Black-Scholes como
	  ecuación de difusión.
	- *Laboratorio*: simular GBM y GARCH(1,1); ¿reproduce GARCH el
	  agrupamiento de volatilidad de la sesión 2?
	- Lectura: RC 2.3.5. Complementaria: SL 2.2.
	  **[Avanzado]** MS 10; RC 3.3, 4.3; Gardiner (consulta).
- **Sesión 6: Crecimiento multiplicativo**
	- *Fenómeno*: ¿crecen igual las empresas grandes y las pequeñas?
	- *Repaso*: el logaritmo de un producto es una suma de logaritmos
	  (de ahí la lognormal).
	- Ley de Gibrat y lognormal.
	- Proceso de Kesten: un mecanismo genérico para leyes de potencias.
	- Distribución de tamaños de empresas y ciudades.
	- **[Avanzado]** Derivación del exponente de Kesten.
	- *Laboratorio*: Gibrat puro frente a Gibrat con barrera inferior; cómo
	  cambia el exponente con los parámetros.
	- Lectura: SL 7.1. Complementaria: SL 7.3.
- **Sesión 7: Correlaciones y matrices aleatorias**
	- *Fenómeno*: las acciones se mueven juntas.
	- *Repaso*: matrices; eigenvalores como "direcciones principales".
	- Matriz de correlación de N activos. Barajar los datos destruye las
	  correlaciones: lo que sobresale del espectro barajado es señal.
	- Modo de mercado y sectores.
	- Distancia entre activos y árbol de expansión mínima: puente hacia redes.
	- **[Avanzado]** Ley de Marchenko-Pastur.
	- *Laboratorio*: ~100 acciones; espectro real frente a barajado; árbol
	  coloreado por sector.
	- Lectura: MS 12–13. Complementaria: SL 1.4; Mantegna (1999).
	  **[Avanzado]** RC 5.2; Laloux et al. (1999).
	- **Entrega**: taller para casa del Módulo I.

## Módulo II: Distribución, redes y transiciones de fase

**Objetivo específico**: Usar máxima entropía, campo medio y teoría de redes
para explicar la desigualdad, el contagio y los cambios abruptos colectivos.

- **Sesión 8: Desigualdad y modelos de intercambio**
	- *Fenómeno*: ¿cómo se reparte el ingreso en México?
	- Curva de Lorenz e índice de Gini.
	- Dinero conservado → distribución exponencial (argumento verbal y numérico).
	- Modelo de Dragulescu-Yakovenko; ahorro uniforme (Chakraborti-Chakrabarti);
	  ahorro heterogéneo (cola de Pareto).
	- **[Avanzado]** Máxima entropía con multiplicadores de Lagrange →
	  Boltzmann-Gibbs; aproximación gamma del modelo con ahorro.
	- *Laboratorio*: simular los tres modelos y compararlos con la ENIGH.
	- Lectura: Yakovenko & Rosser (2009). Complementaria: RC 8.3; SL 7.2.
	  **[Avanzado]** CC 2–5; Garibaldi & Scalas (consulta).
- **Sesión 9: Condensación de riqueza y redistribución**
	- *Fenómeno*: inversión sin intercambio; por qué unos pocos acumulan casi
	  todo.
	- Modelo de Bouchaud-Mézard: intercambio + crecimiento multiplicativo.
	- Transición de condensación explorada con un deslizador.
	- Impuestos y redistribución como experimento numérico.
	- **[Avanzado]** Solución estacionaria de Fokker-Planck → exponente de Pareto.
	- *Laboratorio*: transición de condensación frente a la intensidad de
	  intercambio; efecto de un impuesto sobre el Gini.
	- Lectura: CC 4.5. **[Avanzado]** Bouchaud & Mézard (2000);
	  Ribeiro (consulta).
- **Sesión 10: Redes: medidas y modelos**
	- *Fenómeno*: ¿por qué ganaron los Médici? Los matrimonios de las familias
	  de Florencia.
	- *Repaso*: matriz de adyacencia.
	- Medidas: grado, centralidades, agrupamiento.
	- Erdős-Rényi, Watts-Strogatz, Barabási-Albert.
	- *Laboratorio*: medir la red de Florencia; comparar los tres modelos.
	- Lectura: NW 7. Complementaria: NW 6, 11; RC 8.1.
- **Sesión 11: Redes: el mercado y la robustez**
	- *Fenómeno*: ¿qué pasa si cae el nodo más conectado?
	- La red de correlaciones del S&P 100 (datos de la sesión 7) frente a los
	  modelos.
	- Robustez: fallas aleatorias frente a ataques dirigidos (por simulación).
	- **[Avanzado]** Laplaciano y espectro; asortatividad; modelo de
	  configuración; umbral de percolación.
	- *Laboratorio*: la red del mercado; ataque dirigido frente a fallas
	  aleatorias.
	- Lectura: SL 6. **[Avanzado]** NW 12–13, 15; JK 2–5 (consulta).
- **Sesión 12: Contagio y riesgo sistémico**
	- *Fenómeno*: de una quiebra a una crisis.
	- *Repaso*: un choque que se repite atenuándose.
	- Cascadas de umbral (Watts).
	- Redes interbancarias: contagio de pérdidas (DebtRank simplificado).
	- "Robusto pero frágil": diversificación frente a contagio.
	- **[Avanzado]** Compensación de Eisenberg-Noe como problema de punto fijo.
	- *Laboratorio*: cascadas en ER frente a BA; red interbancaria sintética.
	- Lectura: Battiston et al. (2016). Complementaria: Watts (2002).
	  **[Avanzado]** Eisenberg & Noe (2001); Battiston et al. (2012);
	  Acemoglu, Ozdaglar & Tahbaz-Salehi (2015) y Glasserman & Young (2016)
	  (consulta).
- **Sesión 13: Logit, Boltzmann e Ising**
	- *Fenómeno*: mismos incentivos, resultados distintos.
	- *Repaso*: función logística.
	- Logit ↔ distribución de Boltzmann: el modelo de elección discreta que
	  usan los economistas es física estadística.
	- Elección con interacciones sociales: el modelo de Ising.
	- *Laboratorio*: el algoritmo de Metropolis; transición de fase con un
	  deslizador.
	- Lectura: SC 2. **[Avanzado]** Brock & Durlauf (2001).
- **Sesión 14: Campo medio, histéresis y crisis**
	- *Fenómeno*: cambios abruptos que no se revierten al deshacer la causa.
	- Campo medio resuelto gráficamente: intersección de m y tanh(βJm).
	- Histéresis: el estado depende de la historia.
	- **[Avanzado]** Derivación del campo medio; Brock-Durlauf formal; campo
	  aleatorio, avalanchas y crisis.
	- *Laboratorio*: simulación frente a solución gráfica; ciclo de histéresis.
	- Lectura: RC 6. **[Avanzado]** Bouchaud, Marsili & Nadal (2023).
	- **Entrega**: taller para casa del Módulo II.

## Módulo III: Agentes, estrategias y mercados artificiales

**Objetivo específico**: Construir modelos de agentes con racionalidad acotada
y aprendizaje, y usarlos para reproducir patrones empíricos.

- **Sesión 15: Votante y Kirman: imitación y herding**
	- *Fenómeno*: dos restaurantes idénticos; uno se llena y el otro no.
	- *Repaso*: estados absorbentes frente a ruido.
	- Modelo de votante: consenso, probabilidad de victoria, efecto de la red.
	- Modelo de hormigas de Kirman: herding y bimodalidad.
	- **[Avanzado]** Tiempo de consenso del votante; ecuación maestra del
	  modelo de Kirman.
	- *Laboratorio*: votante en la red completa y en redes; Kirman con un
	  deslizador.
	- Lectura: SC 3. Complementaria: KI 2–3.
	  **[Avanzado]** SL 3; Castellano, Fortunato & Loreto (2009); Kirman (1993).
- **Sesión 16: Schelling: segregación**
	- *Fenómeno*: barrios segregados sin intolerancia.
	- Schelling: preferencias leves → segregación fuerte.
	- Analogía con el modelo de Ising y sus límites.
	- *Laboratorio*: Schelling completo (implementación propia en Julia, sin
	  paquetes de agentes); barrido del umbral de tolerancia.
	- Lectura: KI 6. Complementaria: SC 4; WR (consulta).
	  **[Avanzado]** SL 8.2; SC 6.
- **Sesión 17: Replicador y ciclo de Goodwin**
	- *Fenómeno*: ¿por qué no todos son agresivos? ¿Y por qué la economía
	  oscila sin choques?
	- *Repaso*: una ecuación diferencial como regla de actualización.
	- Dinámica del replicador; el juego Halcón-Paloma.
	- Lotka-Volterra ↔ ciclo de Goodwin.
	- **[Avanzado]** Estabilidad de puntos fijos con el jacobiano; cantidad
	  conservada de Goodwin.
	- *Laboratorio*: Halcón-Paloma por EDO y por agentes; el ciclo de Goodwin.
	- Lectura: RC 7. **[Avanzado]** Goodwin (1967); Strogatz (consulta).
- **Sesión 18: Cooperación en redes**
	- *Fenómeno*: cooperación entre egoístas cuando traicionar paga más.
	- Dilema del prisionero; por qué la cooperación desaparece en una población
	  bien mezclada.
	- Reciprocidad espacial en la retícula y papel de los concentradores en
	  redes de Barabási-Albert.
	- **[Avanzado]** Regla b/c > k de reciprocidad de red (enunciado y
	  simulación; la demostración, como lectura).
	- *Laboratorio*: cooperación en retícula frente a red BA.
	- Lectura: Nowak & May (1992). Complementaria: SL 8.1.
	  **[Avanzado]** Santos & Pacheco (2005); Ohtsuki et al. (2006);
	  Schecter & Gintis (consulta).
- **Sesión 19: Juego de la minoría**
	- *Fenómeno*: el bar El Farol; la trampa de pensar igual.
	- Reglas del juego de la minoría.
	- Volatilidad σ²/N frente a α = 2^M/N: fases hacinada y diluida.
	- El juego de la minoría como modelo de mercado.
	- **[Avanzado]** Predictibilidad y agentes congelados en α_c (intuición; la
	  teoría completa, como lectura).
	- *Laboratorio*: curva σ²/N frente a α para varios N.
	- Lectura: SL 5.1. Complementaria: Challet & Zhang (1997); RC 7.1.1;
	  SL 5.3. **[Avanzado]** SL 5.2; Challet, Marsili & Zhang (2005, consulta).
- **Sesión 20: Mercados artificiales y expectativas heterogéneas**
	- *Fenómeno*: calma y tormenta generadas desde dentro del mercado.
	- Brock-Hommes: selección evolutiva de reglas de predicción.
	- Fundamentalistas frente a chartistas (Lux-Marchesi, Giardina-Bouchaud);
	  en el laboratorio, la versión estimada de Franke & Westerhoff (2012).
	- **[Avanzado]** Estabilidad del equilibrio fundamental y diagrama de
	  bifurcación de Brock-Hommes.
	- *Laboratorio*: burbujas en Brock-Hommes; volatilidad en
	  Franke-Westerhoff.
	- Lectura: LeBaron (2006). Complementaria: HM 1; Lux & Marchesi (1999).
	  **[Avanzado]** HM 2, 6; SL 3.2; Giardina & Bouchaud (2003); Brock &
	  Hommes (1998).
	- **Entrega**: taller para casa del Módulo III. **Cierre del ciclo**:
	  aplicar el análisis de las sesiones 1 y 2 a series simuladas. ¿Qué
	  hechos estilizados reproduce cada modelo?

---

## Calendarios

El contenido es el mismo; cambia el tiempo para practicar.

**Semestre (15 semanas, 30 sesiones).** Las 10 sesiones adicionales no
agregan temas: son talleres de repaso, asesorías del proyecto y
presentaciones.

| Semana | Sesión A | Sesión B |
|---|---|---|
| 1 | 1 | 2 |
| 2 | 3 | 4 |
| 3 | 5 | 6 |
| 4 | 7 | Taller del Módulo I en clase |
| 5 | 8 | 9 |
| 6 | 10 | 11 |
| 7 | 12 | 13 |
| 8 | 14 | Taller del Módulo II en clase |
| 9 | 15 | 16 |
| 10 | 17 | 18 |
| 11 | 19 | 20 |
| 12 | Taller del Módulo III en clase | Asesoría del proyecto final |
| 13 | Presentaciones de artículos | Asesoría del proyecto final |
| 14 | Presentaciones finales | Presentaciones finales |
| 15 | Presentaciones finales | Margen para imprevistos |

**Trimestre (11 semanas, 22 sesiones).** Las 20 sesiones de contenido van
seguidas, dos por semana, y la semana 11 queda para las presentaciones
finales. Los talleres se resuelven en casa. Si se pierde una sesión, se
compactan la 10 con la 11 o la 15 con la 16.

---

## Evaluación

| Componente | Peso |
|---|---|
| Laboratorios, ruta básica. Cada uno cierra con una pregunta de medición o interpretación | 25% |
| Talleres para casa, uno por módulo (10% cada uno) | 30% |
| Presentación de un artículo de investigación (en el trimestre, reseña escrita de dos páginas) | 10% |
| Proyecto final: réplica y extensión de un artículo (escrito + presentación) | 35% |
| **Extra**: ejercicios de la ruta avanzada | hasta +10% |

**Los talleres y los laboratorios evalúan interpretación, no derivaciones.**
Ejemplos de preguntas:
- "Este histograma, ¿es exponencial o ley de potencias? ¿Cómo lo sabes?"
- "Si aumenta el ahorro de los agentes, ¿qué le pasa al Gini? Justifica con tu
  simulación."
- "¿Qué hecho estilizado reproduce este modelo y cuál no?"

**Talleres para casa**:
- **Módulo I** (después de la sesión 7): hechos estilizados y colas de un
  activo que no se vio en clase.
- **Módulo II** (después de la sesión 14): por definir.
- **Módulo III** (después de la sesión 20): cierre del ciclo, con los modelos
  de agentes frente a los datos.

Cada taller incluye una sección opcional de ruta avanzada que cuenta como
extra.

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
- Wilensky, U., & Rand, W. (2015). *An Introduction to Agent-Based Modeling:
  Modeling Natural, Social, and Engineered Complex Systems with NetLogo*.
  MIT Press.
- Railsback, S. F., & Grimm, V. (2019). *Agent-Based and Individual-Based
  Modeling: A Practical Introduction* (2.ª ed.). Princeton University Press.
- Strogatz, S. (2019). *Nonlinear Dynamics and Chaos* (2.ª ed.). CRC Press.
- Schecter, S., & Gintis, H. (2016). *Game Theory in Action*. Princeton
  University Press.
- Challet, D., Marsili, M., & Zhang, Y.-C. (2005). *Minority Games*. Oxford
  University Press.

### Artículos
- Acemoglu, D., Ozdaglar, A., & Tahbaz-Salehi, A. (2015). Systemic risk and
  stability in financial networks. *American Economic Review*, 105(2), 564–608.
  https://doi.org/10.1257/aer.20130456
- Arthur, W. B. (1994). Inductive reasoning and bounded rationality.
  *American Economic Review*, 84(2), 406–411. https://www.jstor.org/stable/2117868
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
- Brock, W. A., & Hommes, C. H. (1998). Heterogeneous beliefs and routes to
  chaos in a simple asset pricing model. *Journal of Economic Dynamics and
  Control*, 22(8–9), 1235–1274. https://doi.org/10.1016/S0165-1889(98)00011-6
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
- Franke, R., & Westerhoff, F. (2012). Structural stochastic volatility in asset
  pricing dynamics: Estimation and model contest. *Journal of Economic Dynamics
  and Control*, 36(8), 1193–1211. https://doi.org/10.1016/j.jedc.2011.10.004
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
- Ohtsuki, H., Hauert, C., Lieberman, E., & Nowak, M. A. (2006). A simple rule
  for the evolution of cooperation on graphs and social networks. *Nature*,
  441(7092), 502–505. https://doi.org/10.1038/nature04605
- Santos, F. C., & Pacheco, J. M. (2005). Scale-free networks provide a
  unifying framework for the emergence of cooperation. *Physical Review
  Letters*, 95(9), 098104. https://doi.org/10.1103/PhysRevLett.95.098104
- Watts, D. J. (2002). A simple model of global cascades on random networks.
  *PNAS*, 99(9), 5766–5771. https://doi.org/10.1073/pnas.082090499
- Yakovenko, V. M., & Rosser, J. B. (2009). Colloquium: Statistical mechanics
  of money, wealth, and income. *Reviews of Modern Physics*, 81(4), 1703–1725.
  https://doi.org/10.1103/RevModPhys.81.1703

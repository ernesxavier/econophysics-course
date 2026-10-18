# Diagnóstico de matemáticas (sesión 1)

Diez preguntas cortas para saber qué conviene repasar antes de empezar. **No
cuenta para la calificación**: sirve para que tú (y el profesor) sepan qué
claves de repaso necesitas. Resuélvelo sin calculadora en unos 25 minutos;
donde haga falta un número, basta una aproximación. Las respuestas están
plegadas debajo de cada pregunta.

Todas las preguntas usan sólo la ruta básica: las claves B de las [notas de
repaso](README.md).

---

**1. Retornos y logaritmos.** Un precio pasa de 100 a 110 y luego a 99.

- (a) ¿Cuáles son los dos retornos simples, en porcentaje? ¿Cuánto suman?
- (b) ¿Cuánto suman los dos retornos logarítmicos $r_t = \ln(P_t/P_{t-1})$?
  ¿Por qué esa suma sí coincide con el cambio total del precio?

<details>
<summary>Respuesta</summary>

(a) $+10\%$ y $-10\%$; suman cero, aunque el precio bajó 1%.
(b) $\ln(110/100) + \ln(99/110) = \ln(99/100) \approx -0.01$. El logaritmo
convierte productos en sumas: $\ln(a\thinspace b) = \ln a + \ln b$, así que la
suma de retornos logarítmicos es el retorno logarítmico total.
</details>

**2. Media y desviación estándar.** Los datos son $2, 4, 6, 8$.

- (a) Calcula la media y la varianza muestral (dividiendo entre $n - 1$).
- (b) Si multiplicas todos los datos por 10, ¿qué pasa con la media y con la
  desviación estándar? ¿Y si en vez de eso les sumas 10?

<details>
<summary>Respuesta</summary>

(a) Media 5. Desviaciones $-3, -1, 1, 3$; suma de cuadrados 20; varianza
$20/3 \approx 6.7$ y desviación estándar $\approx 2.6$.
(b) Al multiplicar por 10, la media y la desviación estándar se multiplican
por 10 (la varianza, por 100). Al sumar 10, la media aumenta 10 y la
desviación estándar no cambia.
</details>

**3. La distribución normal.** Los retornos diarios de un índice tienen media
0 y desviación estándar de 1%.

- (a) Si fueran normales, ¿qué fracción de los días tendría un retorno mayor
  que 2% en valor absoluto?
- (b) Un día el índice cae 5%. ¿Cuántas desviaciones estándar es eso? Si los
  retornos fueran normales, ¿qué tan seguido esperarías un día así?

<details>
<summary>Respuesta</summary>

(a) Alrededor del 5%: con $z = r/\sigma$, $\mathbb P(\lvert Z\rvert > 2)
\approx 0.046$.
(b) Es $z = -5$. Para una normal, $\mathbb P(Z < -5) \approx 3 \times
10^{-7}$: con unos 250 días de operación al año, una vez cada 10,000 años o
más. En los mercados reales pasa varias veces por década; ése es el tema de
la sesión 1.
</details>

**4. Escalas logarítmicas.** Considera $y = 5\thinspace x^{-2}$ y $y = 5\thinspace
e^{-2x}$.

- (a) ¿Cuál de las dos se ve como una recta en una gráfica log-log (logaritmo
  en los dos ejes)? ¿Qué pendiente tiene?
- (b) ¿En qué escala se ve recta la otra?

<details>
<summary>Respuesta</summary>

(a) La ley de potencias: $\ln y = \ln 5 - 2 \ln x$, una recta de pendiente
$-2$ en log-log.
(b) La exponencial: $\ln y = \ln 5 - 2x$ es una recta en escala semilog
(logaritmo sólo en el eje vertical).
</details>

**5. Probabilidad de exceder un valor.** Una variable $X \ge 1$ cumple
$\mathbb P(X > x) = x^{-1.5}$.

- (a) Calcula $\mathbb P(X > 4)$ y $\mathbb P(X > 100)$.
- (b) ¿Qué fracción de los valores cae entre 1 y 4?

<details>
<summary>Respuesta</summary>

(a) $4^{-1.5} = 1/8 = 0.125$ y $100^{-1.5} = 10^{-3}$.
(b) $1 - 0.125 = 0.875$: la función $\mathbb P(X > x)$ es la CCDF, y uno menos
la CCDF es la función de distribución acumulada.
</details>

**6. La derivada como tasa de cambio.** La fracción $x$ de adoptantes de una
tecnología cumple $dx/dt = x(1 - x)$, con $0 \le x \le 1$.

- (a) Si hoy $x = 0.3$, ¿la fracción sube o baja? ¿A qué ritmo?
- (b) ¿En qué valores de $x$ deja de cambiar? ¿Hacia cuál se mueve si empieza
  cerca de cero?

<details>
<summary>Respuesta</summary>

(a) Sube, a ritmo $0.3 \times 0.7 = 0.21$ por unidad de tiempo.
(b) En $x = 0$ y $x = 1$, donde la derivada vale cero (los puntos fijos).
Entre los dos la derivada es positiva, así que desde cerca de cero la
fracción crece hacia 1.
</details>

**7. La serie geométrica.** Un banco pierde 1 millón; cada ronda de contagio
le hace perder a la red el 60% de lo que se perdió en la ronda anterior.

- (a) ¿Cuánto se pierde en total, sumando todas las rondas?
- (b) ¿Qué pasa si el porcentaje es mayor que 100%?

<details>
<summary>Respuesta</summary>

(a) $1 + 0.6 + 0.6^2 + \dots = 1/(1 - 0.6) = 2.5$ millones.
(b) La suma no converge: cada ronda pierde más que la anterior y el choque
crece sin límite. La sesión 12 usa justamente esta frontera.
</details>

**8. Matrices y redes.** En una red de cuatro bancos, $A_{ij} = 1$ si los
bancos $i$ y $j$ tienen un préstamo entre sí:

$$
A = \begin{pmatrix} 0 & 1 & 1 & 0 \cr 1 & 0 & 1 & 0 \cr 1 & 1 & 0 & 1 \cr 0 & 0 & 1 & 0 \end{pmatrix} .
$$

- (a) Calcula $A\thinspace\mathbf 1$, donde $\mathbf 1 = (1, 1, 1, 1)^\top$.
  ¿Qué significa cada entrada?
- (b) ¿Cuántos préstamos (aristas) hay?

<details>
<summary>Respuesta</summary>

(a) $A\thinspace\mathbf 1 = (2, 2, 3, 1)^\top$: cada entrada es la suma de
una fila, el número de vecinos (grado) de ese banco.
(b) Los grados suman 8 y cada arista se cuenta dos veces: hay 4.
</details>

**9. La función logística y el ruido.** Sea $\sigma(u) = 1/(1 + e^{-u})$.
Un agente elige la opción A con probabilidad $\sigma(\beta\thinspace\Delta)$,
donde $\Delta$ es cuánto más gana con A que con B y $\beta \ge 0$.

- (a) Calcula $\sigma(0)$ y comprueba que $\sigma(u) + \sigma(-u) = 1$.
- (b) Con $\Delta = 1$, ¿qué probabilidad tiene A si $\beta = 0$? ¿Y si
  $\beta$ es enorme? ¿Qué papel juega $\beta$?

<details>
<summary>Respuesta</summary>

(a) $\sigma(0) = 1/2$. Además $\sigma(-u) = e^{-u}/(1 + e^{-u})$, y sumado a
$\sigma(u) = 1/(1 + e^{-u})$ da 1.
(b) Con $\beta = 0$, $1/2$: elige al azar. Con $\beta$ enorme, casi 1: siempre
elige lo mejor. $\beta$ mide qué tan poco ruido hay en la decisión; en física
es el inverso de la temperatura, $\beta = 1/T$.
</details>

**10. Caminatas y promedios de simulaciones.** Una caminata da 100 pasos de
$+1$ o $-1$ con probabilidad $1/2$ cada uno.

- (a) ¿Cuál es la posición media al final? ¿Y la desviación estándar?
- (b) Estimas la posición media promediando 100 simulaciones. ¿Cuántas
  necesitas para que el error típico de ese promedio baje a la mitad?

<details>
<summary>Respuesta</summary>

(a) Media 0. Cada paso tiene varianza 1 y los pasos son independientes, así
que la varianza es 100 y la desviación estándar 10, es decir, $\sqrt n$.
(b) El error del promedio de $M$ simulaciones baja como $1/\sqrt M$: hacen
falta 400, cuatro veces más.
</details>

---

## Qué repasar según tus respuestas

| Si fallaste o dudaste en la pregunta | Repasa | Se usa en las sesiones |
|---|---|---|
| 1 | [C.B1](calculo.md#cb1) (exponenciales y logaritmos) | 1, 3, 6, 8 |
| 2 | [P.B2](probabilidad.md#pb2) (media, varianza y desviación estándar) | Todas |
| 3 | [P.B4](probabilidad.md#pb4) (normal y estandarización), [P.B3](probabilidad.md#pb3) (histogramas) | 1, 2, 3, 5, 8, 14 |
| 4 | [C.B2](calculo.md#cb2) (gráficas en escala logarítmica) | 1, 3, 6, 8, 10 |
| 5 | [P.B1](probabilidad.md#pb1) (variable aleatoria y su distribución), [C.B6](calculo.md#cb6) (la integral como área) | 1, 3, 8, 10 |
| 6 | [C.B3](calculo.md#cb3) (la derivada como tasa de cambio) | 5, 8, 13, 17 |
| 7 | [C.B5](calculo.md#cb5) (sumatorias y la serie geométrica) | 4, 10, 12 |
| 8 | [L.B1–B2](algebra_lineal.md#lb1) (vectores y matrices), [L.B4](algebra_lineal.md#lb4) (matriz de adyacencia) | 4, 7, 10, 12, 15, 18 |
| 9 | [C.B4](calculo.md#cb4) (funciones que reaparecen), [F.B3](fisica_estadistica.md#fb3) (temperatura: ruido contra orden) | 8, 13, 15, 18 |
| 10 | [E.B2](procesos_estocasticos.md#eb2) (caminata aleatoria), [E.B3](procesos_estocasticos.md#eb3) (promediar simulaciones), [P.B5](probabilidad.md#pb5) (ley de los grandes números) | Todas |

Las preguntas no cubren la estadística inferencial ni la física estadística
básica. Si nunca has visto intervalos de confianza o p-valores, lee
[S.B1–B3](estadistica.md#sb1) antes de la sesión 1. Si nunca has visto un
modelo de muchas partículas, lee [F.B1](fisica_estadistica.md#fb1) y
[F.B4](fisica_estadistica.md#fb4) antes de la sesión 8.

## Orientación según tu carrera

Es una guía de lo que suele pasar, no una regla; manda el resultado del
diagnóstico.

- **Economía.** Las preguntas 1, 2 y 3 suelen salir bien. Conviene revisar
  las escalas logarítmicas (C.B2), la serie geométrica (C.B5), las matrices
  (L.B1–B4) y, antes del Módulo II, la física estadística básica (F.B1–B4).
- **Física.** Las preguntas 6, 9 y 10 suelen salir bien. Conviene revisar la
  estimación y las pruebas (S.B1–B5) y la lectura de distribuciones empíricas
  (P.B3, P.B5), que se usan desde la sesión 1.
- **Matemáticas.** Las preguntas de cálculo y álgebra suelen salir bien.
  Conviene revisar la parte con datos (P.B3, P.B5, S.B1–B5) y la intuición
  física de la temperatura y las transiciones de fase (F.B1–B4).

Si fallaste tres preguntas o más, no te preocupes: el curso está pensado para
eso. Lee las claves de la tabla poco a poco, en la semana anterior a la
sesión que las usa, y usa el autodiagnóstico al final de cada área de
repaso.

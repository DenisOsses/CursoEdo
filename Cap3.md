---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

<!--######################################################################################################################################################################################################################################################################################################################################################
-->

+++

(U3)=
# Unidad 3

+++

## Sistemas de Ecuaciones Diferenciales

En este capítulo estudiaremos EDO que modelan situaciones donde aparecen varias variables dependientes interrelacionadas entre sí.

+++

### Introducción 

Considere los dos estanques de la figura: 

```{figure} SistEDO1.png
---
height: 200px
name: TEU
---
Sistema de estanques
```

Suponga que el estanque $A$ contiene un volumen de $50~lts.$ de agua en los que hay disueltos $25~grs.$ de sal y que el estanque $B$ contiene $50~lts$ de agua pura. Hacia los estanques entra y sale líquido (como se indica en la figura); se supone que tanto la mezcla intercambiada entre
los dos estanques como el líquido bombeado hacia fuera del tanque $B$ están bien mezclados. 

Construyamos un PVI que describa la cantidad de gramos $x_1(t)$ y $x_2(t)$ de sal en los estanques $A$ y $B$, respectivamente, en el tiempo $t$. Para ello, procedemos como en el modelo de Mezclas que estudiamos la clase [](EDOLineal), al considerar que el ritmo de cambio de la sal en cada estanque es la diferencia entre la razón de entrada y de salida de la sal. Procedemos a plantear una ecuación para la cantidad de sal en cada estanque y obtenemos:

$$
\begin{eqnarray}
\frac{dx_1}{dt}&=&-\frac{2}{25}x_1+\frac{1}{50}x_2\\ 
\frac{dx_2}{dt}&=&\frac{2}{25}x_1-\frac{2}{25}x_2\\
x_1(0)&=&25~,~x_2(0)=0
\end{eqnarray}
$$ (SistMezclas)

Las ecuaciones [](SistMezclas) se conocen como **Sistema de Ecuaciones Diferenciales Lineales de Primer Orden**. ¿Cómo encontrar simultáneamente, si es posible, su solución $x_1(t)$, $x_2(t)$?  

Otro problema relacionado a sistemas es el de **resortes acoplados**, como vemos en la figura:

```{figure} SistEDO2.png
---
height: 250px
name: TEU
---
Sistema de resortes acoplados
```

Dos masas $m_1$ y $m_2$ están conectadas a dos resortes $A$ y $B$ (de masa despreciable) con constantes de resorte $k_1$ y $k_2$, respectivamente. Los dos resortes están unidos. Sean $x_1($t) y $x_2(t)$ los desplazamientos verticales de las masas desde sus posiciones de equilibrio. Cuando el sistema está en movimiento, el resorte $B$ está sujeto a elongación y compresión: su elongación neta es $x_2-x_1$. Por tanto, se deduce de la ley de Hooke que los resortes $A$ y $B$ ejercen fuerzas $-k_1x_1$ y $k_2(x_2-x_1)$ respectivamente, sobre $m_1$. Si ninguna fuerza externa se aplica al sistema y si ninguna fuerza de amortiguamiento está presente, entonces la fuerza neta en $m_1$ es $-k_1x_1+k_2(x_2-x_1)$. Análogamente, la fuerza neta ejercida en la masa $m_2$ se debe sólo a la elongación neta de $B$; es decir, $-k_2(x_2-x_1)$. Por tanto, por la segunda Ley de Newton [](eq1.1), se tiene que 

$$
\begin{eqnarray}
m_1x_1''&=&-k_1x_1+k_2(x_2-x_1)\\ 
m_2x_2''&=&-k_2(x_2-x_1)
\end{eqnarray}
$$ (SistResortes)

Las ecuaciones [](SistResortes)se conocen como **Sistema de Ecuaciones Diferenciales Lineales de Segundo Orden**. Nuevamente nos preguntamos cómo encontrar las soluciones $x_1($t) y $x_2(t)$ para los desplazamientos verticales de las masas. Además, ¿qué ocurre si se aplica una fuerza externa $f(t)$ sobre $m_2$?

+++

### Teoría Preliminar: Sistemas Lineales

Un **Sistema de Ecuaciones Lineales de Primer Orden** tiene la estructura 

$$
\begin{eqnarray}
\frac{dx_1}{dt}&=&a_{11}(t)x_1+a_{12}(t)x_2+\cdots+a_{1n}(t)x_n+f_1(t)\\ 
\frac{dx_2}{dt}&=&a_{21}(t)x_1+a_{22}(t)x_2+\cdots+a_{2n}(t)x_n+f_2(t)\\
&\vdots&\\
\frac{dx_n}{dt}&=&a_{n1}(t)x_1+a_{n2}(t)x_2+\cdots+a_{nn}(t)x_n+f_n(t)
\end{eqnarray}
$$ (SistEDOLineal)

Diremos que [](SistEDOLineal) es simplemente un **sistema lineal**. Los coeficientes $a_ij(t)$ así como las funciones $f_i(t)$ son continuas en un intervalo común $I$. Cuando $f_i(t)=0$, $\forall~i=1, 2,\ldots, n$, se dice que el sistema lineal es **homogéneo**; de otro modo es **no homogéneo**.

Todo sistema lineal se puede representar matricialmente si hacemos 

$$
\mathbf{X}=\begin{pmatrix}x_1(t)\\ x_2(t)\\ \vdots \\ x_n(t)\end{pmatrix}~,~\mathbf{A}=\begin{pmatrix} a_{11}(t)&a_{12}(t)&\cdots&a_{1n}(t)\\ a_{21}(t)&a_{22}(t)&\cdots&a_{2n}(t)\\ \vdots&\vdots&\ddots&\vdots\\ a_{n1}(t)&a_{n2}(t)&\cdots&a_{nn}(t)\end{pmatrix}~,~ \mathbf{F}=\begin{pmatrix}f_1(t)\\ f_2(t)\\ \vdots \\ f_n(t)\end{pmatrix}
$$

y escribimos $\mathbf{X}'=\mathbf{A}\mathbf{X}+\mathbf{F}$. El sistema homogéneo es $\mathbf{X}'=\mathbf{A}\mathbf{X}$. Un **vector solución** $\mathbf{X}$ en un intervalo $I$ es un vector cuyos elementos son funciones derivables que satisfacen todas las ecuaciones del sistema $\mathbf{X}'=\mathbf{A}\mathbf{X}+\mathbf{F}$ en $I$.


Por ejemplo, el sistema [](SistMezclas) se puede escribir como 

$$
\mathbf{X}'=\begin{pmatrix}-\frac{2}{25}&\frac{1}{50}\\ \frac{2}{25}&-\frac{2}{25}\end{pmatrix}\mathbf{X} ~\Leftrightarrow~\mathbf{X}'=\mathbf{A}\mathbf{X}
$$

Para resolver el sistema matricialmente, suponemos una solución de la forma $\mathbf{X}=\begin{pmatrix}k_1\\ k_2\end{pmatrix}e^{\lambda t}=\mathbf{K}e^{\lambda t}$ para $k_1,k_2,\lambda\in\mathbb{R}$. Entonces $\mathbf{X}'=\mathbf{K}\lambda e^{\lambda t}$ y si reemplazamos en $\mathbf{X}'=\mathbf{A}\mathbf{X}$ obtenemos

$$
\mathbf{X}'=\mathbf{A}\mathbf{X}~\Rightarrow~\mathbf{K}\lambda e^{\lambda t}=\mathbf{A}\mathbf{K}e^{\lambda t}~\Rightarrow~\mathbf{A}\mathbf{K}=\lambda\mathbf{K}
$$ 

¿Reconoce esta última igualdad? Es justamente la ecuación que se plantea para obtener los **vectores propios** $\mathbf{K}$ y sus **valores propios** asociados $\lambda$ correspondientes a la matriz $\mathbf{A}$. Recordamos que podemos encontrar tales valores y vectores propios al resolver 

$$
\det(\mathbf{A}-\lambda\mathbf{I})=0
$$

Así, obtenemos los vectores propios $\mathbf{K}_1=\begin{pmatrix}-\frac{1}{2}\\ 1\end{pmatrix}$ con valor propio asociado $\lambda_1=-\frac{3}{25}$, y $\mathbf{K}_2=\begin{pmatrix}\frac{1}{2}\\ 1\end{pmatrix}$ con valor propio asociado $\lambda_2=-\frac{1}{25}$. Por lo tanto, la solución general del sistema es 

$$
\mathbf{X}=c_1\begin{pmatrix}-\frac{1}{2}\\ 1\end{pmatrix}e^{-\frac{3}{25}t}+c_2\begin{pmatrix}\frac{1}{2}\\ 1\end{pmatrix}e^{-\frac{1}{25}t}
$$

Finalmente, si reemplazamos las condiciones iniciales $x_1(0)=25$ y $x_2(0)=0$, obtenemos que la solución del sistema [](SistMezclas) es 

$$
\mathbf{X}=\begin{pmatrix}x_1(t)\\ x_2(t)\end{pmatrix}=-25\begin{pmatrix}-\frac{1}{2}\\ 1\end{pmatrix}e^{-\frac{3}{25}t}+25\begin{pmatrix}\frac{1}{2}\\ 1\end{pmatrix}e^{-\frac{1}{25}t}
$$

```{code-cell}
:tags: [SistMezclas]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
import numpy as np
from sympy import Matrix, Rational, symbols, exp, simplify
import sympy as sp

# Activar impresión bonita
sp.init_printing()

# Definir la matriz A usando fracciones exactas ## SE PUEDE MODIFICAR
A = Matrix([[Rational(-2, 25), Rational(1, 50)],
            [Rational(2, 25), Rational(-2, 25)]])

# Calcular los vectores y valores propios con fracciones
val_propios = A.eigenvals()  # Valores propios
vec_propios = A.eigenvects()  # Vectores propios

# Extraer los valores y vectores propios en listas asegurando el orden correcto
val_propios_lista = []
vectores_propios_lista = []

for val_propios, mult, vectores in vec_propios:
    val_propios_lista.append(val_propios)
    vectores_propios_lista.append(vectores[0])  # Cada valor propio puede tener múltiples vectores asociados, elegimos el primero

# Asignar los vectores propios en el orden de los valores propios
v1 = vectores_propios_lista[0]  # Primer vector propio
v2 = vectores_propios_lista[1]  # Segundo vector propio

# Condiciones iniciales ## SE PUEDE MODIFICAR
x0 = Matrix([25, 0])

# Resolver para las constantes c1 y c2 usando SymPy
C = Matrix([[v1[0], v2[0]],
            [v1[1], v2[1]]]).inv() * x0

# Expresiones de x1(t) y x2(t)
t = symbols('t')

def x1(t):
    return simplify(C[0] * v1[0] * exp(val_propios_lista[0] * t) + C[1] * v2[0] * exp(val_propios_lista[1] * t))

def x2(t):
    return simplify(C[0] * v1[1] * exp(val_propios_lista[0] * t) + C[1] * v2[1] * exp(val_propios_lista[1] * t))

# Expresión general para x1(t) y x2(t) en formato más legible
x1_exp = x1(t)
x2_exp = x2(t)

# Imprimir las expresiones en un formato más adecuado
print("Solución x1(t) y x2(t):")
display(x1_exp)
display(x2_exp)
```
+++

### Teoría General

La teoría de sistemas de $n$ ecuaciones diferenciales de primer orden es similar a la de las ecuaciones diferenciales de $n$-ésimo orden: Si $t_0\in I$ entonces el PVI asociado es 

$$
\mathbf{PVI}~~~~\left\{\begin{array}{ccc}\mathbf{X}'(t)&=&\mathbf{A}(t)\mathbf{X}(t)+\mathbf{F}(t)\\&&\\ \mathbf{X}(t_0)&=&\mathbf{X}_0\end{array}\right.
$$ (PVIsist)

**Teorema**: Si las funciones componentes de $\mathbf{A}(t)$ y $\mathbf{F}(t)$ son continuas en $I$ y $t_0\in I$, entonces el PVI [](PVIsist) tiene una solución única en $I$.

**Teorema**: Si el sistema homogéneo $\mathbf{X}'(t)=\mathbf{A}(t)\mathbf{X}(t)$ tiene $n$ soluciones LI: $\mathbf{X}_1, \mathbf{X}_2,\ldots, \mathbf{X}_n$ (es decir su wronskiano $W(\mathbf{X}_1, \mathbf{X}_2,\ldots, \mathbf{X}_n)\neq0$) entonces éstas forman un conjunto fundamental de soluciones y 

$$
\mathbf{X}_H=c_1\mathbf{X}_1+c_2\mathbf{X}_2+\cdots+c_n\mathbf{X}_n~~,~~c_i\in\mathbb{R}
$$

es la solución general del sistema homogéneo. Para resolver $\mathbf{X}'(t)=\mathbf{A}(t)\mathbf{X}(t)+\mathbf{F}(t)$ completamente, se debe determinar una solución particular $\mathbf{X}_P$, con lo cual 

$$
\mathbf{X}=\mathbf{X}_H+\mathbf{X}_P
$$

es la **solución general** del sistema no homogéneo.

+++

## Sistemas Homogéneos

Dado el sistema homogéneo $\mathbf{X}'(t)=\mathbf{A}\mathbf{X}(t)$, donde $\mathbf{A}$ es una <u>**matriz de constantes**</u>, ¿es posible determinar una solución de la forma $\mathbf{X}=\mathbf{K}e^{\lambda t}$? Sabemos que esto lleva a calcular $\det(\mathbf{A}-\lambda\mathbf{I})=0$ para determinar los valores propios de $\mathbf{A}$ y luego obtener los vectores propios $\mathbf{K}$ tales que $(\mathbf{A}-\lambda\mathbf{I})\mathbf{K}=\mathbf{0}$. Estudiaremos 3 casos:

### Valores Propios Reales y Distintos: 

Si $\mathbf{A}$ (matriz de $n\times n$) tiene $n$ valores propios distintos $\lambda_1,\lambda_2,\ldots,\lambda_n$, entonces siempre se puede encontrar un conjunto de $n$ vectores propios LI $\mathbf{K}_1,\mathbf{K}_2,\ldots,\mathbf{K}_n$ tales que la solución del sistema homogéneo es

$$
\mathbf{X}_H=c_1\mathbf{K}_1e^{\lambda_1 t}+c_2\mathbf{K}_2e^{\lambda_2 t}+\cdots+c_n\mathbf{K}_ne^{\lambda_n t}~~,~~c_i\in\mathbb{R}~,~t\in\mathbb{R}.
$$

**Plano de Fase**: Las funciones componentes $\mathbf{x}_1(t), \mathbf{x}_2(t),\ldots,\mathbf{x}_n(t)$, $t\in I$ de la solución del sistema pueden interpretarse como las ecuaciones paramétricas de curvas en $\mathbb{R}^n$. Para cada condición inicial, tenemos una trayectoria particular que indica el flujo de los puntos sobre la curva. El conjunto de las trayectorias se denomina **plano o diagrama de fase** para el sistema.

```{admonition} Ejercicio Aplicado
En Python, esboce el plano de fase del sistema de mezclas [](SistMezclas). ¿Qué observa para la condición inicial dada?, ¿qué ocurre con las curvas solución en el largo plazo?
```
```{code-cell}
:tags: [PlanoFase1]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Definir la matriz del sistema
A = np.array([[-2/25, 1/50],
              [2/25, -2/25]])

# Calcular los valores propios y vectores propios
eigenvalues, eigenvectors = np.linalg.eig(A)

# Crear un grid de puntos
x = np.linspace(-30, 30, 10)
y = np.linspace(-30, 30, 10)
X, Y = np.meshgrid(x, y)

# Definir las ecuaciones diferenciales del sistema
U = (-2/25) * X + (1/50) * Y
V = (2/25) * X + (-2/25) * Y

# Configurar la gráfica del plano de fase
plt.figure(figsize=(8, 8))
plt.streamplot(X, Y, U, V, color='b')  # Graficar el campo vectorial

# Condiciones iniciales para diferentes trayectorias
initial_conditions = [(5, 25), (30, 30), (25, 0), (10, 30), (20, 20), (10, 10)]

t = np.linspace(0, 10, 200)

# Definimos el sistema de ecuaciones diferenciales
def system(X, t):
    x, y = X
    dxdt = (-2/25)*x + (1/50)*y
    dydt = (2/25)*x + (-2/25)*y
    return [dxdt, dydt]

# Graficamos las trayectorias continuas
for ic in initial_conditions:
    traj = odeint(system, ic, t)
    plt.plot(traj[:, 0], traj[:, 1], lw=2)

# Graficar los ejes formados por los vectores propios. Se multiplica por +-40 para mejor visualización
for eigenvector in eigenvectors.T:
    plt.plot([-40*eigenvector[0], 40*eigenvector[0]],
             [-40*eigenvector[1], 40*eigenvector[1]], 'r--', lw=2)

# Añadir detalles a la gráfica
plt.xlim([0, 30])
plt.ylim([0, 30])
plt.axhline(0, color='black',linewidth=1)
plt.axvline(0, color='black',linewidth=1)
plt.title('Plano de Fase del Sistema')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()
```

### Valores Propios Reales Repetidos

Si $\lambda$ es un valor propio con multiplicidad $m$ (repetido $m$ veces), podemos tener algunos casos como los siguientes (no cubren todas las opciones):

1. Para algunas matrices $\mathbf{A}$ es posible encontrar $m$ vectores propios $\mathbf{K}_1, \mathbf{K}_2,\ldots,\mathbf{K}_m$ LI. Para este caso, la solución contiene la combinación lineal 

$$
c_1\mathbf{K}_1e^{\lambda t}+c_2\mathbf{K}_2e^{\lambda t}+\cdots+c_m\mathbf{K}_me^{\lambda t}
$$

2. Si sólo hay un único vector propio asociado al valor propio $\lambda$ con multiplicidad $m$, entonces siempre se pueden encontrar $m$ soluciones LI de la forma:

$$
\begin{eqnarray*}
\mathbf{X}_1&=&\mathbf{K}_{11}e^{\lambda t}\\
\mathbf{X}_2&=&\mathbf{K}_{21}te^{\lambda t}+\mathbf{K}_{22}e^{\lambda t}\\
&\vdots&\\
\mathbf{X}_m&=&\mathbf{K}_{m1}\frac{t^{m-1}}{(m-1)!}e^{\lambda t}+\mathbf{K}_{m2}\frac{t^{m-2}}{(m-2)!}e^{\lambda t}+\cdots+ \mathbf{K}_{mm}e^{\lambda t}
\end{eqnarray*}
$$

En particular, si $\mathbf{A}$ es de $2\times2$, entonces $\lambda$ tiene multiplicidad 2 y las soluciones LI son: 

$$
\begin{eqnarray*}
\mathbf{X}_1&=&\mathbf{K}_1e^{\lambda t}\\
\mathbf{X}_2&=&\mathbf{K}_1te^{\lambda t}+\mathbf{K}_2e^{\lambda t}
\end{eqnarray*}
$$

donde $\mathbf{K}_2$ se obtiene al resolver: $(\mathbf{A}-\lambda\mathbf{I})\mathbf{K}_2=\mathbf{K}_1$.

También, si $\mathbf{A}$ es de $3\times3$ y $\lambda$ tiene multiplicidad 3, las soluciones LI son:

$$
\begin{eqnarray*}
\mathbf{X}_1&=&\mathbf{K}_1e^{\lambda t}\\
\mathbf{X}_2&=&\mathbf{K}_1te^{\lambda t}+\mathbf{K}_2e^{\lambda t}\\
\mathbf{X}_3&=&\mathbf{K}_1\frac{t^2}{2}e^{\lambda t}+\mathbf{K}_2te^{\lambda t}+\mathbf{K}_3e^{\lambda t}
\end{eqnarray*}
$$

donde $\mathbf{K}_3$ se obtiene al resolver: $(\mathbf{A}-\lambda\mathbf{I})\mathbf{K}_3=\mathbf{K}_2$.

```{admonition} Ejercicio Teórico
Determine la solución del sistema 

$$
\mathbf{X}'=\begin{pmatrix}3&-18\\2&-9\end{pmatrix}\mathbf{X}
$$

Luego esboce su plano de fase
```

```{code-cell}
:tags: [PlanoFase2]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
from sympy import Matrix, Rational
import sympy as sp

# Activar impresión bonita
sp.init_printing()

# Definir la matriz A usando fracciones exactas ## SE PUEDE MODIFICAR
A = Matrix([[Rational(3,1), Rational(-18,1)],
            [Rational(2, 1), Rational(-9, 1)]])

# Calcular los vectores y valores propios con fracciones
val_propios = A.eigenvals()  # Valores propios
vec_propios = A.eigenvects()  # Vectores propios

vec_propios #El output tiene el siguiente orden: valor propio, multiplicidad, vector propio
```

```{code-cell}
:tags: [PlanoFase22]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
import numpy as np
import pylab as plt

# Definimos el sistema de ecuaciones diferenciales
def system(X, t):
    x, y = X
    dxdt = 3*x -18*y
    dydt = 2*x -9*y
    return [dxdt, dydt]

# Definimos el espacio de puntos para las condiciones iniciales
x_vals = np.linspace(-20, 20, 20)
y_vals = np.linspace(-20, 20, 20)

# Creamos una malla de puntos
X, Y = np.meshgrid(x_vals, y_vals)

# Calculamos las derivadas para cada punto de la malla
u = 3*X -18*Y
v = 2*X -9*Y

# Graficamos el campo vectorial usando streamplot
plt.figure(figsize=(8, 8))
plt.streamplot(X, Y, u, v, color='b', linewidth=1)

# Etiquetas de los ejes
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('Plano de Fase')

# Ajustamos los límites del gráfico
plt.xlim([-20, 20])
plt.ylim([-20, 20])

# Mostramos el gráfico
plt.grid()
plt.show()
```

### Valores Propios Complejos

Si $\lambda=\alpha+\beta i$, $\alpha,\beta\in\mathbb{R}$, es un valor propio complejo de la matriz $\mathbf{A}$ y $\mathbf{K}$ es su vector propio correspondiente entonces

$$
\begin{eqnarray*}
\mathbf{X}_1&=&\left(\mathbf{B}_1\cos(\beta t)-\mathbf{B}_2\sin(\beta t)\right)e^{\alpha t}\\
\mathbf{X}_2&=&\left(\mathbf{B}_2\cos(\beta t)+\mathbf{B}_1\sin(\beta t)\right)e^{\alpha t}
\end{eqnarray*}
$$

son soluciones LI del sistema homogéneo para $t\in\mathbb{R}$, donde $\mathbf{B}_1=\frac{1}{2}\left(\mathbf{K}_1+ \overline{\mathbf{K}}_1\right)=Re(\mathbf{B}_1)$ y $\mathbf{B}_2=\frac{i}{2}\left(-\mathbf{K}_1+\overline{\mathbf{K}}_1\right)=Im(\mathbf{B}_1)$.

```{admonition} Ejercicio Teórico
Determine la solución del sistema 

$$
\mathbf{X}'=\begin{pmatrix}2&8\\-1&-2\end{pmatrix}\mathbf{X}
$$

Luego esboce su plano de fase
```

```{code-cell}
:tags: [PlanoFase3]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
import numpy as np
import pylab as plt

# Definimos el sistema de ecuaciones diferenciales
def system(X, t):
    x, y = X
    dxdt = 2*x +8*y
    dydt = -1*x -2*y
    return [dxdt, dydt]

# Definimos el espacio de puntos para las condiciones iniciales
x_vals = np.linspace(-20, 20, 20)
y_vals = np.linspace(-20, 20, 20)

# Creamos una malla de puntos
X, Y = np.meshgrid(x_vals, y_vals)

# Calculamos las derivadas para cada punto de la malla
u = 2*X +8*Y
v = -1*X -2*Y

# Graficamos el campo vectorial usando streamplot
plt.figure(figsize=(8, 8))
plt.streamplot(X, Y, u, v, color='b', linewidth=1)

# Etiquetas de los ejes
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('Plano de Fase')

# Ajustamos los límites del gráfico
plt.xlim([-20, 20])
plt.ylim([-20, 20])

# Mostramos el gráfico
plt.grid()
plt.show()
```

**Nota**: Si consideramos el sistema homogéneo 

$$
\mathbf{X}'=\begin{pmatrix}a&b\\c&d\end{pmatrix}\mathbf{X}
$$


notamos que $\mathbf{A}\mathbf{X}=\mathbf{0}$ si $\mathbf{X}=(0,0)$ (para $\det(\mathbf{A})=0$) y este es el único **punto crítico** del sistema. Podemos establecer su **estabilidad** en términos de la traza $T=a+d$ y el determinante $D=ad-bc$ de $\mathbf{A}$. Así, la ecuación característica asociada a los valores propios de $\mathbf{A}$ es $\lambda^2-T\lambda+D=0$. Esta cuadrática en $\lambda$ es importante ya que sus raíces son los valores propios de $\mathbf{A}$ y su carácter -raíces reales distintas, repetidas o complejas- depende del signo de $\Delta=T^2-4D$. En resumen:

1. **Nodo** si $D>0$ y $\Delta\geq0$.
2. **Punto Silla** si $D<0$.
3. **Punto Espiral** si $T\neq0$ y $\Delta<0$.
4. **Centro** si $T=0$ y $D>0$.

Si ponemos en el eje horizontal a $T$ y en el vertical a $D$, obtenemos el **plano Traza-Determinante** que nos permite estudiar la estabilidad del sistema en el largo plazo. En resumen:

1. **Atractor** (asintóticamente estable) si $D>0$ y $T<0$.
2. **Estable** si $D>0$ y $T=0$.
3. **Repulsor** si $D<0$ o $T>0$.

Gráficamente:

```{figure} PlanoTrazaDet.png
---
height: 500px
name: TrazaDet
---
Estabilidad de un Sistema de Ecuaciones Diferenciales de Primer Orden
```

```{admonition} Ejercicio Teórico
Determine el carácter y estabilidad del punto crítico $(0,0)$ para los sistemas

$$
\mathbf{X}'=\begin{pmatrix}3&-18\\2&-9\end{pmatrix}\mathbf{X}\quad\text{y}\quad\mathbf{X}'=\begin{pmatrix}2&8\\-1&-2\end{pmatrix}\mathbf{X}
$$
```
+++

## Sistemas Autónomos

Consideremos el sistema de ecuaciones diferenciales **autónomas** de orden 2:

$$
\begin{eqnarray*}
\frac{dx}{dt}&=&F(x,y)\\
\frac{dy}{dt}&=&G(x,y)\quad\Leftrightarrow\quad\mathbf{X}'(t)=\mathbf{f}(\mathbf{x})~,~\mathbf{x}(t_0)=\mathbf{x}^0\\
x(t_0)=x_0~&,&y(t_0)=y_0
\end{eqnarray*}
$$ (SistAutonomo)

donde $F$ y $G$ son continuas en $D\in\mathbb{R}$ y $(x_0,y_0)\in D$. Esto implica que el sistema [](SistAutonomo) tiene solución única (una curva $(x(t),y(t))$ para $t\in I$ que contiene a $t_0$). 

Si existen puntos $(x,y)$ tales que $F(x,y)=0$ y $G(x,y)=0$ o, equivalentemente $\mathbf{f}(\mathbf{x})=\mathbf{0}$, los denominamos **puntos críticos** del sistema. También en estos puntos $\mathbf{X}'=\mathbf{0}$ por lo que corresponden a soluciones constantes o de **equilibrio** del sistema. 

```{admonition} Ejercicio Teórico
Considere el sistema de ecuaciones 

$$
\begin{eqnarray*}
\frac{dx}{dt}&=&-(x-y)(1-x-y)\\
\frac{dy}{dt}&=&x(2+y)
\end{eqnarray*}
$$

Determine sus puntos críticos y esboce su plano de fase. A partir del plano fase, establezca el carácter y estabilidad de los puntos críticos.
```

Código para determinar los puntos críticos:

```{code-cell}
:tags: [PuntosCriticos]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
from sympy import symbols, Eq, solve

# Definimos las variables simbólicas
x, y = symbols('x y')

# Definimos las ecuaciones del sistema
eq1 = Eq(-(x - y) * (1 - x - y), 0)  # x' = 0
eq2 = Eq(x * (2 + y), 0)             # y' = 0

# Resolvemos el sistema de ecuaciones
puntos_criticos = solve([eq1, eq2], (x, y))
puntos_criticos
```

Código para determinar los puntos críticos:

```{code-cell}
:tags: [PlanoFase4]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
import numpy as np
import matplotlib.pyplot as plt

# Definimos las funciones del sistema de ecuaciones diferenciales
def dx_dt(x, y):
    return -(x - y) * (1 - x - y)

def dy_dt(x, y):
    return x * (2 + y)

# Definimos el rango de valores para x e y
x_values = np.linspace(-5, 5, 20)
y_values = np.linspace(-5, 5, 20)

# Creamos una malla de puntos
X, Y = np.meshgrid(x_values, y_values)

# Calculamos los valores del campo vectorial
U = dx_dt(X, Y)
V = dy_dt(X, Y)

# Graficamos el plano de fase
plt.figure(figsize=(8, 8))
plt.streamplot(X, Y, U, V, color='b')
plt.title("Retrato de fase del sistema")
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(True)
plt.show()
```

Nos interesa determinar el comportamiento de las soluciones curvas solución del sistema [](SistAutonomo) en la cercanía de sus puntos críticos $(x_0,y_0)$. Para ello, la idea es **linealizar** (mediante el plano tangente) $F(x,y)$ y $G(x,y)$ en torno a tales puntos:

$$
\begin{eqnarray*}
F(x,y)&=&\cancelto{0}{F(x_0,y_0)}+F_x(x_0,y_0)(x-x_0)+F_y(x_0,y_0)(y-y_0)+h_1(x,y)\\
G(x,y)&=&\cancelto{0}{G(x_0,y_0)}+G_x(x_0,y_0)(x-x_0)+G_y(x_0,y_0)(y-y_0)+h_2(x,y)
\end{eqnarray*}
$$

Así, podemos reescribir el sistema como

$$
\begin{eqnarray*}
\frac{d}{dt}\begin{pmatrix}x-x_0\\y-y_0\end{pmatrix}=\begin{pmatrix}F_x(x_0,y_0)& F_y(x_0,y_0)\\ G_x(x_0,y_0)& G_y(x_0,y_0)\end{pmatrix}\begin{pmatrix}x-x_0\\y-y_0\end{pmatrix}+\begin{pmatrix}h_1(x,y)\\h_2(x,y)\end{pmatrix}
\end{eqnarray*}
$$

o vectorialmente $\mathbf{U}'=\mathbf{J}(x_0,y_0)\mathbf{U}+\mathbf{H}(\mathbf{x})$, donde $\mathbf{J}$ es la **matriz jacobiana** de las funciones $F$ y $G$, y asumimos que $\frac{||\mathbf{H}||}{||\mathbf{x}-\mathbf{x}^0||}\to0$ cuando $\mathbf{x}\to\mathbf{x}^0$-(el sistema de EDO es **localmente lineal**-. De este modo, nos quedamos con el sistema linealizado $\mathbf{U}'=\mathbf{J}(x_0,y_0)\mathbf{U}$ para analizar localmente los puntos críticos. 

Una tabla resumen comparativa del tipo de puntos críticos y su estabilidad se presenta a continuación para sistemas lineales y localmente lineales:

|                           | Sistema Lineal         |                         | Sistema Localmente Lineal |                         |
|---------------------------|------------------------|-------------------------|---------------------------|-------------------------|
| Valores Propios           | Tipo                   | Estabilidad             | Tipo                      | Estabilidad             |
| $r_1>r_2>0$               | Nodo                   | Inestable               | Nodo                      | Inestable               |
| $r_1<r_2<0$               | Nodo                   | Asintóticamente Estable | Nodo                      | Asintóticamente Estable |
| $r_2<0<r_1$               | Punto Silla            | Inestable               | Punto Silla               | Inestable               |
| $r_1=r_2>0$               | Nodo Propio o Impropio | Inestable               | Nodo o Espiral            | Inestable               |
| $r_1=r_2<0$               | Nodo Propio o Impropio | Asintóticamente Estable | Nodo o Espiral            | Asintóticamente Estable |
| $r_1,r_2=\lambda\pm\mu i$ |                        |                         |                           |                         |
| $\lambda>0$               | Espiral                | Inestable               | Espiral                   | Inestable               |
| $\lambda<0$               | Espiral                | Asintóticamente Estable | Espiral                   | Asintóticamente Estable |
| $\lambda=0$               | Centro                 | Estable                 | Centro o Espiral          | Indeterminado           |

```{admonition} Ejercicio Teórico
Para el sistema de ecuaciones del ejercicio anterior

$$
\begin{eqnarray*}
\frac{dx}{dt}&=&-(x-y)(1-x-y)\\
\frac{dy}{dt}&=&x(2+y)
\end{eqnarray*}
$$

Establezca el carácter y estabilidad de los puntos críticos.
```

**Nota**: La linealización es un proceso muy útil para estudiar el carácter de los puntos críticos en sistemas de ecuaciones diferenciales de orden $2$, $3$ o superior. Existen muchos sistemas de ecuaciones no lineales que modelan diversas situaciones interesantes: Las [Ecuaciones de Lorenz](https://en.wikipedia.org/wiki/Lorenz_system) usadas inicialmente en Meteorología, el [modelo de Hindmarsh–Rose](https://en.wikipedia.org/wiki/Hindmarsh%E2%80%93Rose_model) para la actividad neuronal o el [modelo de Lotka-Volterra](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations) que describe sistemas biológicos donde interactúan especies. 

+++

## Actividad de Cierre Unidad 3

+++

### Problema 1: Modelamiento Sistema de Mezcla

Considere el sistema de dos tanques de la figura:

```{figure} tanques.png
---
height: 120px
name: TEU
---
Sistema de tanques
```

(a) Si la tasa de entrada y salida de ambos tanques es $r$, el tanque 1 contiene un volumen $V_1$ de agua (con $x_0$ gramos de sal disuelta), el tanque 2 contiene un volumen $V_2$ de agua (con $y_0$ gramos de sal disuelta), establezca un sistema de ecuaciones lineales que modele la cantidad de sal $x$ e $y$ en cada tanque, respectivamente, en el instante $t$.

(b) Considere que $V_1=50$ y $V_2=25$, $r=10$, $x_0=15$ e $y_0=0$, determine explícitamente las soluciones $x(t)$ e $y(t)$, y luego realice su gráfico.

(c) Ahora suponga que agua pura entra al tanque 1 a una tasa un porcentaje menor a $r$ ($kr$ con $0\leq k<1$) y que sale mezcla del tanque 2 a la misma tasa. Para los valores del inciso antrior, estudie el carácter y estabilidad del punto crítico $(0,0)$. Esboce el plano de fase en este caso.

+++

### Problema 2: Modelamiento Competencia entre 2 Especies

Imagina un ecosistema en el que coexisten dos especies de plantas $A$ y $B$ que compiten por un recurso limitado, como los nutrientes del suelo. Ambas especies dependen de los nutrientes para crecer, pero el recurso no es infinito, y a medida que la población de una especie crece, afecta la disponibilidad de nutrientes para la otra especie. Queremos modelar este sistema usando un conjunto de dos ecuaciones diferenciales no lineales que describan cómo varían las poblaciones de plantas con el tiempo bajo competencia.

Sean $x(t)$ e $y(t)$ la poblaciones de las especies de plantas $A$ y $B$ en el instante, respectivamente. Suponga que:

1. Cada especie, cuando no está en competencia con la otra, crece según una tasa de crecimiento natural limitada por la cantidad de nutrientes en el suelo. Para la especie $A$, su tasa de crecimiento en condiciones ideales es $r_Ax$ donde $r_A$  es la tasa de crecimiento natural de $A$. Análogamente $r_By$ para $B$.
2. Dado que el recurso es limitado, el crecimiento de cada población estará sujeto a competencia con la otra. La competencia entre especies puede modelarse utilizando términos que describan la reducción en la tasa de crecimiento debido a la presencia de la otra especie:

(i) **Competencia intraespecífica**: Cada especie afecta su propio crecimiento debido a la sobrepoblación. A medida que aumenta la densidad de una especie, los nutrientes se vuelven menos accesibles para sus propios individuos. Este efecto se puede modelar con un término cuadrático de la forma $a_{AA}x^2$, donde $a_{AA}$ representa la competencia dentro de la especie $A$. Similarmente, la especie $B$ está sujeta a $a_{BB}y^2$.

(ii) **Competencia interespecífica**: La especie $A$ también afecta el crecimiento de la especie $B$ y viceversa -hay una interacción entre ellas-. Este efecto puede modelarse con un término que dependa de ambas poblaciones como $a_{AB}xy$ para la influencia de $A$ sobre $B$, y $a_{BA}xy$ para la influencia de $B$ sobre $A$. 

A partir de las consideraciones anteriores:

(a) Establezca un sistema de ecuaciones diferenciales no lineal que modele esta competencia entre 2 especies de plantas.

(b) Considere los valores: $r_A=0.2$, $r_B=0.1$, $a_{AA}=0.1$, $a_{BB}=0.05$, $a_{AB}=0.01$ y $a_{BA}=0.02$. Determine los puntos críticos del sistema.

(c) Linealice el sistema y estudie la estabilidad de sus puntos críticos.

(d) Esboce el planos de fase correspondiente. Interprete sus resultados en relación a la extinción o subsistencia de las plantas $A$ y $B$ en el largo plazo.



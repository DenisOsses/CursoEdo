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
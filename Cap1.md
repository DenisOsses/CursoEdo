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

(U1)=
# Unidad 1

(Intro)=
## Clase 1: Introducción

**Segunda ley de Newton**: Mutationem motus proportionalem esse vi motrici impressae, \& fieri secundum lineam rectam qua vis illa imprimitur.

"El cambio de movimiento es directamente proporcional a la fuerza motriz impresa y ocurre según la línea recta a lo largo de la cual aquella fuerza se imprime."

+++

Actualmente decimos: La aceleración $a$ de un cuerpo de masa $m$ es proporcional a la fuerza total $F$ ejercida sobre él: 

```{math}
:label: eq1.1
F=ma
```

Supongamos que un cuerpo de masa $m$ cae bajo la única influencia de la gravitación. La única fuerza que actúa sobre él es $mg$, donde $g$ es la aceleración de gravedad. Si $y$ es la altura medida hacia abajo desde una posición fija, entonces $v=\frac{dy}{dt}$ es el ritmo de cambio de su posición y su aceleración $a=\frac{dv}{dt}=\frac{d^2y}{dt^2}$ es el ritmo de cambio de la velocidad. Reemplazando en [](eq1.1) obtenemos 

```{math}
:label: eq1.2
m\frac{d^2y}{dt^2}=mg
```

Esta es una **Ecuación Diferencial Ordinaria (EDO)** de **Segundo Orden**.

La ecuación [](eq1.2) puede ser reescrita como una EDO de **Primer Orden**:

```{math}
:label: eq1.3
m\frac{dv}{dt}=mg
```

Si admitimos que el aire ejerce una fuerza de resistencia proporcional a la velocidad, la fuerza total que actúa sobre el cuerpo es $mg-kv$ y la ecuación [](eq1.3) queda  

\begin{equation*}
    m\frac{dv}{dt}=mg-kv 
\end{equation*}

¿Es posible determinar la velocidad $v(t)$ del cuerpo en cualquier instante $t$ si su velocidad inicial es $v(0)=v_0$?

Si consideramos esta condición, tenemos el **Problema de Valor Inicial (PVI)**

$$\mathbf{PVI}~~~~\left\{\begin{array}{ccc}m\dfrac{dv}{dt}&=&mg-kv\\&&\\ v(0)&=&v_0\end{array}\right.$$

Sin resolver explícitamente el PVI, podemos visualizar esta solución mediante el siguiente código de Python:

```{code-cell}
:tags: [Ejer1]
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# función que entrega dv/dt
def modelo(v,t,k):
    dvdt = g - k * v/m
    return dvdt

# condición inicial
v0 = 0

#parámetros
g = 9.8
m = 1 

# valores de tiempo
t = np.linspace(0,30)

# soluciones EDO para distintos k
k = 0.1
v1 = odeint(modelo,v0,t,args=(k,))
k = 0.2
v2 = odeint(modelo,v0,t,args=(k,))
k = 1
v3 = odeint(modelo,v0,t,args=(k,))

# gráfico
plt.plot(t,v1,'r',linewidth=1,label='k=0.1')
plt.plot(t,v2,'b',linewidth=1,label='k=0.2')
plt.plot(t,v3,'g',linewidth=1,label='k=1')
plt.xlabel('Tiempo')
plt.ylabel('v(t)')
plt.legend(shadow=True)
plt.show()
```

Para determinar una solcuón del PVI usamos el **método de separación de variables**, obteniendo 

$$
v(t)=\frac{mg}{k}+\left(v_0-\frac{mg}{k}\right)e^{-\frac{k}{m}t}
$$

Esta es una **solución explícita** del PVI y su **intervalo de definición** es $t\in\mathbb{R}_0^+$.

Python también nos permite encontrar una solución general del PVI:

```{code-cell}
:tags: [Ejer1-1]
from sympy import symbols, Function, dsolve

# Símbolos
t, g, k , m, v0= symbols("t, g, k, m, v0")

# v es función del tiempo
v = Function('v')(t)

# Derivada de v
dvdt = v.diff(t)

# Ecuación diferencial
eq = dvdt - g + k*v/m

# Solución de la EDO con la condición inicial v(0)=v0
dsolve(eq, v, ics={v.subs(t, 0): v0})
```

```{admonition} Ejercicio
Establezca un PVI que exprese la posición $y(t)$ del cuerpo en cualquier instante $t$ si su posición inicial es $y(0)=y_0$. Resuelva el PVI anterior. ¿Qué ocurre cuando $t\to\infty$? 
```
+++

## Clase 2
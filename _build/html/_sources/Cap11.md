# Unidad 1

(Intro)=
## Introducción

**Segunda ley de Newton**: Mutationem motus proportionalem esse vi motrici impressae, \& fieri secundum lineam rectam qua vis illa imprimitur.

"El cambio de movimiento es directamente proporcional a la fuerza motriz impresa y ocurre según la línea recta a lo largo de la cual aquella fuerza se imprime."

+++

Actualmente decimos: La aceleración $a$ de un cuerpo de masa $m$ es proporcional a la fuerza total $F$ ejercida sobre él: 

$$
F=ma
$$ 

Supongamos que un cuerpo de masa $m$ cae bajo la única influencia de la gravitación. La única fuerza que actúa sobre él es $mg$, donde $g$ es la aceleración de gravedad. Si $y$ es la altura medida hacia abajo desde una posición fija, entonces $v=\frac{dy}{dt}$ es el ritmo de cambio de su posición y su aceleración $a=\frac{dv}{dt}=\frac{d^2y}{dt^2}$ es el ritmo de cambio de la velocidad. Reemplazando en (1.1) obtenemos 

$$
\tag{1.2}
    m\frac{d^2y}{dt^2}=mg
$$

Esta es una **Ecuación Diferencial Ordinaria (EDO)** de **Segundo Orden**.

La ecuación (1.2) puede ser reescrita como una EDO de **Primer Orden**:

$$
\tag{1.3}
    m\frac{dv}{dt}=mg
$$

Si admitimos que el aire ejerce una fuerza de resistencia proporcional a la velocidad, la fuerza total que actúa sobre el cuerpo es $mg-kv$ y la ecuación (1.3) queda  

\begin{equation*}
    m\frac{dv}{dt}=mg-kv 
\end{equation*}

¿Es posible determinar la velocidad $v(t)$ del cuerpo en cualquier instante $t$ si su velocidad inicial es $v(0)=v_0$?

Si consideramos esta condición, tenemos el **Problema de Valor Inicial (PVI)**

$$\mathbf{PVI}~~~~\left\{\begin{array}{ccc}m\dfrac{dv}{dt}&=&mg-kv\\&&\\ v(0)&=&v_0\end{array}\right.$$

Para resolver el PVI usamos el **método de separación de variables**, obteniendo 

$$
v(t)=\frac{mg}{k}+\left(v_0-\frac{mg}{k}\right)e^{-\frac{k}{m}t}
$$

Esta es una **solución explícita** del PVI y su **intervalo de definición** es $t\in\mathbb{R}_0^+$. Algunas preguntas son:

```{admonition} Ejercicio
Establezca un PVI que exprese la posición $y(t)$ del cuerpo en cualquier instante $t$ si su posiciión inicial es $y(0)=y_0$. Resuelva el PVI anterior. ¿Qué ocurre cuando $t\to\infty$? 
```

Podemos visualizar esta solución mediante el siguiente código de Python:
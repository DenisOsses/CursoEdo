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

(Clase1.2)=
# Clase 1.2 

### Método de Variables Separables. Ecuaciones Lineales de Primer Orden. 

**Ley de Enfriamiento de Newton** Esta ley establece que la rapidez con que cambia la temperatura $T(t)$ de un cuerpo en el instante $t$ es proporcional a la diferencia entre la temperatura de dicho cuerpo y la del medio $T_m$ que lo rodea.

Así, obtenemos la EDO

```{math}
:label: eq2.1
\frac{dT}{dt}=k(T-T_m)
```

Si $T_m$ es constante, la ecuación [](eq2.1) es de variables separables; en cambio, si $T_m=f(t)$ (varía en el tiempo), es una EDO **Lineal de Primer Orden** 
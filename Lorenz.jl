__precompile__(true)
"""
LZ

Módulo que resuelve las ecuaciones de Lorenz clásicas

dx/dt = sigma * (y-x),
dy/dt = r * x - y - x * z,
dz/dt = x * y -b * z.

Donde r, sigma, b son parámetros reales positivos.

La solución se hace por medio del método de Taylor con polinomios de Taylor
de grado indicado por el usuario.

Véase las funciones

lorenz(t0,tf,x0,y0,z0,r,sigma,b,p),
diflorenz(t0,tf,x01,y01,z01,x02,y02,z02,r,sigma,b,p)

"""
module LZ
	export lorenz
	export diflorenz
using TS

"""
lorenz(t0,tf,x0,y0,z0,r,sigma,b,p)

Resuelve las ecuaciones de Lorenz con parámetros r, sigma, b y condiciones
iniciales (x0,y0,z0) del tiempo inicial t0 al tiempo final tf. Las solución
se hace mediante el método de Taylor con polinomios de Taylor hasta orden p.

Devuelve cuatro listas

t, x, y, z

con las soluciones.
"""
function lorenz(t0::Real, tf::Real, x0::Real, y0::Real, z0::Real, r::Real, sigma::Real, b::Real, p::Int)
    # Inicializamos la lista de respuestas.
    tl = [t0]
    xl = [x0]
    yl = [y0]
    zl = [z0]
    
    # Ejecutamos mientras que el tiempo inicial sea menor al final.
    while t0 <= tf
        
        # Empezamos los Taylors con la condición inicial.
        x = Taylor(x0)
        y = Taylor(y0)
        z = Taylor(z0)
        
        # A continuación se calcula la serie de Taylor hasta orden p
        # usando las relaciones de recurrencia entre las x,y,z y sus derivadas.
        for i in range(1,p)
           # En cada paso se vuelve a calcular la serie de Taylor de dx/dt = f(x,y,z),
           # para dx/dt, dy/dt, dz/dt cada vez a mayor orden.
            dx = sigma*(y-x)
            dy = r*x - y - x*z
            dz = x*y - b*z
           # De ésta se extrae el nuevo coeficiente de x(t), y(t), z(t) que se anexa.
            x = Taylor(push!(x.taylor_vec,dx.taylor_vec[i]/i))
            y = Taylor(push!(y.taylor_vec,dy.taylor_vec[i]/i))
            z = Taylor(push!(z.taylor_vec,dz.taylor_vec[i]/i))
        end
    
        # Ahora se escoge un paso. Como se recomienda, se toman los dos 
        # últimos términos de la serie para calcular h1 y h2, mismo que
	# se hace para x,y,z.
        hx1 = (1/2)*(eps(1.0)/abs(x.taylor_vec[p+1]))^(1/p)
        hx2 = (1/2)*(eps(1.0)/abs(x.taylor_vec[p]))^(1/(p-1))
        hy1 = (1/2)*(eps(1.0)/abs(y.taylor_vec[p+1]))^(1/p)
        hy2 = (1/2)*(eps(1.0)/abs(y.taylor_vec[p]))^(1/(p-1))
        hz1 = (1/2)*(eps(1.0)/abs(z.taylor_vec[p+1]))^(1/p)
        hz2 = (1/2)*(eps(1.0)/abs(z.taylor_vec[p]))^(1/(p-1))
        # Luego se toma h como el mínimo de las anteriores y se suma
        # al tiempo.
        h = minimum([hx1, hx2, hy1, hy2, hz1, hz2])
        t0 += h
    
        # Ahora sumamos las series de acuerdo al método de Horner:
        sx = x.taylor_vec[p+1]
        sy = y.taylor_vec[p+1]
        sz = z.taylor_vec[p+1]
        for k in range(1,p)
            sx = x.taylor_vec[p+1-k] + h*sx
            sy = y.taylor_vec[p+1-k] + h*sy
            sz = z.taylor_vec[p+1-k] + h*sz
        end
    
        # x0 es ahora la suma de la serie: x(t0+h), y análogo para y,z.
        x0 = sx
        y0 = sy
        z0 = sz
        
        # Vamos a poner una condición por si la solución diverge.
        if (isnan(x0) | (abs(x0) == Inf)) | (isnan(y0) | (abs(y0) == Inf)) | (isnan(z0) | (abs(z0) == Inf))
            println("La solución diverge en t = ", t0)
            println("No pude llegar al tiempo final t = ", tf)
            # Que devuelva las listas antes de la divergencia.
            return tl, xl, yl, zl
        end
        
        # Se anexa a la lista de respuestas.
        push!(tl,t0)
        push!(xl,x0)
        push!(yl,y0)
        push!(zl,z0)
    
    end
    # Se devuelven las lista.
    return tl, xl, yl, zl
end

"""
lorenzdif(t0,tf,x01,y01,z01,x02,y02,z02,r,sigma,b,p)

Resuelve las ecuaciones de Lorenz con parámetros r, sigma, b para dos condiciones
iniciales diferentes. La primera es un sistema con (x01,y01,z01) y el segundo
sistema es (x02,y02,z02). Los dos sistemas se resuelven a los mismos intervalos
de tiempo dados en la lista t.

Devuelve siete listas

t, x1, y1, z1, x2, y2, z2

con las soluciones de los dos sistemas y la lista de tiempos.
"""
function diflorenz(t0::Real, tf::Real, x01::Real, y01::Real, z01::Real, x02::Real, y02::Real, z02::Real, r::Real, sigma::Real, b::Real, p::Int)
# Esta función es prácticamente igual a la anterior, lorentz que calcula dos soluciones.
# La diferencia central es que el paso de tiempo se selcciona como el mínimo del paso de
# cada solución.

    tl = [t0]
    x1l = [x01]
    y1l = [y01]
    z1l = [z01]
    x2l = [x02]
    y2l = [y02]
    z2l = [z02]
    
    while t0 <= tf
        
        x1 = Taylor(x01)
        y1 = Taylor(y01)
        z1 = Taylor(z01)
        
        x2 = Taylor(x02)
        y2 = Taylor(y02)
        z2 = Taylor(z02)
        
        for i in range(1,p)
            dx1 = sigma*(y1-x1)
            dy1 = r*x1-y1-x1*z1
            dz1 = x1*y1 - b*z1
            dx2 = sigma*(y2-x2)
            dy2 = r*x2-y2-x2*z2
            dz2 = x2*y2 - b*z2

            x1 = Taylor(push!(x1.taylor_vec,dx1.taylor_vec[i]/i))
            y1 = Taylor(push!(y1.taylor_vec,dy1.taylor_vec[i]/i))
            z1 = Taylor(push!(z1.taylor_vec,dz1.taylor_vec[i]/i))
            x2 = Taylor(push!(x2.taylor_vec,dx2.taylor_vec[i]/i))
            y2 = Taylor(push!(y2.taylor_vec,dy2.taylor_vec[i]/i))
            z2 = Taylor(push!(z2.taylor_vec,dz2.taylor_vec[i]/i))
        end
    
        hx11 = (1/2)*(eps(1.0)/abs(x1.taylor_vec[p+1]))^(1/p)
        hx12 = (1/2)*(eps(1.0)/abs(x1.taylor_vec[p]))^(1/(p-1))
        hy11 = (1/2)*(eps(1.0)/abs(y1.taylor_vec[p+1]))^(1/p)
        hy12 = (1/2)*(eps(1.0)/abs(y1.taylor_vec[p]))^(1/(p-1))
        hz11 = (1/2)*(eps(1.0)/abs(z1.taylor_vec[p+1]))^(1/p)
        hz12 = (1/2)*(eps(1.0)/abs(z1.taylor_vec[p]))^(1/(p-1))
        hx21 = (1/2)*(eps(1.0)/abs(x2.taylor_vec[p+1]))^(1/p)
        hx22 = (1/2)*(eps(1.0)/abs(x2.taylor_vec[p]))^(1/(p-1))
        hy21 = (1/2)*(eps(1.0)/abs(y2.taylor_vec[p+1]))^(1/p)
        hy22 = (1/2)*(eps(1.0)/abs(y2.taylor_vec[p]))^(1/(p-1))
        hz21 = (1/2)*(eps(1.0)/abs(z2.taylor_vec[p+1]))^(1/p)
        hz22 = (1/2)*(eps(1.0)/abs(z2.taylor_vec[p]))^(1/(p-1))

	# Esta es la diferencia central a la función anterior, se selecciona
	# la misma h considerando las seis distintas. Así los pasos de tiempo
	# son iguales en las dos soluciones.
        h = minimum([hx11, hx12, hy11, hy12, hz11, hz12, hx21, hx22, hy21, hy22, hz21, hz22])
        t0 += h
    
        sx1 = x1.taylor_vec[p+1]
        sy1 = y1.taylor_vec[p+1]
        sz1 = z1.taylor_vec[p+1]
        sx2 = x2.taylor_vec[p+1]
        sy2 = y2.taylor_vec[p+1]
        sz2 = z2.taylor_vec[p+1]
        for k in range(1,p)
            sx1 = x1.taylor_vec[p+1-k] + h*sx1
            sy1 = y1.taylor_vec[p+1-k] + h*sy1
            sz1 = z1.taylor_vec[p+1-k] + h*sz1
            sx2 = x2.taylor_vec[p+1-k] + h*sx2
            sy2 = y2.taylor_vec[p+1-k] + h*sy2
            sz2 = z2.taylor_vec[p+1-k] + h*sz2
        end
    
        x01 = sx1
        y01 = sy1
        z01 = sz1
        x02 = sx2
        y02 = sy2
        z02 = sz2
        
        if (isnan(x01) | (abs(x01) == Inf)) | (isnan(y01) | (abs(y01) == Inf)) | (isnan(z01) | (abs(z01) == Inf))
            println("La solución diverge en t = ", t0)
            println("No pude llegar al tiempo final t = ", tf)
            # Que devuelva las listas antes de la divergencia.
            return tl, x1l, y1l, z1l, x2l, y2l, z2l
        end
        
        if (isnan(x02) | (abs(x02) == Inf)) | (isnan(y02) | (abs(y02) == Inf)) | (isnan(z02) | (abs(z02) == Inf))
            println("La solución diverge en t = ", t0)
            println("No pude llegar al tiempo final t = ", tf)
            # Que devuelva las listas antes de la divergencia.
            return tl, x1l, y1l, z1l, x2l, y2l, z2l
        end
        
        push!(tl,t0)
        push!(x1l,x01)
        push!(y1l,y01)
        push!(z1l,z01)
        push!(x2l,x02)
        push!(y2l,y02)
        push!(z2l,z02)
    
    end
    return tl, x1l, y1l, z1l, x2l, y2l, z2l
end

# Termina el módulo.
end

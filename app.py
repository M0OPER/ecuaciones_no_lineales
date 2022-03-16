#Algoritmo de Bisección y Newton
#Author; Edwin M., José Jimenez

# -*- coding: utf-8 -*-

#import tkinter
import tkinter as tk
from tkinter import * 
from tkinter import ttk

#import sympy
from sympy import *
from sympy import sympify
from sympy.parsing.sympy_parser import parse_expr

#import matplotlib
import matplotlib.pyplot as plt

#import numpy
import numpy as np

def resolver():
  def mostrarBiseccion():
    for i in tblBiseccion.get_children():
      tblBiseccion.delete(i)

    i = 0
    for f in resultB[0]:
      tblBiseccion.insert(parent='', index='end', iid=i, text='', values=(f[0], "{0:.5f}".format(f[1]), "{0:.5f}".format(f[2]), "{0:.5f}".format(f[3]), f[4]))
      i += 1
    graficar(resultB[1])

  def mostrarNewton():
    for i in tblNewton.get_children():
      tblNewton.delete(i)

    j = 0
    for ff in resultN[0]:
      tblNewton.insert(parent='', index='end', iid=j, text='', values=(ff[0], "{0:.5f}".format(ff[1]), "{0:.5f}".format(ff[2]), ff[3]))
      j += 1
    graficar(resultB[1])

  try:
    limite1 = int(tryLimite1.get())
    limite2 = int(tryLimite2.get())
    toleran = float(tryTolerancia.get())
    fun     = parse_expr(tryFuncion.get())

    bise = chkBiseccion.instate(['selected'])
    newt = chkNewton.instate(['selected'])

    resultB = biseccion(limite1, limite2, fun, toleran)
    resultN = newton(limite1, limite2, fun, toleran)

    if bise and newt:
      mostrarBiseccion()
      mostrarNewton()
      btnGrafica["state"] = "active"
    elif bise:
      mostrarBiseccion()
      btnGrafica["state"] = "active"
    elif newt :
      mostrarNewton()
      btnGrafica["state"] = "active"
    else:
      pass

    lblMensaje.config(text = 'Ok ::: => Datos procesados correctamente' )
  except ValueError:
    lblMensaje.config(text = 'Error ::: Introduce un numero')

def limpiar():
  try:
    tryLimite1.delete(0, END)
    tryLimite2.delete(0, END)
    tryFuncion.delete(0, END)
    for i in tblNewton.get_children():
      tblNewton.delete(i)

    for i in tblBiseccion.get_children():
      tblBiseccion.delete(i)
    lblMensaje.config(text = 'Ok ::: => Datos limpiados' )
  except Exception as e:
    lblMensaje.config(text = 'Error ::: => ' + e )  

def biseccion(limI, limF, funcion, toleran):

  res     = []
  tabla   = []
  grafica = []

  # INGRESO
  tolerancia = toleran

  # PROCEDIMIENTO
  tramo = limF-limI

  fa = resolverFx(funcion, limI)
  fb = resolverFx(funcion, limF)
  i  = 1

  while (tramo > tolerancia):
    c  = (limI + limF) / 2
    fc = resolverFx(funcion, c)
    tabla.append([i, limI, limF, c, tramo])
    grafica.append([i,limI,c,limF,fa,fc,fb,tramo])
    i = i + 1

    cambia = fa * fc
    if(cambia < 0):
      limF = c
      fb   = fc
    else:
      limI = c
      fa   = fc
    tramo = limF-limI

  c  = (limI+limF) / 2
  fc = resolverFx(funcion, c)
  tabla.append([i, limI,limF, c, tramo])
  grafica.append([i,limI,c,limF,fa,fc,fb,tramo])


  raiz = c

  # Tabla con formato
  res.append(tabla)
  res.append(grafica)
  return res

def newton(limI, limF, funcion, toleran):

  res      = []
  tabla    = []
  grafica  = []

  # INGRESO
  tolerancia = toleran

  derivada = derivadaFx(funcion)

  tramo = abs(2 * tolerancia)

  i = 1
  while (tramo >= tolerancia):
    derivada = derivadaFx(funcion)
    xnuevo   = limF - (float(resolverFx(funcion, limF))/ float(resolverDfx(derivada, limF)))
    tramo  = abs(xnuevo - limF)
    tabla.append([i, limF, xnuevo, tramo])
    limF = xnuevo
    i += 1

  res.append(tabla)
  res.append(derivadaFx(funcion))

  return res

def graficar(grafica):

  grafica = np.array(grafica)
  fig = plt.figure()

  xi = grafica[:, 2]
  yi = grafica[:, 5]

  # ordena los puntos para la grafica
  orden = np.argsort(xi)
  xi = xi[orden]
  yi = yi[orden]

  plt.plot(xi,yi)
  plt.plot(xi,yi,'o')
  plt.axhline(0, color="black")

  plt.xlabel('X')
  plt.ylabel('Y')
  plt.title('Grafica en f(x)')
  plt.grid()
  
  fig.set_size_inches(4.3, 4)
  plt.savefig('grafica.png')
  img = PhotoImage(file = 'grafica.png')
  imgBiseccion.config(image = img, width=500, height=500)
  imgBiseccion.image = img

def mostrarGrafica():
  try:
    limite1 = int(tryLimite1.get())
    limite2 = int(tryLimite2.get())
    toleran = float(tryTolerancia.get())
    fun     = parse_expr(tryFuncion.get())
    res     = biseccion(limite1, limite2, fun, toleran)
    grafica = np.array(res[1])
    fig = plt.figure()

    xi = grafica[:, 2]
    yi = grafica[:, 5]

    # ordena los puntos para la grafica
    orden = np.argsort(xi)
    xi = xi[orden]
    yi = yi[orden]

    plt.plot(xi,yi)
    plt.plot(xi,yi,'o')
    plt.axhline(0, color="black")

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Grafica en f(x)')
    plt.grid()
    
    fig.set_size_inches(7, 7)
    plt.savefig('big_grafica.png')
    plt.show()
  except Exception as e:
    raise

def resolverFx(f, sustitucion):
  try:
    x = Symbol('x')
    e = Symbol('e')
    res = f.replace(e, 2.71828)
    res = res.replace(x, sustitucion)
    return res
  except Exception as e:
    raise e
  
def resolverDfx(f, sustitucion):
  x = Symbol('x')
  e = Symbol('e')
  res = f.replace(e, 2.71828)
  res = res.replace(x, sustitucion)
  return res

def derivadaFx(f):
  x = Symbol('x')
  res = diff(f, x)
  return res

app = Tk()
app.geometry('2000x2000')
app.title("Ecuaciones no lineales - Biseccion y Newton")

#VENTANA PRINCIPAL
vp     = Frame(app)

canvas = Canvas(vp, bg='white', scrollregion=(0,0,500,500))

hbar = Scrollbar(vp, orient = HORIZONTAL)
hbar.config(command = canvas.xview)

vbar = Scrollbar(vp, orient = VERTICAL)
vbar.config(command = canvas.yview)

#CONFIGURANDO EL ROW
vp.columnconfigure(0, weight = 1)
vp.columnconfigure(1, weight = 1)
vp.columnconfigure(2, weight = 3)
vp.columnconfigure(3, weight = 3)

#COLUMN 1
lblMensaje = Label(vp, text="")
lblMensaje.grid(column=0, row=0, columnspan=3)
Label(vp, text="Limite inferior:").grid(column=0, row=1)
Label(vp, text="Limite superior:").grid(column=0, row=2)
Label(vp, text="Funcion:").grid(column=0, row=3)
btnLimpiar = Button(vp, text = "LIMPIAR", command = limpiar, bg='yellow').grid(column=0, row=4)
Label(vp, text="Autor: Edwin M., José Jimenez").grid(column=0, row=7, columnspan=3, sticky = W)
btnGrafica = Button(vp, text = "INTERACTUAR CON LA GRAFICA", command = mostrarGrafica, bg='red', fg="white", state = DISABLED)
btnGrafica.grid(column=3, row=7)

#COLUMN 2
tryLimite1  = ttk.Entry(vp)
tryLimite1.grid(column=1, row=1, sticky = W)
tryLimite2  = ttk.Entry(vp)
tryLimite2.grid(column=1, row=2, sticky = W)
tryFuncion  = ttk.Entry(vp)
tryFuncion.grid(column=1, row=3, sticky = W)
btnResolver = Button(vp, text = "RESOLVER", command = resolver, bg='green', fg="white")
btnResolver.grid(column=1, row=4, sticky = W)

#COLUMN 3
checkBiseccion = IntVar(value=1)
chkBiseccion = ttk.Checkbutton(vp, text="BISECCION", variable = checkBiseccion)
chkBiseccion.grid(column=2, row=1, sticky = W)
checkNewton = IntVar(value=0)
chkNewton    = ttk.Checkbutton(vp, text="NEWTON", variable = checkNewton)
chkNewton.grid(column=2, row=2, sticky = W)
Label(vp, text='TOLERANCIA').grid(column=2, row=3, sticky = W)
tryTolerancia = ttk.Entry(vp)
tryTolerancia.grid(column=2, row=4, sticky = W)
tryTolerancia.insert(0, "0.001")

#COLUMN 4
Label(vp, text="TABLA DE OPERADORES: \n |                 Suma = + | Resta = - | \n | Multiplicacion = * | División = / | \n | Potencia = ** |").grid(column=3, row=1, rowspan=3, sticky = W)
Label(vp, text="COMO ESCRIBIR LAS FUNCONES: \n x**3 + 4*x**2 - 10 \n e**x - x**2 + 3*x -2 ").grid(column=3, row=4, sticky = W)

#IMAGENES
imgB         = PhotoImage(file = 'default.png')
imgBiseccion = Label(vp, image = imgB)
imgBiseccion.grid(column=2, row=5, columnspan=2, rowspan=2)

#TABLA BISECCION
tblBiseccion = ttk.Treeview(vp)
tblBiseccion['columns'] = ('contador', 'limiteInferior', 'limiteSuperior', 'raiz', 'tramo')
tblBiseccion.column("contador", anchor=CENTER, width=150)
tblBiseccion.column("limiteInferior", anchor=CENTER, width=140)
tblBiseccion.column("limiteSuperior", anchor=CENTER, width=140)
tblBiseccion.column("raiz", anchor=CENTER, width=100)
tblBiseccion.column("tramo", anchor=CENTER, width=200)

tblBiseccion.heading("contador",text="# BISECCION", anchor=CENTER)
tblBiseccion.heading("limiteInferior", text="Limite inferior",anchor=CENTER)
tblBiseccion.heading("limiteSuperior", text="Limite superior",anchor=CENTER)
tblBiseccion.heading("raiz", text="Raiz",anchor=CENTER)
tblBiseccion.heading("tramo", text="Tramo",anchor=CENTER)
tblBiseccion.grid(column=0, row=5, columnspan=2)

#TABLA NEWTON
tblNewton = ttk.Treeview(vp)
tblNewton['columns'] = ('contador', 'limiteSuperior', 'raiz', 'tramo')
tblNewton.column("contador", anchor=CENTER, width=200)
tblNewton.column("limiteSuperior", anchor=CENTER, width=140)
tblNewton.column("raiz", anchor=CENTER, width=100)
tblNewton.column("tramo", anchor=CENTER, width=300)

tblNewton.heading("contador",text="# NEWTON", anchor=CENTER)
tblNewton.heading("limiteSuperior", text="Limite superior",anchor=CENTER)
tblNewton.heading("raiz", text="Raiz",anchor=CENTER)
tblNewton.heading("tramo", text="Tramo",anchor=CENTER)
tblNewton.grid(column=0, row=6, columnspan=2)

vp.pack(side="right", fill="both", expand = True)

app.mainloop()
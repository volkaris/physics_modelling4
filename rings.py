import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from tkinter import *
import tkinter as Tk
from LightPipes import *
import math

root = Tk.Tk()

labda = 600 * nm  # длина волны света
size = 5 * mm  # размер области моделирования
N = 300  # 300 на 300 пикселей

R = 100 * cm  # радиус кривизны линзы
nfilm = 1.0  # показатель преломления пленки (для воздуха 1)

LABDA = DoubleVar()
LABDA.set(labda / nm)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

canvas = FigureCanvasTkAgg(fig, master=root)
canvas._tkcanvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
v = StringVar()

d = np.ndarray((N, N))
X = np.arange(-size / 2, size / 2, size / N)
Y = X

for i in range(0, N - 1):
    for j in range(0, N - 1):
        r2 = X[i] * X[i] + Y[j] * Y[j]
        d[i][j] = r2 / 2 / R  # толщина пленки в каждой точке сетки


def calculate_intensity(labda):
    F = Begin(size, labda, N)
    F1 = F
    Phi = Phase(F)
    k = 2 * math.pi / labda
    p2 = 2 * nfilm * k

    for i in range(1, N):
        for j in range(1, N):
            Phi[i][j] = p2 * d[i][j] + math.pi

    F = SubPhase(Phi, F)
    F = BeamMix(F1, F)
    I = Intensity(1, F)
    return I


def radial_intensity(I):
    center = N // 2
    radial_intensity = np.zeros(center)
    for r in range(center):
        radial_intensity[r] = I[center, center + r]
    return radial_intensity


def get_colormap(labda):
    wavelength = labda / nm
    if wavelength < 380 or wavelength > 780:
        raise ValueError("Wavelength out of visible range")
    color = wavelength_to_rgb(wavelength)
    return mcolors.LinearSegmentedColormap.from_list("custom", [(0, "black"), (1, color)])


def wavelength_to_rgb(wavelength):
    gamma = 0.8
    intensity_max = 255
    factor = 0.0
    R, G, B = 0, 0, 0

    if 380 <= wavelength <= 440:
        R = -(wavelength - 440) / (440 - 380)
        G = 0.0
        B = 1.0
    elif 440 <= wavelength <= 490:
        R = 0.0
        G = (wavelength - 440) / (490 - 440)
        B = 1.0
    elif 490 <= wavelength <= 510:
        R = 0.0
        G = 1.0
        B = -(wavelength - 510) / (510 - 490)
    elif 510 <= wavelength <= 580:
        R = (wavelength - 510) / (580 - 510)
        G = 1.0
        B = 0.0
    elif 580 <= wavelength <= 645:
        R = 1.0
        G = -(wavelength - 645) / (645 - 580)
        B = 0.0
    elif 645 <= wavelength <= 780:
        R = 1.0
        G = 0.0
        B = 0.0

    if 380 <= wavelength <= 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif 420 <= wavelength <= 645:
        factor = 1.0
    elif 645 <= wavelength <= 780:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 645)

    R = round(intensity_max * ((R * factor) ** gamma))
    G = round(intensity_max * ((G * factor) ** gamma))
    B = round(intensity_max * ((B * factor) ** gamma))

    return (R / 255.0, G / 255.0, B / 255.0)


def TheExample(kostil):
    global I
    labda = LABDA.get() * nm


    I = calculate_intensity(labda)
    radial_intensity_mono = radial_intensity(I)

    ax1.clear()
    colormap = get_colormap(labda)
    ax1.imshow(I, cmap=colormap)
    ax1.axis('off')
    ax1.axis('equal')
    ax1.set_title('Кольца Ньютона для монохроматического света')

    ax2.clear()
    ax2.plot(radial_intensity_mono)
    ax2.set_title('График зависимости интенсивности от радиальной координаты')
    ax2.set_xlabel('Радиальная координата')
    ax2.set_ylabel('Интенсивность')

    canvas.draw()


Scale(root,
      takefocus=1,
      orient='horizontal',
      label='длина волны [нанометры]',
      length=200, from_=380.0, to=780.0,
      resolution=0.1,
      variable=LABDA,
      cursor="hand2",
      command=TheExample).pack()



Label(root, textvariable=v).pack(pady=50)

TheExample(0)
root.mainloop()
root.destroy()

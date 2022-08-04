import csv
import math as m

# First we need to import matplotlib module using "pip install matplotlib" in the terminal
import matplotlib.pyplot as plt

# input signal
x = []

# output signal
y = []

# given kernel
h = [1 / 16, 4 / 16, 6 / 16, 4 / 16, 1 / 16]
file = open('data.csv')
reader = csv.reader(file)
next(reader)

# just to append the data in the lists x and y
for r in reader:
    x.append(float(r[0]))
    y.append(float(r[1]))
print("All the values of the given input signal:", x)
print()

# plotting the graph of x
x0 = 0
x193 = 193
L1 = []
for i in range(x0, x193):
    L1.append(i)
plt.plot(L1, x, label="X")

"""Part 1 of the question"""

# Denoising then Deblurring:

# Denoise

denoise1 = []

# for denoising we have selected size 3

for i in range(len(x)):

    if i == 0:

        denoise1.append((y[0] + y[0] + y[1]) / 3)

    elif i == len(x) - 1:

        denoise1.append((y[191] + y[192] + y[192]) / 3)

    else:

        denoise1.append((y[i] + y[i - 1] + y[i + 1]) / 3)

# this is how we denoise the signal by taking the average of size 3

# Deblurring

# For deblurring and getting x[n], we need inverse fourier transform of ((fourier transform[y[n]])/(fourier transform[h[n]])) and that we will do in seperate functions.
# We will use complex numbers of python then evaluate


def F_Dn(w):
    Summation = complex(0, 0)
    j = complex(0, 1)
    for i in range(0, 193):

        # As we have only 193 inputs so the summation of -infinity to +infinity is 0 to 192 terms
        Summation = Summation + (denoise1[i] * (m.e ** ((-j) * w * i)))

    return Summation


# For the fourier transform of h[n], Here we need the mid element of h[n] as the first element , so we need to do some changes to the loop calculations according to that


def F_kernal(w):
    Summation = complex(0, 0)
    j = complex(0, 1)
    for i in range(-2, 3):
       #When i is -2, we need zeroth element from list h[n]

        Summation = Summation + (h[i + 2] * (m.e ** ((-j) * w * i)))
    return Summation


#So we have the fourier transform now, next we need inverse fourier transform as mentioned above to approximate x[n], so moving ahead with that

#We will use summation property of integration to solve the integration of Inverse Fourier transform by taking small parts of omega then summing


def Inv1():
    X = []
    j = complex(0, 1)
    for n in range(193):
        sum = 0

        #We are sampling here by 1000 elements
        for samp in range(0, 1000):
            w = (2 * m.pi * samp) / 1000
            mod = (m.sqrt((F_kernal(w).real) ** 2 + (F_kernal(w).imag) ** 2))

            # To avoid values close to infinity, we have kept a threshold value of 0.314 in the denominator and the use the max the function
            r = (((F_kernal(w).real) * (F_Dn(w).real)) - ((F_Dn(w).imag) * (F_kernal(w).imag))) / max(0.314, mod)
            i = (((F_kernal(w).imag) * (F_Dn(w).real)) + ((F_Dn(w).imag) * (F_kernal(w).real))) / max(0.314, mod)

            # Real and imaginary parts
            v = complex(r, i)
            val = (v * (m.e ** (j * w * n)))
            sum = sum + val.real
        X.append(sum / 1000)
    return X

print()
print("X1: ", Inv1())
L5 = Inv1()
plt.plot(L1, L5, label="X1")
print()

"""Part 2 of the question"""

"""Deblurring then Denoising"""

#Deblurring


def F_y(w):
    Summation = complex(0, 0)
    j = complex(0, 1)
    for i in range(0,193):

        # As we have only 193 inputs so the summation of -infinity to +infinity is 0 to 192 terms
        Summation = Summation + (y[i] * ((m.e) ** ((-j) * w * i)))

    return Summation


# Now fourier transform is known of both y[n],h[n], now calculating x[n], through it by similar procedure as done before:

def Inv2():
    X = []
    j = complex(0, 1)
    for n in range(193):
        sum = 0

        # We are sampling here by 1000 elements
        for samp in range(0, 1000):
            w = (2 * m.pi * samp) / 1000
            mod = (m.sqrt((F_kernal(w).real) ** 2 + (F_kernal(w).imag) ** 2))

            # To avoid values close to infinity, we have kept a threshold value of 0.314 in the denominator and the use the max the function
            r = (((F_kernal(w).real) * (F_y(w).real)) - ((F_y(w).imag) * (F_kernal(w).imag))) / max(0.314, mod)
            i = (((F_kernal(w).imag) * (F_y(w).real)) + ((F_y(w).imag) * (F_kernal(w).real))) / max(0.314, mod)

            # Real and imaginary parts
            v = complex(r, i)
            val = (v * (m.e ** (j * w * n)))
            sum += val.real
        X.append(sum / 1000)
    return X

# Denoising

y1 = Inv2()

# This is the deblurred signal
denoise2 = []

# for denoising we have selected size 3

for i in range(len(x)):

    if i == 0:

        denoise2.append((y1[0] + y1[0] + y1[1]) / 3)

    elif i == len(x) - 1:

        denoise2.append((y1[191] + y1[192] + y1[192]) / 3)

    else:

        denoise2.append((y1[i] + y1[i - 1] + y1[i + 1]) / 3)

print("X2:  ", denoise2)

L6 = Inv2()

plt.plot(L1, denoise2, label="X2")

plt.legend()

plt.show()


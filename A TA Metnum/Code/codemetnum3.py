#!/usr/bin/env python
# coding: utf-8

# In[22]:


import sys
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
from numpy import float32, single, double, array, zeros, diag, diagflat, dot
from pprint import pprint
from scipy.linalg import solve
from math import sin

def modul_2():
    def MetodeSetengahInterval(X1,X2):
        X1=X1
        X2=X2
        error = 1
        iterasi = 0
        while(error > 0.0001):
            iterasi +=1
            #-0.0923x^3 + 5.179x^2 - 95.099x + 571.92
            FXi = (-0.0923*(float(X1)**3))+(5.179*(float(X1)**2))-(95.099*(X1))+571.92
            FXii = (-0.0923*(float(X2)**3))+(5.179*(float(X2)**2))-(95.099*(X2))+571.92

            #===Interval Metode Setengah Interval===
            Xt = (X1+X2)/2

            #===Memeriksa Persamaan Baru===
            FXt = (-0.0923*(float(Xt)**3))+(5.179*(float(Xt)**2))-(95.099*(Xt))+571.92

            #===Update Interval===
            if FXi * FXt > 0:
                X1 = Xt
            elif FXi * FXt < 0:
                X2 = Xt
            else:
                print ("Akar Penyelesaisan:", Xt)

            #===Check Konvergensi
            if FXt < 0:
                error = FXt * (-1)
            else:
                error = FXt
            if iterasi > 100:
                print("angka tak hingga")
                break
            print(iterasi, "|", FXi, "|", FXii, "|", Xt, "|", FXt, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar Persamaan: ", Xt)
        print("Toleransi Error: ", error)

    def MetodeInterpolasiLinear(X1):
        X1 = X1
        X2 = X1 + 1
        error = 1
        iterasi = 0
        while(error > 0.0001):
            iterasi +=1
            #-0.0923x^3 + 5.179x^2 - 95.099x + 571.92
            FX1 = (-0.0923*(float(X1)**3))+(5.179*(float(X1)**2))-(95.099*(X1))+571.92
            FX2 = (-0.0923*(float(X2)**3))+(5.179*(float(X2)**2))-(95.099*(X2))+571.92

            #===Metode Interpolasi Linear===
            Xt = X2 - ((FX2/(FX2-FX1)))*(X2-X1)

            #===Memeriksa Persamaan Baru===
            FXt = (-0.0923*(float(Xt)**3))+(5.179*(float(Xt)**2))-(95.099*(Xt))+571.92

            #===Update Interval===
            if FXt * FX1 > 0:
                X2 = Xt
                FX2 = FXt
            else: 
                X1 = Xt
                FX1 = FXt

            #===Check Konvergensi
            if FXt < 0:
                error = FXt * (-1)
            else:
                error = FXt
            if iterasi > 500:
                print("angka tak hingga")
                break
            print(iterasi, "|", FX1, "|", FX2, "|", Xt, "|", FXt, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar Persamaan: ", Xt)
        print("Toleransi Error: ", error)

    def MetodeSecant (X1):
        X1 = X1
        X2 = X1 - 1
        error = 1
        iterasi = 0
        while(error > 0.0001):
            iterasi +=1
            #-0.0923x^3 + 5.179x^2 - 95.099x + 571.92
            FX1 = (-0.0923*(float(X1)**3))+(5.179*(float(X1)**2))-(95.099*(X1))+571.92
            FXmin = (-0.0923*(float(X2)**3))+(5.179*(float(X2)**2))-(95.099*(X2))+571.92

            #===Metode Secant untuk Akar Baru===
            X3 = X1 - ((FX1)*(X1-(X2)))/((FX1)-(FXmin))

            #===Memeriksa Persamaan Baru===
            FXplus = (-0.0923*(float(X3)**3))+(5.179*(float(X3)**2))-(95.099*(X3))+571.92

            #===Update Interval===
            if FXplus < 0:
                error = FXplus * (-1)
            else: 
                error = FXplus                       
            #===Check Konvergensi
            if error > 0.0001:
                X2 = X1
                X1 = X3
            else:
                print("Selesai")
            if iterasi > 500:
                print("angka tak hingga")
                break
            print(iterasi, "|", FX1, "|", FXmin, "|", X3, "|", FXplus, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar Persamaan: ", X3)
        print("Toleransi Error: ", error)

    def MetodeNewtonRapson (X1):
        X1 = X1
        iterasi = 0
        akar = 1

        while (akar > 0.0001):
            iterasi +=1
            #-0.0923x^3 + 5.179x^2 - 95.099x + 571.92
            Fxn = (-0.0923*(float(X1)**3))+(5.179*(float(X1)**2))-(95.099*X1)+571.92
            Fxxn = (-0.0923*(float(3*X1)**2))+(5.179*(float(2*X1)))-(95.099)

            #===Memeriksa Hasil Akar Persamaan===
            xnp1 = (X1-(Fxn/Fxxn))
            fxnp1 = (-0.0923*(xnp1**3))+(5.179*(xnp1**2))-(95.099*xnp1)+571.92

            #===Analisa Konvergensi===
            Ea = ((xnp1-X1)/xnp1)*100
            if Ea < 0.0001:
                X1 = xnp1
                akar = Ea*(-1)
            else:
                akar = xnp1
                print("Nilai akar adalah: ", akar)
                print("Nilai error adalah: ", Ea)
            if iterasi > 100:
                break
            print(iterasi, "|", X1, "|", xnp1, "|", akar, "|", Ea)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", xnp1)
        print("Toleransi Error: ", akar)

    def MetodeIterasi (X1):
        X1 = X1
        error = 1
        iterasi = 0
        while (error > 0.0001):
            iterasi += 1
            Fxn = (float(0.0923*X1)**3)-(float(5.179*X1)**2)+(95.099*X1)-571.92
            X2 = ((-X1**2)+(3*X1+3)**0.333334)
            Ea = ((X2-X1)/(X2))*100
            if Ea < error:
                X1 = X2
                if Ea > 0:
                    error = Ea
                else:
                    error = Ea*(-1)
            else:
                error = Ea
            if iterasi > 100:
                print("Angka tak hingga")
                break
            print(iterasi, "|", X1, "|", X2, "|", Ea, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", X2)
        print("Toleransi Error: ", error)


    print("code akar persamaan:\n",
         "1.Metode Setengah Interval\n",
         "2.Metode Interpolasi Linear\n",
         "3.Metode Secant\n",
         "4.Metode NewtonRapson\n",
         "5.Metode Iterasi\n")
    setting =int(input("masukan code akar persamaan:"))
    if (setting ==1):
        X1 = float(input("Masukkan Nilai Pertama: "))
        X2 = float(input("Masukkan Nilai Kedua: "))
        X = MetodeSetengahInterval(X1,X2)
        print (X)
    elif (setting ==2):
        X1 = float(input("Masukkan Nilai Pertama: "))
        X = MetodeInterpolasiLinear(X1)
        print (X)
    elif (setting ==3):
        X1 = float(input("Masukkan Nilai Pertama: "))
        X = MetodeSecant(X1)
        print (X)
    elif (setting ==4):
        X1 = float(input("Masukkan Nilai Pertama: "))
        X = MetodeNewtonRapson(X1)
        print (X)
    else:
        X1 = float(input("Masukkan Nilai Pertama: "))
        X = MetodeIterasi(X1)
        print (X)
        
def modul_3():
     # Metode Gauss
    def Gauss(A, f):
        A = np.array((A), dtype=float)
        f = np.array((f), dtype=float)
        n = len(f)
        for i in range(0, n - 1):  # Looping untuk kolom matriks
            if np.abs(A[i, i]) == 0:
                for k in range(i + 1, n):
                    if np.abs(A[k, i]) > np.abs(A[i, i]):
                        A[[i, k]] = A[[k, i]]  # Tukar antara baris i dan k
                        f[[i, k]] = f[[k, i]]
                        break
            for j in range(i + 1, n):  # Ulangi baris yang ada di bawah diagonal untuk setiap kolom
                m = A[j, i] / A[i, i]
                A[j, :] = A[j, :] - m * A[i, :]
                f[j] = f[j] - m * f[i]
        return A, f
    # Metode Gauss Jordan
    def GaussJordan(a,n):
        #Step1 ===> Looping untuk pengolahan metode Gauss Jordan
        print('==============Mulai Iterasi===============')
        for i in range(n):
            if a[i][i]==0:
                sys.exit('Dibagi dengan angka nol (proses tidak dapat dilanjutkan)')
            for j in range(n):
                if i !=j:
                    ratio=a[j][i]/a[i][i]
                    #print('posisi nol di:[',j,i,']', 'nilai rasio:',ratio)
                    for k in range(n+1):
                        a[j,k]=a[j][k]-ratio*a[i][k]
                    print(a)
                    print(f'============================================')
        # Step 2 ====> Membuat semua variabel(x,y,z,...)==1
        ax=np.zeros((n,n+1))
        for i in range(n):
            for j in range(n+1):
                ax[i,j]=a[i][j]/a[i][i]
        print('===================Akhir Iterasi===================')
        return ax
    # Metode Jacobi
    def Jacobi(A,b, N=10, x=None, info=True):
        print("Hasil perhitungan : ")
        if x is None:
            x = zeros(len(A[0]))
        D = diag(A)
        R = A - diagflat(D)
        for i in range(N):
            x = (b - dot(R,x))/D
            if info:
                pprint(x)
        return x
    # Metode Gauss Seidel
    def GaussSeidel(A,B,x,n):
        L= np.tril(A)
        U=A-L
        for i in range(n):
            x = np.dot(np.linalg.inv(L), B-np.dot(U, x))
            print(x)
        return x
    print("Kode untuk penyelesaian sistem persamaan linier dan matriks: \n",
         "1. Metode Gauss \n",
         "2. Metode Gauss Jordan \n",
         "3. Metode Jacobi \n", 
         "4. Metode Gauss Seidel \n")
    setting = int(input("Masukkan kode penyelesaian sistem persamaan linier dan matriks: "))
    if (setting == 1):
        A = np.array([[4, 6, 10, 1], [2, 6, 3, 8], [2, 6, 4, 2], [1, 2, 3, 4]])
        f = np.array([5.099, 6.099, 2.099, 12.099])
        print('A = \n%s and f = %s \n' % (A, f))
        Cel = Gauss(A, f)
        x = np.linalg.solve(A, f)
        print('Hasil perhitungan dengan metode Gauss adalah x = \n %s \n' % x)
    elif (setting == 2):
        m = np.array([[4, 6, 10, 1, 5.099],
                      [2, 6, 3, 8, 6.099],
                      [2, 6, 4, 2, 2.099],
                      [1, 2, 3, 4, 12.099]],dtype=float)
        n = 4
        print('Matriks Persamaan') # Menampilkan matrix awal
        print(m)
        m = GaussJordan(m,n) # Menampilkan Hasil
        print(f"""Hasil Pengolahan menggunakan metode Gauss Jordan didapatkan hasil sebagai berikut:
    {m} \n""")
    elif (setting == 3):
        A = np.array([[4, 6, 10, 1], [2, 6, 3, 8], [2, 6, 4, 2], [1, 2, 3, 4]])
        b = np.array([3.5, -0.54, -1, 4.8])
        guess = np.array([0,0,0,0])
        Solusi = Jacobi(A,b, N=5, x=guess)
        print("A :")
        pprint(A)
        print("b :")
        pprint(b)
        print("x :")
        pprint(Solusi)
    elif (setting == 4):
        A = np.array([[4, 6, 10, 1], [2, 6, 3, 8], [2, 6, 4, 2], [1, 2, 3, 4]]) # Untuk masuk ke script input, rumus A cara menulisnya : np.array([[2, -6, -7, 4], [6, -7, -7, -1], [-6, 4, 5, -2], [-3, 2, 4, 6]])
        B = np.array([5.099, 6.099, 2.099, 12.099]) # Untuk masuk ke script input, rumus B cara menulisnya : [hasil nilai D pers.1, hasil nilai D pers.2, hasil nilai D pers.3, hasil nilai D pers.4]
        x = np.array([0,0,0,0]) # Disarankan untuk memasukan nilai 0,0,0,0 untuk nilai awal x
        n = eval(input('Masukan Batas Iterasi:')) # Disarankan batas iterasi dimulai dari 1 sampai 15 agar tidak terlalu banyak
        Cel = GaussSeidel (A,B,x,n)
        print('Hasil Perhitungan Gauss-Seidel: \n', solve(A,B))
    else:
        print("Periksa Kembali Kode yang Diminta!")

def modul_4():
    #Metode Trapesium Banyak Pias
    import numpy as np
    import matplotlib.pyplot as plt

    def trapesium(f,a,b,N):
        x = np.linspace(a,b,N+1)
        y = f(x)
        y_right = y[1:]
        y_left = y[:-1]
        dx = (b-a)/N
        T = (dx/2) * np.sum(y_right + y_left)
        return T 

    f = lambda x : 4*x*3 + 4*x*2 + 10
    a = 0
    b = 10
    N = 10

    x = np.linspace(a,b,N+1)
    y = f(x)

    X = np.linspace(a,b+1,N)
    Y = f(X)

    import math

    def simpson1per3(x0,xn,n):
      f = lambda x : 4*x*3+4*x*2+10
      h = (xn - x0) / n 

      integral = f(x0) + f(xn)

      for i in range(1,n):
        k = x0 + 1*h

        if 1%2 == 0:
          integral = integral + 2 * f(k)
        else:
          integral = integral + 4 * f(k)

      integral = integral * h/3

      return integral


    print("Kode untuk Integrasi Numerik: \n",
         "1. Metode Trapesium Banyak Pias \n",
         "2. Metode Simpson 1/3 \n",)
    setting = int(input("Masukkan kode Integrasi Numerik: "))

    if (setting == 1):
        plt.plot(X,Y)
        for i in range(N):
            xs = [x[i],x[i],x[i+1],x[i+1]]
            ys = [0,f(x[i]),f(x[i+1]),0]
            plt.fill(xs,ys,'b',edgecolor='b',alpha=0.2)
        plt.title('Trapesium Banyak Pias, N = ()'.format(N))
        plt.savefig("C:/A TA Metnum/Image/Banyak Pias.png")
        L = trapesium(f,a,b,N)
        print(L)
    elif (setting == 2):
        x1 = float(input("batas bawah (x1: "))
        x2 = float(input("batas atas kudu nulis kelipatan 3.14: "))
        hasil = simpson1per3(x1, x2, 2)
        print("nilai integral metode Simpson 1/3", hasil )
    else:
        print("Periksa Kembali Kode yang Diminta!")

def modul_5():
    def euler():

        #PENDEKATAN PERSAMAAN DIFFERENSIAL BIASA DENGAN METODE EULER
        Nama = input("Masukkan Nama: ")
        NIM = (input("Masukkan NIM: "))
        Kelas = input("Masukkan Kelas: ")

        #=== Import Library ===#
        import numpy as np
        import matplotlib.pyplot as plt
        from IPython import get_ipython

        plt.style.use('seaborn-poster')
        ipy = get_ipython()
        if ipy is not None: 
            ipy.run_line_magic('matplotlib', 'inline')

        #==== Mendefinisikan parameter dan fungsi ===#
        h = float(input("Masukkan nilai h: ")) # Step size
        x0 = float(input("Masukkan nilai x awal: "))
        xn = float(input("Masukkan nilai x akhir: "))
        x = np.arange(x0, xn + h, h) # Numerical grid
        y0 = float(input("Masukkan nilai y awal: ")) # Initial Condition
        G = (2*x**3)+(9*x**2)+(4*x)+12
        f = lambda x, y: (2*x**3)+(9*x**2)+(4*x)+12 # Persamaan Differensial Biasa

        #=== Metode Euler Eksplisit ===#
        y = np.zeros(len(x))
        y[0] = y0
        for i in range(0, len(x) - 1):
            y[i + 1] = y[i] + h*f(x[i], y[i])

        #=== Menghitung Error ===#
        Galat = G-y
        print (Galat)

        #=== Menampilkan Grafik ===#
        Judul = ("\n Grafik Pendekatan Persamaan Differensial Biasa Dengan Metode Euler \n\n %s_%s_%s \n" % (Nama,NIM,Kelas))
        plt.figure(figsize = (12, 12))
        plt.plot(x, y, 'b--', label='Hasil Pendekatan') 
        plt.plot(x, G, '-g', label='Hasil Analitik')
        plt.title(Judul) # Judul plot
        plt.xlabel('x') # Label sumbu x
        plt.ylabel('F(x)') # Label sumbu y
        plt.grid() # Menampilkan grid
        plt.legend(loc='lower right')
        plt.savefig("C:/A TA Metnum/Image/euler.png")

    def Heun():
        #PENDEKATAN PERSAMAAN DIFFERENSIAL BIASA DENGAN METODE HEUN
        Nama = input("Masukkan Nama: ")
        NIM = (input("Masukkan NIM: "))
        Kelas = input("Masukkan Kelas: ")

        #=== Import Library ===#
        import numpy as np
        import matplotlib.pyplot as plt
        from IPython import get_ipython

        plt.style.use('seaborn-poster')
        ipy = get_ipython()
        if ipy is not None:
            ipy.run_line_magic('matplotlib', 'inline')

        #==== Mendefinisikan parameter dan fungsi ===#
        h = float(input("Masukkan nilai h: ")) # Step size
        x0 = float(input("Masukkan nilai x awal: "))
        xn = float(input("Masukkan nilai x akhir: "))
        x = np.arange(x0, xn + h, h) # Numerical grid
        y0 = float(input("Masukkan nilai y awal: ")) # Initial Condition
        G = (2*x**3)+(9*x**2)+(4*x)+12
        f = lambda x, y: (2*x**3)+(9*x**2)+(4*x)+12 # Persamaan Differensial Biasa

        #=== Metode Heun / Runge Kutta Orde 2 ===#
        y = np.zeros(len(x))
        y[0] = y0
        for i in range(0, len(x) - 1):
            k1 = f(x[i], y[i])
            k2 = f((x[i]+h), (y[i]+(h*k1)))
            y[i+1] = y[i]+(0.5*h*(k1+k2))

        #=== Menghitung Error ===#
        Galat = G-y
        print (Galat)

        #=== Menampilkan Grafik ===#
        Judul = ("\n Grafik Pendekatan Persamaan Differensial Biasa Dengan Metode Heun \n\n %s_%s_%s \n" % (Nama,NIM,Kelas))
        plt.figure(figsize = (12, 12))
        plt.plot(x, y, 'b--', label='Hasil Pendekatan') 
        plt.plot(x, G, '-g', label='Hasil Analitik')
        plt.title(Judul) # Judul plot
        plt.xlabel('x') # Label sumbu x
        plt.ylabel('F(x)') # Label sumbu y
        plt.grid() # Menampilkan grid
        plt.legend(loc='lower right')
        plt.savefig("C:/A TA Metnum/Image/heun.png")

    print("Kode untuk Persamaan Diferensial Biasa: \n",
         "1. Metode euler \n",
         "2. Metode Heun 1/3 \n",)
    setting = int(input("Masukkan kode Persamaan Diferensial Biasa: "))

    if (setting == 1):
        Hasil = euler()
    elif (setting == 2):
        Hasil = Heun()
    else:
        print("Periksa Kembali Kode yang Diminta!")
print("Kode untuk materi metode numerik: \n",
      "1. Modul 2 Akar-akar Persamaan \n",
      "2. Modul 3 Sistem Persamaan Linier dan Matriks \n",
      "3. Modul 4 Integrasi Numerik \n",
      "4. Modul 5 Persamaan Persamaan Differensial Biasa \n")
setting = int(input("Masukkan kode untuk materi metode numerik: "))
if (setting == 1):
    Cel = modul_2()
elif (setting == 2):
    Cel = modul_3()
elif (setting == 3):
    Cel = modul_4()
elif (setting == 4):
    Cel = modul_5()
else:
    print("Periksa Kembali Kode yang Diminta!")


# In[ ]:





# In[ ]:





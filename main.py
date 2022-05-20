import math
import time
import copy
import matplotlib.pyplot as plt
from tabulate import tabulate


class Matrix:
    def __init__(self, rows,columnes):
        self.rows = rows
        self.columnes = columnes
        self.value = [[0.0 for x in range(columnes)] for y in range(rows)]

    def __str__(self):
        my_word = ""
        for i in range(self.rows):
            for j in range(self.columnes):
                my_word = my_word + str(self.value[i][j]) + " "
            my_word = my_word + "\n"
        return my_word

    def __mul__(self, other):
        if self.columnes == other.rows:
            new_matrix = Matrix(self.rows, other.columnes)
            for i in range(new_matrix.rows):
                for j in range(new_matrix.columnes):
                    for k in range(self.columnes):
                        new_matrix.value[i][j] += self.value[i][k] * other.value[k][j]
            return new_matrix
        else:
            raise Exception("Matrixes have wrong size")

    def __add__(self, other):
        if self.rows == other.rows and self.columnes == other.columnes:
            new_matrix = Matrix(self.rows, self.columnes)
            for i in range(self.rows):
                for j in range(self.columnes):
                    new_matrix.value[i][j] = self.value[i][j] + other.value[i][j]
            return new_matrix
        else:
            raise Exception("Matrixes have diffrent size")

    def __sub__(self, other):
        if self.rows == other.rows and self.columnes == other.columnes:
            new_matrix = Matrix(self.rows, self.columnes)
            for i in range(self.rows):
                for j in range(self.columnes):
                    new_matrix.value[i][j] = self.value[i][j] - other.value[i][j]
            return new_matrix
        else:
            raise Exception("Matrixes have diffrent size")

def Jacobi_method(A, b):
    x_pre = Matrix(A.rows,1)
    x = 0
    while True:
        x_nxt = Matrix(A.rows, 1)
        for i in range(x_nxt.rows):
            sum_one = 0.0
            if i > 0:
                for j in range(0, i - 1 + 1):
                    sum_one += A.value[i][j] * x_pre.value[j][0]
            sum_two = 0.0
            for j in range(i + 1, A.rows):
                sum_one += A.value[i][j]*x_pre.value[j][0]
            x_nxt.value[i][0] =(b.value[i][0] - sum_one - sum_two)/(A.value[i][i])
        x_pre = x_nxt
        res = A*x_pre - b
        norma = math.sqrt(sum([(res.value[j][0])**2 for j in range(res.rows)]))
        if norma <= 10.0**(-9):
            break
        x += 1
        if x == 100000:
            print("nie zbiega sie")
            break
    return x_pre, x


def Gauss_Seidel_method(A, b):
    x_pre = Matrix(A.rows, 1)
    x = 0
    while True:
        x_nxt = Matrix(A.rows, 1)
        for i in range(x_nxt.rows):
            sum_one = 0.0
            if i > 0:
                for j in range(0, i - 1 + 1):
                    sum_one += A.value[i][j] * x_nxt.value[j][0]
            sum_two = 0.0
            for j in range(i + 1, A.rows):
                sum_one += A.value[i][j] * x_pre.value[j][0]
            x_nxt.value[i][0] = (b.value[i][0] - sum_one - sum_two) / (A.value[i][i])
        x_pre = x_nxt
        res = A * x_pre - b
        norma = math.sqrt(sum([(res.value[j][0]) ** 2 for j in range(res.rows)]))
        if norma <= 10.0 ** (-9):
            break
        x += 1
        if x == 100000:
            print("nie zbiega sie")
            break
    return x_pre, x

def prepre_A_matrix(N,a1,a2,a3):
    A = Matrix(N, N)
    for i in range(N):
        A.value[i][i] = a1
        if i - 1 > -1:
            A.value[i][i - 1] = a2
        if i + 1 < N:
            A.value[i][i + 1] = a2
        if i - 2 > -1:
            A.value[i][i - 2] = a3
        if i + 2 < N:
            A.value[i][i + 2] = a3
    return A

def prepre_b_vector(N):
    b = Matrix(A.rows, 1)
    for i in range(N):
        b.value[i][0] = math.sin((i + 1) * (4 + 1))
    return b

def LU_decomposition(A,b):
    U = copy.deepcopy(A)
    L = Matrix(A.rows,A.columnes)
    for i in range(A.rows):
        L.value[i][i] = 1
    for k in range(A.rows - 1):
        for j in range(k+1,A.rows):
            #if U.value[j][k] == 0: # mozliwosc optymalizacji i przyspieszenia o rzedy wielkosci
            #    continue
            L.value[j][k] = U.value[j][k]/U.value[k][k]
            for i in range(k, A.rows):
                U.value[j][i] = U.value[j][i] - L.value[j][k] * U.value[k][i]
    y = Matrix(A.rows, 1)
    for i in range(A.rows):
        suma = 0.0
        for j in range(i):
            suma += L.value[i][j] * y.value[j][0]
        y.value[i][0] = (b.value[i][0] - suma)/L.value[i][i]
    x = Matrix(A.rows, 1)
    for i in range(A.rows - 1,-1,-1):
        suma = 0.0
        for j in range(i+1,A.rows):
            suma += U.value[i][j] * x.value[j][0]
        x.value[i][0] = (y.value[i][0] - suma)/U.value[i][i]
    return x


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    N = 968
    a1 = 5 + 5
    a2 = -1
    a3 = -1
    A = prepre_A_matrix(N,a1,a2,a3)
    b = prepre_b_vector(N)
    tic = time.perf_counter()
    _, x =Jacobi_method(A,b)
    toc = time.perf_counter()
    print("Jacobi "+str(x)+" iteracji")
    print("Jacobi " + str(toc- tic) + " s")
    tic = time.perf_counter()
    _, x = Gauss_Seidel_method(A, b)
    toc = time.perf_counter()
    print("Gauss-Seidel " + str(x) + " iteracji")
    print("Gauss-Seidel " + str(toc - tic) + " s")
    N = 968
    a1 = 3
    a2 = -1
    a3 = -1
    A = prepre_A_matrix(N, a1, a2, a3)
    # nie zbiega sie
    #Jacobi_method(A, b)
    #Gauss_Seidel_method(A, b)
    x = LU_decomposition(A, b)
    res = A * x - b
    norma = math.sqrt(sum([(res.value[j][0]) ** 2 for j in range(res.rows)]))
    print("Norma z residuum dla LU faktoryzacji: " + str(norma))
    a1 = 5 + 5
    a2 = -1
    a3 = -1
    Jacobi_times = []
    Gauss_Seidel_times = []
    LU_decomposition_times = []
    Ns = [100,500,1000,2000,3000,4000]
    for N in Ns:
        print("rozpoczynam dla N = "+str(N))
        A = prepre_A_matrix(N, a1, a2, a3)
        b = prepre_b_vector(N)
        time_for_Jacobi = 0.0
        time_for_Gauss_Seidel = 0.0
        time_for_LU_decomposition = 0.0
        tic = time.perf_counter()
        Jacobi_method(A, b)
        toc = time.perf_counter()
        time_for_Jacobi = toc - tic
        print("dla Jacobi " + str(time_for_Jacobi) + " s")
        tic = time.perf_counter()
        Gauss_Seidel_method(A, b)
        toc = time.perf_counter()
        time_for_Gauss_Seidel = toc - tic
        print("dla Gauss-Seidel " + str(time_for_Gauss_Seidel) + " s")
        tic = time.perf_counter()
        LU_decomposition(A, b)
        toc = time.perf_counter()
        time_for_LU_decomposition = toc - tic
        print("dla faktoryzacji LU " + str(time_for_LU_decomposition) + " s")
        Jacobi_times.append(time_for_Jacobi)
        Gauss_Seidel_times.append(time_for_Gauss_Seidel)
        LU_decomposition_times.append(time_for_LU_decomposition)
    plt.figure()
    plt.plot(Ns, Jacobi_times , 'r-', Ns, Gauss_Seidel_times, 'b-',Ns, LU_decomposition_times, 'y-')
    plt.legend(['Jacobi', 'Gauss-Seidel', 'faktoryzacja LU'])
    plt.yscale('log')
    plt.title('Czas od N dla roznych metod')
    plt.xlabel('N')
    plt.ylabel('t [s]')
    plt.grid(True)
    plt.savefig('wykres.png')
    print("Dane: ")
    Jacobi_times.insert(0,"Jacobi: ")
    Gauss_Seidel_times.insert(0, "Gauss-Seidel: ")
    LU_decomposition_times.insert(0, "faktoryzacja LU: ")
    print(tabulate([["N: ",100,500,1000,2000,3000,4000],Jacobi_times,Gauss_Seidel_times,LU_decomposition_times],))










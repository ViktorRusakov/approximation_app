import time
import pickle
import sympy as sym
from api.approximation_tools import ControlSystem
# from calc_indeces import calculate_indeces
sym.init_printing()

x1, x2, x3, x4, t = sym.symbols('x_1 x_2 x_3 x_4 t')

# Пример 1

# a = sym.Matrix([0, x1, x2 + x1**2])
# b = sym.Matrix([1, 0, 0])

a = sym.Matrix([x1, x1 + x2**2, x2*x1])
b = sym.Matrix([1, 0, 0])


system = ControlSystem(a, b)
approx_series = system.calc_series(6)
system.calc_approx_system()
system.generate_pdf('test_1_new')

# Пример 2
a = sym.Matrix([0, -sym.Rational(1, 2)*x1**2 - 4*t*x1 - 3*t**2*x1, -x1**2 - 2*t*x1 - 3*t**2*x1])
b = sym.Matrix([-1, 0, 0])

system = ControlSystem(a, b)
approx_series = system.approximate(6)
system.calc_approx_system()
system.generate_pdf()

# Пример 3

a = sym.Matrix([0, -x1, -2*x2, 6*x1, -3*x3, -2*x4 + x1*x2, 0])
b = sym.Matrix([-1, -5, 0, x2, x1*x2, x3, sym.Rational(3, 2)*x4 + sym.Rational(1, 2)*x1*x2])

system = ControlSystem(a, b)
approx_series = system.approximate(6)
system.calc_approx_system()
system.generate_pdf('test_2')


# Пример 4

a = sym.Matrix([t + x2, x1**2, t*x1**2 + x3 + x1])
b = sym.Matrix([-1, -t**2, -(t**3 + x2)])

system = ControlSystem(a, b)
approx_series = system.calc_series(7)
system.calc_approx_system()
system.generate_pdf('test_4')

# start = time.time()
# calculate_indeces(15)
# end = time.time()
#file = open('coeffs.pickle', 'rb')
#res = pickle.load(file)
#for order, indeces in res.items():
#    res[order] = ['.'.join(str(index) for index in index_set) for index_set in indeces]
#file = open('coeffs.pickle', 'wb')
#pickle.dump(res, file)

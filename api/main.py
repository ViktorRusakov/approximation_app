# import time
# import pickle
import sympy as sym
from api.approximation_tools import ControlSystem
# from calc_indeces import calculate_indeces
sym.init_printing()

x1, x2, x3, x4, x5, x6, x7, x8, c6, c7, c8, t = sym.symbols(
    'x_{1} x_{2} x_{3} x_{4} x_{5} x_{6} x_{7} x_{8} c_{6} c_{7} c_{8} t')

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
approx_series = system.calc_series(7, fliess=True)
system.calc_approx_system(fliess=True)
system.generate_pdf('test_12345678', fliess=True)


a = sym.Matrix([0, 0, 0, 0, 0, 0, 0, 1])
b = sym.Matrix([x5, x3*x5, x4*x5, 1, c6 + x6, c7 + x7, c8 + x8,  0])
# b = sym.Matrix([1, x3, x4, x5, c6 + x6, c7 + x7, c8 + x8,  0])
system = ControlSystem(a, b)
res = system.calculate_lie_projections(True)
res = system.sort_lie_elements(True)

# start = time.time()
# calculate_indeces(15)
# end = time.time()
#file = open('coeffs.pickle', 'rb')
#res = pickle.load(file)
#for order, indeces in res.items():
#    res[order] = ['.'.join(str(index) for index in index_set) for index_set in indeces]
#file = open('coeffs.pickle', 'wb')
#pickle.dump(res, file)

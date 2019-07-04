import os
import subprocess
import pickle
import sympy as sym
import numpy as np
from functools import reduce
from itertools import groupby
from api.Lie_tools import unfold_lie_bracket
from api.shuffle_tools import calc_shuffle_lin_comb


class ControlSystem:
    
    t = sym.Symbol('t')

    def __init__(self, a, b):
        assert len(a) == len(b), 'Размерности векторов a и b не совпадают'
        self.a = a
        self.b = b
        self.dim = len(a)
        self.X = self.get_variables()
        self.E = sym.Matrix(self.X)
        self.text = ''
        self.series_text = ''
        self.basis_lie_elements = {}
        self.main_part_lie = {}
        self.right_ideal = {}
        self.projections = {}
        self.approx_system = sym.zeros(self.dim, 1)
        
    def get_variables(self):
        n = self.dim
        variables = sym.symbols('x_{{1:{}}}'.format(n + 1))
        return list(variables)

    def get_init_text(self):
        res = r'''
        \documentclass{article}
        \usepackage{amsmath}
        \usepackage{amsthm}
        \usepackage{amssymb}
        \usepackage[english]{babel}
        \textheight=235mm\textwidth=170mm \oddsidemargin=-5mm
        \topmargin=-20mm\hoffset=8mm%\voffset=mm
        \renewcommand{\baselinestretch}{1.3}
        \begin{document}
        \centerline{\Large\bf Initial system}
        \[
        \left\{
        \begin{aligned}'''
        a, b = self.a, self.b
        n = self.dim
        x_dots = sym.symbols(r'\dot{{x}}_1:{}'.format(n + 1))
        u = sym.Symbol('u', commutative=False)
        for i in range(n):
            res += r'''{} &= {}\\'''.format(sym.latex(x_dots[i]), sym.latex(a[i] + b[i]*u))
        res += r'''\end{aligned}
        \right.
        \]
        \\
        '''
        return res

    def get_lie_basis_text(self):
        lie_elems = self.basis_lie_elements
        text = r'''\centerline{\Large\bf Basis Lie elements and their coefficients}
        \begin{align*}
        '''
        j = 1
        for order, brackets_info in lie_elems.items():
            text += r'''\mathcal{{A}}_{{{}}}: '''.format(order)
            for bracket, v in brackets_info.items():
                if ']' not in bracket:
                    if '.' not in bracket:
                        lie_for_latex = r'''b_{{{}}}&=\xi_{{{}}} '''.format(j, bracket)
                    else:
                        bracket = bracket.split('.')
                        lie_for_latex = r'''b_{{{}}}&=[\xi_{{{}}}, \xi_{{{}}}] '''.format(j, bracket[0], bracket[1])
                else:
                    lie_split = bracket.split(']')
                    lie_for_latex = r'''b_{{{}}}&=[[\xi_{{{}}}, \xi_{{{}}}],'''.format(j, lie_split[0][0],
                                                                                       lie_split[0][-1])
                    for a in lie_split[1:]:
                        lie_for_latex += r'''\xi_{{{}}}],'''.format(a)
                lie_for_latex = lie_for_latex[:-1] + r' &' + r'''v(b_{{{}}})&={}\\'''.format(j, sym.latex(v.T))
                text += lie_for_latex
                j += 1
        text += r'''\end{align*}'''
        return text

    def get_lie_main_text(self):
        main_part = self.main_part_lie
        text = r'''\centerline{\Large\bf Lie elements for principal part}
        \begin{align*}
        '''
        j = 1
        for order, brackets_info in main_part.items():
            text += r'''\mathcal{{A}}_{{{}}}: '''.format(order)
            for bracket in brackets_info:
                if ']' not in bracket:
                    if '.' not in bracket:
                        lie_for_latex = r'''\ell_{{{}}}&=\xi_{{{}}} '''.format(j, bracket)
                    else:
                        bracket = bracket.split('.')
                        lie_for_latex = r'''\ell_{{{}}}&=[\xi_{{{}}}, \xi_{{{}}}] '''.format(j, bracket[0], bracket[1])
                else:
                    lie_split = bracket.split(']')
                    first_elem = lie_split[0].split('.')
                    lie_for_latex = r'''\ell_{{{}}}&=[[\xi_{{{}}}, \xi_{{{}}}],'''.format(j, first_elem[0],
                                                                                          first_elem[1])
                    for a in lie_split[1:]:
                        lie_for_latex += r'''\xi_{{{}}}],'''.format(a)
                lie_for_latex = lie_for_latex[:-1] + r'''\\'''
                text += lie_for_latex
                j += 1
        text += r'''\end{align*}'''
        return text

    def get_right_ideal_text(self):
        ideal = self.right_ideal
        text = r'''\centerline{\Large\bf Right ideal elements}
        \begin{align*}
        '''
        j = 1
        for order, brackets_info in ideal.items():
            if brackets_info:
                text += r'''\mathcal{{A}}_{{{}}}: '''.format(order)
                for bracket, v in brackets_info.items():
                    split_to_addends = bracket.split('|')
                    latex_text = r'''z_{{{}}}&='''.format(j)
                    value = r''''''
                    for addend in split_to_addends:
                        if 'x' in addend:
                            coef, element = addend.split('x')
                        else:
                            coef, element = '1', addend
                        if ']' not in element:
                            if '.' not in element:
                                element_text = r'''\xi_{{{}}} '''.format(element)
                            else:
                                element = element.split('.')
                                element_text = r'''[\xi_{{{}}}, \xi_{{{}}}] '''.format(element[0], element[1])
                        else:
                            lie_split = element.split(']')
                            first_elem = lie_split[0].split('.')
                            element_text = r'''[[\xi_{{{}}}, \xi_{{{}}}],'''.format(first_elem[0], first_elem[1])
                            for a in lie_split[1:]:
                                element_text += r'''\xi_{{{}}}],'''.format(a)
                        if not coef.startswith('-'):
                            value += '+' + sym.latex(sym.Rational(coef) * sym.Symbol(element_text[:-1]))
                        else:
                            value += sym.latex(sym.Rational(coef) * sym.Symbol(element_text[:-1]))
                    if value.startswith('+'):
                        latex_text += value[1:]
                    else:
                        latex_text += value
                    text += latex_text + r' &' + r'''v(z_{{{}}})&={}\\'''.format(j, sym.latex(v.T))
                    j += 1
        text += r'''\end{align*}'''
        return text

    def get_projections_text(self):
        projections = self.projections
        with open('api/lie_basis.pickle', 'rb') as lb:
            lie_basis = pickle.load(lb)
        text = r'''\centerline{\Large\bf Principal part}
        \begin{align*}
        '''
        for order, proj_info in projections.items():
            text += r'''\mathcal{{A}}_{{{}}}: '''.format(order)
            order_matrix = []
            order_lie_basis = []
            for lie_element, r in lie_basis[order].items():
                order_matrix.append(r['repr'])
                order_lie_basis.append(lie_element)
            order_matrix = np.concatenate(order_matrix).reshape((-1, len(order_matrix)), order='F')
            order_matrix = sym.Matrix(order_matrix)
            for elem in proj_info:
                latex_text = r'''\tilde\ell_{{{}}}&='''.format(elem['index'])
                algebra_repr = sym.Matrix(elem['vector_repr'])
                value = r''''''
                try:
                    lie_repr = order_matrix.LUsolve(algebra_repr)
                except ValueError:
                    split_to_addends = elem['algebra_repr'].split('|')
                    value = r''''''
                    for addend in split_to_addends:
                        if 'x' in addend:
                            coef, element = addend.split('x')
                        else:
                            coef, element = '1', addend
                        element_text = r'''\xi_{{{}}}'''.format(element.replace('.', ''))
                        if not coef.startswith('-'):
                            value += '+' + sym.latex(sym.Rational(coef) * sym.Symbol(element_text))
                        else:
                            value += sym.latex(sym.Rational(coef) * sym.Symbol(element_text))
                else:
                    for i in range(len(lie_repr)):
                        if lie_repr[i]:
                            coef = str(lie_repr[i])
                            basis_elem = order_lie_basis[i]
                            if ']' not in basis_elem:
                                if '.' not in basis_elem:
                                    element_text = r'''\xi_{{{}}} '''.format(basis_elem)
                                else:
                                    basis_elem = basis_elem.split('.')
                                    element_text = r'''[\xi_{{{}}}, \xi_{{{}}}] '''.format(basis_elem[0], basis_elem[1])
                            else:
                                lie_split = basis_elem.split(']')
                                first_elem = lie_split[0].split('.')
                                element_text = r'''[[\xi_{{{}}}, \xi_{{{}}}],'''.format(first_elem[0], first_elem[1])
                                for a in lie_split[1:]:
                                    element_text += r'''\xi_{{{}}}],'''.format(a)
                            if not coef.startswith('-'):
                                value += '+' + sym.latex(sym.Rational(coef) * sym.Symbol(element_text[:-1]))
                            else:
                                value += sym.latex(sym.Rational(coef) * sym.Symbol(element_text[:-1]))
                if value.startswith('+'):
                    latex_text += value[1:]
                else:
                    latex_text += value
                text += latex_text + r'''\\'''
        text += r'''\end{align*}'''
        return text

    def get_approx_system_text(self):
        b = self.approx_system
        text = r'''\centerline{\Large\bf Approximating system}
        \[
        \left\{
        \begin{aligned}
        '''
        n = self.dim
        x_dots = sym.symbols(r'\dot{{x}}_1:{}'.format(n + 1))
        u = sym.Symbol('u', commutative=False)
        for i in range(n):
            text += r'''{} &= {}\\'''.format(sym.latex(x_dots[i]), sym.latex(b[i] * u))
        text += r'''\end{aligned}
        \right.
        \]
        '''
        return text

    def generate_pdf(self, name='test'):
        max_order = max(self.main_part_lie)
        self.calc_series(max_order)
        init_text = self.get_init_text()
        series_text = self.series_text
        lie_basis_text = self.get_lie_basis_text()
        lie_main_part_text = self.get_lie_main_text()
        right_ideal_text = self.get_right_ideal_text()
        projections_text = self.get_projections_text()
        approx_system_text = self.get_approx_system_text()
        text = init_text + series_text + lie_basis_text + lie_main_part_text + right_ideal_text +\
               projections_text + approx_system_text
        text += r'''\end{document}'''
        # os.chdir(path=r'C:\Users\User\Desktop\Python\approximation_web_app\pdfs')
        with open('{}.tex'.format(name), 'w') as f:
            f.write(text)

        cmd = ['pdflatex', '-interaction', 'nonstopmode', '{}.tex'.format(name)]
        proc = subprocess.Popen(cmd)
        proc.communicate()

        retcode = proc.returncode
        if not retcode == 0:
            os.unlink('{}.pdf'.format(name))
            raise ValueError('Error {} executing command: {}'.format(retcode, ' '.join(cmd)))

        os.unlink('{}.tex'.format(name))
        os.unlink('{}.log'.format(name))

    def R_0(self, func):
        return func.diff(self.t) + func.jacobian(self.X) * self.a
    
    def R_1(self, func):
        return func.jacobian(self.X) * self.b
    
    def ad(self, power, func):
        if power == 0:
            return self.R_1(func)
        else:
            return self.R_0(self.ad(power - 1, func)) - self.ad(power - 1, self.R_0(func))

    def coefficient(self, indeces):
        """
        Функция, считающая определенный коэффициент ряда.
        На вход подается строка из индексов, разделитель между индексами - '.'
        (разделитель необходим так как индексы бывают двузначные)
        
        Пример: indeces = '1.2.3' => coefficient(indeces) = v_123
        """
        
        func = self.E
        indeces = list(map(int, indeces.split('.')))
        for index in reversed(indeces):
            func = self.ad(index, func)
        res = func.subs(list(zip(self.X + [self.t], [0] * (self.dim + 1))))
        zero_vect = sym.zeros(self.dim, 1)
        # Если ады занулились, то коэф перед ними можно не считать
        if res == zero_vect:
            return zero_vect
        else:
            k = len(indeces)
            if k == 1:
                factor = sym.Rational(-1, sym.factorial(indeces[0]))
            else:
                numerator = (-1)**k
                factorials = map(sym.factorial, indeces)
                denominator = reduce(lambda x, y: x*y, factorials)
                factor = sym.Rational(numerator, denominator)
            return factor * res

    def calculate_coefficients(self, max_order):
        """
        Подсчет коэффициентов ряда до заданного порядка
        включительно. На выходе получаем словарь, в котором ключ обозначает
        порядок, а значение ключа - соответствующие этому порядку элементы
        (также задаются с помощью словаря, ключ - индекс коэффициента,
        значение - посчитанное значение коэффициента)
        
        Пример: 
            max_order = 3 => calculate_coefficients(max_order) =
            {
                1: {'0': v_0},
                2: {'0.0': v_00, '1': v_1},
                3: {'2': v_2, '0.1': v_01, '1.0': v_10, '0.0.0': v_000}
            }
        """

        res = {}
        with open('api/moments_grading.pickle', 'rb') as f:
            indeces = pickle.load(f)
        for i in range(1, max_order + 1):
            order_res = {}
            for index_set in indeces[i]:
                order_res[index_set] = self.coefficient(index_set)
            res[i] = order_res
        return res

    def calc_series(self, max_order):
        """
            Подсчет ряда системы до заданного порядка.
        """
        
        res = sym.zeros(self.dim, 1)
        zero_vect = sym.zeros(self.dim, 1)
        text = r'''
        \centerline{\Large\bf Series}
        '''
        eq_text = r'''x^0 = '''
        coeffs = self.calculate_coefficients(max_order)
        multline = False
        i = 1
        for key, value in coeffs.items():
            for index_set, coef in coeffs[key].items():
                if coef != zero_vect:
                    res += coef * sym.Symbol(r'\xi_{{{}}}'.format(
                            index_set.replace('.', '')))
                    if i == 9:
                        eq_text += r'''\\ + {}{} +'''.format(sym.latex(coef), sym.latex(sym.Symbol(r'\xi_{{{}}}'.format(
                            index_set.replace('.', '')))))
                        i = 0
                        multline = True
                    else:
                        eq_text += r'''{}{} +'''.format(sym.latex(coef), sym.latex(sym.Symbol(r'\xi_{{{}}}'.format(
                                    index_set.replace('.', '')))))
                    i += 1
        if multline:
            eq_text = r'''\begin{multline*}\\''' + eq_text[:-1] + r'''\end{multline*}'''
        else:
            eq_text = r'''\[\\''' + eq_text[:-1] + r'''\]'''
        text += eq_text
        self.series_text = text + r'''\\'''
        return res

    def sort_lie_elements(self):
        """
            Находим первые n ЛНЗ Ли-элементов (из которых будет построена главная часть системы),
            остальные забрасываем в правый идеал.
        """

        with open('api/lie_basis.pickle', 'rb') as lb:
            lie_basis = pickle.load(lb)
        lie_coef = {}
        main_part = {}
        right_ideal = {}
        count = 0                           # счетчик кол-ва эл-ов в главной части
        independent_coeffs = sym.Matrix()   # матрица из ЛНЗ коэф-ов
        rank = 0                            # ранг матрицы 
        max_order = len(lie_basis)          # макс допустимый порядок
        zero_vect = sym.zeros(self.dim, 1)  # нулевой вектор, нужен для проверки
        for order in range(1, max_order + 1):

            # Если нашли n ЛНЗ Ли-элементов, то возвращаем результат
            if count == self.dim:
                self.basis_lie_elements = lie_coef
                self.main_part_lie = main_part
                self.right_ideal = right_ideal
                return main_part, right_ideal

            # Все Ли-элементы текущего порядка
            lie_elements = lie_basis[order].keys()
            lie_coef_order = {}
            
            # Независимые Ли-элементы текущего порядка
            independent = lie_basis[order].copy()
            
            # Правый идеал текущего порядка
            ideal = {}
            
            # Проходим по каждому элементу
            for bracket in lie_elements:
                
                # Коэффициент при элементе
                v = sym.zeros(self.dim, 1)

                # Раскрываем скобку и считаем коэф при ней
                unfolded = unfold_lie_bracket(bracket).split('|')
                for element in unfolded:
                    if element.startswith('-'):
                        v -= self.coefficient(element[1:])
                    else:
                        v += self.coefficient(element)
                lie_coef_order[bracket] = v
                        
                # Если нулевой, то сразу забрасываем в идеал
                if v == zero_vect:
                    ideal[bracket] = v
                    independent.pop(bracket)
                    continue

                # Добавляем коэф в матрицу всех ЛНЗ из предыдущих порядков
                # для проверки на ЛЗ, если ранг м-цы не поменялся - забрасываем
                # в идеал
                new_mat = sym.Matrix.hstack(independent_coeffs, v)
                new_rank = new_mat.rank()
                if new_rank == rank:
                    ideal[bracket] = independent[bracket]
                    independent.pop(bracket)
                else:
                    independent[bracket]['coef'] = v

            lie_coef[order] = lie_coef_order
            # Если нашлись ЛНЗ элементы, то нужно проверить, что их лин. комб.
            # тоже является ЛНЗ. Если существуют зависимые лин. комб. каких-то
            # элементов, то забрасываем эти комбинации в идеал.
            if independent:
                brackets, data = zip(*independent.items())
                rep = [p['coef'] for p in data]
                
                # Добавляем все найденные незав. коэф. в матрицу коэф-ов
                # предыдущих порядков
                new_coef_mat = np.concatenate(rep).reshape((-1, len(rep)), order='F')
                new_mat = sym.Matrix.hstack(independent_coeffs, sym.Matrix(new_coef_mat))
                new_rank = new_mat.rank()
                
                # Если ранг м-цы повысился на добавленное кол-во элементов,
                # то забрасываем все в главную часть, а в идеал - ничего,
                # если нет, то ищем какие лин. комб. пойдут в идеал
                if new_rank == rank + len(independent):
                    independent_coeffs = new_mat
                    rank = new_rank
                    main_part[order] = independent
                    count += len(independent)
                else:
                    kernel = new_mat.nullspace()
                    nonzero = list(range(len(independent)))  # для нахождения эл-ов из главной части
                    for basis_vector in kernel:
                        coefs = basis_vector[-len(independent):]
                        bad_elem = ''
                        value = sym.zeros(self.dim, 1)
                        for i, c in enumerate(coefs):
                            if c:
                                bad_elem += str(c) + 'x' + brackets[i] + '|'
                                value += c * rep[i]
                            else:
                                try:
                                    nonzero.remove(i)
                                except ValueError:
                                    pass
                        ideal[bad_elem[:-1]] = value
                    if len(kernel) == 1:
                        nonzero = nonzero[:-1]
                    main_part[order] = {
                            brackets[i]: data[i] for i in nonzero
                    }
                    for j in nonzero:
                        independent_coeffs = sym.Matrix.hstack(independent_coeffs, rep[j])
                    count += len(nonzero)
            right_ideal[order] = ideal
        self.basis_lie_elements = lie_coef
        self.main_part_lie = main_part
        self.right_ideal = right_ideal
        return main_part, right_ideal

    def calculate_lie_projections(self):
        """
        Подсчет проекций первых self.dim ЛНЗ Ли-элементов на ортогональное дополнение к правому идеалу
        """

        main_part, right_ideal = self.sort_lie_elements()
        projections = {}
        with open('api/moments_grading.pickle', 'rb') as f:
            indeces = pickle.load(f)
        main_orders, bad_orders = main_part.keys(), right_ideal.keys()
        i = 1

        # Проходим по всем порядкам элементов для главной части
        for order in main_orders:
            order_matrix = sym.Matrix()
            order_projections = []
            moments = indeces[order]
            # moments_repr = indeces[order]
            dim = len(moments)
            basis = np.eye(dim, dtype=int)
            moments_repr = {moments[i]: basis[:, [i]] for i in range(dim)}

            # По всем порядкам элементов правого идеала (дополняем правый идеал домножением элементов из идеала
            # на подалгебру нужного порядка)
            for bad_order in bad_orders:

                # Если дошли до текущего порядка элементов главной части, то домножать ни на что не надо
                if bad_order == order:
                    for lie_elem in right_ideal[bad_order]:
                        elem = lie_elem.split('|')
                        new_elem_repr = np.zeros((dim, 1), dtype=int)
                        for element in elem:
                            if 'x' not in element:
                                unfolded = unfold_lie_bracket(element).split('|')
                                for elem in unfolded:
                                    if elem.startswith('-'):
                                        new_elem_repr -= moments_repr[elem[1:]]
                                    else:
                                        new_elem_repr += moments_repr[elem]
                            else:
                                element = element.split('x')
                                unfolded = unfold_lie_bracket(element[1]).split('|')
                                for elem in unfolded:
                                    if elem.startswith('-'):
                                        new_elem_repr -= moments_repr[elem[1:]]
                                    else:
                                        new_elem_repr += moments_repr[elem]
                                new_elem_repr *= int(element[0])
                        order_matrix = order_matrix.row_join(sym.Matrix(new_elem_repr))
                    break
                else:
                    # Элементы подалгебры, которыми будем дополнять идеал
                    algebra_elems = indeces[order - bad_order]
                    for lie_elem in right_ideal[bad_order]:
                        unfolded = lie_elem.split('|')
                        for alg_elem in algebra_elems:
                            new_elem_repr = np.zeros((dim, 1))
                            for element in unfolded:
                                if 'x' not in element:
                                    new_unfolded = unfold_lie_bracket(element).split('|')
                                    for elem in new_unfolded:
                                        if elem.startswith('-'):
                                            new_elem_repr -= moments_repr[elem[1:] + '.' + alg_elem]
                                        else:
                                            new_elem_repr += moments_repr[elem + '.' + alg_elem]
                                else:
                                    element = element.split('x')
                                    new_unfolded = unfold_lie_bracket(element[1]).split('|')
                                    for elem in new_unfolded:
                                        if elem.startswith('-'):
                                            new_elem_repr -= moments_repr[elem[1:] + '.' + alg_elem]
                                        else:
                                            new_elem_repr += moments_repr[elem + '.' + alg_elem]
                                    new_elem_repr *= int(element[0])
                            order_matrix = order_matrix.row_join(sym.Matrix(new_elem_repr))

            # если идеал не пустой, то ищем проекции Ли-элементов на его ортогональное дополнение
            if order_matrix:
                transpose = order_matrix.T
                multi = transpose * order_matrix
                for a, v in main_part[order].items():
                    a_2 = v['repr'] - order_matrix * multi.LUsolve(transpose * v['repr'])
                    a2_algebra = ''
                    for j in range(len(a_2)):
                        if a_2[j]:
                            a2_algebra += str(a_2[j]) + 'x' + moments[j] + '|'
                    order_projections.append({
                        'index': i,
                        'vector_repr': a_2,
                        'algebra_repr': a2_algebra[:-1]
                    })
                    i += 1

            # если пустой, то проекциями являются сами Ли-элементы
            else:
                for a, v in main_part[order].items():
                    unfolded = unfold_lie_bracket(a)
                    order_projections.append({
                            'index': i,
                            'vector_repr': v['repr'],
                            'algebra_repr': unfolded
                    })
                    i += 1
            projections[order] = order_projections
        self.projections = projections
        return projections
    
    def calc_approx_system(self):
        """
        Подсчет аппроксимирующей системы.
        """
        projections = self.calculate_lie_projections()
        with open('api/moments_grading.pickle', 'rb') as f:
            indeces = pickle.load(f)
        orders = list(projections.keys())
        b_final = sym.zeros(self.dim, 1)
        for order in orders:
            order_projections = projections[order]
            for projection in order_projections:
                # Если проекция это однородный момент, то можно сразу записывать результат
                if '.' not in projection['algebra_repr']:
                    value = projection['algebra_repr'].split('x')
                    if len(value) == 2:
                        coef = -sym.Rational(value[0])
                    else:
                        coef = sym.Rational('-1')
                    b_final[projection['index'] - 1] = coef * self.t**int(value[-1])
                else:
                    # группировка моментов по одинаковому последнему индексу и сортировка по нему
                    grouped = {k: [m for m in g] for k, g in groupby(sorted(projection['algebra_repr'].split('|'),
                                                                            key=lambda x: x.split('.')[-1]),
                                                                     key=lambda x: x[-1])}
                    b = 0
                    for decomp_index, value in grouped.items():
                        cur_order = order - 1 - int(decomp_index)  # порядок элементов без последнего индекса
                        if cur_order == 0:
                            # если момент был однородным, то сразу пишем результат
                            split_value = value[0].split('x')
                            if len(split_value) == 2:
                                coef = -sym.Rational(split_value[0])
                            else:
                                coef = sym.Rational('-1')
                            b += coef * self.t**(order - 1)
                            continue
                        moments = indeces[cur_order]
                        # moments_repr = indeces[cur_order]
                        dim = len(moments)
                        basis = np.eye(dim, dtype=int)
                        moments_repr = {moments[i]: basis[:, [i]] for i in range(dim)}
                        order_matrix = sym.Matrix()    # матрица представлений шафл-элементов текущего порядка
                        order_shuffle_basis = []      # шафл-элементы текущего порядка

                        # Проходим по всем предыдущим порядкам, чтобы посчитать всевозможные шафл-произведения
                        # предыдущих элементов так, чтобы получить текущий порядок
                        for i, order_1 in enumerate(orders):
                            order_1_elems = projections[order_1]

                            # Если дошли до текущего порядка, то останавливаемся
                            if order_1 == cur_order:
                                for element in order_1_elems:
                                    order_shuffle_basis.append(str(element['index']))
                                    order_matrix = order_matrix.row_join(sym.Matrix(element['vector_repr']))
                                break

                            # Если попался порядок, который делит текущий, то значит элементы этого порядка можно
                            # повысить до текущего с помошью самоумножения (шафлом)
                            if cur_order % order_1 == 0:
                                count = cur_order // order_1
                                for elem in order_1_elems:
                                    index = elem['index']
                                    order_shuffle_basis.append(str(index) + (count - 1) * '*{}'.format(index))
                                    shuffle_res = calc_shuffle_lin_comb(elem['algebra_repr'], x_count=count)
                                    shuffle_res_repr = sym.zeros(dim, 1)
                                    shuffle_split = shuffle_res.split('|')
                                    for shuffle_elem in shuffle_split:
                                        if 'x' not in shuffle_elem:
                                            coef = sym.Rational('1')
                                            moment = shuffle_split
                                        else:
                                            split = shuffle_elem.split('x')
                                            coef, moment = sym.Rational(split[0]), split[1]
                                        shuffle_res_repr += coef * sym.Matrix(moments_repr[moment])
                                    order_matrix = order_matrix.row_join(shuffle_res_repr)

                            # Проходимся по порядком старше и проверяем можно ли составить шафл-комбинации
                            for order_2 in orders[i+1:]:
                                if order_2 + order_1 > cur_order:
                                    break
                                limit_1 = cur_order // order_1
                                limit_2 = cur_order // order_2
                                order_2_elems = projections[order_2]
                                for j in range(1, limit_1 + 1):
                                    for k in range(1, limit_2 + 1):
                                        if j * order_1 + k * order_2 == cur_order:
                                            for elem_1 in order_1_elems:
                                                index_1 = elem_1['index']
                                                addend_1 = str(index_1) + (j - 1) * '*{}'.format(index_1) + '*'
                                                for elem_2 in order_2_elems:
                                                    index_2 = elem_2['index']
                                                    addend_2 = str(index_2) + (k - 1) * '*{}'.format(index_2)
                                                    order_shuffle_basis.append(addend_1 + addend_2)
                                                    shuffle_res = calc_shuffle_lin_comb(elem_1['algebra_repr'],
                                                                                        elem_2['algebra_repr'],
                                                                                        x_count=j,
                                                                                        y_count=k)
                                                    shuffle_res_repr = sym.zeros(dim, 1)
                                                    shuffle_split = shuffle_res.split('|')
                                                    for shuffle_elem in shuffle_split:
                                                        if 'x' not in shuffle_elem:
                                                            coef = sym.Rational('1')
                                                            moment = shuffle_split
                                                        else:
                                                            split = shuffle_elem.split('x')
                                                            coef, moment = sym.Rational(split[0]), split[1]
                                                        shuffle_res_repr += coef * sym.Matrix(moments_repr[moment])
                                                    order_matrix = order_matrix.row_join(shuffle_res_repr)
                        algebra_repr = sym.zeros(dim, 1)
                        for elem in value:
                            element = '.'.join(elem.split('.')[:-1])
                            element = element.split('x')
                            if len(element) == 1:
                                if elem.startswith('-'):
                                    algebra_repr -= sym.Matrix(moments_repr[element[0][1:]])
                                else:
                                    algebra_repr += sym.Matrix(moments_repr[element[0]])
                            else:
                                coef = element[0]
                                algebra_repr += sym.Rational(coef) * sym.Matrix(moments_repr[element[1]])
                        shuffle_repr = order_matrix.LUsolve(algebra_repr)
                        for i in range(len(shuffle_repr)):
                            if shuffle_repr[i]:
                                basis_elem = order_shuffle_basis[i]
                                basis_elem = basis_elem.split('*')
                                elem_res = 1
                                for e in basis_elem:
                                    elem_res *= sym.Symbol('x{}'.format(e))
                                b += -shuffle_repr[i] * self.t**int(decomp_index) * elem_res
                    b_final[projection['index'] - 1] = b
        self.approx_system = b_final

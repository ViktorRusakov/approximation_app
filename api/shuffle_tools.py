import sympy as sym


def shuffle_product(x, y):
    """
    Вычисление шафл-произведения двух элементов.
    """
    if '.' not in x:
        if '.' not in y:
            return x + '.' + y + '|' + y + '.' + x
        else:
            y_split = y.split('.')
            res_1 = x + '.' + y
            res_2 = '|'.join(map(lambda z: y_split[0] + '.' + z, shuffle_product('.'.join(y_split[1:]), x).split('|')))
            return res_1 + '|' + res_2
    elif '.' not in y:
        x_split = x.split('.')
        res_1 = '|'.join(map(lambda z: x_split[0] + '.' + z, shuffle_product('.'.join(x_split[1:]), y).split('|')))
        res_2 = y + '.' + x
        return res_1 + '|' + res_2
    else:
        x_split = x.split('.')
        y_split = y.split('.')
        res_1 = '|'.join(map(lambda z: x_split[0] + '.' + z, shuffle_product('.'.join(x_split[1:]), y).split('|')))
        res_2 = '|'.join(map(lambda z: y_split[0] + '.' + z, shuffle_product('.'.join(y_split[1:]), x).split('|')))
        return res_1 + '|' + res_2


def calc_shuffle_lin_comb(x='', y='', x_count=0, y_count=0):
    """
        Шафл-произведение линейной комбинации моментов. x_count, y_count - количество раз соответствующий
        элемент попадет в произведение.
    """
    x_times = x_count - 1
    y_times = y_count - 1
    if x_times >= 0:
        x_split = x.split('|')
        res = [x]
        times = 1
        while times <= x_times:
            new_res = []
            for element in res:
                temp = element.split('|')
                for moment_1_repr in temp:
                    if 'x' not in moment_1_repr:
                        moment_1, coef_1 = moment_1_repr, sym.Rational('1')
                    else:
                        split_moment_1 = moment_1_repr.split('x')
                        moment_1, coef_1 = split_moment_1[1], sym.Rational(split_moment_1[0])
                    for moment_2_repr in x_split:
                        if 'x' not in moment_2_repr:
                            moment_2, coef_2 = moment_2_repr, sym.Rational('1')
                        else:
                            split_moment_2 = moment_2_repr.split('x')
                            moment_2, coef_2 = split_moment_2[1], sym.Rational(split_moment_2[0])
                        coef = str(coef_1 * coef_2)
                        shuffle_prod = coef + 'x' + shuffle_product(moment_1, moment_2)
                        new_res.append(shuffle_prod.replace('|', '|' + coef + 'x'))
            times += 1
            res = new_res.copy()
        res = '|'.join(res)
        res = res.split('|')
        times = 0
    else:
        res = [y]
        times = 1
    y_split = y.split('|')
    while times <= y_times:
        new_res = []
        for element in res:
            temp = element.split('|')
            for moment_1_repr in temp:
                if 'x' not in moment_1_repr:
                    moment_1, coef_1 = moment_1_repr, sym.Rational('1')
                else:
                    split_moment_1 = moment_1_repr.split('x')
                    moment_1, coef_1 = split_moment_1[1], sym.Rational(split_moment_1[0])
                for moment_2_repr in y_split:
                    if 'x' not in moment_2_repr:
                        moment_2, coef_2 = moment_2_repr, sym.Rational('1')
                    else:
                        split_moment_2 = moment_2_repr.split('x')
                        moment_2, coef_2 = split_moment_2[1], sym.Rational(split_moment_2[0])
                    coef = str(coef_1 * coef_2)
                    shuffle_prod = coef + 'x' + shuffle_product(moment_1, moment_2)
                    new_res.append(shuffle_prod.replace('|', '|' + coef + 'x'))
        times += 1
        res = new_res.copy()
    return '|'.join(res)

from django import forms
import sympy as sym
from sympy.parsing.latex import parse_latex


def latex_validator(value):
    valid_chars = ['f', 'r', 'a', '\\', 'x', '_', '{', '}', 't', '+', '-', '^', ' ', 's', 'i', 'n', 'c', 'o', 'l', 'g',
                   'e', 'x', 'p', '*', ')', '(']
    for s in value:
        if s not in valid_chars:
            try:
                test = int(s)
            except ValueError:
                raise forms.ValidationError('Invalid LaTeX expression.')
    return value


def latex_number_validator(value):
    valid_chars = [r'\frac', '\\', '{', '}', ' ']
    for s in value:
        if s not in valid_chars:
            try:
                test = int(s)
            except ValueError:
                raise forms.ValidationError('Invalid LaTeX expression. Must be only numbers or LaTeX fractions.')
    return value


def equilibrium_validator(value):
    parsed = parse_latex(value)
    variables = list(sym.symbols('x_{1:10}'))
    res = parsed.subs(list(zip(variables, [0] * 10)))
    if res != 0:
        raise forms.ValidationError('The point x=0 should be an equilibrium point of the function.')
    return value


class InputForm(forms.Form):

    def __init__(self, *args, **kwargs):
        self.n = kwargs.pop('n')
        self.algebra = kwargs.pop('algebra')
        super().__init__(*args, **kwargs)
        if self.algebra == 'fliess':
            a_validators = [latex_validator]
        else:
            a_validators = [latex_validator, equilibrium_validator]
        for i in range(1, self.n + 1):
            self.fields['a{}'.format(i)] = forms.CharField(max_length=100, widget=forms.TextInput,
                                                           validators=a_validators,
                                                           label='a<sub>{}</sub>(t, x)'.format(i))
            self.fields['b{}'.format(i)] = forms.CharField(max_length=100, widget=forms.TextInput,
                                                           validators=[latex_validator],
                                                           label='b<sub>{}</sub>(t, x)'.format(i))
            if self.algebra == 'fliess':
                self.fields['x{}'.format(i)] = forms.CharField(max_length=100, widget=forms.TextInput,
                                                               validators=[latex_validator],
                                                               label="x<sub>{}</sub><sup style='position: relative; left: -.5em; top: -.6em'>0</sup>".format(i),
                                                               initial='0', required=True)


class DimensionForm(forms.Form):
    ALGEBRA_TYPES = (
        ('fliess', 'Fliess algebra'),
        ('nonlinear', 'Algebra of nonlinear moments')
    )
    n = forms.IntegerField(max_value=10, min_value=1, help_text="Select dimension of the system",
                           widget=forms.NumberInput, label="Dimension", initial=3)
    algebra = forms.ChoiceField(choices=ALGEBRA_TYPES, widget=forms.Select, help_text='Choose algebra',
                                initial='fliess')

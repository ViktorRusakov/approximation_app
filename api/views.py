import sympy as sym
from django.shortcuts import render, redirect
from .forms import InputForm, DimensionForm
from api.approximation_tools import ControlSystem
from sympy.parsing.latex import parse_latex
from django.http import FileResponse
from django.contrib import messages
from .Lie_tools import LieElementsNotFound

# Create your views here.


def index(request):
    if request.method == 'POST':
        form = DimensionForm(request.POST)
        if form.is_valid():
            n = form.cleaned_data['n']
            request.session['dimension'] = n
            return redirect('approximate')
    else:
        form = DimensionForm()
    return render(request, 'index.html', {'form': form})


def approximate(request):
    if request.method == 'POST':
        dimension = len(request.POST) // 2
        form = InputForm(dimension, request.POST)
        if form.is_valid():
            clean = form.cleaned_data
            n = len(clean) // 2
            a, b = [], []
            for i in range(n):
                a.append(parse_latex(clean['a{}'.format(i + 1)]))
                b.append(parse_latex(clean['b{}'.format(i + 1)]))
            a = sym.Matrix(a)
            b = sym.Matrix(b)
            system = ControlSystem(a, b)
            try:
                system.calc_approx_system()
            except LieElementsNotFound:
                messages.error(request, 'Could not find necessary amount of linearly independent Lie elements '
                                        'among currently calculated Lie algebra grading (up to 10th order).'
                                        ' Try another system.')
                return render(request, 'approximate.html', {'form': form, 'dimension': dimension})
            else:
                messages.success(request, 'The system has been successfully approximated!')
                system.generate_pdf('test_app')
                return FileResponse(open(r'test_app.pdf', 'rb'), as_attachment=True, filename='test.pdf')
    else:
        dimension = request.session['dimension']
        form = InputForm(dimension)
    return render(request, 'approximate.html', {'form': form, 'dimension': dimension})

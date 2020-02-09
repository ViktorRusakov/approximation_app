import sympy as sym
import os
from api.models import Report
from django.shortcuts import render, redirect
from .forms import InputForm, DimensionForm
from api.approximation_tools import ControlSystem
from sympy.parsing.latex import parse_latex
from django.http import FileResponse
from django.contrib import messages
from .Lie_tools import LieElementsNotFound
from approximation_web_app.settings import BASE_DIR


def index(request):
    if request.method == 'POST':
        form = DimensionForm(request.POST)
        if form.is_valid():
            n = form.cleaned_data['n']
            algebra = form.cleaned_data['algebra']
            request.session['dimension'] = n
            request.session['algebra'] = algebra
            return redirect('approximate')
    else:
        form = DimensionForm()
    return render(request, 'index.html', {'form': form})


def approximate(request):
    if request.method == 'POST':
        dimension = request.session['dimension']
        algebra = request.session['algebra']
        form = InputForm(request.POST, **{'algebra': algebra, 'n': dimension})
        if form.is_valid():
            clean = form.cleaned_data
            fliess = algebra == 'fliess'
            a, b, init_point = [], [], []
            for i in range(dimension):
                a.append(parse_latex(clean['a{}'.format(i + 1)]))
                b.append(parse_latex(clean['b{}'.format(i + 1)]))
                if fliess:
                    init_point.append(parse_latex(clean['x{}'.format(i + 1)]))
            a = sym.Matrix(a)
            b = sym.Matrix(b)
            init_point = sym.Matrix(init_point)
            system = ControlSystem(a, b, init_point)
            try:
                system.calc_approx_system(fliess)
            except LieElementsNotFound:
                messages.error(request, 'Could not find necessary amount of linearly independent Lie elements '
                                        'among currently calculated Lie algebra grading (up to 10th order).'
                                        ' Try another system.')
                return render(request, 'approximate.html', {
                    'form': form,
                    'dimension': dimension,
                    'algebra': algebra})
            else:
                result = Report.objects.create()
                system.generate_pdf('report_{}'.format(result.id), fliess)
                messages.success(request, 'The system has been successfully approximated! You can now download'
                                          ' the resulting file below.')
                return render(request, 'success.html', {'report_id': result.id})
        else:
            return render(request, 'approximate.html', {
                'form': form,
                'dimension': dimension,
                'algebra': algebra})
    else:
        dimension = request.session['dimension']
        algebra = request.session['algebra']
        form = InputForm(**{'algebra': algebra, 'n': dimension})
        return render(request, 'approximate.html', {'form': form, 'dimension': dimension, 'algebra': algebra})


def download(request, file_id):
    file = open(os.path.join(BASE_DIR, 'pdfs/report_{}.pdf').format(file_id), 'rb')
    return FileResponse(file, filename='report.pdf')

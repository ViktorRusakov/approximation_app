from django import forms
from django.core import validators


def latex_validator(value):
    valid_chars = [r'\frac', '\\', 'x', '_', '{', '}', 't', '+', '-', '^', ' ']
    for s in value:
        if s not in valid_chars:
            try:
                test = int(s)
            except ValueError:
                raise forms.ValidationError('Invalid LaTeX expression.')
    return value


# class LatexValidator(object):
#     def __init__(self, dimension):
#         self.n = dimension
#
#     def __call__(self, value):
#         valid_chars = [r'\frac', '\\', 'x', '_', '{', '}', 't', '+', '-', '^', ' ']
#         for i, s in enumerate(value):
#             if s not in valid_chars:
#                 try:
#                     test = int(s)
#                     if test < 0 or test > self.n and value[i-2] == '_':
#                         raise forms.ValidationError('Invalid index of a variable')
#                 except ValueError:
#                     raise forms.ValidationError('Invalid LaTeX expression')
#         return value


class InputForm(forms.Form):

    def __init__(self, n,  *args, **kwargs):
        super().__init__(*args, **kwargs)
        # fields = self.fields
        # for i in range(len(fields) + 1):
        #     field_name = 'interest_%s' % (i,)
        #     self.fields[field_name] = forms.CharField(required=False)
        #     try:
        #         self.initial[field_name] = interests[i].interest
        #     except IndexError:
        #         self.initial[field_name] = “”
        #         # create an extra blank field
        #         field_name = 'interest_%s' % (i + 1,)
        #         self.fields[field_name] = forms.CharField(required=False)
        for i in range(1, n + 1):
            self.fields['a{}'.format(i)] = forms.CharField(max_length=100, widget=forms.TextInput,
                                                           validators=[latex_validator])
            self.fields['b{}'.format(i)] = forms.CharField(max_length=100, widget=forms.TextInput,
                                                           validators=[latex_validator])

    # def clean(self):
    #     cleaned_data = super().clean()


class DimensionForm(forms.Form):
    n = forms.IntegerField(max_value=10, min_value=1, help_text="Select dimension of the system",
                           widget=forms.NumberInput, label="Dimension", initial=3)

from django.urls import path

from . import views

urlpatterns = [
    path('approximate', views.approximate, name='approximate'),
    path('', views.index, name='index'),
]
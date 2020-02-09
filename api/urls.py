from django.urls import path

from . import views

urlpatterns = [
    path('approximate', views.approximate, name='approximate'),
    path('download/<int:file_id>', views.download, name='download'),
    path('', views.index, name='index'),
]
from django.db import models


class Report(models.Model):
    report = models.FileField(default=None, blank=True, null=True)

from .base import *

import os

DEBUG = True

ROOT_URLCONF = 'NGPhylogeny_fr.urls.local'

# INSTALLED_APPS += (
#     'debug_toolbar',
# )

# MIDDLEWARE += (
#     'debug_toolbar.middleware.DebugToolbarMiddleware',
# )

if os.environ.get('NGPHYLO_DATABASE_HOST') is not None:
    DATABASES = {
        'default': {
            'ENGINE': os.environ.get('NGPHYLO_DATABASE_ENGINE'),
            'NAME': os.environ.get('NGPHYLO_DATABASE_NAME'),
            'USER': os.environ.get('NGPHYLO_DATABASE_USER'),
            'PASSWORD': os.environ.get('NGPHYLO_DATABASE_PASSWORD'),
            'HOST': os.environ.get('NGPHYLO_DATABASE_HOST'),
            'PORT': os.environ.get('NGPHYLO_DATABASE_PORT'),
        }
    }
else:
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
        }
    }


from .base import *

DEBUG = True

ROOT_URLCONF = 'NGPhylogeny_fr.urls.local'

# INSTALLED_APPS += (
#     'debug_toolbar',
# )

# MIDDLEWARE += (
#     'debug_toolbar.middleware.DebugToolbarMiddleware',
# )

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}

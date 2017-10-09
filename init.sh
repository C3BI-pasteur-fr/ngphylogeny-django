#!/bin/bash
pip install -r requirement.txt
python manage.py makemigrations
python manage.py migrate
python manage.py loaddata tool_flags
python manage.py importtools
python manage.py compute_tools_links
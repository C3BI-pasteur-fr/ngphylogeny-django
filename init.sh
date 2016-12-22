#!/bin/bash
pip install -r requirement
python manage.py loaddata tool_flags
python manage.py importtools
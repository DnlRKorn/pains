#!/bin/bash
source /home/dkorn/miniconda2/bin/activate
#uwsgi --http-socket 127.0.0.1:8001 --wsgi-file main.py --callable app
uwsgi --ini pains_uwsgi.ini 

#!/usr/bin/python
#coding=utf-8
"""
Python relative imports are strange creatures.
In order to keep a relatively sane package structure,
we place SiGC inside the SiC³ module,
and then run the server from here. If you don't do this,
Python will barf with a ValueError. Beware Python relative imports.
"""
from sigc.server import app

#NOTE: Update the below to use host=0.0.0.0, port=80 and debug=False for production deployments
app.run(host="0.0.0.0",debug=True)

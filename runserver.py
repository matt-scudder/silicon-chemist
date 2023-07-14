#!/usr/bin/env python3
"""
Python imports are strange creatures.
This is in the root of the repository so that this file can be executed from any working directory
This file is not part of a module, but it's location here allows 'sic' to exist as the top level package. 
"""

from sic.sigc.server import app

#NOTE: Update the below to use host=0.0.0.0, port=80 and debug=False for production deployments
app.run(host="0.0.0.0",debug=True)

# -*- coding: utf-8 -*-
#
# sphinx configuration file
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import os
import subprocess

# -- General configuration ------------------------------------------------
extensions = ['sphinx.ext.mathjax']
master_doc = 'index'
source_suffix = '.rst'
rst_prolog = '''
.. |s| replace:: *shark*
.. |ss| replace:: ``shark-submit``
'''

# General information about the project.
project = u'shark'
author = u'Claudia Lagos, Rodrigo Tobar'
copyright = u"""UWA (in the framework of ICRAR).
ICRAR - International Centre for Radio Astronomy Research.
UWA - The University of Western Australia, 2018."""
with open('../VERSION') as f:
    version = f.read()

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

html_theme = 'sphinx_rtd_theme'
latex_documents = [
    (master_doc, 'shark.tex', u'shark Documentation',
     author, 'manual'),
]

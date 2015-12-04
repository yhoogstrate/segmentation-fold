#!/usr/bin/env python3

import sys
import os
import re
import shutil
from subprocess import call

def get_git_revision():
    with open("../.git/HEAD","r") as fh:
        target_file = fh.read().strip()
        if target_file[0:5] == "ref: ":
            with open("../.git/"+target_file[5:],"r") as fh2:
                return fh2.read().strip()
        else:
            return "?"

def get_code_version_from_cmake():
    with open("../CMakeLists.txt","r") as fh:
        content = fh.read()
    return re.search("set\(PROJECT_VERSION[ ]+['\"]([^'\"]+)['\"]\)",content).group(1)

def get_project_name_from_cmake():
    with open("../CMakeLists.txt","r") as fh:
        content = fh.read()
    return re.search("project\(([^\)]+)\)",content).group(1)



# General information about the project.
project = get_project_name_from_cmake()
copyright = ''

version = get_code_version_from_cmake()
release = version+"-"+get_git_revision()


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'breathe',
]

# Breathe extension variables
breathe_projects = { "segmentation-fold": "../doc/xml/" }
breathe_default_project = "segmentation-fold"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'


# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = 'cpp'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Output file base name for HTML help builder.
htmlhelp_basename = 'test'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    ('index', 'index.tex', 'segmentation-fold', 'me','manual'),
]

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'index', 'Documentation', ['me'], 1)
]


# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'segmentation-fold', 'segmentation-fold',
   'me', 'segmentation-fold', 'descr.',
   'Miscellaneous'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
#man_pages = [
#    ('index', 'test', 'test',
#     [''], 1)
#]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)


# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#texinfo_no_detailmenu = False


on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    # Do stuff that cmake is supposed to do...
    fh_out = open('../doc/Doxyfile', "w")
    with open('../doc/Doxyfile.in', "r") as fh:
        for line in fh:
            line = line.replace("@CMAKE_PROJECT_NAME@",project)
            line = line.replace("@PROJECT_VERSION@",version)
            
            
            line = line.replace("CLASS_DIAGRAMS         = YES","CLASS_DIAGRAMS         = NO")
            line = line.replace("GENERATE_XML           = NO", "GENERATE_XML           = YES")
            line = line.replace("GENERATE_HTML          = YES","GENERATE_HTML          = NO")
            line = line.replace("GENERATE_LATEX         = YES","GENERATE_LATEX         = NO")
            
            line = line.replace("EXTRACT_ALL            = YES","EXTRACT_ALL            = NO")
            line = line.replace("EXTRACT_PRIVATE        = YES","EXTRACT_PRIVATE        = NO")
            #line = line.replace("                         include \\","                         include")
            #line = line.replace("                         test", "")
            line = line.replace("HAVE_DOT               = YES","HAVE_DOT               = NO")
            
            fh_out.write(line)
    fh_out.close()
    
    call('cd .. ; nice doxygen doc/Doxyfile',shell=True)
    
#    print "\n\n\n---\n\n\n"
#    call('ls -als ..',shell=True)
#    print "\n\n\n---\n\n\n"
#    call('ls -als ../doc',shell=True)
#    print "\n\n\n---\n\n\n"
#    call('ls -als ../doc/latex',shell=True)
#    print "\n\n\n---\n\n\n"

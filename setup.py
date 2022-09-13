long_description = """\
Code accompanying the paper "Counting essential surfaces in
3-manifolds" by Nathan M. Dunfield, Stavros Garoufalidis, and Hyam
Rubinstein.
"""

from setuptools import setup

setup(
    name = 'nscomplex',
    version = '1.0.1',
    description = 'Counting essential surfaces for fun and profit',
    long_description = long_description,
    url = 'http://t3m.computop.org',
    author = 'Nathan M. Dunfield, Stavros Garoufalidis, and Hyam Rubinstein',
    author_email = 'nathan@dunfield.info',
    license='GPLv2+',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],
    packages = ['nscomplex'],
    package_dir = {'nscomplex':'nscomplex'}, 
)


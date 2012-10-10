#!/usr/bin/env python
try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("error: install setuptools")

setup(
    name='blast2cap3',
    version=0.91,
    author='Vince Buffalo',
    author_email='vsbuffalo@gmail.com',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points = {
        'console_scripts': [
            'blast2cap3 = blast2cap3.blast2cap3:main']
        },
    url='http://github.com/DubcovskyLab/blast2cap3',
    license='GPL 2.0',
    description="""Try to join transcriptome sequencing contigs by leveraging BLASTX
queries in protein space, and joining in using CAP3 in nucleotide
space.""",
    requires=["BioPython"]
    )   

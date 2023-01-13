from setuptools import setup
from setuptools import Extension

setup(
	name='multiplacement',
	version='1.0.1',
    author='Quinn Mood',
    author_email='qmood1@umbc.edu',
    description='C extension to expedite calculations for complex transcription factor binding',
	ext_modules=[
        Extension("_multiplacement",
        ['src/_interface.c', 'src/_multiplacement.c', 'src/_aux.c'],
        include_dirs = ['src'],
        )
    ],
)

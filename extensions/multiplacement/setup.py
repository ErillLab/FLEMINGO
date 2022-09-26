from setuptools import setup
from setuptools import Extension

setup(
	name='multiplacement',
	version='1.0.0',
    author='Quinn Mood',
    author_email='qmood1@umbc.edu',
    description='C extension to expedite calculations for complex transcription factor binding',
	ext_modules=[Extension("_multiplacement",['src/_multiplacement.c'])],
)

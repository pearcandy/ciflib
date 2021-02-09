''' 
               setup.py

          author : PearCandy
          date   : 2021/1/26               
                                           '''

#coding:utf-8

from glob import glob
from os.path import basename
from os.path import splitext

from setuptools import setup
from setuptools import find_packages


def _requires_from_file(filename):
    return open(filename).read().splitlines()


setup(
    name="ciflib",
    version="0.1",
    license="GNU lv3",
    description="Tools for quantum chemistry with quantum computer",
    long_description="README.md",
    long_description_content_type="text/markdown",
    author="PearCandy",
    url="https://github.com/pearcandy",
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    install_requires=_requires_from_file('requirements.txt'),
    setup_requires=["pytest-runner"],
    tests_require=["pytest", "pytest-cov"],
    python_requires='>=3.7',
)

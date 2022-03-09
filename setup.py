from setuptools import setup, Extension
import setuptools_scm  # noqa  Ensure it’s installed

setup(
    ext_modules=[
        Extension("dnaio._core", sources=["src/dnaio/_core.pyx"],
                  extra_compile_args=["-msse4.1"]),
    ],
)

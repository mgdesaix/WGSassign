from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [Extension(
				"WGSassign.reader_cy",
				["WGSassign/reader_cy.pyx"],
				extra_compile_args=['-Xpreprocessor', '-fopenmp', '-g0'],
				extra_link_args=['-Xpreprocessor', '-fopenmp'],
				include_dirs=[numpy.get_include()],
				language="c++"
			),
			Extension(
				"WGSassign.emMAF_cy",
				["WGSassign/emMAF_cy.pyx"],
				extra_compile_args=['-Xpreprocessor', '-fopenmp', '-g0'],
				extra_link_args=['-Xpreprocessor', '-fopenmp'],
				include_dirs=[numpy.get_include()]
			),
			Extension(
				"WGSassign.glassy_cy",
				["WGSassign/glassy_cy.pyx"],
				extra_compile_args=['-Xpreprocessor', '-fopenmp', '-g0'],
				extra_link_args=['-Xpreprocessor', '-fopenmp'],
				include_dirs=[numpy.get_include()]
			)]

setup(
	name="WGSassign",
	version="0.00",
	author="Matt DeSaix",
	description="Population assignment methods for whole-genome sequence and genotype likelihood data",
	packages=["WGSassign"],
	entry_points={
		"console_scripts": ["WGSassign=WGSassign.WGSassign:main"]
	},
	python_requires=">=3.6",
	install_requires=[
		'numpy',
		'cython',
		'scipy'
    ],
    ext_modules=cythonize(extensions, compiler_directives={'language_level':'3'}),
    include_dirs=[numpy.get_include()]
)

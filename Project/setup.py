from distutils.core import setup, Extension

module1 = Extension(
	'methods', # module name in interpreter
	sources = ['test.c']
)

setup(
	name = 'methods',
	version = '1.1',
	description = 'Simple module',
	ext_modules= [module1]
)

try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

config = {
	'description' : 'Visual Statistics',
	'author' : 'Ryan L. Melvin',
	'url' : 'http://salsburygroup.squarespace.com/',
	'download_url' : 'http://salsburygroup.squarespace.com/',
	'author_email' : 'melvrl13@wfu.edu',
	'version' : '0.1',
	'install_requirements': ['nose'],
	'packages' : ['vstats'],
	'scripts' : [],
	'name' : 'projectname'
}

setup(**config)

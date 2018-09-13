from setuptools import setup
setup(name                 ='SIAC',
      version              ='2.0',
      description          = 'A sensor invariant Atmospheric Correction (SIAC)',
      author               = 'Feng Yin',
      author_email         = 'ucfafyi@ucl.ac.uk',
      classifiers          = ['Development Status :: 2 - BETA',
                              'Programming Language :: Python :: 2.7, 3.6',
                              'Topic :: Atmospheric Correction'],
      install_requires     = ['gdal>=2.1', 'numpy>=1.13', 'scipy>=1.0', 'lightgbm>=2.1.0',\
                              'requests', 'scikit-learn', 'scikit-image', 'pandas', 'psutil'],
      url                  = 'https://github.com/MarcYin/SIAC',
      license              = "GNU Affero General Public License v3.0",
      include_package_data = True,
      packages             = ['SIAC'],
     )

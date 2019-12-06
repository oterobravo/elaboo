from setuptools import setup

setup(name='elaboorate',
      version='0.1',
      description='Branch Extract-&-Leave-All-But-One-Out Reconstruction',
      long_description='See README',
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.0',
      ],
      author='Alejandro Otero Bravo',
      author_email='alejandro.otero.b@gmail.com',
      license='MIT',
      packages=['elaboorate'],
      install_requires=[
          'logging','markdown',
      ],
      zip_safe=False)

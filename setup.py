from setuptools import setup, find_packages

setup(
    name='scRNA-seq-Integrated-Cancer-Attractor-Analysis',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'scipy',
        'scikit-learn',
        'networkx',
        'sympy',
        'tqdm',
        'joblib',
        'openpyxl',
        'pytest',
        'pytest-cov',
        'coveralls'
    ],
)

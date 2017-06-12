from setuptools import setup, find_packages

__version__ = "0.0.0"
exec(open('brainmap/_version.py').read())

setup(
    name="brainmap",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'numpy',
        'allensdk',
        'matplotlib',
        'colormap',
        'ipywidgets',
        'scikit-image'
    ],
    # scripts=[],
    author="Gioele La Manno",
    author_email="gioelelamanno@gmail.com",
    description="AllenBrainAtlas utils module",
    license="MIT",
    url="https://github.com/linnarsson-lab/brainmap",
)

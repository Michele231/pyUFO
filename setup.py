from setuptools import find_packages, setup

setup(
    name='pyUFO',
    packages=find_packages(include=['pyUFO','pyUFO.optic']),
    version='0.1',
    setup_requires=['setuptools','numpy'],
    install_requires=['numpy'],
    description='This library is used to compute the FOV lat-lon coordinate given the satellite view geometry.',
    author='Michele Martinazzo',
    author_email='michele.martinazzo2@unibo.it',
    license='MIT',
)

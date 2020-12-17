from setuptools import setup, find_packages

requirements = [
    "qutip",
    "qiskit>=0.13.0",
    "picos>=2.0.0"
]

MAJOR = 1
MINOR = 5
MICRO = 6
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

setup(
    name="QuAwesome",
    version=VERSION,
    packages=find_packages(),

    author="Yi-Te Huang",
    author_email="y.t.d.huang@phys.ncku.edu.tw",
    description="Software for dealing with some Quantum Problems",
    url="https://gitlab.com/research_by_y.t.d.huang/quawesome",

    classifiers=[
        "Environment :: Console",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering",
    ],
    keywords="Quantum",
    install_requires=requirements,
    setup_requires=['Cython>=0.29.14', 'pyscf>=1.6.5'],
    include_package_data=True,
    python_requires=">=3.6",
)
from setuptools import setup, find_packages

setup(name="pyvam",
      version="1.0.2",
      description='A Python package for comparative visualization of animal mitochondrial genomes',
      url='https://github.com/thecgs/PyVAM',
      author='Guisen Chen',
      author_email='thecgs001@foxmail.com',
      long_description="A Python package for visualizing animal mitochondrial",
      license='MIT License',
      packages=find_packages(exclude=["doc", ".github"]),
      keywords=['bioinformatics', 'mitochondrial', 'genome', 'rearrangement'],
      install_requires=["biopython", "plotly", "matplotlib"]
     )

from setuptools import setup, find_packages

def parse_requirements(filename):
    ''' Load requirements from a pip requirements file '''
    with open(filename, 'r') as fd:
        lines = []
        for line in fd:
            line.strip()
            if line and not line.startswith("#"):
                lines.append(line)
    return lines

requirements = parse_requirements('requirements.txt')


if __name__ == '__main__':
    setup(
        name='vaster',
        version='0.2.0',
        description = "Pipeline for image-plane fast transients search with ASKAP data", 
        author='Yuanming Wang',
        author_email='yuanmingwang@swin.edu.au',
        maintainer='Yuanming Wang',
        maintainer_email='yuanmingwang@swin.edu.au',
        install_requires=requirements,
        keywords=['vaster'],
        packages=find_packages(where='src'),
        package_dir={'': 'src'},
        entry_points={
            'console_scripts': [
                'askapsoft_rescale=vaster.askapsoft_rescale:main',
                'fix_dir=vaster.fix_dir:main', 
                'prepare_scripts=vaster.prepare_scripts:_main',
                'select_candidates=vaster.select_candidates:_main', 
                'inspect_usage=vaster.inspect_usage:_main', 
                'submit_slurm_jobs=vaster.submit_slurm_jobs:_main', 
                'prepare_data=vaster.prepare_data:_main',
                'plot_candidates=vaster.plot_candidates:_main', 
                'check_complete=vaster.check_complete:_main', 
            ],
        },
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ]
    )

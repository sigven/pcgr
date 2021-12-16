from setuptools import setup, find_packages

setup(
    name='pcgr',
    version='0.9.4',
    license='MIT',
    author='Sigve Nakken',
    author_email='sigven@gmail.com',
    description='Personal Cancer Genome Reporter (PCGR) - variant interpretation for precision cancer medicine',
    url='https://github.com/sigven/pcgr',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'pcgr = pcgr.main:cli'
        ]
    },
    scripts=[
        'scripts/annoutils.py',
        'scripts/cpsr_validate_input.py',
        'scripts/pcgr_summarise.py',
        'scripts/pcgr_validate_input.py',
        'scripts/pcgr_vcfanno.py',
        'scripts/vcf2tsv.py',
        'scripts/pcgrr.R'
    ]
)

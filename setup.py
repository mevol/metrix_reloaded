import os
from setuptools import setup

version = os.environ.get("MODULEVER", "0.0")

setup(
    #    install_requires = ['cothread'], # require statements go here
    name="metrix_reloaded",
    version=version,
    description="Module",
    author="Melanie Vollmar",
    author_email="melanie.vollmar@diamond.ac.uk",
    packages=["modules",
              "modules.metrix_integrate",
              "modules.metrix_db",
              "modules.metrix_ml"],
    install_requires=[
         "procrunner",
#        "tensorflow-gpu<=1.13",
#        "Keras",
#        "Pillow",
#        "procrunner",
#        "PyYaml",
#        "scikit-learn",
#        "mrcfile",
#        "pandas",
#        "logconfig",
#        "matplotlib",
#        "configargparse",
    ],
    scripts=[
            "bin/metrix_integrate/run_my_xia2"],
    
#    entry_points={
#        "console_scripts": [
#            "topaz3.prepare = topaz3.command_line_preparation:main",
#            "topaz3.test_split = topaz3.train_test_split:command_line",
#            "topaz3.predict_from_maps = topaz3.predictions:command_line",
#            "topaz3.filter = topaz3.filters:filter_command_line",
#        ]
#    },
    #    entry_points = {'console_scripts': ['test-python-hello-world = topaz3.topaz3:main']}, # this makes a script
    #    include_package_data = True, # use this to include non python files
    license="BSD3 license",
    zip_safe=False,
)

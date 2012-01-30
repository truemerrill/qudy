from distutils.core import setup
import qudylib
import os


def check_dir( target ):
    # Recursively check a target directory for python scripts.  This
    # function makes it easy to automatically package scripts and
    # input files in the binary distribution.
    
    # Save current path
    root = os.getcwd()
    
    # Change path to target
    path = os.path.abspath( target )
    os.chdir( path )
    
    files = [o for o in os.listdir('.') if os.path.isfile(o)]
    subdirs = [o for o in os.listdir('.') if os.path.isdir(o)]
    
    SCRIPTS = []
    for item in files:
        extension = os.path.splitext( item )[1]
        if extension == '.py' or extension == '.inp':
            SCRIPTS.append(path + '/' + item)
    
    for subdir in subdirs:
        M = check_dir( subdir )
        SCRIPTS = SCRIPTS + M
        # Return to parent directory
        os.chdir( path )
        
    # Return to original directory
    os.chdir( root )
    return SCRIPTS


# Grab all scripts in the examples directory
SCRIPTS = check_dir('examples')

# Add the "qudy" driver script
SCRIPTS.append( 'qudy' )

setup (
    name = 'qudy',
    version = qudylib.release,
    author = 'True Merrill, John Patrick Addison, Chingiz K, Kenneth Brown',
    author_email = 'true.merrill@gatech.edu',
    url = 'http://ww2.chemistry.gatech.edu/~brown/',
    packages = ['qudylib','qudylib/error_models','qudylib/hardware_models'],
    scripts = SCRIPTS,
    license = 'LICENSE.txt',
    description = 'A qubit dynamics simulator',
    long_description = open('README.txt').read()
)

#

from distutils.core import setup
import qudy
import os


def check_dir( target ):
    # Recursively check a target directory for python scripts.  This
    # function makes it easy to automatically package scripts and
    # input files in the source distribution.
    
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

# Try to build documentation using sphinx.  First check to see if the
# current tree is a development version. Query the user for their
# build preference.

root = os.getcwd()
docs = root + '/docs'    
if os.path.exists( docs + '/source' ):
    
    query = raw_input("Build the documentation (y/n) : ")
    if query.lower() == "y":
    
        print "Attempting to build documentation."
        try:
            os.chdir( docs )
            os.system("make html")
            os.system("make latexpdf")
            os.chdir( root )
            
        except:
            print "Could not build documentation."
    
    else:
        
        print "Documentation will not be built"


setup (
    name = 'qudy',
    version = qudy.release,
    author = 'True Merrill and collaborators',
    author_email = 'true.merrill@gatech.edu',
    url = 'http://ww2.chemistry.gatech.edu/~brown/',
    packages = ['qudy','qudy/error_models','qudy/hardware_models', \
                'qudy/optimal_control'],
    license = 'LICENSE.txt',
    description = 'A qubit dynamics simulator',
    long_description = open('README.txt').read(),
    requires = [
        "numpy (>= 1.4.1)",
        "scipy (>= 0.7.2)",
    ]
)

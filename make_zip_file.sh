mkdir AlphaImpute

# Assumes that the program and manual have both been built.

# To build the program run:
# NOTE: Binaries should be moved to the "binaries" folder and uploaded to bitbucket after builds.
#cmake . ; make

# to build the manual:
# NOTE: Building the docs requires two additional libraries, pandoc-tablenos and pandoc-fignos. See the readme in alphaimpute-docs for more information.
#( cd alphaimpute-docs ; ./make-complete.sh )

cp -r example AlphaImpute

# Copy in the documentation.
cp alphaimpute-docs/complete.pdf AlphaImpute/AlphaImputeUserManual.pdf

if [ $? != 0 ]; then                   # last command: echo
    echo "The manual needs to be built." # last command: [
    exit 1
fi


# Copy in the binaries
cp binaries/current_release/* AlphaImpute

# Create a version file

version=`git describe --tags --abbrev=0`
commit=`git rev-parse --short HEAD`

echo Version: $version > AlphaImpute/version.txt
echo Commit: $commit >> AlphaImpute/version.txt

zip -r AlphaImpute.zip AlphaImpute

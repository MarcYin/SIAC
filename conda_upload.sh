conda install conda-build  
conda install anaconda-client 
PKG_NAME=siac && USER=f0xy  
OS=$TRAVIS_OS_NAME-64  
mkdir ~/conda-bld 
conda config --set anaconda_upload no  
export CONDA_BLD_PATH=~/conda-bld 
export VERSION=$SIAC_VERSION
conda build .
export CONDA_PACKAGE=`conda build --output . | grep bz2`
echo $CONDA_PACKAGE
ls -lah $CONDA_BLD_PATH/$OS 
ls -lah $CONDA_BLD_PATH/
ls -lah $CONDA_BLD_PATH/noarch
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $CONDA_PACKAGE --force

# Only need to change these two variables
#PKG_NAME=siac
#USER=f0xy
#OS=$TRAVIS_OS_NAME-64
#mkdir ~/conda-bld
#conda config --set anaconda_upload no
#export CONDA_BLD_PATH=~/conda-bld
#export VERSION=`date +%Y.%m.%d`
#conda build . -c conda-forge
#ls $CONDA_BLD_PATH/$OS
#anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $(ls $CONDA_BLD_PATH/$OS/$PKG_NAME-$VERSION*.tar.bz2) --force

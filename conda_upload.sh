conda install conda-build  
conda install anaconda-client 
PKG_NAME=siac && USER=f0xy  
OS=$TRAVIS_OS_NAME-64  
mkdir ~/conda-bld 
conda config --set anaconda_upload no  
export CONDA_BLD_PATH=~/conda-bld 
conda build . -c conda-forge
ls -lah $CONDA_BLD_PATH/$OS 
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $(ls $CONDA_BLD_PATH/$OS/$PKG_NAME-$SIAC_VERSION*.tar.bz2) --force

# Only need to change these two variables
#PKG_NAME=siac
#USER=f0xy
#OS=$TRAVIS_OS_NAME-64
#mkdir ~/conda-bld
#conda config --set anaconda_upload no
#export CONDA_BLD_PATH=~/conda-bld
#export VERSION=`date +%Y.%m.%d`
#conda build . -c conda-forge
#ls $CONDA_BLD_PATH
#ls $CONDA_BLD_PATH/$OS
#anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $(ls $CONDA_BLD_PATH/$OS/$PKG_NAME-$VERSION*.tar.bz2) --force

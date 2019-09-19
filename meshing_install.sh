sudo apt-get -y install libcurl4-openssl-dev freeglut3-dev freeglut3

BASE_DIR=~/software/binaries

# MESHTOOL
MESHTOOL_SOURCE_DIR=${BASE_DIR}/MESHTOOL_source
MESHTOOL_BUILD_DIR=${BASE_DIR}/MESHTOOL_build
MESHTOOL_INSTALL_DIR=${BASE_DIR}/MESHTOOL_install
git clone https://c4science.ch/diffusion/9312/meshtool.git ${MESHTOOL_SOURCE_DIR}
cd ${MESHTOOL_SOURCE_DIR}
git checkout master

# VTK
VTK_SOURCE_DIR=${BASE_DIR}/VTK_source
VTK_BUILD_DIR=${BASE_DIR}/VTK_build
VTK_INSTALL_DIR=${BASE_DIR}/VTK_install
git clone https://github.com/Kitware/VTK ${VTK_SOURCE_DIR}
cd ${VTK_SOURCE_DIR}
git checkout v8.1.0
mkdir ${VTK_BUILD_DIR}
cd ${VTK_BUILD_DIR}
cmake -DCMAKE_INSTALL_PREFIX=${VTK_INSTALL_DIR} -DBUILD_EXAMPLES:BOOL=OFF -DBUILD_TESTING:BOOL=OFF ${VTK_SOURCE_DIR}
make
make install

# XERCES
XERCES_SOURCE_DIR=${BASE_DIR}/XERCES_source
XERCES_BUILD_DIR=${BASE_DIR}/XERCES_build
XERCES_INSTALL_DIR=${BASE_DIR}/XERCES_install
git clone https://github.com/apache/xerces-c ${XERCES_SOURCE_DIR}
cd ${XERCES_SOURCE_DIR}
git checkout Xerces-C_3_2_0
mkdir ${XERCES_BUILD_DIR}
cd ${XERCES_BUILD_DIR}
cmake -DCMAKE_INSTALL_PREFIX=${XERCES_INSTALL_DIR} ${XERCES_SOURCE_DIR}
make
make install

# # ZLIB
# ZLIB_SOURCE_DIR=${BASE_DIR}/ZLIB_source
# ZLIB_BUILD_DIR=${BASE_DIR}/ZLIB_build
# ZLIB_INSTALL_DIR=${BASE_DIR}/ZLIB_install
# git clone https://github.com/madler/zlib.git ${ZLIB_SOURCE_DIR}
# cd ${ZLIB_SOURCE_DIR}
# git checkout v1.2.11
# mkdir ${ZLIB_BUILD_DIR}
# cd ${ZLIB_BUILD_DIR}
# cmake -DCMAKE_INSTALL_PREFIX=${ZLIB_INSTALL_DIR} ${ZLIB_SOURCE_DIR}
# make
# make install

sudo apt-get update && \
     apt-get -y install libgmp3-dev libmpfr-dev libboost-dev


# BOOST
BOOST_SOURCE_DIR=${BASE_DIR}/BOOST_source
BOOST_TAR_FILE=${BASE_DIR}/BOOST.tar.bz2
wget -O ${BOOST_TAR_FILE} https://sourceforge.net/projects/boost/files/boost/1.66.0/boost_1_66_0.tar.bz2/download
mkdir ${BOOST_SOURCE_DIR}
tar -xjf ${BOOST_TAR_FILE} -C ${BOOST_SOURCE_DIR} --strip-components=1



# CGAL
CGAL_SOURCE_DIR=${BASE_DIR}/CGAL_source
CGAL_BUILD_DIR=${BASE_DIR}/CGAL_build
CGAL_INSTALL_DIR=${BASE_DIR}/CGAL_install
git clone https://github.com/CGAL/cgal.git ${CGAL_SOURCE_DIR}
cd ${CGAL_SOURCE_DIR}
git checkout releases/CGAL-4.8.1
mkdir ${CGAL_BUILD_DIR}
cd ${CGAL_BUILD_DIR}
cmake -DCMAKE_INSTALL_PREFIX=${CGAL_INSTALL_DIR} -DBUILD_SHARED_LIBS=OFF ${CGAL_SOURCE_DIR} -DBOOST_ROOT=${BOOST_SOURCE_DIR}
make
make install

# # XSD
XSD_TAR_FILE=${BASE_DIR}/XSD.tar.bz2
XSD_INSTALL_DIR=${BASE_DIR}/XSD_install
wget -O ${XSD_TAR_FILE} https://www.codesynthesis.com/download/xsd/4.0/linux-gnu/x86_64/xsd-4.0.0-x86_64-linux-gnu.tar.bz2
tar -xjf ${XSD_TAR_FILE} -C ${XSD_INSTALL_DIR} --strip-components=1

# # MESHTOOL
cd ${MESHTOOL_SOURCE_DIR}
mkdir ${MESHTOOL_BUILD_DIR}
cd ${MESHTOOL_BUILD_DIR}
cmake -DCMAKE_INSTALL_PREFIX=${MESHTOOL_INSTALL_DIR} -DCGAL_DIR=${CGAL_BUILD_DIR} -DVTK_DIR=${VTK_BUILD_DIR} -DXERCESC_ROOT_DIR=${XERCES_INSTALL_DIR} -DXSD_DIR=${XSD_INSTALL_DIR}/libxsd ${MESHTOOL_SOURCE_DIR}
make
make install

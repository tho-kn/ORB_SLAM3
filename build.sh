# Your package name, while not written in document, this packages are needed for compilation
PACKAGE="libssl-dev libboost-dev libboost-serialization-dev"

# Check if we are root
if [ "$EUID" -ne 0 ]
  then echo "Not root user, using sudo for installation"
  # Not root user, use sudo
  sudo apt-get install -y $PACKAGE
else
  echo "Root user, installing without sudo"
  # Root user, no need for sudo
  apt-get install -y $PACKAGE
fi

echo "Configuring and building Thirdparty/DBoW2 ..."

cd Thirdparty/DBoW2
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

cd ../../g2o

echo "Configuring and building Thirdparty/g2o ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

cd ../../Sophus

echo "Configuring and building Thirdparty/Sophus ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

cd ../../../

echo "Uncompress vocabulary ..."

cd Vocabulary
tar -xf ORBvoc.txt.tar.gz
cd ..

echo "Configuring and building ORB_SLAM3 ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4


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

<<<<<<< HEAD
#mkdir build
#cd build
#cmake .. -DCMAKE_BUILD_TYPE=Release
#make -j
=======
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
>>>>>>> 4452a3c4ab75b1cde34e5505a36ec3f9edcdc4c4

# Retrieve biscuit
git clone git@github.com:huishenlab/biscuit.git biscuit_master

# Move to biscuit directory and build
cd biscuit_master
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../ ../
make install

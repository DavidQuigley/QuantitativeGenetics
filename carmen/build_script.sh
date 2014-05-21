wget http://downloads.sourceforge.net/project/boost/boost/1.54.0/boost_1_54_0.tar.gz
gunzip boost_1_54_0.tar.gz
tar -xf boost_1_54_0.tar
cd boost_1_54_0
mkdir build
./bootstrap.sh --with-libraries=system,filesystem,regex,thread --prefix=/home/dquigley/boost_1_54_0/build
./b2
rm -f /home/dquigley/boost_1_54_0/stage/lib/*.dylib
rm -f /home/dquigley/boost_1_54_0/stage/lib/*.so*
    
# compiler include path: /home/dquigley/boost_1_54_0
# linker path: /home/dquigley/boost_1_54_0/stage/lib

mkdir /home/dquigley/notebook
mkdir /home/dquigley/notebook/code
cd /home/dquigley/notebook/code
git config --global user.name "David Quigley"
git config --global user.email "dquigley@cc.ucsf.edu"
git clone https://github.com/DavidQuigley/QuantitativeGenetics.git src
git pull https://github.com/DavidQuigley/QuantitativeGenetics.git

cd /home/dquigley/notebook/code/src/carmen/build

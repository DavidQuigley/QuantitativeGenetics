wget http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.gz
gunzip boost_1_57_0.tar.gz
tar -xf boost_1_57_0.tar
cd boost_1_57_0
sudo ./bootstrap.sh --with-libraries=system,filesystem,regex,thread
sudo ./b2 install link=static runtime-link=static
rm -f /usr/local/lib/libboost_filesystem.so*
rm -f /usr/local/lib/libboost_regex.so*
rm -f /usr/local/lib/libboost_system.so*
rm -f /usr/local/lib/libboost_thead.so*

# compiler include path: /usr/local/lib/boost_1_57_0
# linker path: /usr/local/lib

cd /notebook/code
git config --global user.name "David Quigley"
git config --global user.email "dquigley@cc.ucsf.edu"
git clone https://github.com/DavidQuigley/QuantitativeGenetics.git src
cd /notebook/code/src/carmen/build
mkdir obj
make
sudo mv libCarmen.a /usr/local/lib
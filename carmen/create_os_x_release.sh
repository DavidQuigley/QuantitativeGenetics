# copy OS X binaries into Carmen.app
# Note that we rename the jasper executable to carmen
cp /notebook/code/src/jasper/build/Release/jasper /notebook/code/release/Carmen.app/Contents/MacOS/Carmen
cp /notebook/code/src/spear/build/Release/spear /notebook/code/release/Carmen.app/Contents/MacOS/spear
cp /notebook/code/src/CorrelationNetwork/build/Release/CorrelationNetwork /notebook/code/release/Carmen.app/Contents/MacOS/CorrelationNetwork
cp /notebook/code/src/miner/build/Release/miner /notebook/code/release/Carmen.app/Contents/MacOS/miner
cp /notebook/code/src/coreminer/build/Release/coreminer /notebook/code/release/Carmen.app/Contents/MacOS/coreminer
cp /notebook/code/src/difference/build/Release/difference /notebook/code/release/Current/Carmen.app/Contents/MacOS/differences
cp /notebook/code/release/Current/mouse_annotation.txt /notebook/code/release/Current/Carmen.app/Contents/MacOS
cp /notebook/code/release/Current/human_annotation.txt /notebook/code/release/Current/Carmen.app/Contents/MacOS
cp /notebook/code/release/Current/gene_ontology.1_2.obo.compressed /notebook/code/release/Current/Carmen.app/Contents/MacOS
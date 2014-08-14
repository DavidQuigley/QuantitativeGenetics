# copy OS X binaries into Carmen.app
# Note that we rename the jasper executable to carmen
copy \notebook\code\src\jasper\Release\jasper.exe \notebook\code\release\Current\Carmen.exe /Y
copy \notebook\code\src\spear\Release\spear.exe \notebook\code\release\Current /Y
copy \notebook\code\src\CorrelationNetwork\Release\CorrelationNetwork.exe \notebook\code\release\Current /Y
copy \notebook\code\src\miner\Release\miner.exe \notebook\code\release\Current /Y
copy \notebook\code\src\coreminer\Release\coreminer.exe \notebook\code\release\Current /Y
copy \notebook\code\src\difference\Release\difference.exe \notebook\code\release\Current /Y
copy \notebook\code\src\annotations\mouse_annotation.txt \notebook\code\release\Current /Y
copy \notebook\code\src\annotations\human_annotation.txt \notebook\code\release\Current /Y
copy \notebook\code\src\annotations\gene_ontology.obo.compressed \notebook\code\release\Current /Y
copy \notebook\code\src\annotations\CARMEN_documentation.pdf \notebook\code\release\Current /Y
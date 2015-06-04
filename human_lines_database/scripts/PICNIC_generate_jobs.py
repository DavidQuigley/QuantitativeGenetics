# files in outdir/raw should have been processed by celConverter.jar
# for each file in outdir/raw, generate a shell script that will call processing and HMM
# then generate a shell script that will start all of these jobs when run

import os
dir_HOME = "/home/quigleyd"
dir_PICNIC_outdir = '/mnt/speed/quigleyd/CN/PICNIC_outdir'
dir_PICNIC_software = '/mnt/speed/quigleyd/software/PICNIC'
existing = os.listdir(dir_PICNIC_outdir + '/raw')

f = open(dir_HOME + '/scripts/PICNIC_LIBRARY_HEADER')
LD_header = ''.join( f.readlines() )
f.close()

for filename in existing:
    root_filename = filename.replace(".feature_intensity", "").replace("CGP_", "")
    f = open( dir_HOME + '/scripts/process_PICNIC/' + filename + '.sh', 'w')
    f.write( "#!/bin/bash\n" )
    f.write( "ulimit -c 0\n" )
    f.write( 'FILENAME="CGP_' + root_filename + '"\n' )
    f.write( LD_header )
    f.write( dir_PICNIC_software + "/preprocessing " )
    f.write( "CGP_" + root_filename + ".feature_intensity " )
    f.write( dir_PICNIC_software + "/info/ " ) # intentional trailing /
    f.write( dir_PICNIC_outdir + "/raw/ " ) # intentional trailing /
    f.write( dir_PICNIC_outdir + "/output/ " ) # intentional trailing /
    f.write( dir_PICNIC_outdir + "/\n") # intentional trailing /
    f.write( "\n")
    f.write( "IN_PI=$(cat " + dir_PICNIC_outdir + "/output2/" + root_filename + "_feature.TXT/ploidy_" + root_filename + "_feature.TXT.csv | cut -f 1 -d ',')\n")
    f.write( "PLOIDY=$(cat " + dir_PICNIC_outdir + "/output2/" + root_filename + "_feature.TXT/ploidy_" + root_filename + "_feature.TXT.csv | cut -f 2 -d ',')\n")
    f.write( "ALPHA=$(cat " + dir_PICNIC_outdir + "/output2/" + root_filename + "_feature.TXT/ploidy_" + root_filename + "_feature.TXT.csv | cut -f 3 -d ',')\n")
    f.write( 'echo "Starting HMM"\n')
    f.write( 'echo $IN_PI\n')
    f.write( 'echo $PLOIDY\n')
    f.write( 'echo $ALPHA\n')
    f.write( "/mnt/speed/quigleyd/software/PICNIC/HMM ")
    f.write( "CGP_" + root_filename + ".feature_intensity " )
    f.write( "/mnt/speed/quigleyd/software/PICNIC/info/ ") # intentional trailing /
    f.write( dir_PICNIC_outdir + "/output/ " ) # intentional trailing /
    f.write( dir_PICNIC_outdir + "/ 8 $IN_PI $PLOIDY $ALPHA" )
    f.close()

fo = open( dir_HOME + '/scripts/run_PICNIC_jobs.sh', 'w')
for filename in sorted(existing):
    root_filename = filename.replace(".feature_intensity", "").replace("CGP_", "").split("_")
    nickname = root_filename[0] + "_" + root_filename[-1]
    fo.write( "qsub -N " + nickname + " -o " + dir_HOME + "/logs/out_" + nickname )
    fo.write( " -e " + dir_HOME + "/logs/err_" + nickname + " -S /bin/bash ")
    fo.write( dir_HOME + "/scripts/process_PICNIC/" + filename + ".sh\n")

fo.close()



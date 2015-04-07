import argparse
parser = argparse.ArgumentParser(description='Lift over probe reference file')

parser.add_argument('original', help='original Proberef file')
parser.add_argument('output', help='new Proberef file')
parser.add_argument('liftover_perfect', help='Perfect matches for liftover results')
parser.add_argument('liftover_best_guess', help='best fit for failed liftover results')

args = parser.parse_args()

fn_pos_new = args.liftover_perfect
id2chr = {}
id2pos = {}

print( "********************************************************************************" )
print( "Generating new ProbeRef file" )
print "Reading from perfect match liftover file:" + args.liftover_perfect
f = open(fn_pos_new)
for line in f:
    a = line.rstrip('\r\n').split('\t')
    id2pos[ a[3] ] = a[1]
    id2chr[ a[3] ] = a[0]

f.close()
print( "Reading from best guess liftover file:" + args.liftover_best_guess )
f = open(args.liftover_best_guess)
for line in f:
    a = line.rstrip('\r\n').split('\t')
    id2pos[ a[3] ] = a[1]
    id2chr[ a[3] ] = a[0]

f.close()

print( "Reading from original ProbeRef file:" + args.original )
f = open(args.original)
fo = open(args.output, 'w')
fo.write( "# Lifted over to HG19\n")
fo.write(f.readline())
for line in f:
    a = line.rstrip('\r\n').split(',')
    probeset = a[0]
    a[2] = id2chr[probeset]
    a[3] = id2pos[probeset]
    fo.write( ','.join(a) + '\n' )

fo.close()
f.close()

print "wrote to " + args.output
import sys
import argparse
parser = argparse.ArgumentParser(description='Round any matrix by element')
parser.add_argument('--signif', dest="signif", help='number of significant digits', default=3, type=int)

args = parser.parse_args()
nsig = args.signif

for line in sys.stdin:
    a=line.rstrip('\r\n').split('\t')
    for i, val in enumerate(a):
        try:
            a[i] = round( float(val), nsig)
        except ValueError:
            pass
    print '\t'.join([str(x) for x in a])


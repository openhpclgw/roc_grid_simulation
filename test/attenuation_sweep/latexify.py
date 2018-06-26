import sys
"""
Creates the core of a latex table from the errors file in the given path
"""
case_dir = sys.argv[1]
f1 = open(case_dir+'/size8/errors', 'r')
f2 = open(case_dir+'/size16/errors', 'r')
f3 = open(case_dir+'/size32/errors', 'r')
line_frmt = '& {l1_vals}\t&\t{l2_vals}\t&\t{l3_vals} \\\\'
for l1, l2, l3 in zip(f1, f2, f3):
    print(line_frmt.format(
            l1_vals=' & '.join('{:0.2f}'.format(float(v)) for v in l1.split()),
            l2_vals=' & '.join('{:0.2f}'.format(float(v)) for v in l2.split()),
            l3_vals=' & '.join('{:0.2f}'.format(float(v)) for v in l3.split())),
          end='\n\n')

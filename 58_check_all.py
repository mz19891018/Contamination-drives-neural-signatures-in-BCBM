import gzip
lines=[]
with gzip.open(r'D:\BCBM_Project\data\raw\GSE2603_series_matrix.txt.gz','rt',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_characteristics_ch1'):
            vals=[v.strip().strip('"') for v in line.split('\t')[1:]]
            keys=set()
            for v in vals:
                if ':' in v: keys.add(v.split(':')[0].strip())
            print('Line keys:', sorted(keys))
            # if bm event in this line
            for v in vals:
                if 'bm event' in v.lower() or 'bmfs' in v.lower():
                    print('  sample:', v)
            lines.append(vals)
print('total char lines:', len(lines))

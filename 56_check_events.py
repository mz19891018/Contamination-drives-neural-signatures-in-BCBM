import gzip
# GSE2603
with gzip.open(r'D:\BCBM_Project\data\raw\GSE2603_series_matrix.txt.gz','rt',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_characteristics_ch1'):
            vals=[v.strip().strip('"') for v in line.split('\t')[1:]]
            bm_events=[v for v in vals if v.startswith('bm event')]
            bmfs=[v for v in vals if v.startswith('bmfs')]
            tissues=[v for v in vals if v.startswith('tissue type')]
            cell_lines=[v for v in tissues if 'cell line' in v.lower() or 'MDA' in v or 'subpopulation' in v.lower()]
            print('GSE2603 bm event values:', set(bm_events))
            print('GSE2603 bmfs count:', len(bmfs))
            print('GSE2603 cell line samples:', len(cell_lines))
            # print first few bm event + bmfs pairs
            break
# GSE12276 keys
keys=set()
with gzip.open(r'D:\BCBM_Project\data\raw\GSE12276_series_matrix.txt.gz','rt',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_characteristics_ch1'):
            vals=[v.strip().strip('"') for v in line.split('\t')[1:]]
            for v in vals:
                if ':' in v: keys.add(v.split(':')[0].strip())
print('GSE12276 keys:', sorted(keys))

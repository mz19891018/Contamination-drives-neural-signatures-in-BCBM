import gzip
with gzip.open(r'D:\BCBM_Project\data\raw\GSE2603_series_matrix.txt.gz','rt',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_characteristics_ch1'):
            vals=[v.strip().strip('"') for v in line.split('\t')[1:]]
            bm_vals=[v for v in vals if 'bm' in v.lower()]
            print('GSE2603 bm-related values (first 10):', bm_vals[:10])
            print('total bm-related:', len(bm_vals))
            # event-like
            event_vals=[v for v in vals if 'event' in v.lower()]
            print('event-related:', set(event_vals))
            # tissue
            tissue_vals=[v for v in vals if 'tissue' in v.lower()]
            print('tissue types:', set(tissue_vals))
            break

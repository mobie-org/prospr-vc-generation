import csv
import json
import sys
import os


if __name__ == '__main__':
    expr_table_file = sys.argv[1]
    dict_file = sys.argv[2]
    new_expr_table_file = os.path.splitext(expr_table_file)[0] + '_upd' + os.path.splitext(expr_table_file)[1]

    with open(dict_file) as json_file:
        names_dict = json.load(json_file)

    short_dict = {'-'.join(key.split('-')[4:]): '-'.join(value.split('-')[4:])
                  for key, value in names_dict.items() if key.startswith('prospr')}
    short_dict.pop('virtual-cells-labels', None)
    filt_dict = {(key.split('-')[-1] if key.startswith('segmented') else key.split('-')[0]):
                 (value.split('-')[-1] if value.startswith('segmented') else value)
                 for key, value in short_dict.items()}

    with open(expr_table_file) as csv_read:
        csv_reader = csv.reader(csv_read, delimiter='\t')
        gene_names = csv_reader.__next__()
        gene_names = [name.split('--')[0] for name in gene_names]
        new_gene_names = [filt_dict[name] if name in filt_dict else name for name in gene_names]
        new_gene_names[0] = 'VCName'
        with open(new_expr_table_file, 'w') as csv_write:
            csv_writer = csv.writer(csv_write, delimiter='\t')
            csv_writer.writerow(['label_id', ] + new_gene_names)
            for i, row in enumerate(csv_reader):
                _ = csv_writer.writerow([i + 1, ] + row)

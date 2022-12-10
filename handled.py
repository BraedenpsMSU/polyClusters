import os
from typing import Optional


def enforce_str_dict(inp):
    for item in inp:
        if inp[item] is None:
            pass
        else:
            inp[item] = str(inp[item])
    return inp


def check_memo_asm(inp_key, store=[]):
    out_file = 'assembly.csv'
    primary_key = 'accession'
    if not store:
        saved = load_to_dict(out_file)
        saved_dict = dict()
        for entry in saved:
            saved_dict[entry[primary_key]] = entry
        store.append(saved_dict)
    return store[0].get(inp_key)


def check_memo_gene(inp_key, store=[]):
    out_file = 'protein.csv'
    primary_key = 'gene_id'
    if not store:
        saved = load_to_dict(out_file)
        saved_dict = dict()
        for entry in saved:
            saved_dict[entry[primary_key]] = entry
        store.append(saved_dict)
    return store[0].get(inp_key)


def check_memo_protein(inp_key, store=[]):
    out_file = 'protein.csv'
    primary_key = 'protein_id'
    if not store:
        saved = load_to_dict(out_file)
        saved_dict = dict()
        for entry in saved:
            saved_dict[entry[primary_key]] = entry
        store.append(saved_dict)
    return store[0].get(inp_key)


def extract_primary(input_list, inpfile, key):
    header = False
    outp = []
    primary_key_index: Optional[int] = None
    try:
        with open(inpfile, 'r') as f:
            for line in f.readlines():
                templine = line.replace('\n', '')
                if not header:
                    headers = templine.split(',')
                    header = True
                    primary_key_index = headers.index(key)
                else:
                    values = templine.split(',')
                    assert values[primary_key_index] != ''
                    outp.append(
                        values[primary_key_index]
                    )
                    input_list.append(
                        values[primary_key_index]
                    )
    except FileNotFoundError:
        pass
    return outp


def load_to_dict(inpfile: str):
    header = False
    headers = []
    outp = []
    with open(inpfile, 'r') as f:
        for line in f.readlines():
            templine = line.replace('\n', '')
            if not header:
                headers = templine.split(',')
                header = True
            else:
                values = templine.split(',')
                outp.append(
                    {headers[i]: str(values[i]) if values[i] != '' else None for i in range(0, len(headers))}
                )
    return outp


def write_asm(inp, store=[]):
    out_file = 'assembly.csv'
    primary_key = 'accession'
    if not store:
        extract_primary(store, out_file, primary_key)
    colm = ['taxon_id', 'organism', 'accession', 'gc_count', 'seq_length', 'chromosomes',
            'gene_total', 'gene_num_protein_coding', 'gene_num_non_coding', 'gene_num_pseudogene']
    return write_to_table(inp, out_file, columns=colm, store=store, primary_key=primary_key)


def write_gene(inp, store=[]):
    out_file = 'gene.csv'
    primary_key = 'gene_id'
    if not store:
        extract_primary(store, out_file, primary_key)
    colm = ['gene_id', 'gene_symbol', 'proteins_for_gene', 'assembly', 'taxon_id', 'description']
    return write_to_table(inp, out_file, columns=colm, store=store, primary_key=primary_key)


def write_protein(inp, store=[]):
    out_file = 'protein.csv'
    primary_key = 'protein_id'
    if not store:
        extract_primary(store, out_file, primary_key)
    colm = ['protein_id', 'protein_title', 'cdd_id', 'gene_id', 'taxon_id', 'poly_type']
    return write_to_table(inp, out_file, columns=colm, store=store, primary_key=primary_key)


def write_taxon(inp, store=[]):
    out_file = 'taxon.csv'
    primary_key = 'taxon_id'
    colm = ['taxon_id', 'name',
            'species_id', 'species',
            'genus_id', 'genus',
            'family_id', 'family',
            'order_id', 'order',
            'class_id', 'class',
            'phylum_id', 'phylum']
    if not store:
        extract_primary(store, out_file, primary_key)
    return write_to_table(inp, out_file, columns=colm, store=store, primary_key=primary_key)


def write_to_table(inp, out_file, columns, store: Optional[list] = None, primary_key: Optional[str] = None):
    assert (store is None) == (primary_key is None)
    try:
        if os.path.getsize(out_file) == 0:
            with open(out_file, 'w') as out:
                out.write(f"{','.join(columns)}\n")
                new_line = [str(inp[column]).replace(',', ' ') for column in columns]
                out.write(f"{','.join(new_line)}\n")
        else:
            with open(out_file, 'a') as out:
                if primary_key:
                    if inp[primary_key] not in store:
                        new_line = [str(inp[column]).replace(',', ' ') for column in columns]
                        out.write(f"{','.join(new_line)}\n")
                        store.append(primary_key)
                else:
                    new_line = [str(inp[column]).replace(',', ' ') for column in columns]
                    out.write(f"{','.join(new_line)}\n")

    except FileNotFoundError:
        with open(out_file, 'w') as out:
            out.write(f"{','.join(columns)}\n")
            new_line = [str(inp[column]).replace(',', ' ') for column in columns]
            out.write(f"{','.join(new_line)}\n")

if __name__ == "__main__":
    pass

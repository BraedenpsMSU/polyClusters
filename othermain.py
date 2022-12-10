from collections import OrderedDict
from time import sleep

import pandas
from requests import post, get


import NCBIQueries
from QueryDriver import QueryDriver
from pandas import DataFrame

polymerasesTable = pandas.read_csv('pssmIds_csv.csv')
polymerases= {str(row.cdd_uid): row.Type for row in polymerasesTable.itertuples()}

# polymerases = {
#   '176474': 'A',
#   '176476': 'A',
#   '176477': 'A',
#   '176479': 'A',
#   '176480': 'A',
#   '176454': 'Y',
#   '176458': 'Y',
#   '176459': 'Y',
#   '99914': 'B',
#   '99920': 'B',
#   '99921': 'B',
#   '234767': 'C',
#   '273161': 'C',
#   '225087': 'C',
#   '273601': 'C',
#   '213988': 'C',
#   '213989': 'C',
#   '213990': 'C',
#   '213997': 'C',
#   '239931': 'C',
#   '239930': 'C',
#   '213991': 'X',
#   '431250': 'X',
#   '143386': 'X'
# }

qd = QueryDriver()

def write_out_dict(value: OrderedDict | list, name=None):
    output = [] if name is None else [name]
    if isinstance(value, OrderedDict):
        for key in value:
            if not isinstance(value[key], OrderedDict) and not isinstance(value[key], list):
                yield output.copy() + [key, value[key]]
            else:
                for item in write_out_dict(value[key], name=key):
                    yield output.copy() + item
    if isinstance(value, list):
        for index, value in enumerate(value):
            if not isinstance(value, OrderedDict) and not isinstance(value, list):
                yield output.copy() + [index, value]
            else:
                for item in write_out_dict(value, name=index):
                    yield output.copy() + item


def list_paths(value: OrderedDict, as_list=False):
    for item in write_out_dict(value):
        if as_list:
            print(item)
        print(''.join(f'[{repr(key)}]' for key in item))
    print('')


def iterate_if_list(inp: list | OrderedDict, key: str|list[str]|None):
    def extractor(_inp: OrderedDict):
        if isinstance(key, list):
            place = _inp
            for subkey in key:
                place = place[subkey]
            return place
        elif isinstance(key, str):
            return _inp[key]
        else:
            return _inp
    if isinstance(inp, list):
        for item in inp:
            yield extractor(item)
    else:
        yield extractor(inp)


def run_taxon_assembly(taxon_id):
    qd('taxonomy', 'assembly', [str(taxon_id)], linkname="taxonomy_assembly")
    outp = []
    for results in qd.run_query():
        for _id in iterate_if_list(results['eLinkResult']['LinkSet']['LinkSetDb']['Link'], 'Id'):
            outp.append(_id)
    return outp


def run_gene_protein(gene_ids: list[str]):
    chunks = qd.factory.list_chunker(gene_ids, chunk_size=149)
    qd('gene', 'protein', chunks.pop(), linkname="gene_protein_refseq")
    outp = []
    while True:
        for results in qd.run_query():
            for _id in iterate_if_list(results['eLinkResult']['LinkSet']['LinkSetDb']['Link'], 'Id'):
                outp.append(_id)
        if chunks:
            qd('gene', 'protein', chunks.pop(), linkname="gene_protein_refseq")
        else:
            break
    return outp


def run_protein_cdd(gene_ids: list[str]):
    chunks = qd.factory.list_chunker(gene_ids, chunk_size=149)
    qd('protein', 'cdd', chunks.pop(), linkname="protein_cdd" )
    outp = []
    assm = {'A': 0, 'Y': 0, 'B': 0, 'C': 0, 'X': 0}
    while True:
        for results in qd.run_query():
            if 'LinkSetDb' in results['eLinkResult']['LinkSet']:
                for _id in iterate_if_list(results['eLinkResult']['LinkSet']['LinkSetDb']['Link'], 'Id'):
                    outp.append(_id)
                    if _id in polymerases.keys():
                        assm[polymerases[_id]] += 1
            else:
                print(results)
        if chunks:
            qd('protein', 'cdd', chunks.pop(), linkname="protein_cdd")
        else:
            break
    return assm, outp


def run_taxon_genes(taxon, reference):
    qd(taxon=taxon)
    outp = []
    while True:
        page_token = None
        for results in qd.run_query():
            list_paths(results)
            for gene in iterate_if_list(results['reports'], None):
                for accession in iterate_if_list(gene['gene']['annotations'],['assembly_accession']):
                    if accession == reference:
                        outp.append(gene['gene']['gene_id'])
            page_token = results.get("next_page_token")
        if page_token:
            qd(taxon=taxon, page_token=page_token)
        else:
            break
    return outp


if __name__ == "__main__":
    cmd = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report"

    json_arg = {"filters":
                    {"reference_only": True,  # only return references
                     "assembly_level": ["complete_genome"],  # only entire complete genomes
                     "exclude_paired_reports": True,
                     "assembly_version": "current"},
                "page_size": 2, "page_token": None,
                "returned_content": "COMPLETE",
                "taxons": ["2"]  # bacteria (eubacteria) id
                }

    resp = post(cmd, json=json_arg).json(object_pairs_hook=OrderedDict)
    # print(resp["total_count"])
    values = [(item['organism']['tax_id'], item['accession']) for item in resp['reports']]
    # print(values)
    run_taxon_genes(*values[0])
    print(*values)



    sleep(.4)

    # for taxon_id in values:
    #     if taxon_id[0] in result_df[['TaxonID']].values:
    #         continue
    #     results_list = dict()
    #     results_list['TaxonID'] = taxon_id[0]
    #     results_list['Assembly'] = taxon_id[1]
    #     genes = run_taxon_genes(*taxon_id)
    #     proteins = run_gene_protein(genes)
    #     assm = run_protein_cdd(proteins)
    #     results_list.update(assm[0])
    #     result_df = pandas.concat([
    #         result_df,
    #         pandas.DataFrame(results_list, index=[0])
    #         ],
    #         ignore_index=True
    #     )
    #     result_df.to_csv('Results.csv', sep='\t')


        # print(genes)
        # print(len(genes))
        # print(proteins)
        # print(len(proteins))
        # print(assm)
        # input()


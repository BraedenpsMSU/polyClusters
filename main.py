import json
import os
from collections import OrderedDict

import pandas
from requests import post, get
from othermain import *

import NCBIQueries
import QueryDriver
from pandas import DataFrame
import handled


qd = QueryDriver.QueryDriver()

Error_file_path = 'ErrorFile.txt'

polymerasesTable = pandas.read_csv('pssmIds_csv.csv')
polymerases = {str(row.cdd_uid): row.Type for row in polymerasesTable.itertuples()}
accepted_domains = [str(row.cdd_uid) for row in polymerasesTable.itertuples()]


def store_bad(results):
    try:
        with open(Error_file_path, 'r') as fp:
            next_line = len(fp.readlines())
        print(f'writing bad results starting at line {next_line} of error file')
        with open(Error_file_path, 'a') as fp:
            fp.write('---------------------------------------\n')
            fp.write(
                json.dumps(results, indent=4)
            )
            fp.write('\n')
    except Exception as e:
        print('Error writing error messsage - not written')


def get_Taxonomy(taxon_ids):
    chunks = qd.factory.list_chunker(taxon_ids, chunk_size=20)
    output = []
    for chunk in chunks:
        qd(lineage=chunk)
        for results in qd.run_query():
            for ind, taxon_id in enumerate(chunk):
                try:
                    result = dict()
                    lineage = results['taxonomy_nodes'][ind]['taxonomy']['lineage']
                    result['name'] = results['taxonomy_nodes'][ind]['taxonomy']['organism_name']
                    result['taxon_id'] = results['taxonomy_nodes'][ind]['taxonomy']['tax_id']
                    lineage = [str(taxon_id) for taxon_id in lineage]
                    result.update(get_name_rank(lineage))
                    output.append(result)
                except (ValueError,  KeyError):
                    print(f"Missing value or server response was bad - skipping {taxon_id=}")
                    store_bad(results)
    return output


def get_name_rank(taxon_ids):
    _levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    qd('taxonomy', taxon_ids, summary=True)
    outp = {lvl: None for lvl in _levels}
    outp.update({lvl + '_id': None for lvl in _levels})
    for result in qd.run_query():
        for taxon_id in taxon_ids:
            if result['result'][taxon_id]['rank'] in _levels:
                outp[result['result'][taxon_id]['rank']] = result['result'][taxon_id]['scientificname']
                outp[result['result'][taxon_id]['rank'] + '_id'] = str(result['result'][taxon_id]['uid'])
    return outp


def gene_to_protein(genes: list) -> list[str]:
    chunks = qd.factory.list_chunker(genes, chunk_size=125)
    outp = []
    for chunk in chunks:
        qd('gene', 'protein', chunk, linkname="gene_protein_refseq")
        for results in qd.run_query():
            try:
                for linksetdbs in iterate_if_list(results['linksets'], None):
                    for link in iterate_if_list(linksetdbs['linksetdbs'], None):
                        for protein_id in iterate_if_list(link['links'], None):
                            outp.append(str(protein_id))
            except (ValueError,  KeyError):
                print(f"Missing value or server response was bad - skipping {chunk=}")
                store_bad(results)
    return outp


def protein_to_cdd(genes: list) -> list[str]:
    chunks = qd.factory.list_chunker(genes, chunk_size=125)
    outp = []
    for chunk in chunks:
        qd('protein', 'cdd', chunk, linkname="protein_cdd")
        for results in qd.run_query():
            try:
                # list_paths(results)
                for cdd_id in results['linksets'][0]['linksetdbs'][0]['links']:
                    # print(type(cdd_id), f"{cdd_id=}")
                    if polymerases.get(cdd_id):
                        outp.append(str(cdd_id))
            except (ValueError, KeyError):
                print(f"Missing value or server response was bad - skipping {chunk=}")
                store_bad(results)
    return outp


def taxon_to_genes(taxon, reference) -> list[dict]:
    qd(taxon=taxon)
    outp = []
    while True:
        page_token = None
        for results in qd.run_query():
            # list_paths(results)
            try:
                for gene in iterate_if_list(results['reports'], None):
                    try:
                        for accession in iterate_if_list(gene['gene']['annotations'], 'assembly_accession'):
                            if accession == reference:
                                gene_data = dict()
                                gene_data['assembly'] = accession
                                gene_data['gene_id'] = gene['gene']['gene_id']
                                gene_data['taxon_id'] = gene['gene']['tax_id']
                                gene_data['gene_symbol'] = gene['gene'].get('symbol', None)
                                gene_data['proteins_for_gene'] = gene['gene'].get('protein_count', 0)
                                desc = gene['gene'].get('description', 'No description')[:100].replace(',', ' ')
                                gene_data['description'] = desc
                                outp.append(
                                    handled.enforce_str_dict(
                                        gene_data
                                    )
                                )
                    except (ValueError,  KeyError):
                        print(f"Missing value or server response was bad - skipping this {gene}")
                        store_bad(results)
            except (ValueError,  KeyError):
                print(f"Missing value or server response was bad - skipping this page with {page_token=}")
                store_bad(results)
            page_token = results.get("next_page_token")
        if page_token:
            qd(taxon=taxon, page_token=page_token)
        else:
            break
    return outp


def match_protein_cdd(proteins: list[dict]) -> list[dict]:
    def match_helper(helper_proteins: list):
        if not helper_proteins:
            return ()
        output = []
        cdds = protein_to_cdd(
            helper_proteins
        )
        if cdds:
            if len(helper_proteins) == 1:
                for cdd in cdds:
                    entry = dict()
                    entry['protein_id'] = helper_proteins[0]
                    entry['cdd_id'] = cdd
                    output.append(handled.enforce_str_dict(entry))
                return tuple(output)
            else:
                even_results = match_helper(helper_proteins[::2])
                odd_results = match_helper(helper_proteins[1::2])
                output += list(even_results)
                output += list(odd_results)
                return tuple(output)
        else:
            return tuple[output]
    return match_helper(proteins)


def match_gene_to_cdd(genes: list[dict]) -> list[dict]:
    def match_helper(helper_genes: list):
        if not helper_genes:
            return ()
        output = []
        gproteins = gene_to_protein(
            helper_genes
        )
        cdds = protein_to_cdd(
            gproteins
        )
        if cdds:
            if len(helper_genes) == 1:
                matches = match_protein_cdd(gproteins)
                for match in matches:
                    entry = match.copy()
                    entry['gene_id'] = helper_genes[0]
                    output.append(
                        handled.enforce_str_dict(
                            entry
                        )
                    )
                return tuple(output)
            else:
                even_results = match_helper(helper_genes[::2])
                odd_results = match_helper(helper_genes[1::2])
                output += list(even_results)
                output += list(odd_results)
                return tuple(output)
        else:
            return tuple(output)
    output = []
    genes_ids = [gene['gene_id'] for gene in genes]
    for key in genes_ids:
        if handled.check_memo_gene(key):
            value = handled.check_memo_gene(key)
            entry = dict()
            entry['gene_id'] = key
            entry['protein_id'] = value['protein_id']
            entry['cdd_id'] = value['cdd_id']
            output.append(handled.enforce_str_dict(entry))
    chunks = qd.factory.list_chunker(genes_ids, chunk_size=250)
    for chunk in chunks:
        output += list(match_helper(chunk))
    return output


def get_protein(proteins) -> list[dict]:
    chunks = qd.factory.list_chunker(proteins, chunk_size=75)
    new_data = dict()
    for chunk in chunks:
        protein_ids = [protein['protein_id'] for protein in chunk]
        qd('protein', protein_ids, summary='True')
        for results in qd.run_query():
            try:
                for protein_id in results['result']:
                    if protein_id == 'uids':
                        continue
                    try:
                        new_data[protein_id] = dict()
                        new_data[protein_id]['taxon_id'] = results['result'][protein_id]['taxid']
                        new_data[protein_id]['protein_title'] = results['result'][protein_id]['title'].replace(',', ' ')
                    except (ValueError,  KeyError):
                        print(f"Missing value or server response was bad - skipping {chunk=}")
                        store_bad(results)
            except (ValueError,  KeyError):
                print(f"Missing value or server response was bad - skipping {chunk=}")
                store_bad(results)
    for value in proteins:
        value.update(new_data[value['protein_id']])
        value['poly_type'] = polymerases[value['cdd_id']]
        handled.enforce_str_dict(
            value
        )
    return proteins


def build_assembly_data(page_token=None):
    cmd = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report"
    json_arg = {"filters":
                    {"reference_only": True,  # only return references
                     "assembly_level": ["complete_genome"],  # only entire complete genomes
                     "exclude_paired_reports": True,
                     "assembly_version": "current"},
                "page_size": 700, "page_token": page_token,
                "returned_content": "COMPLETE",
                "taxons": ["2"]  # bacteria (eubacteria) id
                }

    resp = post(cmd, json=json_arg).json(object_pairs_hook=OrderedDict)
    outp = []
    page_token = resp.get('next_page_token', None)
    for item in resp['reports']:
        asm_data = dict()
        try:
            asm_data['accession'] = item['accession']
            asm_data['taxon_id'] = item['organism']['tax_id']
            asm_data['organism'] = item['organism'].get('organism_name')
            if item.get('assembly_stats'):
                asm_data['seq_length'] = item['assembly_stats'].get('total_sequence_length')
                asm_data['gc_count'] = item['assembly_stats'].get('gc_count', '')
                asm_data['chromosomes'] = item['assembly_stats'].get('total_number_of_chromosomes')
            try:
                asm_data['gene_total'] = item['annotation_info']['stats']['gene_counts']['total']
                asm_data['gene_num_protein_coding'] = item['annotation_info']['stats']['gene_counts']['protein_coding']
                asm_data['gene_num_non_coding'] = item['annotation_info']['stats']['gene_counts']['non_coding']
                asm_data['gene_num_pseudogene'] = item['annotation_info']['stats']['gene_counts']['pseudogene']
            except (ValueError,  KeyError):
                print(f"No gene count data for {item['organism']['tax_id']} - omitting it")
                asm_data['gene_total'] = None
                asm_data['gene_num_protein_coding'] = None
                asm_data['gene_num_non_coding'] = None
                asm_data['gene_num_pseudogene'] = None
            outp.append(asm_data)
            handled.write_asm(
                handled.enforce_str_dict(asm_data)
            )
        except (ValueError,  KeyError):
            print('Skipping this block')
    if page_token:
        sleep(.4)
        return outp + build_assembly_data(page_token=page_token)
    else:
        return outp



def get_assembly_data():
    try:
        if os.path.getsize('assembly.csv') == 0:
            return build_assembly_data()
        else:
            return handled.load_to_dict('assembly.csv')
    except FileNotFoundError:
        return build_assembly_data()


def run(id=None):
    assemblies = get_assembly_data()
    _values = [(asm['taxon_id'], asm['accession']) for asm in assemblies]
    sleep(.36)

    for taxon in _values:
        if id is not None:
            if taxon[0] not in id:
                continue
        full_taxon = get_Taxonomy([str(taxon[0])])
        for item in iterate_if_list(full_taxon, None):
            handled.write_taxon(item)
        print(f"{full_taxon=}")
        final_genes_all = taxon_to_genes(*taxon)
        print(f"{final_genes_all=}")
        # get proteins that are polymerases
        final_proteins = match_gene_to_cdd(final_genes_all)
        print(f"{final_proteins=}")
        # finish protein entries
        final_proteins = get_protein(final_proteins)
        for item in iterate_if_list(final_proteins, None):
            handled.write_protein(item)
        # extract gids that correspond to final protein
        gids = [item['gene_id'] for item in final_proteins]
        genes_with_poly = []
        # write relevent genes to file
        for gene in final_genes_all:
            if gene["gene_id"] in gids:
                genes_with_poly.append(gene)
                handled.write_gene(
                    handled.enforce_str_dict(gene)
                )



if __name__ == "__main__":
    # taxon =  ('208964','GCF_000006765.1')
    # final_genes = taxon_to_genes(*taxon)
    # print(json.dumps(get_Taxonomy(['208964']), indent=4))
    # print(repr(get_Taxonomy(['208964'])))
    # # test_value = ['452909', '425885', '179607', '176459', '176458', '176457', '176456',
    # #               '176455', '176454', '176453', '882393']
    # test_genes = ['883126', '883094', '882789', '882652', '882473', '882393',
    #                '882249', '882125', '882038', '881987','880683']
    # test_value3 = ['15595867']
    # t2 = gene_to_protein(test_genes)
    # t4 = protein_to_cdd(t2)
    #
    # t5 = match_gene_to_cdd(genes=test_genes)
    # t6 = get_protein(t5)
    run(test=True)

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

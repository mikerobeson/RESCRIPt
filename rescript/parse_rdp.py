import re
import pandas as pd
from collections import OrderedDict

ALLOWED_CHARS = set(''.join(['0123456789',
                             'abcdefghijklmnopqrstuvwxyz',
                             'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
                             '_-[]()/.\\']))

ALLOWED_RDP_RANKS =  OrderedDict({'rootrank': 'root__', 'domain': 'd__',
                                  'phylum': 'p__', 'class': 'c__',
                                  'subclass': 'cs__', 'order': 'o__',
                                  'suborder': 'os__', 'family': 'f__',
                                  'genus': 'g__'})

DEFAULT_RDP_RANKS = ['domain', 'phylum', 'class', 'subclass', 'order',
                     'suborder', 'family', 'genus']



# def _make_upper_case_sequence(sequences: skbio.DNA) -> DNAFASTAFormat:
#     result = DNAFASTAFormat()
#     with result.open() as out_fasta:
#         for seq in skbio.read(sequences, format='fasta', lowercase=True, constructor=skbio.DNA):
#             seq.write(out_fasta)
#     return result


def parse_rdp(rdp_ref_seqs: skbio.DNA,
              rank_propagation: bool = True,
              ranks: list = None) -> (TSVTaxonomyFormat, DNAFASTAFormat):

        taxd = defaultdict()

        with TSVTaxonomyFormat.open() as rdp_taxonomy,
             DNAFASTAFormat.open() as rdp_seqs:
             # add ability to import lower-case this to q2-types
             # then we can remove the import command below
             # see:
             # https://github.com/qiime2/q2-types/issues/91
             # i.e. just make a transformer to convert lower-case to upper-case
             for seq in skbio.io.read(rdp_ref_seqs, format='fasta',
                                      lowercase=True, constructor=skbio.DNA):
                # write sequence to file
                seq.write(rdp_seqs)
                # parse taxonomy
                seq_id = seq.metadata['id']
                tax_data = seq.metadata['description'].split(
                                                'Lineage=')[1].replace('"','')
                known_lineage = tax_data.rsplit(';')

                if tax_data.endswith(';'):
                    known_lineage = known_lineage[:-1]

                taxd[seq_id] = dict([(x[1], x[0])
                    for x in zip(known_lineage[0::2], known_lineage[1::2])])

            # make dataframe with ranks as column labels.
            df = pd.DataFrame.from_dict(taxd, orient='index',
                                        columns=DEFAULT_RDP_RANKS))
            # propagate ranks (or not), fillna, and append prefixes
            if rank_propagation:
                df.ffill(axis=1, inplace=True)

            df.rename(columns=ALLOWED_RDP_RANKS, inplace=True)

            # In case we decide NOT to ffill. We'll need an empty string, so
            # that we can have empty prefixes (e.g. 'p__') appear in our
            # taxonomy file.Otherwise these will be empty cells in the
            # dataframe.
            df.fillna('', inplace=True

            # use default ranks if no rank list supplied
            if ranks is None:
                ranks = DEFAULT_RDP_RANKS
            else:
                # Sort user ranks relative to allowed_ranks
                # this ensures that the taxonomic ranks supplied by the user
                # are in order
                ranks = [p for r, p in ALLOWED_RDP_RANKS.items()
                                if r in ranks]

            # subset dataframe based on user-selected ranks, and add prefixes
            df.loc[:, ranks.values()] = \
                df.loc[:, ranks.values()].apply(lambda x: x.name + x)

            #df.to_csv(rdp_taxonomy, sep='\t')
            #return series, e.g.
            taxonomy = df.loc[:, ranks].agg('; '.join, axis=1)
            taxonomy.rename('Taxon', inplace=True)

    return (rdp_taxonomy, rdp_seqs)

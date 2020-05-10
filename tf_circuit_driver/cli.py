"""Console script for tf_circuit_driver."""
import click
import sys
import pandas as pd
import numpy as np
import os


class DGE:
    df = None

    def __init__(self, dge_filepath):
        self.df = pd.read_csv(dge_filepath, sep=None, engine="python")

    def gene_score(self, symbol):
        rows = self.df[self.df['symbol'] == symbol].to_dict('records')
        if len(rows) > 0:
            gene = rows[0]
            defaultP = 1e-314
            if gene['pvalue'] == 0:
                raw_score = gene['fc'] * (-1 * np.log10(defaultP))
                return abs(raw_score)
            else:
                raw_score = gene['fc'] * (-1 * np.log10(defaultP))
                return abs(raw_score)
        else:
            return None


class MARA:
    target_dir = None
    df = None
    _tf_symbols = []

    def __init__(self, target_dir):
        self.target_dir = target_dir
        self._build()

    def _build(self):
        file_list = os.listdir(self.target_dir)
        collection = []
        for file in file_list:
            with open(os.path.join(self.target_dir, file), 'r') as f:
                lines = [line.rstrip() for line in f]

                for line in lines:
                    output = self.mara_line_parse(line)

                    for obj in output:
                        collection.append(obj)

        df = pd.DataFrame(collection)
        self.df = df.groupby(['motif', 'target']).sum().reset_index()

    def mara_line_parse(self, line):
        output = []

        contents = line.split('\t')
        targets = [hit.split('|')[1] for hit in contents[3:]]
        target_set = set(targets)
        motifs = contents[2].split('_')
        for motif in motifs:
            self._tf_symbols.append(motif)
            for target in target_set:
                res = dict(score=float(
                    contents[1]), motif=motif, target=target)
                output.append(res)
        return output

    def search(self, symbol):
        return self.df[self.df['motif'] == symbol].to_dict('records')

    def get_tf_symbols(self):
        return list(set(self._tf_symbols))


class PPI:
    mappings = None
    actions = None
    df = None

    def __init__(self, mappings_filepath, actions_filepath):
        self.mappings = pd.read_csv(mappings_filepath, sep="\t")
        actions = pd.read_csv(actions_filepath, sep="\t")
        self.actions = actions[actions['a_is_acting'] == 't']

    def search(self, symbol):
        rows = self.mappings[self.mappings['preferred_name']
                             == symbol].to_dict('records')

        if len(rows) == 0:
            return []
        stringID = rows[0]['protein_external_id']

        ppi_actions = self.actions[self.actions['item_id_a'] == stringID].to_dict(
            'records')
        target_ids = [row['item_id_b'] for row in ppi_actions]
        return self._ids_to_symbol(target_ids)

    def _ids_to_symbol(self, ids):
        collection = []
        for ID in ids:
            rows = self.mappings[self.mappings['protein_external_id'] == ID].to_dict(
                'records')
            for row in rows:
                collection.append(row['preferred_name'])
        return list(set(collection))


class TF_CIRCUIT:
    mara = None
    ppi = None
    dge = None
    tfs = []
    tfs_calculated = {}

    def __init__(self, mappings, actions, targets, dge):
        self.mara = MARA(targets)
        self.ppi = PPI(mappings_filepath=mappings, actions_filepath=actions)
        self.dge = DGE(dge)
        self.tfs = self.mara.get_tf_symbols()

    def _compute_mara_score(self, symbol, include="differential_only"):
        mara_hits = self.mara.search(symbol)
        mara_hit_symbols = [hit['target'] for hit in mara_hits]
        scores = []
        for mara_hit in mara_hits:
            gene_score = self.dge.gene_score(mara_hit['target'])
            if gene_score != None:
                scores.append(gene_score * mara_hit['score'])
            else:
                scores.append(gene_score)

        valid_scores = [score for score in scores if score != None]
        total = sum(valid_scores)
        if len(valid_scores) == 0:
            return dict(score=0, hits=mara_hit_symbols)
        if include == 'all':
            num_hits = len(scores)
            computed = total / num_hits
            return dict(score=computed, hits=mara_hit_symbols)
        else:
            num_hits = len(valid_scores)
            computed = total / num_hits
            return dict(score=computed, hits=mara_hit_symbols)

    def _compute_ppi_score(self, symbol, include="differential_only"):
        ppi_hits = self.ppi.search(symbol)
        scores = []
        for ppi_hit in ppi_hits:
            gene_score = self.dge.gene_score(ppi_hit)
            scores.append(gene_score)
        valid_scores = [score for score in scores if score != None]
        total = sum(valid_scores)
        if len(valid_scores) == 0:
            return dict(score=0, hits=ppi_hits)
        if include == 'all':
            num_hits = len(scores)
            computed = total / num_hits
            return dict(score=computed, hits=ppi_hits)
        else:
            num_hits = len(valid_scores)
            computed = total / num_hits
            return dict(score=computed, hits=ppi_hits)

    def _calculate_level_one(self):
        for tf in self.tfs:
            print("Computing level 1 for {}".format(tf))
            m_1 = self._compute_mara_score(tf)
            p_1 = self._compute_ppi_score(tf)

            self.tfs_calculated[tf] = dict(mara=m_1, ppi=p_1)

    def _calculate_level_two(self):
        for tf in self.tfs:
            print("Computing level 2 for {}".format(tf))
            m_2 = []
            p_2 = []
            for m_hit in self.tfs_calculated[tf]['mara']['hits']:
                if m_hit in self.tfs:
                    m_2.append(self.tfs_calculated[m_hit]['mara'])
            for p_hit in self.tfs_calculated[tf]['ppi']['hits']:
                if p_hit in self.tfs:
                    p_2.append(self.tfs_calculated[p_hit]['ppi'])
            self.tfs_calculated[tf]['level_2'] = dict(mara=m_2, ppi=p_2)

    def _calculate_level_three(self):
        for tf in self.tfs:
            print("Computing level 3 for {}".format(tf))
            m_3 = []
            p_3 = []
            mara_3_hitlist = []
            ppi_3_hitlist = []
            for m_hit in self.tfs_calculated[tf]['level_2']['mara']:
                mara_3_hitlist = mara_3_hitlist + m_hit['hits']
            for p_hit in self.tfs_calculated[tf]['level_2']['ppi']:
                ppi_3_hitlist = ppi_3_hitlist + p_hit['hits']

            for l3_m_hit in list(set(mara_3_hitlist)):
                if l3_m_hit in self.tfs:
                    m_3.append(self.tfs_calculated[l3_m_hit]['mara'])
            for l3_p_hit in list(set(ppi_3_hitlist)):
                if l3_p_hit in self.tfs:
                    p_3.append(self.tfs_calculated[l3_p_hit]['ppi'])
            self.tfs_calculated[tf]['level_3'] = dict(mara=m_3, ppi=p_3)

    def calculate(self):
        self._calculate_level_one()
        self._calculate_level_two()
        self._calculate_level_three()

    def export(self, export_filepath):
        collection = []
        for tf in self.tfs:
            m1 = self.tfs_calculated[tf]['mara']['score']
            m2 = sum([x['score']
                      for x in self.tfs_calculated[tf]['level_2']['mara']]) / 2
            m3 = sum([x['score']
                      for x in self.tfs_calculated[tf]['level_3']['mara']])/3
            p1 = self.tfs_calculated[tf]['ppi']['score']
            p2 = sum([x['score']
                      for x in self.tfs_calculated[tf]['level_2']['ppi']]) / 2
            p3 = sum([x['score']
                      for x in self.tfs_calculated[tf]['level_3']['ppi']])/3
            mara_aggregate = m1 + m2 + m3
            ppi_aggregate = p1 + p2 + p3

            collection.append(dict(tf=tf, mara_aggregate=mara_aggregate,
                                   ppi_aggregate=ppi_aggregate, gene_score=self.dge.gene_score(tf)))
        df = pd.DataFrame(collection)
        df.sort_values(by="mara_aggregate", ascending=False, inplace=True)
        df["MARA_RANK"] = list(range(1, df.shape[0] + 1))
        df.sort_values(by="ppi_aggregate", ascending=False, inplace=True)
        df["PPI_RANK"] = list(range(1, df.shape[0] + 1))
        df.sort_values(by="gene_score", ascending=False, inplace=True)
        df["GENE_RANK"] = list(range(1, df.shape[0] + 1))
        df["AGGREGATE_RANK_SCORE"] = df["MARA_RANK"] + \
            df["PPI_RANK"] + df["GENE_RANK"]
        df.sort_values(by="AGGREGATE_RANK_SCORE", inplace=True)
        df.to_csv(export_filepath, index=None)

    def run(self, export_filepath):
        self.calculate()
        self.export(export_filepath)


@click.command()
@click.option('--actions', help="filepath to the StringDB proteins actions table")
@click.option('--mappings', help="filepath to the StringDB protein mapping")
@click.option('--targets', help="directory path to your unzipped targets directory form the MARA report")
@click.option('--dge', help="filepath to your differential expression file")
@click.option('--output', help="filepath to write output results")
def main(actions, mappings, targets, dge, output):
    circuit = TF_CIRCUIT(mappings=mappings, targets=targets,
                         dge=dge, actions=actions)
    circuit.run(output)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

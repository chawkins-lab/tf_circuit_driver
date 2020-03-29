import pandas as pd
import numpy as np
import os
import click
from tf_circuit_driver.client import driver


def parse_tf_to_target(target_dir):
    click.echo('PARSING TRANSCRIPTION FACTORS TO TARGETS')
    file_list = os.listdir(target_dir)
    items = []
    for file in file_list:
        with open(os.path.join(target_dir, file), 'r') as f:
            for line in f:
                line_contents = line.strip().split('\t')
                motif = line_contents[2]
                tf_array = motif.split('_')
                targets = line_contents[3:]
                target_array = []
                for target in targets:
                    target_array.append(target.split('|')[1])

                for _tf in tf_array:
                    for _target in np.unique(target_array):
                        items.append({
                            "tf": _tf,
                            "target": _target,
                            "score": line_contents[1]
                        })
    counter = 0
    for i in range(round(len(items) / 10000)):
        start = counter
        end = counter + 10000
        batch = items[start:end]
        with driver.session() as session:
            session.run(
                """
                UNWIND $batch as x
                MATCH (p:Protein {displayName: x.tf})
                SET p.is_tf = true
                WITH p, x
                MATCH (t:Protein {displayName: x.target})
                MERGE (p)-[:ACTIVATES_TRANSCRIPTION {score: x.score}]->(t)
                """,
                batch=batch
            )
        counter = counter + 10000

    with driver.session() as session:
        session.run(
            """
            UNWIND $batch as x
            MATCH (p:Protein {displayName: x.tf})
            SET p.is_tf = true
            WITH p, x
            MATCH (t:Protein {displayName: x.target})
            MERGE (p)-[:ACTIVATES_TRANSCRIPTION {score: x.score}]->(t)
            """,
            batch=items[counter:]
        )


def build_protein_nodes(mapping_filepath):
    click.echo('BUILDING PROTEIN NODES')
    mapping = pd.read_table(mapping_filepath, sep="\t", header=0)[
        ['protein_external_id', 'preferred_name']].to_dict('records')
    with driver.session() as session:
        session.run(
            "CREATE CONSTRAINT ON (p:Protein) ASSERT p.stringDbId IS UNIQUE")
        session.run("CREATE INDEX ON :Protein(displayName)")
        session.run(
            "UNWIND $proteins as protein MERGE (p:Protein {stringDbId: protein.protein_external_id}) ON CREATE SET p.displayName = protein.preferred_name",
            proteins=mapping
        )


def parse_dge(dge_filepath):
    dge = pd.read_table(dge_filepath, sep="\t", header=0)
    batch = []
    click.echo('PARSING DIFFERENTIAL EXPRESSION DATA')
    for index, row in dge.iterrows():
        data = row.to_dict()
        absLogFC = abs(data['log_FC'])
        log10pValue = np.log10(data['pvalue']) * -1
        dgeScore = absLogFC * log10pValue
        batch.append({
            "gene_id": data['gene_id'],
            "dgeScore": dgeScore if dgeScore else 0
        })
    with driver.session() as session:
        session.run(
            """
            UNWIND $batch as x
            MATCH (p:Protein {displayName: x.gene_id})
            SET p.dgeScore = x.dgeScore
            """,
            batch=batch
        )
        session.run(
            "MATCH (p:Protein) WHERE p.dgeScore IS NULL SET p.dgeScore = 0")


def build_undirected_ppi(actions_filepath):
    click.echo('BUILDING UNDIRECTED PPIs')
    ppi = pd.read_table(actions_filepath, sep="\t", header=0)
    undirected = ppi[ppi['mode'] == 'binding'][[
        'item_id_a', 'item_id_b', 'mode']].to_dict('records')
    click.echo('Identified {length} directed ppis in manifest'.format(
        length=len(undirected)))
    counter = 0
    for i in range(round(len(undirected) / 10000)):
        start = counter
        end = counter + 10000
        batch = undirected[start:end]
        click.echo('parsing {start} to {end} of {total}'.format(
            start=start, end=end, total=len(undirected)))
        with driver.session() as session:
            session.run("UNWIND $undirected as x MATCH (p1:Protein {stringDbId: x.item_id_a}) MATCH (p2:Protein {stringDbId: x.item_id_b}) MERGE (p1)-[:INTERACTS_WITH {mode: x.mode}]->(p2) MERGE (p2)-[:INTERACTS_WITH {mode: x.mode}]->(p1)",
                        undirected=batch)
        counter = counter + 10000
    with driver.session() as session:
        session.run("UNWIND $undirected as x MATCH (p1:Protein {stringDbId: x.item_id_a}) MATCH (p2:Protein {stringDbId: x.item_id_b}) MERGE (p1)-[:INTERACTS_WITH {mode: x.mode}]->(p2) MERGE (p2)-[:INTERACTS_WITH {mode: x.mode}]->(p1)",
                    undirected=undirected[counter:])
        click.echo('undirected ppi generation completed')


def build_directed_ppi(actions_filepath):
    click.echo('BUILDING DIRECTED PPIs')
    ppi = pd.read_table(actions_filepath, sep="\t", header=0)
    directed = ppi[ppi['a_is_acting'] == 't'][[
        'item_id_a', 'item_id_b', 'mode']].to_dict('records')
    click.echo('Identified {length} directed ppis in manifest'.format(
        length=len(directed)))
    counter = 0
    for i in range(round(len(directed) / 10000)):
        start = counter
        end = counter + 10000
        click.echo('parsing {start} to {end} of {total}'.format(
            start=start, end=end, total=len(directed)))
        batch = directed[start:end]
        with driver.session() as session:
            session.run("UNWIND $directed as x MATCH (p1:Protein {stringDbId: x.item_id_a}) MATCH (p2:Protein {stringDbId: x.item_id_b}) MERGE (p1)-[:ACTS_ON]->(p2)",
                        directed=batch)
        counter = counter + 10000
    with driver.session() as session:
        session.run("UNWIND $directed as x MATCH (p1:Protein {stringDbId: x.item_id_a}) MATCH (p2:Protein {stringDbId: x.item_id_b}) MERGE (p1)-[:ACTS_ON]->(p2)",
                    directed=directed[counter:])
        click.echo('directed ppi generation completed')


def bootstrap(actions, mapping, targets, dge):
    click.echo('INITIALIZING DB BOOTSTRAPPING')
    build_protein_nodes(mapping)
    build_undirected_ppi(actions)
    build_directed_ppi(actions)
    parse_tf_to_target(targets)
    parse_dge(dge)

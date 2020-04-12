"""CALCULATE.py"""
from tf_circuit_driver.client import driver
from tf_circuit_driver.models.tf import TranscriptionFactor
import pandas as pd
import click


def calculate(output):
    tf_table_rows = []
    with driver.session() as session:
        result = session.run(
            "MATCH (p:Protein) WHERE p.is_tf = true AND p.dgeScore > 0 RETURN properties(p) as tf")

        for record in result:
            tf = TranscriptionFactor(record['tf'])

            tf.calc_self_mara()
            tf.calc_self_string()
    with driver.session() as session:
        result = session.run(
            "MATCH (p:Protein) WHERE p.is_tf = true AND p.dgeScore > 0 RETURN properties(p) as tf")
        for record in result:
            click.echo("Constructing {} circuit".format(
                tf.properties['displayName']))
            tf = TranscriptionFactor(record['tf'])
            click.echo("Calculating Protein to DNA network")
            mara = tf.calc_mara()
            click.echo("Calculating Protein to Protein network")
            string = tf.calc_string()

            row = {
                "tf": tf.properties['displayName'],
                "geneScore": tf.properties['dgeScore']
            }
            row.update(mara)
            row.update(string)
            tf_table_rows.append(row)
    df = pd.DataFrame(data=tf_table_rows)
    df['mara_aggregate'] = df['level1_mara'] + \
        df['level2_mara'] + df['level3_mara']
    df['string_aggregate'] = df['level1_string'] + \
        df['level2_string'] + df['level3_string']
    df['mara_rank'] = pd.Series(index=df.index)
    df['string_rank'] = pd.Series(index=df.index)
    df['dge_rank'] = pd.Series(index=df.index)

    mara_rank = list(df['mara_aggregate'].sort_values(ascending=False).index)
    floor_mara_rank = len(df[df['mara_aggregate'] > 0].index)
    for index, position in enumerate(mara_rank):
        if index == floor_mara_rank:
            df['mara_rank'][position] = floor_mara_rank
        else:
            df['mara_rank'][position] = index + 1
    string_rank = list(
        df['string_aggregate'].sort_values(ascending=False).index)
    floor_string_rank = len(df[df['string_aggregate'] > 0].index)
    for index, position in enumerate(string_rank):
        if index == floor_string_rank:
            df['string_rank'][position] = floor_string_rank
        else:
            df['string_rank'][position] = index + 1

    dge_rank = list(df['geneScore'].sort_values(ascending=False).index)
    floor_dge_rank = len(df[df['geneScore'] > 0].index)
    for index, position in enumerate(dge_rank):
        if index == floor_dge_rank:
            df['dge_rank'][position] = floor_dge_rank
        else:
            df['dge_rank'][position] = index + 1

    df['aggregate_rank'] = df['mara_rank'] + df['string_rank'] + df['dge_rank']

    df['RANK'] = pd.Series(index=df.index)
    RANK = list(df['aggregate_rank'].sort_values(ascending=True).index)
    for index, position in enumerate(RANK):
        df['RANK'][position] = index + 1
    df.set_index('RANK', inplace=True)
    df.to_csv(output)

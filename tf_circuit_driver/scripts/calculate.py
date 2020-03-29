"""CALCULATE.py"""
from tf_circuit_driver.client import driver
from tf_circuit_driver.models.tf import TranscriptionFactor
import pandas as pd


def calculate(output):
    tf_table_rows = []
    with driver.session() as session:
        result = session.run(
            "MATCH (p:Protein) WHERE p.is_tf = true RETURN properties(p) as tf")
        for record in result:
            tf = TranscriptionFactor(record['tf'])
            tf.calc_self_mara()
            tf.calc_self_string()
    with driver.session() as session:
        result = session.run(
            "MATCH (p:Protein) WHERE p.is_tf = true RETURN properties(p) as tf")
        for record in result:
            tf = TranscriptionFactor(record['tf'])
            mara = tf.calc_mara()
            string = tf.calc_string()

            row = {
                "tf": tf.properties['displayName'],
                "geneScore": tf.properties['dgeScore']
            }
            row.update(mara)
            row.update(string)
            tf_table_rows.append(row)
    df = pd.DataFrame(data=tf_table_rows)
    df.to_csv(output)

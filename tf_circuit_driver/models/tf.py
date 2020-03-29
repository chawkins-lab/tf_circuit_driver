from tf_circuit_driver.client import driver


class TranscriptionFactor:
    driver = driver
    properties = {}

    def __init__(self, properties):
        self.properties = properties

    def calc_self_string(self):
        with self.driver.session() as session:
            level1 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                CALL apoc.neighbors.athop(p, "ACTS_ON>|INTERACTS_WITH>", 1)
                YIELD node
                WHERE node.dgeScore > 0
                RETURN node

                """,
                stringDbId=self.properties['stringDbId']
            )

            values = [record['node']['dgeScore'] for record in level1]
            stringScore = 0
            if len(values) > 0:
                stringScore = sum(values)/len(values)
            session.run(
                "MATCH (p:Protein {stringDbId: $stringDbId}) SET p.stringScore = $stringScore",
                stringDbId=self.properties['stringDbId'], stringScore=stringScore
            )

    def calc_self_mara(self):
        with self.driver.session() as session:
            level1 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                CALL apoc.neighbors.athop(p, "ACTIVATES_TRANSCRIPTION>", 1)
                YIELD node
                WHERE node.dgeScore > 0
                RETURN node
                """,
                stringDbId=self.properties['stringDbId']
            )

            values = [record['node']['dgeScore'] for record in level1]
            maraScore = 0
            if len(values) > 0:
                maraScore = sum(values)/len(values)

            session.run(
                "MATCH (p:Protein {stringDbId: $stringDbId}) SET p.maraScore = $maraScore",
                stringDbId=self.properties['stringDbId'], maraScore=maraScore
            )

    def calc_mara(self):
        output = {
            'level1_mara': 0,
            'level2_mara': 0,
            'level3_mara': 0
        }
        with self.driver.session() as session:
            level1 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                RETURN properties(p) as node
                """,
                stringDbId=self.properties['stringDbId']
            )
            level1Values = [record['node']['maraScore'] for record in level1]
            output['level1_mara'] = 0
            if len(level1Values) > 0:
                output['level1_mara'] = sum(level1Values) / len(level1Values)
            level2 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                MATCH (p)-[:ACTIVATES_TRANSCRIPTION]->(p2:Protein)
                WHERE p2.is_tf = true
                RETURN properties(p2) as node
                """,
                stringDbId=self.properties['stringDbId']
            )

            level2Values = [record['node']['maraScore'] for record in level2]
            output['level2_mara'] = 0
            if len(level2Values) > 0:
                output['level2_mara'] = (
                    sum(level2Values) / len(level2Values)) / 2
            level3 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                MATCH (p)-[:ACTIVATES_TRANSCRIPTION]->(p2:Protein)
                WHERE p2.is_tf = true
                MATCH (p2)-[:ACTIVATES_TRANSCRIPTION]->(p3:Protein)
                WHERE p3.is_tf = true
                RETURN properties(p3) as node
                """,
                stringDbId=self.properties['stringDbId']
            )
            level3Values = [record['node']['maraScore'] for record in level3]
            output['level3_mara'] = 0
            if len(level3Values) > 0:
                output['level3_mara'] = (
                    sum(level3Values) / len(level3Values)) / 3
        return output

    def calc_string(self):
        output = {
            'level1_string': 0,
            'level2_string': 0,
            'level3_string': 0
        }
        with self.driver.session() as session:
            level1 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                RETURN p as node
                """,
                stringDbId=self.properties['stringDbId']
            )
            level1Values = [record['node']['stringScore'] for record in level1]
            output['level1_string'] = sum(level1Values) / len(level1Values)
            level2 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                MATCH (p)-[:ACTS_ON|:INTERACTS_WITH]->(p2:Protein)
                WHERE p2.is_tf = true
                RETURN properties(p2) as node
                """,
                stringDbId=self.properties['stringDbId']
            )
            level2Values = [record['node']['stringScore'] for record in level2]
            output['level2_string'] = 0
            if len(level2Values) > 0:
                output['level2_string'] = (
                    sum(level2Values) / len(level2Values)) / 2
            level3 = session.run(
                """
                MATCH (p:Protein {stringDbId: $stringDbId})
                MATCH (p)-[:ACTS_ON|:INTERACTS_WITH]->(p2:Protein)
                WHERE p2.is_tf = true
                MATCH (p2)-[:ACTS_ON|:INTERACTS_WITH]->(p3:Protein)
                WHERE p3.is_tf = true
                RETURN properties(p3) as node
                """,
                stringDbId=self.properties['stringDbId']
            )

            level3Values = [record['node']['stringScore'] for record in level3]
            output['level3_string'] = 0
            if len(level3Values) > 0:
                output['level3_string'] = (
                    sum(level3Values) / len(level3Values)) / 3
        return output

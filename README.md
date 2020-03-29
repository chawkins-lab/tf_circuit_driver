# tf_circuit_driver

```{sh}
docker-compose up -d neo4j
pip install --editable .
tf_circuit_driver --mode init \
--actions path/to/actions/file \
--mappings path/to/mappings/file \
--dge path/to/dge \
--targets path/to/targets/dir/

tf_circuit_driver --mode calculate --output ~/test_out.csv
```

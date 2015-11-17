spark-submit --executor-cores 5 --num-executors 100 --executor-memory 32g --driver-memory 6g --conf spark.yarn.executor.memoryOverhead=8000  --conf "spark.executor.extraJavaOptions=-XX:-UseGCOverheadLimit"  dense-subgraph-computation_2.10-1.0.jar $1 $2 $3
Rscript UndirectedDense.R $3

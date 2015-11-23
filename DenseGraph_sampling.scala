import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import org.apache.spark.HashPartitioner
import org.apache.spark.rdd.RDD
import org.apache.spark.storage.StorageLevel
//
import org.apache.spark.graphx.Edge
import org.apache.spark.graphx.EdgeRDD
import org.apache.spark.graphx.Graph
import org.apache.spark.graphx.GraphLoader
import org.apache.spark.graphx.PartitionStrategy
import org.apache.spark.graphx.VertexRDD
import math.{Pi,log1p}
import java.io._

object denseGraph {
  def main(args: Array[String]) {
  val conf = new SparkConf().setAppName("Compute Dense Graph stats")
  if (args.length < 1) {
        System.err.println("Usage:  <epsilon> <delta> <pathtographfie>")
        System.exit(1)
      }
  	
  val cacheType = StorageLevel.MEMORY_AND_DISK 
  val parts = 1000
  val epsilon = args(0).toDouble
  val delta = args(1).toDouble
  val fileloc = args(2)
  val sc = new SparkContext(conf)
val sqlContext = new org.apache.spark.sql.hive.HiveContext(sc)
val edgeqry = sqlContext.sql("SELECT src,dest from bitcoin_graph")
  val edges = edgeqry.map(ids => (ids(0).toString.trim.toLong, ids(1).toString.trim.toLong)).map(e => Edge(e._1, e._2, 0))
  val edgeRDD = EdgeRDD.fromEdges(edges)
  val graph = Graph.fromEdges(edgeRDD, defaultValue=0, cacheType, cacheType).partitionBy(PartitionStrategy.RandomVertexCut, parts)
  val n = graph.vertices.count
  val m = graph.edges.count

  var C = (12.0*n*(4.0+delta)*log1p(m)/(epsilon*epsilon)).toInt

  if(C >10000)C=10000
  val samplededges = edges.takeSample(true,C)
  val file = new File(fileloc)
  val writer = new BufferedWriter(new FileWriter(file))
  val edrdd:RDD[String] = edgeqry.map(x=> "%s %s\n".format(x(0).toString  , x(1).toString ))
  edrdd.takeSample(true,C).foreach(writer.write)
  writer.close
  
  
  
}  
  
  }

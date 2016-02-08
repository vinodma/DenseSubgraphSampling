//Implements weighted node sampling using page rank score as the weights

import scala.util.Random
import org.apache.spark.graphx._
import org.apache.spark.graphx.util._

//Example usage

val g = GraphGenerators.logNormalGraph(sc,1000)

val g1=sampleWeightedGraph(g,.3)
/////////
def sampleWeightedGraph(g:Graph[Long,Int], samplerate:Double):Graph[Double,Double] = {

val samplesize:Int = (samplerate*g.vertices.count).toInt
val pg = g.pageRank(.01)

val ki = pg.mapVertices((id,attr)=>Math.pow(Random.nextDouble(),1.0/attr))

val cutoffval = ki.vertices.map(x=>(x._1,x._2)).takeOrdered(samplesize)(Ordering[Double].reverse.on(x=>x._2)).map(a=>a._2).reduce((a,b)=>if(a < b ) a else b )

return(ki.subgraph(vpred=(id,attr)=> attr > cutoffval))

}

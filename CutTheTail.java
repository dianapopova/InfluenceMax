import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.Random;
import java.util.Arrays;

import java.io.IOException;

import it.unimi.dsi.webgraph.ImmutableGraph;
import it.unimi.dsi.webgraph.IncrementalImmutableSequentialGraph;
import it.unimi.dsi.webgraph.BVGraph;

import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.IntStream;
import java.util.concurrent.Callable;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.BufferedOutputStream;

import java.util.Iterator;


public class CutTheTail {
	ImmutableGraph G, GO;
    int n, k;
	long m, sketch_num;
	double eps, p, c; // given error allowed
    
    long[] node_infl;
    
    String basename;
    
    int out80; // outdegree in 80th percentile in original (not transpose) graph GO
    
    int node_tail; // the max value of total_degree of the dropped sketches
	
	public CutTheTail(String basename, double  p, double eps, int k) throws Exception {
		long readFile = System.currentTimeMillis();
        G = ImmutableGraph.load(basename+"t");  //  Transpose graph
        
        double fileTime = (System.currentTimeMillis() - readFile)/1000.0;
        
        System.out.println();
        System.out.println("Reading transpose graph took " + fileTime + " seconds");
		
		n = G.numNodes();
		m = G.numArcs();
        
        GO = ImmutableGraph.load(basename); //Origal (not transpose) graph
        
        out80 = outdegree80(GO.copy());
        
        c = 4.0*(1 + eps)*(1 + 1/k); //constant for full Borgs et al. formula
        
        node_tail = out80; // nodes in GO - in 80th percentile of outdegree
        
        
        System.out.println();
        System.out.println("node tail = " + node_tail);
        System.out.println();
        
		System.out.println("n="+n + ", m=" +m  + ", R=" + ( c * m * k * (1/Math.pow(eps, 2)) * Math.log(n) ) );
        System.out.println("k = " + k);
        System.out.println("p = " + p);
        System.out.println("eps = " + eps);
        
        this.p = p;
        this.c = c;
        this.eps = eps;
        this.basename = basename;
        this.k = k;
		
    }


    public int outdegree80(ImmutableGraph GO) {
        Iterator<Integer> degIter = GO.outdegrees();
        int[] degArr = new int[n];
        int i = 0;
        
        while (degIter.hasNext()) {
            Integer deg = degIter.next();
            degArr[i] = deg;
            i++;
            if(i == n-1)
                break;
        }
        
        Arrays.sort(degArr);
        
        int nodes80 = (n * 8)/10;

        return degArr[nodes80]; // Value of outdegree in GO that cuts out 80% of nodes
    }

 
    void get_sketch() throws Exception {
        long getHypergraph = System.currentTimeMillis();
        
        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("int cores = " + cores);
        
        
		double R = ( c * m * k * (1/Math.pow(eps, 2)) * Math.log(n) )/cores;
        
        ImmutableGraph G_sk;
  
        long node_infl[][] = new long[cores][];
        for(int i=0; i<cores; i++)
            node_infl[i] = new long[n+1];
        
        IntStream.range(0, cores)
        .parallel()
        .forEach(core -> {
            try {
                node_infl[core] = get_sketch(G.copy(),GO.copy(),core,R);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
        
        
        /////////////////////////////////////////////////////////////////////////////
        
        //Merge the g graphs produced by the cpu cores:
        
        final IncrementalImmutableSequentialGraph g = new IncrementalImmutableSequentialGraph();
        ExecutorService executor = Executors.newSingleThreadExecutor();
        final Future<Void> future = executor.submit( new Callable<Void>() {
            public Void call() throws IOException {
                BVGraph.store( g, G.basename()+"-sketches");
                return null;
            }
        } );
        
        for(int core=0; core<cores; core++) {
            ImmutableGraph h = ImmutableGraph.load(G.basename()+"-sketch"+core);
            
            for(int i=0;i<h.numNodes();i++) {
                int[] sketch_vertex_list = h.successorArray(i);
                g.add(sketch_vertex_list, 0, h.outdegree(i));
            }
//            System.out.println("Core = " + core + "; Number of sketches saved = " + h.numNodes() );
            
        }
        
        g.add( IncrementalImmutableSequentialGraph.END_OF_GRAPH );
        future.get();
        executor.shutdown();
        
        double gTime = (System.currentTimeMillis() - getHypergraph)/1000.0;
        
        System.out.println();
        System.out.println("Creating hypergraph took " + gTime + " seconds");
        System.out.println();
        
        
        System.gc();
        
        /////////////////////////////////////////////////////////////////////////////

        
        // Getting the seeds:
        
        //Calculate the number of sketches and Merge the node_infl[] arrays produced by the cpu cores:
        sketch_num = 0;
        for(int ii = 0; ii<cores;ii++)
            sketch_num = sketch_num + node_infl[ii][0];
        
        System.out.println("Total number of sketches = " + sketch_num );
        
        double coeff = 1.0 * n/sketch_num;
        
        long node_cover[] = new long[n+1];
        
        for(int j=1;j<(n+1);j++) {
            for(int i = 0; i<cores;i++)
                node_cover[j-1] = node_cover[j-1] + node_infl[i][j];
        }
        
        node_cover[n] = sketch_num;
        System.out.println("Writing nodes influence into a binary file...");
        writeLongArrayToBinaryFile(node_cover, G.basename()+"-sketches.cover");
        
        
        G_sk = ImmutableGraph.load(G.basename()+"-sketches");
        
        int hyper_rows = G_sk.numNodes();
        
        System.out.println(
                           "Hypergraph rows = " + hyper_rows);
        System.out.println();
        
        double total_infl = 0.0;
        BitSet sk_gone = new BitSet(hyper_rows);
        
        long startSeeds = System.currentTimeMillis();
        
        get_seeds(G_sk.copy(), node_cover, k, hyper_rows, total_infl, sk_gone, coeff);
        
        double getSeeds = (System.currentTimeMillis() - startSeeds)/1000.0;
        
        System.out.println();
        System.out.println("Calculating seeds took " + getSeeds + " seconds");
        
    }
    
    
    long[] get_sketch(ImmutableGraph G, ImmutableGraph GO, int core, double R) throws Exception {
        
        double weight_of_current_index = 0.0;
            
        Random gen_rnd = new Random();
            
        //Let's create the new sketch "hypergraph" (Borgs et al. term)
        final IncrementalImmutableSequentialGraph g = new IncrementalImmutableSequentialGraph();
        ExecutorService executor = Executors.newSingleThreadExecutor();
        final Future<Void> future = executor.submit( new Callable<Void>() {
            public Void call() throws IOException {
                BVGraph.store( g, G.basename()+"-sketch"+core );
                return null;
            }
        } );
            
        long[] node_infl = new long[n+1];
        
        long num_sketches = 0;
            
        BitSet marked = new BitSet(n);
            
        while(weight_of_current_index < R) {
                
            int v = gen_rnd.nextInt(n);
            marked.clear();
            BFS(G.copy(),v,marked);  // All nodes for the sketch found
            
            num_sketches++;
            
            int[] sketch_vertex_list = new int[marked.cardinality()];
            int j = 0;
            int total_out_degree = 0;  // The number of edges "touched" by BFS in transpose graph
            int sk_outdegree = 0;  // The sum of node degrees in the original graph, for this sketch
                
            for (int u = marked.nextSetBit(0); u >= 0; u = marked.nextSetBit(u+1)) {
                
                sketch_vertex_list[j] = u;
                node_infl[u+1] = node_infl[u+1] + 1;
                    
                j++;
                    
                total_out_degree += G.outdegree(u);
                sk_outdegree += GO.outdegree(u);
            }
            
            weight_of_current_index += total_out_degree;

            // Do not store the sketches with just one node or the sketches with total_out_degree less than node_tail.
            if(marked.cardinality() < 2 || sk_outdegree <= node_tail) {
                continue;
            }
            else {
                g.add(sketch_vertex_list, 0, marked.cardinality());
            }
        }
        
        node_infl[0] = num_sketches;// Real number of sketches computed - including the ones with just one node or total_degree less than tail.
        
        g.add( IncrementalImmutableSequentialGraph.END_OF_GRAPH );
        future.get();
        executor.shutdown();
            
        return (node_infl);
    }
    
    
    void BFS(ImmutableGraph G, int v, BitSet marked) {
        
        Random random = new Random();
        
        Deque<Integer> queue = new ArrayDeque<Integer>();
        
        queue.add(v);
        marked.set(v);
        
        while (!queue.isEmpty()) {
            int u = queue.remove();
            int[] u_neighbors = G.successorArray(u);
            int u_deg = G.outdegree(u);
            
            for (int ni = 0; ni < u_deg; ni++) {
                int uu = u_neighbors[ni];
                double xi = random.nextDouble();
                
                if (!marked.get(uu) && xi < p) {
                    queue.add(uu);
                    marked.set(uu);
                }
            }
        }
    }

    public static void writeLongArrayToBinaryFile(long[] arr, String filename) throws Exception {
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename)));
        
        for (int i=0; i<arr.length; i++)
            dos.writeLong(arr[i]);
        
        dos.close();
    }
    
    void get_seeds(ImmutableGraph G_sk, long[] node_infl, int k_left, int hyper_rows, double total_infl, BitSet sk_gone, double coeff) throws Exception {
        
        // Calculating the node with max influence
        
        int max_node = 0;
        int min_sketch = 1;
        long infl_max = 0;
        
        for(int v=0;v<n;v++)
        {
            if(node_infl[v] < min_sketch)
                continue;
            else
            {
                if(node_infl[v] > infl_max)
                {
                    max_node = v;
                    infl_max = node_infl[v];
                }
            }
        }
        
        //        System.out.println("Maximum Influence = " + infl_max);
        
        total_infl =  total_infl + infl_max * coeff;
        
        System.out.println(max_node);
        
//        System.out.println( "Max Node: " + max_node);
//                           ", Set Influence = " + total_infl );
//        System.out.println();
        
        // Stopping condition: no need to re-calculate the influence, if we already got the k seeds
        if((k_left - 1)==0)
        {
            System.out.println();
            System.out.println( "Set Influence = " + total_infl );
            return;
        }
        
        // Remove max node and the sketches it participated in
        // plus re-calculate the influence of the remaining nodes
        
        //        int jj=0;
        
        for(int s=0;s<hyper_rows;s++)
        {
            if(sk_gone.get(s))
                continue;
            else
            {
                int[] s_nodes = G_sk.successorArray(s);
                int nodes_sk = s_nodes.length;
                
                for(int i=0;i<nodes_sk;i++)
                {
                    if(s_nodes[i] == max_node)
                    {
                        for(int ii=0;ii<nodes_sk;ii++)
                        {
                            if(node_infl[s_nodes[ii]] > 0)
                                node_infl[s_nodes[ii]] = node_infl[s_nodes[ii]] - 1;
                        }
                        sk_gone.set(s);
                        //                                jj++;
                        break;
                    }
                }
            }
        }
        //        }
        //        System.out.println("Sketches with max_node found:" + jj + "; when max_node is " + max_node + " with infl_max = " + infl_max);
        node_infl[max_node] = 0;

        get_seeds(G_sk.copy(), node_infl, k_left-1, hyper_rows, total_infl, sk_gone, coeff);
        
    }


		
	public static void main(String[] args) throws Exception {
		long startTime = System.currentTimeMillis();
		long estimatedTime;
        
        if(args.length < 4) {
            System.out.println("Specify: basename, p, eps, k");
            System.exit(1);
        }
		
        String basename  = args[0];
        double p = Double.valueOf(args[1]);
        double eps = Double.valueOf(args[2]);
        int k = Integer.valueOf(args[3]);
        
        CutTheTail imw = new CutTheTail(basename, p, eps, k);
		imw.get_sketch();
        
        System.gc();
        
//        double appr = 1 - 1/Math.E - eps;
//        System.out.println("Approximation to OPT = " + appr + " with the probability of at least 60%");

			
		estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Total time elapsed = " + estimatedTime / 1000.0 + " seconds");
	}
}

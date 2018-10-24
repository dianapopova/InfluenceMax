import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.Random;
import java.util.Arrays;

import java.lang.Math;

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


public class CutTheTail2 {
	ImmutableGraph G;
    int n, k, sk_limit;
	long m, sketch_num;
	double eps, p, c; // given error allowed and the probability of edge existence
    
    long[] node_infl;
    
    String basename;
    
    int sk_tail; // calculated minimum cardinality for a sketch to be saved
	
	public CutTheTail2(String basename, double  p, double eps, int k) throws Exception {
		long readFile = System.currentTimeMillis();
        G = ImmutableGraph.load(basename+"t");  //  Transpose graph
        
        double fileTime = (System.currentTimeMillis() - readFile)/1000.0;
        
        System.out.println();
        System.out.println("Reading transpose graph took " + fileTime + " seconds");
		
		n = G.numNodes();
		m = G.numArcs();
        sk_limit = Math.max(n/100, 10000); //number of samples to take for cardinality calculation
        c = 4.0*(1 + eps)*(1 + 1/k); //constant for full Borgs et al. formula
        
        Random card_rnd = new Random();
        
        BitSet marked = new BitSet(n);
        int sk_taken = 0;
        
        int[] sk_card = new int[sk_limit];
        
        int j=0;
        while(sk_taken < sk_limit) {
            
            int v = card_rnd.nextInt(n);
            marked.clear();
            
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

            sk_card[j] = marked.cardinality(); //saving cardinality of the sketch in sk_card array
            
            sk_taken++;
            j++;
        }
        
        int max_card = 1;
        for(int i=0;i<sk_limit;i++) {
            if(sk_card[i] > max_card) {
                max_card = sk_card[i];
            }
        }
        
        System.out.println("Maximum cardinality = " + max_card);
        
//        if(max_card > 20) {
            double k_log = Math.log(k);
        
            int card = max_card / (int) Math.round(k_log); //The larger the number of seeds to compute, the more sketches we need to save
        
            int drop_card = Math.min(card/10, 100);
        
            sk_tail = Math.max( drop_card, 2 );
/*        }
        else {
            sk_tail = 1;
        }
*/        
        System.out.println();
        System.out.println("Sketch tail (sketches longer than the tail are saved) = " + sk_tail);
        System.out.println();
        
		System.out.println("n="+n + ", m=" +m  + ", R=" + ( c * m * k * (1/Math.pow(eps, 2)) * Math.log(n) ) );
        System.out.println("k = " + k);
        System.out.println("p = " + p);
        System.out.println("eps = " + eps);
        
        System.gc();
        
        this.p = p;
        this.c = c;
        this.eps = eps;
        this.basename = basename;
        this.k = k;
		
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
                node_infl[core] = get_sketch(G.copy(),core,R);
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
    
    
    long[] get_sketch(ImmutableGraph G, int core, double R) throws Exception {
        
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
            
            for (int u = marked.nextSetBit(0); u >= 0; u = marked.nextSetBit(u+1)) {
                
                sketch_vertex_list[j] = u;
                node_infl[u+1] = node_infl[u+1] + 1;
                    
                j++;
                    
                total_out_degree += G.outdegree(u);
            }
            
            weight_of_current_index += total_out_degree;

            // Do not store the sketches with cardinality less than tail.
            if(marked.cardinality() <= sk_tail) {
                continue;
            }
            else {
                g.add(sketch_vertex_list, 0, marked.cardinality());
            }
        }
        
        node_infl[0] = num_sketches;// Real number of sketches computed - including the ones with number of nodes less than tail.
        
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
        
        CutTheTail2 imw = new CutTheTail2(basename, p, eps, k);
		imw.get_sketch();
        
        System.gc();
        
//        double appr = 1 - 1/Math.E - eps;
//        System.out.println("Approximation to OPT = " + appr + " with the probability of at least 60%");

			
		estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Total time elapsed = " + estimatedTime / 1000.0 + " seconds");
        
	}
}

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Arrays;

import java.util.Comparator;
import java.util.PriorityQueue;

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



public class NoSinglesTopNodes {
	ImmutableGraph G, G_sk;
	Map<Integer, List<Integer>> node_sketches_map = new HashMap<>(); //partial graph
    int n, k;
    long m, sketch_num;
    double p, beta; // given error allowed
    long[] node_infl, node_cover;
    
    String basename;
    int top;
    
    
    public NoSinglesTopNodes(String basename, double  p, int k, double beta) throws Exception {
        long readFile = System.currentTimeMillis();
        G = ImmutableGraph.load(basename+"t"); //  Transpose graph
        
        double fileTime = (System.currentTimeMillis() - readFile)/1000.0;
        
        System.out.println();
        System.out.println("Reading transpose graph took " + fileTime + " seconds");
        
        n = G.numNodes();
        m = G.numArcs();
        top = (int) Math.round(0.2*n);
        
        System.out.println("n="+n + ", m=" +m  + ", top=" + top +
                           ", R=" + (beta*m*k*Math.log(n) ));
        System.out.println("beta = " + beta);
        System.out.println("k = " + k);
        System.out.println("p = " + p);

        this.p = p;
        this.beta = beta;
        this.basename = basename;
        this.k = k;
        this.top = top;
        
    }
    
	
	void get_sketch() throws Exception {
        long getHypergraph = System.currentTimeMillis();
        
        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("int cores = " + cores);
        
        double R = beta*m*k*Math.log(n)/cores;
        
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
        System.out.println("Creating NoSingles hypergraph took " + gTime + " seconds");
        System.out.println();
        
        System.gc();
        
        /////////////////////////////////////////////////////////////////////////////
        
        //Calculate the number of sketches and Merge the node_infl[] arrays produced by the cpu cores:
        
        sketch_num = 0;
        for(int ii = 0; ii<cores;ii++) {
            sketch_num = sketch_num + node_infl[ii][0];
            //            System.out.println("Number of sketches taken by core" + ii + " = " + node_infl[ii][0] );
        }
        
        System.out.println("Total number of sketches = " + sketch_num );
        
        double coeff = 1.0 * n/sketch_num;
        
        node_cover = new long[n+2];
        
        for(int j=1;j<(n+1);j++) {
            for(int i = 0; i<cores;i++)
                node_cover[j-1] = node_cover[j-1] + node_infl[i][j];
        }
        
        node_cover[n] = sketch_num;
        node_cover[n+1] = n;
        
        //        System.out.println("Writing nodes influence into a binary file...");
        writeLongArrayToBinaryFile(node_cover, G.basename()+"-sketches.cover");
        
        G_sk = ImmutableGraph.load(G.basename()+"-sketches"); //Load hypergraph
        
        int hyper_rows = G_sk.numNodes();
        long hyper_size = G_sk.numArcs();
        
        System.out.println("NoSingles Hypergraph rows = " + hyper_rows +
                           "; hypergraph size (number of arcs) = " + hyper_size);
        System.out.println();

        long startTop = System.currentTimeMillis();
        
        //Create a priority queue for the nodes by their influence in node_cover
        
        PriorityQueue<NodeInfluence> pq  = new PriorityQueue<NodeInfluence>(n, new NodeInfluenceComparator());
        PriorityQueue<NodeInfluence> pq2 = new PriorityQueue<NodeInfluence>(n, new NodeInfluenceComparator());
        //Alex: pq2 is just a copy of pq. Needed for creating the partial transpose for the top nodes.
        
        for(int ii=0;ii<n;ii++) {
            pq. add( new NodeInfluence(ii,node_cover[ii]) );
            pq2.add( new NodeInfluence(ii,node_cover[ii]) );
        }
        
        // Calling the get_top method:
        
        node_sketches_map = get_top(pq2, hyper_rows, G_sk.copy());
        
        double getTop = (System.currentTimeMillis() - startTop)/1000.0;
        System.out.println();
        System.out.println("Calculating node_sketches_map took " + getTop + " seconds");
        
        // Calling the get_seeds method:
        
        long startSeeds = System.currentTimeMillis();
        
        double total_infl = 0.0;
        BitSet sk_gone = new BitSet(hyper_rows);
        BitSet reCalc = new BitSet(n);
        
        get_seeds(G_sk.copy(), node_sketches_map, node_cover, k, hyper_rows, total_infl, sk_gone, reCalc, pq, pq2, coeff);
        
        double getSeeds = (System.currentTimeMillis() - startSeeds)/1000.0;
        System.out.println();
        System.out.println("Calculating seeds took " + getSeeds + " seconds");
        
    }
    
    // Scan hypergraph G_sk and put, for each node in  top 20% pq, sketchID into the node_sketches_map:
    
    Map<Integer, List<Integer>> get_top(PriorityQueue<NodeInfluence> pq2, int hyper_rows, ImmutableGraph G_sk) {
        
        //Create a partial hypergraph {nodeID:sketchIDs}, only for the top nodes:
        node_sketches_map.clear();
        
        for(int cnt=0; cnt<top; cnt++) {
        	NodeInfluence ni = pq2.poll();
        	int nodeId = ni.node;
        	node_sketches_map.put(nodeId, new ArrayList<Integer>());
        }
        
        for(int sk=0; sk<hyper_rows; sk++) {
        	int[] sk_nodes = G_sk.successorArray(sk);
        	int sk_deg = G_sk.outdegree(sk);
        	
        	for(int i=0; i<sk_deg; i++) {
        		int nodeId = sk_nodes[i];
        		if(node_sketches_map.containsKey(nodeId))
        			node_sketches_map.get(nodeId).add(sk);
        	}
        }
        return node_sketches_map;
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
            BFS(G.copy(),v,marked);
            
            num_sketches++;
            
            int[] sketch_vertex_list = new int[marked.cardinality()];
            int j = 0;
            int total_out_degree = 0;
            
            for (int u = marked.nextSetBit(0); u >= 0; u = marked.nextSetBit(u+1)) {
                
                sketch_vertex_list[j] = u;
                node_infl[u+1] = node_infl[u+1] + 1;
                
                j++;
                
                total_out_degree += G.outdegree(u);
            }
            
            // Do not store the sketches with just one node.
            if(marked.cardinality() > 1)
                g.add(sketch_vertex_list, 0, marked.cardinality());
            
            weight_of_current_index += total_out_degree;
        }
        
        node_infl[0] = num_sketches;// Real number of sketches computed - including the ones with just one node.
        
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
    
    class NodeInfluence {
        int node;
        long influence;
        NodeInfluence(int node, long influence) {
            this.node = node;
            this.influence = influence;
        }
    }
    
    public class NodeInfluenceComparator implements Comparator<NodeInfluence> {
        public int compare(NodeInfluence x, NodeInfluence y) {
            if (x.influence < y.influence)
                return 1;
            if (x.influence > y.influence)
                return -1;
            return 0;
        }
    }
    
    int[] toArray(List<Integer> arrlist) {
        int[] res = new int[arrlist.size()];
        for(int i=0; i<arrlist.size(); i++)
            res[i] = arrlist.get(i);
        return res;
    }
    
    void get_seeds(ImmutableGraph G_sk, Map<Integer, List<Integer>> n_sk_map, long[] node_infl, int k_left, int hyper_rows, double total_infl, BitSet sk_gone, BitSet reCalc, PriorityQueue<NodeInfluence> pq, PriorityQueue<NodeInfluence> pq2, double coeff) throws Exception {

        // Calculating the node with max influence
       
    	NodeInfluence ni = pq.poll();
        int max_node = ni.node;
        long infl_max = ni.influence;
        
        while(!reCalc.get(max_node)) {

            if(n_sk_map.containsKey(max_node)) {

                int[] node_sketches = toArray(n_sk_map.get(max_node));
                int node_deg = n_sk_map.get(max_node).size();
            
                for (int nii = 0; nii < node_deg; nii++) {
                    int uu = node_sketches[nii];
                    if(sk_gone.get(uu)) {
                        if(node_infl[max_node] > 0)
                            node_infl[max_node] = node_infl[max_node] - 1;
                    }
                }
            
                reCalc.set(max_node); //mark max_node as re-calculated
                pq.add( new NodeInfluence(max_node,node_infl[max_node]) );
            
                ni = pq.poll();
                max_node = ni.node;
                infl_max = ni.influence;
            }
            else {
                System.out.println("Rebuiding nodes in the Map...");
                pq2.clear();
                for(int ii=0;ii<n;ii++) {
                    pq2.add( new NodeInfluence(ii,node_cover[ii]) );
                }
                n_sk_map = get_top(pq2,hyper_rows,G_sk.copy());
            }
        }

        total_infl = total_infl + 1.0*infl_max * coeff;
        
        System.out.println(max_node);
        
        // Stopping condition: no need to re-calculate the influence, if we got k seeds
        if((k_left - 1)==0)
        {
            System.out.println();
            System.out.println( "Set Influence = " + total_infl );
            return;
        }
        
        // Remove max node and the sketches it participated in
        // plus re-calculate the influence of the remaining nodes
        
        node_infl[max_node] = 0;
        
        int[] max_sketches = toArray(n_sk_map.get(max_node));
        int max_deg = n_sk_map.get(max_node).size();
        
        for (int nn = 0; nn < max_deg; nn++) {
            int uu = max_sketches[nn];
            sk_gone.set(uu);
        }
        
        
        reCalc.clear();
            
        get_seeds(G_sk.copy(), n_sk_map, node_infl, k_left-1, hyper_rows, total_infl, sk_gone, reCalc, pq, pq2, coeff);
        
    }
  
  


    public static void main(String[] args) throws Exception {
        long startTime = System.currentTimeMillis();
        long estimatedTime;
        
        
        if(args.length < 4) {
            System.out.println("Specify: basename, p, k, beta");
            System.exit(1);
        }
        
        String basename  = args[0];
        double p = Double.valueOf(args[1]);
        int k = Integer.valueOf(args[2]);
        double beta = Double.valueOf(args[3]);
        
        NoSinglesTopNodes imw = new NoSinglesTopNodes(basename, p, k, beta);
        imw.get_sketch();
        
        System.gc();
        
        
        estimatedTime = System.currentTimeMillis() - startTime;
        System.out.println();
        System.out.println("Total time elapsed = " + estimatedTime / 1000.0 + " seconds");
	}
}

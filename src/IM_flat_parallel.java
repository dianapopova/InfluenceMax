// This implementation of Reversed Influence Sampling algorithm (Borgs et al., Maximizing social influence in nearly optimal time) was written by Diana Popova, University of Victoria, Canada. Data structure for the hypergraph is flat, one-dimensional, arrays. The sampling is done in parallel using Java 8 parallel streams and lambda expressions.

import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.Random;

import java.util.Arrays;

import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.stream.IntStream;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;


public class IM_flat_parallel {
	ImmutableGraph G;
    int n, k;
    long m;
    String basename;
    double p, beta;
    int nMAX = 1400000000; // maximum possible for 16 GB main memory
    
//    BitSet marked, sk_gone, nodes_gone;
    
    // Flat arrays to keep: sketches - all sketch numbers in order, like 0000...1111...2222...;
    //                      nodes - all corresponding node IDs participating in sketch, like 13 17 100 230 ... 13 23 300 10000 ...
    //                      node_infl - counts of sketches the node is participating into;
	
	public IM_flat_parallel(String basename, Double  p, Double beta, int k) throws Exception {
        
        long readFile = System.currentTimeMillis();
        
        G = ImmutableGraph.load(basename+"t"); //  Transpose graph with inverse edge directions
        
        double fileTime = (System.currentTimeMillis() - readFile)/1000.0;
        
        System.out.println();
        System.out.println("Reading transpose graph took " + fileTime + " seconds");
        
		n = G.numNodes();
		m = G.numArcs();

        System.out.println("beta = " + beta);
        System.out.println("k = " + k);

        System.out.println("n="+n + ", m=" +m  + ", hypergraph weight =" + ( beta * k * m * Math.log(n) ) );
		
//		marked = new BitSet(n);
        
        this.basename = basename;
        this.p = p;
        this.beta = beta;
        this.k = k;
	}
    
    void get_sketch() throws Exception {
        long getHypergraph = System.currentTimeMillis();
        
        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("int cores = " + cores);
        
        double R = beta*m*k*Math.log(n)/cores;
        
        int node_infl[][] = new int[cores][];
        for(int i=0; i<cores; i++)
            node_infl[i] = new int[n+2];
        
        IntStream.range(0, cores)
        .parallel()
        .forEach(core -> {
            try {
                node_infl[core] = get_sketch(G.copy(),core,cores,R);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
        /////////////////////////////////////////////////////////////////////////////
        
        //Calculate the number of sketches and the size of sketches array:
        int sketch_num = 0;
        int count_sketches = 0;
        
        for(int ii = 0; ii<cores;ii++) {
            sketch_num = sketch_num + node_infl[ii][n];
            System.out.println("Number of sketches taken by core" + ii + " = " + node_infl[ii][n] );
            count_sketches = count_sketches + node_infl[ii][n+1];
        }
        System.out.println("Total number of sketches = " + sketch_num );
        System.out.println("Size of sketches array = " + count_sketches );

        
        //Calculate the combined node influence (from all cores):
        
        int node_cover[] = new int[n];
        
        for(int j=0;j<n;j++) {
            for(int i = 0; i<cores;i++)
                node_cover[j] = node_cover[j] + node_infl[i][j];
        }
        
        //Merge the sketches and nodes arrays produced by the cpu cores:
       
        
        int[] sketches = new int[count_sketches];
        int[] nodes = new int[count_sketches];
        
        for(int i=0;i<count_sketches;i++)
        {
            sketches[i] = -1;
            nodes[i] = -1;
        }
        
        int start = 0;
        for(int core=0; core<cores; core++) {
            int sketch_size = node_infl[core][n+1];
            int[] sketchArr = new int[sketch_size];
            int[] nodeArr = new int[sketch_size];
            
            readIntArrayFromBinaryFile(sketchArr, G.basename()+"-sketches"+core);
            readIntArrayFromBinaryFile(nodeArr, G.basename()+"-nodes"+core);
            
            for(int i=0;i<sketch_size;i++) {
                sketches[i+start] = sketchArr[i] + start;
                nodes[i+start] = nodeArr[i];
            }
            start = start + sketch_size;
            
        }
        
        /////////////////////////////////////////////////////////////////////////////
        
        
        // Getting the seeds:
        
        double set_infl = 0.0;
        double coeff = 1.0 * n/sketch_num;
        BitSet sk_gone = new BitSet(sketch_num);
        BitSet nodes_gone = new BitSet(count_sketches);
        
        long startSeeds = System.currentTimeMillis();
        
        get_seeds(sketches, nodes, node_cover, k, count_sketches, sketch_num, set_infl, coeff, sk_gone, nodes_gone);
        
        double getSeeds = (System.currentTimeMillis() - startSeeds)/1000.0;
        
        System.out.println();
        System.out.println("Calculating seeds took " + getSeeds + " seconds");

    }
    
	
	int[] get_sketch(ImmutableGraph G, int core, int cores, double R) throws Exception {
        
	    double weight_of_current_index = 0.0;
        int sketch_num = 0;
        int count_sketches = 0;
        Random gen_rnd = new Random();
        
        int[] node_infl = new int[n+2];
        
        BitSet marked = new BitSet(n);
        
        int array_size = nMAX/cores;
        
        int[] sketches = new int[array_size];
        int[] nodes = new int[array_size];
        
        for(int i=0;i<array_size;i++)
        {
            sketches[i] = -1;
            nodes[i] = -1;
        }

	    while(weight_of_current_index < R)
	    {
	    	int v = gen_rnd.nextInt(n);
	        marked.clear();
	        BFS(G.copy(),v,marked);

	        int total_out_degree = 0;
            int iteration = 0;
	        for (int u = marked.nextSetBit(0); u >= 0; u = marked.nextSetBit(u+1))
	        {
                sketches[count_sketches + iteration] = sketch_num;
                nodes[count_sketches + iteration]  = u;
                node_infl[u] = node_infl[u] + 1;
                iteration = iteration + 1;
	            total_out_degree = total_out_degree + G.outdegree(u);
	        }

	        weight_of_current_index = weight_of_current_index + total_out_degree;
            sketch_num = sketch_num + 1;
            count_sketches = count_sketches + marked.cardinality();
	    }
        
        // Cutting off the -1 elements of sketches and nodes arrays, making the arrays shorter
        int[] iSketch = new int[count_sketches];
        System.arraycopy(sketches,0,iSketch,0,count_sketches);
        System.gc();
        
        int[] iNode = new int[count_sketches];
        System.arraycopy(nodes,0,iNode,0,count_sketches);
        System.gc();
        
        //        System.out.println("Writing arrays into binary files...");
        writeIntArrayToBinaryFile(iSketch, G.basename()+"-sketches"+core);
        writeIntArrayToBinaryFile(iNode, G.basename()+"-nodes"+core);
        
        node_infl[n] = sketch_num;
        node_infl[n+1] = count_sketches;
        
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

    public static void writeIntArrayToBinaryFile(int[] arr, String filename) throws Exception {
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename)));

        for (int i=0; i<arr.length; i++)
            dos.writeInt(arr[i]);

        dos.close();
    }


    public static void readIntArrayFromBinaryFile(int[] arr, String filename) throws Exception {

        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(filename)));

        for (int i=0; i<arr.length; i++)
            arr[i] = dis.readInt();

        dis.close();
    }

    void get_seeds(int[] iSketch,int[] iNode, int[] node_infl, int k_left, int count_sketches, int sketch_num, double set_infl, double coeff, BitSet sk_gone, BitSet nodes_gone) {

        // Calculating the node with max influence
        int infl_max = 0;
        int max_node = 0;
        
        for(int v=0;v<n;v++)
        {
               if(node_infl[v] > infl_max)
                {
                    infl_max = node_infl[v];
                    max_node = v;
                }
        }
        
        set_infl = set_infl + infl_max * coeff;

        System.out.println(max_node +
                           ", Its Influence = " + infl_max);
        
        
        // Stopping condition: no need to re-calculate the influence, if we already got the k seeds
        if((k_left - 1)==0) {
            System.out.println("Set Influence = " + set_infl);
            return;
        }
        
        // Re-calculating the influence of the remaining nodes: remove max node and the sketches it participated in
        // plus re-calculate the influence
        node_infl[max_node] = 0;

        for(int j=0;j<count_sketches;j++)
        {
            if(nodes_gone.get(j))
                continue;
            
            else
            {
                if((iNode[j] == max_node) && (!sk_gone.get(iSketch[j])))
                {
                    int redundant_sketch = iSketch[j];
                    sk_gone.set(redundant_sketch);

                // As sketches are added to the array in numerical order, the same redundant sketch can be found before and after the max node
                    if(j < (count_sketches - 1))
                    {
                        int l = j+1;
                        while((l<count_sketches) && (iSketch[l] == redundant_sketch))
                        {
                            if(node_infl[iNode[l]] > 0)
                                node_infl[iNode[l]] = node_infl[iNode[l]] - 1;
                            
                            nodes_gone.set(iNode[l]);
                            l++;
                        }
                    }
            
                    if(j>0)
                    {
                        int ll = j-1;
                        
                        while(((ll+1) > 0) && (iSketch[ll] == redundant_sketch))
                        {
                            if(node_infl[iNode[ll]] > 0)
                                node_infl[iNode[ll]] = node_infl[iNode[ll]] - 1;
                            
                            nodes_gone.set(iNode[ll]);
                            ll--;
                        }
                    }
                }
            }
        }
            
            nodes_gone.set(max_node);
            get_seeds(iSketch, iNode, node_infl, k_left-1, count_sketches, sketch_num, set_infl, coeff, sk_gone, nodes_gone);
    }

		
	public static void main(String[] args) throws Exception {
		long startTime = System.currentTimeMillis();
		long estimatedTime;

        if(args.length < 4) {
            System.out.println("Specify: basename, p, beta, k");
            System.exit(1);
        }
        
        String basename  = args[0];
        double p = Double.valueOf(args[1]);
        double beta = Double.valueOf(args[2]);
        int k = Integer.valueOf(args[3]);
		
        IM_flat_parallel imfl = new IM_flat_parallel(basename, p, beta, k);
		imfl.get_sketch();
			
		estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Time elapsed = " + estimatedTime /(1000.0) + " sec");
	}
}

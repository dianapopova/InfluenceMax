// This implementation of Reversed Influence Sampling algorithm (Borgs et al., Maximizing social influence in nearly optimal time) was written by Diana Popova, University of Victoria, Canada. Data structure for the hypergraph is flat, one-dimensional, arrays.

import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.Random;

import java.util.Arrays;

import it.unimi.dsi.webgraph.ImmutableGraph;

public class IM_flat {
	ImmutableGraph G;
	int n;
	long m;
	double eps = .1;
    double p; // given probability of arc existence
    int k = 5; // given number of maximum influence nodes
    int nMAX = 1400000000; // maximum possible for 16 GB main memory
    
    BitSet marked, sk_gone, nodes_gone;
    
    double beta = 2; // beta parameter
    
    // Flat arrays to keep: sketches - all sketch numbers in order, like 0000...1111...2222...;
    //                      nodes - all corresponding node IDs participating in sketch, like 13 17 100 230 ... 13 23 300 10000 ...
    //                      node_infl - counts of sketches the node is participating into;
    int[] sketches;
    int[] nodes;
    int[] node_infl;

    int count_sketches; // the length of sketches and nodes arrays
	
	public IM_flat(String basename, Double  p) throws Exception {
		G = ImmutableGraph.load(basename);
		
		n = G.numNodes();
		m = G.numArcs();

        System.out.println("beta = " + beta);
        System.out.println("k = " + k);

        System.out.println("n="+n + ", m=" +m  + ", R=" +( beta * k * m*Math.log(n)  ));
		
		marked = new BitSet(n);
        
        sketches = new int[nMAX];
        nodes = new int[nMAX];
        node_infl = new int[n];
        
        for(int i=0;i<nMAX;i++)
        {
            sketches[i] = -1;
            nodes[i] = -1;
        }
        
        this.p = p;
	}
    
	
	void get_sketch() {

        double R = beta * k * m * Math.log(n);
                           
	    double weight_of_current_index = 0.0;
//	    int index = 0;
        int sketch_num = 0;

	    long startTime = System.currentTimeMillis();
        count_sketches = 0;
        Random gen_rnd = new Random();

	    while(weight_of_current_index < R)
	    {
/*	    	if(index % 1000000 == 0) {
	    		System.out.println(
                        "sketch=" + sketch_num + 
	    				",  index=" + index +
                        ", weight_of_current_index=" + weight_of_current_index);
	    	}
*/
	    	int v = gen_rnd.nextInt(n);
	        marked.clear();
	        BFS(v,marked);

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
//	        index = ( index + 1 ) % n;
            sketch_num = sketch_num + 1;
            count_sketches = count_sketches + marked.cardinality();
	    }
        
        System.out.println();
	    System.out.println("Number of Sketches: " + sketch_num);
//                           ", Size of array iSketch: " + count_sketches);
        System.out.println();
        
        // Cutting off the tails of sketches and nodes arrays, making the arrays shorter
//        int[] iSketch = new int[count_sketches + 1];
//        System.arraycopy(sketches,0,iSketch,0,count_sketches);
//        System.gc();

        
//        int[] iNode = new int[count_sketches + 1];
//        System.arraycopy(nodes,0,iNode,0,count_sketches);

/*        System.out.println("First sketch numbers in iSketch = " + iSketch[0] + "; " + iSketch[1] + "; " + iSketch[2] + "; " + iSketch[3] + ".");
        System.out.println("First nodes in iNode = " + iNode[0] + "; " + iNode[1] + "; " + iNode[2]+ "; " + iNode[3] + ".");
*/
        System.gc();
        
        long startSeeds = System.currentTimeMillis();
        
        double set_infl = 0.0;
        double coeff = 1.0 * n/sketch_num;
        sk_gone = new BitSet(sketch_num);
        nodes_gone = new BitSet(count_sketches);

        get_seeds(sketches, nodes, node_infl, k, count_sketches, sketch_num, set_infl, coeff, sk_gone, nodes_gone);
        
        double getSeeds = (System.currentTimeMillis() - startSeeds)/1000.0;
        
        System.out.println();
        System.out.println("Calculating seeds took " + getSeeds + " seconds");
	}
		
	void BFS(int v, BitSet marked) {
		
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
                    int redundant_sketch = sketches[j];
                    sk_gone.set(redundant_sketch);

                // As sketches are added to the array in numerical order, the same redundant sketch can be found before and after the max node
                    if(j < (count_sketches - 1))
                    {
                        int l = j+1;
                        while((l<count_sketches) && (iSketch[l] == redundant_sketch))
                        {
                            if(node_infl[iNode[l]] > 0)
                                node_infl[iNode[l]] = node_infl[iNode[l]] - 1;
                            
                            nodes_gone.set(nodes[l]);
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
                            
                            nodes_gone.set(nodes[ll]);
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

		//args = new String[] {"cnr-nlt", "0.1"};
        args = new String[] {"uk100-nlt", "0.1"};
        //args = new String[] {"eu-nlt", "0.1"};


		
		String basename  = args[0];
		double p = Double.valueOf(args[1]);
		
        IM_flat imfl = new IM_flat(basename, p);
		imfl.get_sketch();
			
		estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Time elapsed = " + estimatedTime /(1000.0) + " sec");
	}
}

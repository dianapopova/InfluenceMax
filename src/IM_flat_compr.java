// This implementation of Reversed Influence Sampling algorithm (Borgs et al., Maximizing social influence in nearly optimal time) was written by Diana Popova, University of Victoria, Canada. Data structure for the hypergraph is flat, one-dimensional, arrays.

import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.Random;

import java.util.Arrays;

import it.unimi.dsi.webgraph.ImmutableGraph;

public class IM_flat_compr {
    ImmutableGraph G;
    int n, k;
    long m;
    String basename;
    double p, beta;
    int nMAX = 2147483645;
    
    BitSet marked, sk_gone, nodes_gone;
    
    // Flat arrays to keep: sketches - number of nodes participating in each sketch, accumulating, like 4 ||  7 || 11 || ...; index of the element is the sketch ID;
    //                      nodes - all corresponding node IDs participating in a sketch, like 13 17 100 230 || 1 23 300 || 13 23 300 10000 || ...
    //                      node_infl - counts of sketches the node is participating in;
    int[] sketches;
    int[] nodes;
    int[] node_infl;

    int count_sketches; // the length of nodes array
	
	public IM_flat_compr(String basename, Double  p, Double beta, int k) throws Exception {
		G = ImmutableGraph.load(basename);
		
		n = G.numNodes();
		m = G.numArcs();

        System.out.println("beta = " + beta);
        System.out.println("k = " + k);

        System.out.println("n="+n + ", m=" +m  + ", R=" +( beta * k * m*Math.log(n)  ));
		
		marked = new BitSet(n);
        
        sketches = new int[nMAX/50];
        nodes = new int[nMAX];
        node_infl = new int[n];
        
        for(int i=0;i<nMAX;i++)
        {
            nodes[i] = -1;
        }
        
        for(int i=0;i<nMAX/50;i++)
        {
            sketches[i] = -1;
        }
        
        this.basename = basename;
        this.p = p;
        this.beta = beta;
        this.k = k;
	}
    
	
	void get_sketch() {

        double R = beta * k * m * Math.log(n);
                           
	    double weight_of_current_index = 0.0;
//	    int index = 0;
        int sketch_num = 0;

	    long startTime = System.currentTimeMillis();
        int accumulated_sketches = 0;

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
                nodes[accumulated_sketches + iteration]  = u;
                node_infl[u] = node_infl[u] + 1;
                iteration = iteration + 1;
	            total_out_degree = total_out_degree + G.outdegree(u);
	        }
            accumulated_sketches += marked.cardinality();
            sketches[sketch_num] = accumulated_sketches;

	        weight_of_current_index = weight_of_current_index + total_out_degree;
//	        index = ( index + 1 ) % n;
            sketch_num = sketch_num + 1;
	    }
        
        System.out.println();
	    System.out.println("Number of Sketches: " + sketch_num);
        System.out.println();
        
        System.gc();
        
        long startSeeds = System.currentTimeMillis();
        
        double set_infl = 0.0;
        double coeff = 1.0 * n/sketch_num;
        sk_gone = new BitSet(sketch_num);
        nodes_gone = new BitSet(count_sketches);

        get_seeds(sketches, nodes, node_infl, k, accumulated_sketches, sketch_num, set_infl, coeff, sk_gone, nodes_gone);
        
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
    
    void get_seeds(int[] sketches,int[] nodes, int[] node_infl, int k_left, int count_sketches, int sketch_num, double set_infl, double coeff, BitSet sk_gone, BitSet nodes_gone) {

        // Calculating the node with max influence
        int infl_max = 0;
        int max_node = -100;
        int redundant_sketch = -100;
        
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

        int start_sk = 0;
        
        for(int j=0;j<count_sketches;j++)
        {
            if(nodes_gone.get(j))
                continue;
            
            else
            {
                if(nodes[j] == max_node)
                {
                    for(int sn = start_sk; sn < sketch_num; sn++)
                    {
                        if(sk_gone.get(sn) || sketches[sn] <= j)
                            continue;
                        else
                        {
                            redundant_sketch = sn;
                            sk_gone.set(redundant_sketch);
                        

                // As nodes are added to the array in numerical order, the affected nodes can be found before and after the max node
                            if(j < (count_sketches - 1))
                            {
                                int l = j+1;
                    
                                while((l<sketches[redundant_sketch]) && (l<count_sketches))
                                {
                                    if(node_infl[nodes[l]] > 0)
                                        node_infl[nodes[l]] = node_infl[nodes[l]] - 1;
                                    
                                    nodes_gone.set(nodes[l]);
                                    l++;
                                }

                            }
            
                            if(j>0)
                            {
                                int ll = j- 1;

                                if(redundant_sketch == 0) {
                                    while((ll+1) > 0)
                                    {
                                        if(node_infl[nodes[ll]] > 0)
                                            node_infl[nodes[ll]] = node_infl[nodes[ll]] - 1;
                                        
                                        nodes_gone.set(nodes[ll]);
                                        ll--;
                                    }
                                    start_sk = redundant_sketch + 1;
                                    break;
                                    }
                    
                                else {
                                    while((ll+1) > sketches[redundant_sketch - 1])
                                    {
                                        if(node_infl[nodes[ll]] > 0)
                                            node_infl[nodes[ll]] = node_infl[nodes[ll]] - 1;
                                        
                                        nodes_gone.set(nodes[ll]);
                                        ll--;
                                    }
                                    start_sk = redundant_sketch + 1;
                                    break;

                                }
                            }
                            start_sk = redundant_sketch + 1;
                            break;
                        }
                    }
                }
            }
        }
        nodes_gone.set(max_node);

        get_seeds(sketches, nodes, node_infl, k_left - 1, count_sketches, sketch_num, set_infl, coeff, sk_gone, nodes_gone);
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

		//args = new String[] {"cnr-nlt", "0.1"};
        //args = new String[] {"uk100-nlt", "0.1"};
        //args = new String[] {"eu-nlt", "0.1"};
		
		//String basename  = args[0];
		//double p = Double.valueOf(args[1]);
		
        IM_flat_compr imfl = new IM_flat_compr(basename, p, beta, k);
		imfl.get_sketch();
			
		estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Time elapsed = " + estimatedTime /(1000.0) + " sec");
	}
}

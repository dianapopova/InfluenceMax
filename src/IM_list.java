// This implementation of Reversed Influence Sampling algorithm (Borgs et al., Maximizing social influence in nearly optimal time) was written by Diana Popova, University of Victoria, Canada. Data structure for the hypergraph is two-dimensional list.

import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import it.unimi.dsi.webgraph.ImmutableGraph;

import java.util.Comparator;
import java.util.PriorityQueue;

import java.util.stream.IntStream;
import java.io.IOException;

public class IM_list {
	ImmutableGraph G;
	int n, k;
    long m;
    String basename;
    double p, beta; // given probability of arc existence
    
    BitSet marked, sk_gone;

	
	public IM_list(String basename, Double  p, Double beta, int k) throws Exception {
		G = ImmutableGraph.load(basename);
		
		n = G.numNodes();
		m = G.numArcs();
        System.out.println("beta = " + beta +
                           "; k = " + k);

        System.out.println("n = "+n + ", m = " +m  + ", R = " + beta*m*k*Math.log(n) );
      
        marked = new BitSet(n);
        
        this.basename = basename;
        this.p = p;
        this.beta = beta;
        this.k = k;

	}
    
	void get_sketch() throws Exception {
        long getHypergraph = System.currentTimeMillis();

        double R = beta*m*k*Math.log(n);
        
        List<List<Integer>> I = new ArrayList<List<Integer>>();
        
        for(int j=0;j<n;j++)
        {
            I.add(new ArrayList<Integer>());
        }
//        System.out.println("index I created.");
                           
	    double weight_of_current_index = 0.0;
//	    int index = 0;
        int sketch_num = 0;
        
        Random gen_rnd = new Random();

	    while(weight_of_current_index < R)
	    {
/*	    	if(index % 100000 == 0) {
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
	        for (int u = marked.nextSetBit(0); u >= 0; u = marked.nextSetBit(u+1))
	        {
                I.get(u).add(sketch_num);
	            total_out_degree += G.outdegree(u);
	        }
	        weight_of_current_index += total_out_degree;
//	        index = ( index + 1 ) % n;
            sketch_num++;
	    }
	    System.out.println(
//                           "Index: " + index +
                           "Number of Sketches: " + sketch_num);
        
        double hypergraphTime = (System.currentTimeMillis() - getHypergraph)/1000.0;
        
        System.out.println();
        System.out.println("Creating Hypergraph took " + hypergraphTime + " seconds");
        System.out.println();
        
        System.gc();
        
        int[] node_degrees = new int[n];
        
        int a = 0;
        while(a < n) {
            node_degrees[a] = I.get(a).size();
            a++;
        }
        long startSeeds = System.currentTimeMillis();
        
        double coeff = 1.0 * n/sketch_num;

        PriorityQueue<NodeInfluence> pq = new PriorityQueue<NodeInfluence>(n, new NodeInfluenceComparator());
        
        for(int ii=0;ii<n;ii++) {
            pq.add( new NodeInfluence(ii,node_degrees[ii]) );
        }
        
        System.gc();
        
/*
         // Getting the spread for the given set of seeds:
         
         int[] seeds = {110604, 87976, 229573, 233431, 262756}; //for cnr2000
         infl_est(I, node_degrees, seeds, coeff);
*/

        //Calculating the seeds:
        
        // Calculating the first node with max influence
        
        int nodeID = pq.poll().node;
        int nodeInfl = pq.poll().influence;
        double total_infl = nodeInfl * coeff;
        

        System.out.println(nodeID + ", its sketch count = " + nodeInfl);

        
        // Stopping condition: no need to re-calculate the influence, if we need only one seed
        if((k - 1)==0) {
            System.out.println("Max Node Influence = " + total_infl );
            System.out.println();

            double getOneSeed = (System.currentTimeMillis() - startSeeds)/1000.0;

            System.out.println("Calculating one seed took " + getOneSeed + " seconds");
            System.out.println();
            return;
        }
        else {
            //Create a priority queue for the nodes by their influence in node_degrees[]
            
            BitSet reCalc = new BitSet(n);
            BitSet sk_gone = new BitSet(sketch_num);
            
            // Calling the get_seeds method:
            
            get_seeds(I, k-1, pq, nodeID, total_infl, sk_gone, reCalc, coeff);
            
            double getSeeds = (System.currentTimeMillis() - startSeeds)/1000.0;
            
            System.out.println();
            System.out.println("Calculating seeds took " + getSeeds + " seconds");

            
        }

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
    

    class NodeInfluence {
        int node;
        int influence;
        NodeInfluence(int node, int influence) {
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

        
    void get_seeds(List<List<Integer>> I, int k_left,  PriorityQueue<NodeInfluence> pq, int max_node, double total_infl, BitSet sk_gone, BitSet reCalc, double coeff) throws Exception {
        
        // Remove max node and the sketches it participated in
        // plus re-calculate the influence of the head node if needed
        
        int[] max_sketches = I.get(max_node).stream().mapToInt(i -> i).toArray();
        int max_deg = max_sketches.length;
        
        for (int ni = 0; ni < max_deg; ni++) {
            int ss = max_sketches[ni];
            sk_gone.set(ss);
        }
        
        int nodeID = pq.poll().node;
        int nodeInfl = pq.poll().influence;
        
        
        while(!reCalc.get(nodeID)) {
            
            nodeInfl = 0;
            int[] node_sketches = I.get(nodeID).stream().mapToInt(i -> i).toArray();
            int node_deg = node_sketches.length;
            
            for (int ni = 0; ni < node_deg; ni++) {
                int ss = node_sketches[ni];
                if(!sk_gone.get(ss))
                    nodeInfl = nodeInfl + 1;
            }
            
            reCalc.set(nodeID);
            pq.add( new NodeInfluence(nodeID,nodeInfl) );
            
            nodeID = pq.poll().node;
            nodeInfl = pq.poll().influence;
            
        }
        
        total_infl = total_infl + nodeInfl * coeff;
        

        System.out.println(nodeID +
                           ", its sketch count = " + nodeInfl);
        
        
        // Stopping condition: no need to re-calculate the influence, if we already got the k seeds
        if((k_left - 1)==0) {
            System.out.println("Seed Set Influence = " + total_infl );
            return;
        }

        else {
            reCalc.clear();
            get_seeds(I, k_left- 1, pq, nodeID, total_infl, sk_gone, reCalc, coeff);
        }
    }

 
    void infl_est(List<List<Integer>> I, int[] node_degrees, int[] seeds, double coeff) throws Exception {
        
        long getSpread = System.currentTimeMillis();
        
        System.out.println("Seeds:");
        System.out.println(Arrays.toString(seeds));
        
        int seeds_num = seeds.length;
        
        int arrays[][] = new int[seeds_num][];
        for(int ii=0; ii<seeds_num; ii++)
            arrays[ii] = I.get(seeds[ii]).stream().mapToInt(i -> i).toArray();
        
        int maxSize = 0;
        int counter = 0;
        
        for(int[] array : arrays) {
            maxSize += array.length;
        }
        
        int[] accumulator = new int[maxSize];
        
        for(int[] array : arrays) {
            for(int i : array) {
                if(!isDuplicated(accumulator, counter, i))
                    accumulator[counter++] = i;
            }
            System.out.println();
            System.out.println("Set influence = " + (counter * coeff));
        }
        long spreadTime = System.currentTimeMillis() - getSpread;
        System.out.println();
        System.out.println("Spread calculation time = " + spreadTime / 1000.0 + " seconds");
        
    }
    
    boolean isDuplicated(int[] array, int counter, int value)
    {
        for(int i = 0; i < counter; i++) {
            if(array[i] == value)
                return true;
        }
        return false;
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
		
		IM_list iml = new IM_list(basename, p, beta, k);
		iml.get_sketch();
			
		estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Time elapsed = " + estimatedTime / 1000.0 + " sec");	
	}
}

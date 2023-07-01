1. run olsr_simulation.cc
    //for original olsr
    $ ./ns3 run scratch/olsr_simulation.cc -- --trName="node20/olsr_ori_20" --numNodes=20 --SimulationTime=150.0
    // Remember to comment the for-loop in sendHello() in olsr-routing-protocol.cc file
    // to simulate the original OSLR
    $ ./ns3 run scratch/olsr_simulation.cc -- --trName="node20/olsr_adj_20" --numNodes=20 --SimulationTime=150.0
    //adj = adjusted
    //the last number of trName represents the # of nodes in this simulation
    
2. generate py-plot
	$ python3 flow.py olsr_ori_20.flowmon olsr_adj_20.flowmon result_20.pdf
	//                [file path 1]       [file path 2]       [output file path]        

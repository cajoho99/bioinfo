CFLAGS= -Wall -g

all:	 make_distance_map domain_partition strudl

make_distance_map: src/make_distance_map.c
	gcc $(CFLAGS) src/make_distance_map.c -o make_distance_map -lm 

domain_partition: src/domain_partition.c
	gcc $(CFLAGS) -D RESULT src/domain_partition.c -o domain_partition  -lm

domain_partition_debug: src/domain_partition.c
	gcc $(CFLAGS) -D DEBUG src/domain_partition.c -o domain_partition  -lm


domain_partition_graph: src/domain_partition.c
	gcc $(CFLAGS) -D GRAPH src/domain_partition.c -o domain_partition  -lm
	./domain_partition ./data/1cdh.pdb > ./out/1cdh.map
	./domain_partition ./data/2csn.pdb > ./out/2csn.map
	python3 ./display_scripts/graph_split.py ./out/1cdh.map &
	python3 ./display_scripts/graph_split.py ./out/2csn.map &

strudl: src/strudl_domain_partition.c
	gcc $(CFLAGS) src/strudl_domain_partition.c -o strudl -lm


clean:
	rm -f domain_partition
	rm -f make_distance_map
	rm -f strudl

# Formats adjacency.tsv files into appropriate format.

# Remove third column (edge weight).
sed -i -r 's/(\s+)?\S+//3' $1

# Sort by source nodes then destination nodes (neighbors).
sort -k1,1n -k2,2n <$1 > $1.tmp
mv -f $1.tmp $1

# Get number of nodes and edges.
NUM_EDGES=$(wc -l $1  | awk '{print $1}')
NUM_NODES_1=$(sort -nrk1,1 $1 | head -1 | awk '{print $1}')
NUM_NODES_2=$(sort -nrk2,2 $1 | head -1 | awk '{print $2}')
NUM_NODES=$NUM_NODES_1

if [ $NUM_NODES_2 -gt $NUM_NODES_1 ]
then
  NUM_NODES=$NUM_NODES_2
fi

# Prepend number of nodes and number of edges.
sed -i "1i$NUM_NODES $NUM_EDGES" $1

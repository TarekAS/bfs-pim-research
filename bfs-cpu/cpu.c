// Contains the CPU version of non-parallel BFS for comparison purposes.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SIZE 100000000
#define PRINT_ERROR(fmt, ...) \
  fprintf(stderr, "\033[0;31mERROR:\033[0m   " fmt "\n", ##__VA_ARGS__)
#define PRINT_WARNING(fmt, ...) \
  fprintf(stderr, "\033[0;35mWARNING:\033[0m " fmt "\n", ##__VA_ARGS__)
#define PRINT_INFO(fmt, ...) \
  fprintf(stderr, "\033[0;32mINFO:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define PRINT_STATUS(status) \
  fprintf(stderr, "Status: %s\n", dpu_api_status_to_string(status))

struct COO {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_nonzeros;
  uint32_t *row_idxs;
  uint32_t *col_idxs;
};

struct CSR {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_nonzeros;
  uint32_t *row_ptrs;
  uint32_t *col_idxs;
};

uint32_t *node_levels;

// Reads a coo-formated file into memory.
struct COO read_coo_matrix(char *file) {

  PRINT_INFO("Loading COO-formated graph from %s.", file);

  struct COO coo;

  // Initialize fields.
  uint32_t num_nodes = 0;
  FILE *fp = fopen(file, "r");
  fscanf(fp, "%u", &num_nodes);
  fscanf(fp, "%u", &coo.num_nonzeros);
  coo.num_rows = num_nodes;
  coo.num_cols = num_nodes;

  if (coo.num_rows % 2 == 1) {
    PRINT_WARNING("Number of rows must be even. Padding with an extra row.");
    coo.num_rows++;
  }
  if (coo.num_cols % 2 == 1) {
    PRINT_WARNING("Number of columns must be even. Padding with an extra column.");
    coo.num_cols++;
  }

  coo.row_idxs = malloc(coo.num_nonzeros * sizeof(uint32_t));
  coo.col_idxs = malloc(coo.num_nonzeros * sizeof(uint32_t));

  PRINT_INFO("Reading COO-formated matrix - %u rows, %u columns, %u nonzeros.", coo.num_rows, coo.num_cols, coo.num_nonzeros);

  // Read nonzeros.
  for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
    uint32_t rowIdx;
    uint32_t colIdx;
    fscanf(fp, "%u", &rowIdx);
    fscanf(fp, "%u", &colIdx);
    coo.row_idxs[i] = rowIdx;
    coo.col_idxs[i] = colIdx;
  }

  return coo;
}

// Converts COO matrix to CSR format.
struct CSR coo_to_csr(struct COO coo) {

  struct CSR csr;

  // Initialize fields
  csr.num_rows = coo.num_rows;
  csr.num_cols = coo.num_cols;
  csr.num_nonzeros = coo.num_nonzeros;
  csr.row_ptrs = malloc((csr.num_rows + 1) * sizeof(uint32_t));
  csr.col_idxs = malloc(csr.num_nonzeros * sizeof(uint32_t));

  // Histogram rowIdxs
  memset(csr.row_ptrs, 0, (csr.num_rows + 1) * sizeof(uint32_t));
  for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
    uint32_t rowIdx = coo.row_idxs[i];
    csr.row_ptrs[rowIdx]++;
  }

  // Prefix sum rowPtrs
  uint32_t sumBeforeNextRow = 0;
  for (uint32_t rowIdx = 0; rowIdx < csr.num_rows; ++rowIdx) {
    uint32_t sumBeforeRow = sumBeforeNextRow;
    sumBeforeNextRow += csr.row_ptrs[rowIdx];
    csr.row_ptrs[rowIdx] = sumBeforeRow;
  }
  csr.row_ptrs[csr.num_rows] = sumBeforeNextRow;

  // Bin the nonzeros
  for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
    uint32_t rowIdx = coo.row_idxs[i];
    uint32_t nnzIdx = csr.row_ptrs[rowIdx]++;
    csr.col_idxs[nnzIdx] = coo.col_idxs[i];
  }

  // Restore rowPtrs
  for (uint32_t rowIdx = csr.num_rows - 1; rowIdx > 0; --rowIdx) {
    csr.row_ptrs[rowIdx] = csr.row_ptrs[rowIdx - 1];
  }
  csr.row_ptrs[0] = 0;

  return csr;
}

void free_coo(struct COO coo) {
  free(coo.row_idxs);
  free(coo.col_idxs);
}

void free_csr(struct CSR csr) {
  free(csr.row_ptrs);
  free(csr.col_idxs);
}

struct queue {
  int items[SIZE];
  int front;
  int rear;
};

struct queue *createQueue();
void enqueue(struct queue *q, int);
int dequeue(struct queue *q);
void display(struct queue *q);
int isEmpty(struct queue *q);
void printQueue(struct queue *q);

struct node {
  int vertex;
  struct node *next;
};

struct node *createNode(int);

struct Graph {
  int numVertices;
  struct node **adjLists;
  int *visited;
};

struct Graph *createGraph(int vertices);
void addEdge(struct Graph *graph, int src, int dest);
void printGraph(struct Graph *graph);
void bfs(struct Graph *graph, int startVertex);

void bfs(struct Graph *graph, int startVertex) {

  struct queue *q = createQueue();

  graph->visited[startVertex] = 1;
  enqueue(q, startVertex);
  node_levels[startVertex] = 0;

  while (!isEmpty(q)) {
    printQueue(q);
    int currentVertex = dequeue(q);
    // printf("Visited %d\n", currentVertex);

    struct node *temp = graph->adjLists[currentVertex];

    while (temp) {
      int adjVertex = temp->vertex;

      if (graph->visited[adjVertex] == 0) {
        graph->visited[adjVertex] = 1;
        node_levels[adjVertex] = node_levels[currentVertex] + 1; // Distance from root.
        if (node_levels)
          enqueue(q, adjVertex);
      }
      temp = temp->next;
    }
  }
}

struct node *createNode(int v) {
  struct node *newNode = malloc(sizeof(struct node));
  newNode->vertex = v;
  newNode->next = NULL;
  return newNode;
}

struct Graph *createGraph(int vertices) {
  struct Graph *graph = malloc(sizeof(struct Graph));
  graph->numVertices = vertices;

  graph->adjLists = malloc(vertices * sizeof(struct node *));
  graph->visited = malloc(vertices * sizeof(int));

  int i;
  for (i = 0; i < vertices; i++) {
    graph->adjLists[i] = NULL;
    graph->visited[i] = 0;
  }

  return graph;
}

void addEdge(struct Graph *graph, int src, int dest) {
  // Add edge from src to dest
  struct node *newNode = createNode(dest);
  newNode->next = graph->adjLists[src];
  graph->adjLists[src] = newNode;

  // Add edge from dest to src
  newNode = createNode(src);
  newNode->next = graph->adjLists[dest];
  graph->adjLists[dest] = newNode;
}

struct queue *createQueue() {
  struct queue *q = malloc(sizeof(struct queue));
  q->front = -1;
  q->rear = -1;
  return q;
}

int isEmpty(struct queue *q) {
  if (q->rear == -1)
    return 1;
  else
    return 0;
}

void enqueue(struct queue *q, int value) {
  if (q->rear == SIZE - 1)
    printf("\nQueue is Full!!");
  else {
    if (q->front == -1)
      q->front = 0;
    q->rear++;
    q->items[q->rear] = value;
  }
}

int dequeue(struct queue *q) {
  int item;
  if (isEmpty(q)) {
    // printf("Queue is empty");
    item = -1;
  } else {
    item = q->items[q->front];
    q->front++;
    if (q->front > q->rear) {
      // printf("Resetting queue");
      q->front = q->rear = -1;
    }
  }
  return item;
}

void printQueue(struct queue *q) {
  int i = q->front;

  if (isEmpty(q)) {
    // printf("Queue is empty");
  } else {
    // printf("\nQueue contains \n");
    for (i = q->front; i < q->rear + 1; i++) {
      // printf("%d ", q->items[i]);
    }
  }
}

int main() {

  // Load coo-matrix from file and convert to csr.
  struct COO coo = read_coo_matrix("data/loc-gowalla_edges.txt");
  struct CSR csr = coo_to_csr(coo);

  node_levels = calloc(csr.num_rows, sizeof(uint32_t));

  struct Graph *graph = createGraph(csr.num_rows);

  for (int node = 0; node < csr.num_rows; ++node) {
    uint32_t from = csr.row_ptrs[node];
    uint32_t to = csr.row_ptrs[node + 1];

    for (int n = from; n < to; ++n) {
      uint32_t nb = csr.col_idxs[n]; // neighbor.
      addEdge(graph, node, nb);
    }
  }

  bfs(graph, 0);

  printf("node\tlevel\n");
  for (uint32_t node = 0; node < csr.num_rows; ++node) {
    uint32_t level = node_levels[node];
    if (node != 0 && level == 0) // Filters out "padded" rows.
      continue;
    printf("%u\t%u\n", node, node_levels[node]);
  }

  free(node_levels);
  free_coo(coo);
  free_csr(csr);

  return 0;
}

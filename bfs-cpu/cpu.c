// Contains the CPU version of BFS for comparison purposes.

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

struct COOMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *rowIdxs;
  uint32_t *colIdxs;
};

struct CSRMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *rowPtrs;
  uint32_t *colIdxs;
};

uint32_t *nodeLevels;

// Reads a coo-formated file into memory.
struct COOMatrix read_coo_matrix(char *file) {

  PRINT_INFO("Loading COO-formated graph from %s.", file);

  struct COOMatrix coo;

  // Initialize fields.
  FILE *fp = fopen(file, "r");
  fscanf(fp, "%u", &coo.numRows);
  fscanf(fp, "%u", &coo.numCols);
  fscanf(fp, "%u", &coo.numNonzeros);

  if (coo.numRows % 2 == 1) {
    PRINT_WARNING("Number of rows must be even. Padding with an extra row.");
    coo.numRows++;
  }
  if (coo.numCols % 2 == 1) {
    PRINT_WARNING("Number of columns must be even. Padding with an extra column.");
    coo.numCols++;
  }

  coo.rowIdxs = malloc(coo.numNonzeros * sizeof(uint32_t));
  coo.colIdxs = malloc(coo.numNonzeros * sizeof(uint32_t));

  PRINT_INFO("Reading COO-formated matrix - %u rows, %u columns, %u nonzeros.", coo.numRows, coo.numCols, coo.numNonzeros);

  // Read nonzeros.
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx;
    uint32_t colIdx;
    fscanf(fp, "%u", &rowIdx);
    fscanf(fp, "%u", &colIdx);
    coo.rowIdxs[i] = rowIdx;
    coo.colIdxs[i] = colIdx;
  }

  return coo;
}

// Converts COO matrix to CSR format.
struct CSRMatrix coo_to_csr(struct COOMatrix coo) {

  struct CSRMatrix csr;

  // Initialize fields
  csr.numRows = coo.numRows;
  csr.numCols = coo.numCols;
  csr.numNonzeros = coo.numNonzeros;
  csr.rowPtrs = malloc((csr.numRows + 1) * sizeof(uint32_t));
  csr.colIdxs = malloc(csr.numNonzeros * sizeof(uint32_t));

  // Histogram rowIdxs
  memset(csr.rowPtrs, 0, (csr.numRows + 1) * sizeof(uint32_t));
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    csr.rowPtrs[rowIdx]++;
  }

  // Prefix sum rowPtrs
  uint32_t sumBeforeNextRow = 0;
  for (uint32_t rowIdx = 0; rowIdx < csr.numRows; ++rowIdx) {
    uint32_t sumBeforeRow = sumBeforeNextRow;
    sumBeforeNextRow += csr.rowPtrs[rowIdx];
    csr.rowPtrs[rowIdx] = sumBeforeRow;
  }
  csr.rowPtrs[csr.numRows] = sumBeforeNextRow;

  // Bin the nonzeros
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    uint32_t nnzIdx = csr.rowPtrs[rowIdx]++;
    csr.colIdxs[nnzIdx] = coo.colIdxs[i];
  }

  // Restore rowPtrs
  for (uint32_t rowIdx = csr.numRows - 1; rowIdx > 0; --rowIdx) {
    csr.rowPtrs[rowIdx] = csr.rowPtrs[rowIdx - 1];
  }
  csr.rowPtrs[0] = 0;

  return csr;
}

void free_coo_matrix(struct COOMatrix coo) {
  free(coo.rowIdxs);
  free(coo.colIdxs);
}

void free_csr_matrix(struct CSRMatrix csr) {
  free(csr.rowPtrs);
  free(csr.colIdxs);
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
  nodeLevels[startVertex] = 0;

  while (!isEmpty(q)) {
    printQueue(q);
    int currentVertex = dequeue(q);
    // printf("Visited %d\n", currentVertex);

    struct node *temp = graph->adjLists[currentVertex];

    while (temp) {
      int adjVertex = temp->vertex;

      if (graph->visited[adjVertex] == 0) {
        graph->visited[adjVertex] = 1;
        nodeLevels[adjVertex] =
            nodeLevels[currentVertex] + 1; // Distance from root.
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
  struct COOMatrix coo = read_coo_matrix("data/loc-gowalla_edges.txt");
  struct CSRMatrix csr = coo_to_csr(coo);

  nodeLevels = calloc(csr.numRows, sizeof(uint32_t));

  struct Graph *graph = createGraph(csr.numRows);

  for (int node = 0; node < csr.numRows; ++node) {
    uint32_t from = csr.rowPtrs[node];
    uint32_t to = csr.rowPtrs[node + 1];

    for (int n = from; n < to; ++n) {
      uint32_t nb = csr.colIdxs[n]; // neighbor.
      addEdge(graph, node, nb);
    }
  }

  bfs(graph, 0);

  for (int node = 0; node < csr.numRows; ++node) {
    printf("nodeLevels[%d]=%d\n", node, nodeLevels[node]);
  }

  free(nodeLevels);
  free_coo_matrix(coo);
  free_csr_matrix(csr);

  return 0;
}

#include <mpi.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define DEPTH 7
#define MAX_COLORS 65536

typedef struct { int w, h; unsigned char *pix; } Image;
typedef struct Node Node;
struct Node {
    uint64_t r, g, b, count;
    unsigned char depth;
    Node *child[8];
};
typedef struct {
    uint32_t prefix;
    uint64_t r, g, b, count;
} PackedLeaf;
typedef struct {
    uint32_t prefix;
    unsigned char depth, r, g, b;
} PaletteEntry;
typedef struct { Node *node; size_t leaves; int found; } Candidate;

static void die(const char *msg, int rank)
{
    fprintf(stderr, "Processo %d: %s\n", rank, msg);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

static void *xmalloc(size_t n, int rank)
{
    void *p = malloc(n ? n : 1);
    if (!p) die("memoria insuficiente", rank);
    return p;
}

static int read_num(FILE *f)
{
    int c, n;
    for (;;) {
        c = fgetc(f);
        if (c == '#') while ((c = fgetc(f)) != '\n' && c != EOF) {}
        else if (c == EOF) return -1;
        else if (!isspace((unsigned char)c)) { ungetc(c, f); break; }
    }
    return fscanf(f, "%d", &n) == 1 ? n : -1;
}

static Image *read_ppm(const char *name)
{
    FILE *f = fopen(name, "rb");
    Image *im = NULL;
    int w, h, maxv, sep;
    size_t bytes;
    if (!f) return NULL;
    if (fgetc(f) != 'P' || fgetc(f) != '6') goto done;
    w = read_num(f); h = read_num(f); maxv = read_num(f); sep = fgetc(f);
    if (w <= 0 || h <= 0 || maxv != 255 || sep == EOF ||
        !isspace((unsigned char)sep) || (size_t)w > SIZE_MAX / (size_t)h / 3)
        goto done;
    bytes = (size_t)w * (size_t)h * 3;
    im = malloc(sizeof(*im));
    if (!im) goto done;
    im->pix = malloc(bytes);
    if (!im->pix || fread(im->pix, 1, bytes, f) != bytes) {
        free(im->pix); free(im); im = NULL; goto done;
    }
    im->w = w; im->h = h;
done:
    fclose(f);
    return im;
}

static int write_ppm(const Image *im, const char *name)
{
    FILE *f = fopen(name, "wb");
    size_t bytes = (size_t)im->w * (size_t)im->h * 3;
    int ok;
    if (!f) return 0;
    ok = fprintf(f, "P6\n%d %d\n255\n", im->w, im->h) >= 0 &&
         fwrite(im->pix, 1, bytes, f) == bytes;
    if (fclose(f)) ok = 0;
    return ok;
}

static void free_image(Image *im)
{
    if (im) { free(im->pix); free(im); }
}

static Node *node_new(unsigned char depth, int rank)
{
    Node *n = calloc(1, sizeof(*n));
    if (!n) die("memoria insuficiente para a octree", rank);
    n->depth = depth;
    return n;
}

static int index_at(const unsigned char *p, int bit)
{
    return (((p[0] >> bit) & 1) << 2) |
           (((p[1] >> bit) & 1) << 1) | ((p[2] >> bit) & 1);
}

static void add(Node *n, uint64_t r, uint64_t g, uint64_t b, uint64_t count)
{
    n->r += r; n->g += g; n->b += b; n->count += count;
}

static void insert_pixel(Node *root, const unsigned char *p, int rank)
{
    Node *n = root;
    add(n, p[0], p[1], p[2], 1);
    for (int depth = 1; depth <= DEPTH; depth++) {
        int idx = index_at(p, 8 - depth);
        if (!n->child[idx]) n->child[idx] = node_new((unsigned char)depth, rank);
        n = n->child[idx];
        add(n, p[0], p[1], p[2], 1);
    }
}

static int has_children(const Node *n)
{
    for (int i = 0; i < 8; i++) if (n->child[i]) return 1;
    return 0;
}

static void tree_free(Node *n)
{
    if (!n) return;
    for (int i = 0; i < 8; i++) tree_free(n->child[i]);
    free(n);
}

static size_t leaf_count(const Node *n)
{
    size_t total = 0;
    if (!n) return 0;
    if (!has_children(n)) return n->count != 0;
    for (int i = 0; i < 8; i++) total += leaf_count(n->child[i]);
    return total;
}

static void pack_leaves(const Node *n, uint32_t prefix, PackedLeaf *out, size_t *at)
{
    if (!has_children(n)) {
        if (n->count) out[(*at)++] = (PackedLeaf){prefix, n->r, n->g, n->b, n->count};
        return;
    }
    for (int i = 0; i < 8; i++) if (n->child[i])
        pack_leaves(n->child[i], (prefix << 3) | (uint32_t)i, out, at);
}

static void merge_leaf(Node *root, const PackedLeaf *leaf, int rank)
{
    Node *n = root;
    add(n, leaf->r, leaf->g, leaf->b, leaf->count);
    for (int depth = 1; depth <= DEPTH; depth++) {
        int idx = (int)((leaf->prefix >> (3 * (DEPTH - depth))) & 7U);
        if (!n->child[idx]) n->child[idx] = node_new((unsigned char)depth, rank);
        n = n->child[idx];
        add(n, leaf->r, leaf->g, leaf->b, leaf->count);
    }
}

static int better(const Node *n, size_t leaves, const Candidate *best, size_t budget)
{
    size_t reduction = leaves - 1;
    if (!best->found) return 1;
    size_t best_reduction = best->leaves - 1;
    int fits = reduction <= budget, best_fits = best_reduction <= budget;
    if (fits != best_fits) return fits;
    if (!fits && reduction != best_reduction) return reduction < best_reduction;
    if (n->depth != best->node->depth) return n->depth > best->node->depth;
    return n->count < best->node->count;
}

static size_t find_candidate(Node *n, size_t budget, Candidate *best)
{
    size_t leaves = 0;

    if (!n)
        return 0;
    if (!has_children(n))
        return n->count != 0;

    for (int i = 0; i < 8; i++)
        leaves += find_candidate(n->child[i], budget, best);
    if (leaves >= 2 && better(n, leaves, best, budget))
        *best = (Candidate){n, leaves, 1};
    return leaves;
}

static size_t reduce_tree(Node *root, size_t max_colors)
{
    size_t leaves = leaf_count(root);
    while (leaves > max_colors) {
        Candidate c = {0};
        find_candidate(root, leaves - max_colors, &c);
        if (!c.found) break;
        for (int i = 0; i < 8; i++) { tree_free(c.node->child[i]); c.node->child[i] = NULL; }
        leaves -= c.leaves - 1;
    }
    return leaves;
}

static void collect_palette(const Node *n, uint32_t prefix,
                            PaletteEntry *out, size_t *at)
{
    if (!has_children(n)) {
        if (n->count) {
            PaletteEntry *e = &out[(*at)++];
            e->prefix = prefix; e->depth = n->depth;
            e->r = (unsigned char)((n->r + n->count / 2) / n->count);
            e->g = (unsigned char)((n->g + n->count / 2) / n->count);
            e->b = (unsigned char)((n->b + n->count / 2) / n->count);
        }
        return;
    }
    for (int i = 0; i < 8; i++) if (n->child[i])
        collect_palette(n->child[i], (prefix << 3) | (uint32_t)i, out, at);
}

static uint32_t pixel_prefix(const unsigned char *p)
{
    uint32_t prefix = 0;
    for (int depth = 1; depth <= DEPTH; depth++)
        prefix = (prefix << 3) | (uint32_t)index_at(p, 8 - depth);
    return prefix;
}

static void quantize(unsigned char *pix, int pixels, const PaletteEntry *pal,
                     int colors, int rank)
{
    for (int i = 0; i < pixels; i++) {
        unsigned char *p = pix + 3 * i;
        uint32_t prefix = pixel_prefix(p);
        int found = 0;
        for (int j = 0; j < colors; j++) {
            int shift = 3 * (DEPTH - pal[j].depth);
            if ((prefix >> shift) == pal[j].prefix) {
                p[0] = pal[j].r; p[1] = pal[j].g; p[2] = pal[j].b;
                found = 1; break;
            }
        }
        if (!found) die("pixel fora das regioes da paleta", rank);
    }
}

/* PCAM: particionamento/aglomeração/mapeamento em blocos balanceados. */
static void distribution(int pixels, int processes, int *counts, int *displs)
{
    int base = pixels / processes, remainder = pixels % processes, at = 0;
    for (int i = 0; i < processes; i++) {
        counts[i] = (base + (i < remainder)) * 3;
        displs[i] = at; at += counts[i];
    }
}

static MPI_Datatype packed_leaf_type(void)
{
    MPI_Datatype raw, result;
    int lengths[2] = {1, 4};
    MPI_Aint displs[2] = {offsetof(PackedLeaf, prefix), offsetof(PackedLeaf, r)};
    MPI_Datatype types[2] = {MPI_UINT32_T, MPI_UINT64_T};
    MPI_Type_create_struct(2, lengths, displs, types, &raw);
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(PackedLeaf), &result);
    MPI_Type_commit(&result); MPI_Type_free(&raw);
    return result;
}

static int process_image(const char *input, const char *output, int max_colors,
                         int rank, int processes, MPI_Datatype leaf_type,
                         double *communication_time, double *computation_time,
                         double *io_time)
{
    Image *im = NULL;
    int dims[2] = {0}, total, local_bytes, local_pixels, status = EXIT_SUCCESS;
    int *counts, *displs, local_n, *leaf_counts = NULL, *leaf_displs = NULL;
    int all_n = 0, palette_n = 0;
    unsigned char *local;
    Node *tree;
    PackedLeaf *local_leaves, *all_leaves = NULL;
    PaletteEntry *palette = NULL;
    size_t n, at = 0;
    double image_start = MPI_Wtime();
    double communication_start = *communication_time;
    double io_start = *io_time;
    double phase_start;
    MPI_Request request;

    if (rank == 0) {
        phase_start = MPI_Wtime();
        im = read_ppm(input);
        *io_time += MPI_Wtime() - phase_start;
        if (im) { dims[0] = im->w; dims[1] = im->h; }
        else { fprintf(stderr, "Erro ao ler: %s\n", input); dims[0] = -1; }
    }
    phase_start = MPI_Wtime();
    MPI_Ibcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    if (dims[0] < 0) { free_image(im); return EXIT_FAILURE; }
    if ((size_t)dims[0] * (size_t)dims[1] > INT_MAX / 3)
        die("imagem grande demais para as coletivas MPI", rank);
    total = dims[0] * dims[1];
    counts = xmalloc((size_t)processes * sizeof(*counts), rank);
    displs = xmalloc((size_t)processes * sizeof(*displs), rank);
    distribution(total, processes, counts, displs);
    local_bytes = counts[rank]; local_pixels = local_bytes / 3;
    local = xmalloc((size_t)local_bytes, rank);

    /* C do PCAM: Scatterv, folhas compactadas, broadcast da paleta e Gatherv. */
    phase_start = MPI_Wtime();
    MPI_Iscatterv(rank == 0 ? im->pix : NULL, counts, displs, MPI_BYTE,
                  local, local_bytes, MPI_BYTE, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    tree = node_new(0, rank);
    for (int i = 0; i < local_pixels; i++) insert_pixel(tree, local + 3 * i, rank);
    n = leaf_count(tree);
    if (n > INT_MAX) die("folhas locais demais", rank);
    local_n = (int)n;
    local_leaves = xmalloc(n * sizeof(*local_leaves), rank);
    pack_leaves(tree, 0, local_leaves, &at); tree_free(tree);

    if (rank == 0) leaf_counts = xmalloc((size_t)processes * sizeof(int), rank);
    phase_start = MPI_Wtime();
    MPI_Igather(&local_n, 1, MPI_INT, leaf_counts, 1, MPI_INT, 0,
                MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    if (rank == 0) {
        leaf_displs = xmalloc((size_t)processes * sizeof(int), rank);
        for (int i = 0; i < processes; i++) {
            leaf_displs[i] = all_n;
            if (leaf_counts[i] > INT_MAX - all_n) die("folhas globais demais", rank);
            all_n += leaf_counts[i];
        }
        all_leaves = xmalloc((size_t)all_n * sizeof(*all_leaves), rank);
    }
    phase_start = MPI_Wtime();
    MPI_Igatherv(local_leaves, local_n, leaf_type, all_leaves, leaf_counts,
                 leaf_displs, leaf_type, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    free(local_leaves);

    if (rank == 0) {
        tree = node_new(0, rank);
        for (int i = 0; i < all_n; i++) merge_leaf(tree, &all_leaves[i], rank);
        /* A octree global e reduzida para n cores antes do broadcast. */
        n = reduce_tree(tree, (size_t)max_colors);
        palette_n = (int)n; palette = xmalloc(n * sizeof(*palette), rank); at = 0;
        collect_palette(tree, 0, palette, &at); tree_free(tree);
    }
    phase_start = MPI_Wtime();
    MPI_Ibcast(&palette_n, 1, MPI_INT, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    if (rank != 0) palette = xmalloc((size_t)palette_n * sizeof(*palette), rank);
    if ((size_t)palette_n > INT_MAX / sizeof(*palette)) die("paleta grande demais", rank);
    phase_start = MPI_Wtime();
    MPI_Ibcast(palette, palette_n * (int)sizeof(*palette), MPI_BYTE, 0,
               MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    quantize(local, local_pixels, palette, palette_n, rank);
    phase_start = MPI_Wtime();
    MPI_Igatherv(local, local_bytes, MPI_BYTE, rank == 0 ? im->pix : NULL,
                 counts, displs, MPI_BYTE, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;

    if (rank == 0) {
        phase_start = MPI_Wtime();
        if (!write_ppm(im, output)) { fprintf(stderr, "Erro ao escrever: %s\n", output); status = EXIT_FAILURE; }
        else printf("%s -> %s: %d cores, %d processos\n", input, output, palette_n, processes);
        *io_time += MPI_Wtime() - phase_start;
    }
    phase_start = MPI_Wtime();
    MPI_Ibcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    *communication_time += MPI_Wtime() - phase_start;
    free(palette); free(all_leaves); free(leaf_displs); free(leaf_counts);
    free(local); free(displs); free(counts); free_image(im);

    {
        double elapsed = MPI_Wtime() - image_start;
        double communication = *communication_time - communication_start;
        double io = *io_time - io_start;
        double computation = elapsed - communication - io;
        if (computation > 0.0) *computation_time += computation;
    }
    return status;
}

static int positive_int(const char *s, int *value)
{
    char *end; long n = strtol(s, &end, 10);
    if (!*s || *end || n <= 0 || n > INT_MAX) return 0;
    *value = (int)n; return 1;
}

static void usage(const char *program)
{
    fprintf(stderr,
            "Uso: mpirun -np P %s n_cores P entrada1.ppm saida1.ppm "
            "[entrada2.ppm saida2.ppm ...]\n",
            program);
}

static void print_timing_stats(const char *metric, const double *values, int count)
{
    double min = values[0], max = values[0], sum = 0.0, variance = 0.0;

    for (int i = 0; i < count; i++) {
        if (values[i] < min) min = values[i];
        if (values[i] > max) max = values[i];
        sum += values[i];
    }

    double mean = sum / count;
    for (int i = 0; i < count; i++) {
        double difference = values[i] - mean;
        variance += difference * difference;
    }

    printf("TIMING_STATS metric=%s min=%.9f max=%.9f mean=%.9f stddev=%.9f\n",
           metric, min, max, mean, sqrt(variance / count));
}

int main(int argc, char **argv)
{
    int rank, processes, colors, requested, result = EXIT_SUCCESS;
    double start, local_total, phase_start;
    double communication_time = 0.0, computation_time = 0.0, io_time = 0.0;
    double local_timings[3], *all_timings = NULL;
    MPI_Datatype leaf_type;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &processes);
    if (argc < 5 || (argc - 3) % 2 || !positive_int(argv[1], &colors) ||
        !positive_int(argv[2], &requested) || colors > MAX_COLORS) {
        if (rank == 0)
            usage(argv[0]);
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    if (requested != processes) {
        if (rank == 0)
            fprintf(stderr, "Solicitados %d processos, mas MPI iniciou %d.\n",
                    requested, processes);
        MPI_Finalize(); return EXIT_FAILURE;
    }
    leaf_type = packed_leaf_type();
    MPI_Ibarrier(MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    start = MPI_Wtime();
    for (int i = 3; i < argc; i += 2)
        if (process_image(argv[i], argv[i + 1], colors, rank, processes, leaf_type,
                          &communication_time, &computation_time, &io_time))
            result = EXIT_FAILURE;

    phase_start = MPI_Wtime();
    MPI_Ibarrier(MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    communication_time += MPI_Wtime() - phase_start;
    local_total = MPI_Wtime() - start;

    for (int current_rank = 0; current_rank < processes; current_rank++) {
        if (rank == current_rank) {
            printf("RANK_TIMING rank=%d total=%.9f comp=%.9f comm=%.9f io=%.9f\n",
                   rank, local_total, computation_time, communication_time, io_time);
            fflush(stdout);
        }
        MPI_Ibarrier(MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    local_timings[0] = local_total;
    local_timings[1] = computation_time;
    local_timings[2] = communication_time;
    if (rank == 0)
        all_timings = xmalloc((size_t)processes * 3 * sizeof(*all_timings), rank);
    MPI_Igather(local_timings, 3, MPI_DOUBLE, all_timings, 3, MPI_DOUBLE,
                0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    if (rank == 0) {
        double *values = xmalloc((size_t)processes * sizeof(*values), rank);
        for (int metric = 0; metric < 3; metric++) {
            for (int process = 0; process < processes; process++)
                values[process] = all_timings[3 * process + metric];
            print_timing_stats(metric == 0 ? "total" : metric == 1 ? "comp" : "comm",
                               values, processes);
        }
        free(values);
    }
    free(all_timings);
    MPI_Type_free(&leaf_type); MPI_Finalize(); return result;
}

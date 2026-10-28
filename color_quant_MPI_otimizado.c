#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

typedef struct {
	int w, h;
	unsigned char *pix;
} image_t, *image;

int write_ppm(image im, char *fn)
{
	FILE *fp = fopen(fn, "wb");
	if (!fp) return 0;
	size_t bytes = (size_t)3 * (size_t)im->w * (size_t)im->h;
	int ok = fprintf(fp, "P6\n%d %d\n255\n", im->w, im->h) >= 0
	      && fwrite(im->pix, 1, bytes, fp) == bytes;
	if (fclose(fp) != 0) ok = 0;
	return ok;
}

image img_new(int w, int h)
{
	image im = malloc(sizeof(image_t) + (size_t)h * (size_t)w * 3);
	if (!im) return NULL;
	im->w = w; im->h = h;
	im->pix = (unsigned char *)(im + 1);
	return im;
}

int read_num(FILE *f)
{
	int n;
	while (!fscanf(f, "%d ", &n)) {
		if ((n = fgetc(f)) == '#') {
			while ((n = fgetc(f)) != '\n')
				if (n == EOF) return 0;
		} else
			return 0;
	}
	return n;
}

image read_ppm(char *fn)
{
	FILE *fp = fopen(fn, "rb");
	int w, h, maxval;
	image im = 0;
	if (!fp) return 0;

	if (fgetc(fp) != 'P' || fgetc(fp) != '6' || !isspace(fgetc(fp)))
		goto bail;

	w = read_num(fp);
	h = read_num(fp);
	maxval = read_num(fp);
	if (w <= 0 || h <= 0 || maxval != 255) goto bail;

	im = img_new(w, h);
	if (!im) goto bail;
	size_t bytes = (size_t)3 * (size_t)w * (size_t)h;
	if (fread(im->pix, 1, bytes, fp) != bytes) {
		free(im);
		im = NULL;
	}
bail:
	if (fp) fclose(fp);
	return im;
}


#define ON_INHEAP	1

typedef struct oct_node_t oct_node_t, *oct_node;
struct oct_node_t{
	int64_t r, g, b; /* sum of all child node colors */
	int count, heap_idx;
	unsigned char n_kids, kid_idx, flags, depth;
	oct_node kids[8], parent;
};

typedef struct {
	int alloc, n;
	oct_node* buf;
} node_heap;

int cmp_node(oct_node a, oct_node b)
{
	if (a->n_kids < b->n_kids) return -1;
	if (a->n_kids > b->n_kids) return 1;

	int ac = a->count >> a->depth;
	int bc = b->count >> b->depth;
	return ac < bc ? -1 : ac > bc;
}

void down_heap(node_heap *h, oct_node p)
{
	int n = p->heap_idx, m;
	while (1) {
		m = n * 2;
		if (m >= h->n) break;
		if (m + 1 < h->n && cmp_node(h->buf[m], h->buf[m + 1]) > 0) m++;

		if (cmp_node(p, h->buf[m]) <= 0) break;

		h->buf[n] = h->buf[m];
		h->buf[n]->heap_idx = n;
		n = m;
	}
	h->buf[n] = p;
	p->heap_idx = n;
}

void up_heap(node_heap *h, oct_node p)
{
	int n = p->heap_idx;
	oct_node prev;

	while (n > 1) {
		prev = h->buf[n / 2];
		if (cmp_node(p, prev) >= 0) break;

		h->buf[n] = prev;
		prev->heap_idx = n;
		n /= 2;
	}
	h->buf[n] = p;
	p->heap_idx = n;
}

void heap_add(node_heap *h, oct_node p)
{
	if ((p->flags & ON_INHEAP)) {
		down_heap(h, p);
		up_heap(h, p);
		return;
	}

	p->flags |= ON_INHEAP;
	if (!h->n) h->n = 1;
	if (h->n >= h->alloc) {
		while (h->n >= h->alloc) h->alloc += 1024;
		h->buf = realloc(h->buf, sizeof(oct_node) * h->alloc);
	}

	p->heap_idx = h->n;
	h->buf[h->n++] = p;
	up_heap(h, p);
}

oct_node pop_heap(node_heap *h)
{
	if (h->n <= 1) return 0;

	oct_node ret = h->buf[1];
	h->buf[1] = h->buf[--h->n];

	h->buf[h->n] = 0;

	h->buf[1]->heap_idx = 1;
	down_heap(h, h->buf[1]);

	return ret;
}

static oct_node pool = 0;
static int pool_len = 0;
oct_node node_new(unsigned char idx, unsigned char depth, oct_node p)
{
	if (pool_len <= 1) {
		oct_node p = calloc(sizeof(oct_node_t), 2048);
		p->parent = pool;
		pool = p;
		pool_len = 2047;
	}

	oct_node x = pool + pool_len--;
	x->kid_idx = idx;
	x->depth = depth;
	x->parent = p;
	if (p) p->n_kids++;
	return x;
}

void node_free()
{
	oct_node p;
	while (pool) {
		p = pool->parent;
		free(pool);
		pool = p;
	}
	pool_len = 0;
}

oct_node node_insert(oct_node root, unsigned char *pix)
{
	unsigned char i, bit, depth = 0;

	for (bit = 1 << 7; ++depth < 8; bit >>= 1) {
		i = !!(pix[1] & bit) * 4 + !!(pix[0] & bit) * 2 + !!(pix[2] & bit);
		if (!root->kids[i])
			root->kids[i] = node_new(i, depth, root);

		root = root->kids[i];
	}

	root->r += pix[0];
	root->g += pix[1];
	root->b += pix[2];
	root->count++;
	return root;
}

oct_node node_fold(oct_node p)
{
	if (p->n_kids) abort();
	oct_node q = p->parent;
	q->count += p->count;

	q->r += p->r;
	q->g += p->g;
	q->b += p->b;
	q->n_kids --;
	q->kids[p->kid_idx] = 0;
	return q;
}

void color_replace(oct_node root, unsigned char *pix)
{
	unsigned char i, bit;

	for (bit = 1 << 7; bit; bit >>= 1) {
		i = !!(pix[1] & bit) * 4 + !!(pix[0] & bit) * 2 + !!(pix[2] & bit);
		if (!root->kids[i]) break;
		root = root->kids[i];
	}

	pix[0] = root->r;
	pix[1] = root->g;
	pix[2] = root->b;
}

void error_diffuse(image im, node_heap *h)
{
	oct_node nearest_color(int *v) {
		int i;
		int diff, max = 100000000;
		oct_node o = 0;
		for (i = 1; i < h->n; i++) {
			diff =	  3 * abs(h->buf[i]->r - v[0])
				+ 5 * abs(h->buf[i]->g - v[1])
				+ 2 * abs(h->buf[i]->b - v[2]);
			if (diff < max) {
				max = diff;
				o = h->buf[i];
			}
		}
		return o;
	}

#	define POS(i, j) (3 * ((i) * im->w + (j)))
	int i, j;
	int *npx = calloc(sizeof(int), im->h * im->w * 3), *px;
	int v[3];
	unsigned char *pix = im->pix;
	oct_node nd;

#define C10 7
#define C01 5
#define C11 2
#define C00 1 
#define CTOTAL (C00 + C11 + C10 + C01)

	for (px = npx, i = 0; i < im->h; i++) {
		for (j = 0; j < im->w; j++, pix += 3, px += 3) {
			px[0] = (int)pix[0] * CTOTAL;
			px[1] = (int)pix[1] * CTOTAL;
			px[2] = (int)pix[2] * CTOTAL;
		}
	}
#define clamp(x, i) if (x[i] > 255) x[i] = 255; if (x[i] < 0) x[i] = 0
	pix = im->pix;
	for (px = npx, i = 0; i < im->h; i++) {
		for (j = 0; j < im->w; j++, pix += 3, px += 3) {
			px[0] /= CTOTAL;
			px[1] /= CTOTAL;
			px[2] /= CTOTAL;
			clamp(px, 0); clamp(px, 1); clamp(px, 2);

			nd = nearest_color(px);

			v[0] = px[0] - nd->r;
			v[1] = px[1] - nd->g;
			v[2] = px[2] - nd->b;

			pix[0] = nd->r; pix[1] = nd->g; pix[2] = nd->b;
			if (j < im->w - 1) {
				npx[POS(i, j+1) + 0] += v[0] * C10;
				npx[POS(i, j+1) + 1] += v[1] * C10;
				npx[POS(i, j+1) + 2] += v[2] * C10;
			}
			if (i >= im->h - 1) continue;

			npx[POS(i+1, j) + 0] += v[0] * C01;
			npx[POS(i+1, j) + 1] += v[1] * C01;
			npx[POS(i+1, j) + 2] += v[2] * C01;

			if (j < im->w - 1) {
				npx[POS(i+1, j+1) + 0] += v[0] * C11;
				npx[POS(i+1, j+1) + 1] += v[1] * C11;
				npx[POS(i+1, j+1) + 2] += v[2] * C11;
			}
			if (j) {
				npx[POS(i+1, j-1) + 0] += v[0] * C00;
				npx[POS(i+1, j-1) + 1] += v[1] * C00;
				npx[POS(i+1, j-1) + 2] += v[2] * C00;
			}
		}
	}
	free(npx);
}

void color_quant(image im, int n_colors, int dither)
{
	int i;
	unsigned char *pix = im->pix;
	node_heap heap = { 0, 0, 0 };

	oct_node root = node_new(0, 0, 0), got;
	for (i = 0; i < im->w * im->h; i++, pix += 3)
		heap_add(&heap, node_insert(root, pix));

	while (heap.n > n_colors + 1)
		heap_add(&heap, node_fold(pop_heap(&heap)));

	double c;
	for (i = 1; i < heap.n; i++) {
		got = heap.buf[i];
		c = got->count;
		got->r = got->r / c + .5;
		got->g = got->g / c + .5;
		got->b = got->b / c + .5;
	}

	if (dither) error_diffuse(im, &heap);
	else
		for (i = 0, pix = im->pix; i < im->w * im->h; i++, pix += 3)
			color_replace(root, pix);

	node_free();
	free(heap.buf);
}



#define PALETTE_DEPTH 7
#define LUT_ENTRIES (1U << (3 * PALETTE_DEPTH))
#define LUT_BYTES (3U * LUT_ENTRIES)

static int palette_index_opt(const unsigned char *pixel, int bit)
{
    return (!!(pixel[1] & (1U << bit)) << 2) |
           (!!(pixel[0] & (1U << bit)) << 1) |
           !!(pixel[2] & (1U << bit));
}

/*
 * Constrói a paleta com o mesmo heap e a mesma ordem de inserção da versão
 * serial. A LUT representa exatamente a caminhada de color_replace() para
 * todos os prefixos RGB de 7 bits, inclusive nós internos parcialmente
 * reduzidos pelo heap serial.
 */
static unsigned char *build_lut_opt(image im, int colors, int *palette_count)
{
    node_heap heap = {0, 0, 0};
    oct_node root = node_new(0, 0, 0), current;
    unsigned char *pixel = im->pix;
    int pixels = im->w * im->h;

    for (int i = 0; i < pixels; i++, pixel += 3)
        heap_add(&heap, node_insert(root, pixel));

    while (heap.n > colors + 1)
        heap_add(&heap, node_fold(pop_heap(&heap)));

    for (int i = 1; i < heap.n; i++) {
        current = heap.buf[i];
        double count = current->count;
        current->r = current->r / count + .5;
        current->g = current->g / count + .5;
        current->b = current->b / count + .5;
    }

    *palette_count = heap.n - 1;
    unsigned char *lut = malloc(LUT_BYTES);
    if (!lut) {
        fprintf(stderr, "Memória insuficiente para a LUT.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (uint32_t prefix = 0; prefix < LUT_ENTRIES; prefix++) {
        oct_node node = root;
        for (int depth = 1; depth <= PALETTE_DEPTH; depth++) {
            int shift = 3 * (PALETTE_DEPTH - depth);
            int child = (int)((prefix >> shift) & 7U);
            if (!node->kids[child]) break;
            node = node->kids[child];
        }
        lut[3U * prefix] = (unsigned char)node->r;
        lut[3U * prefix + 1] = (unsigned char)node->g;
        lut[3U * prefix + 2] = (unsigned char)node->b;
    }

    node_free();
    free(heap.buf);
    return lut;
}

static void replace_pixels_opt(unsigned char *pixels, int pixel_count,
                               const unsigned char *lut)
{
    for (int i = 0; i < pixel_count; i++) {
        unsigned char *pixel = pixels + 3 * i;
        uint32_t prefix = 0;
        for (int depth = 1; depth <= PALETTE_DEPTH; depth++)
            prefix = (prefix << 3) |
                     (uint32_t)palette_index_opt(pixel, 8 - depth);

        pixel[0] = lut[3U * prefix];
        pixel[1] = lut[3U * prefix + 1];
        pixel[2] = lut[3U * prefix + 2];
    }
}

static void distribution_opt(int pixels, int processes, int *counts, int *displs)
{
    int base = pixels / processes, remainder = pixels % processes, offset = 0;
    for (int rank = 0; rank < processes; rank++) {
        counts[rank] = (base + (rank < remainder)) * 3;
        displs[rank] = offset;
        offset += counts[rank];
    }
}

static int positive_int_opt(const char *text, int *value)
{
    char *end;
    long parsed = strtol(text, &end, 10);
    if (!*text || *end || parsed <= 0 || parsed > INT32_MAX) return 0;
    *value = (int)parsed;
    return 1;
}

static void print_stats_opt(const char *metric, const double *values, int count)
{
    double min = values[0], max = values[0], sum = 0.0, variance = 0.0;
    for (int i = 0; i < count; i++) {
        if (values[i] < min) min = values[i];
        if (values[i] > max) max = values[i];
        sum += values[i];
    }
    double mean = sum / count;
    for (int i = 0; i < count; i++) {
        double delta = values[i] - mean;
        variance += delta * delta;
    }
    printf("TIMING_STATS metric=%s min=%.9f max=%.9f mean=%.9f stddev=%.9f\n",
           metric, min, max, mean, sqrt(variance / count));
}

static void usage_opt(const char *program)
{
    fprintf(stderr,
            "Uso: mpirun -np P %s n_cores P entrada1.ppm saida1.ppm "
            "[entrada2.ppm saida2.ppm ...]\n",
            program);
}

static int process_image_opt(const char *input, const char *output, int colors,
                             int rank, int processes, double *comp_time,
                             double *comm_time, double *io_time)
{
    image im = NULL;
    int dims[2] = {0, 0}, status = EXIT_SUCCESS;
    int total_pixels, local_bytes, local_pixels, palette_count = 0;
    int *counts = NULL, *displs = NULL;
    unsigned char *local_pixels_data = NULL;
    unsigned char *lut = NULL;
    double phase_start;

    if (rank == 0) {
        phase_start = MPI_Wtime();
        im = read_ppm((char *)input);
        *io_time += MPI_Wtime() - phase_start;
        if (im) {
            dims[0] = im->w;
            dims[1] = im->h;
        } else {
            fprintf(stderr, "Erro ao ler a imagem: %s\n", input);
            dims[0] = -1;
        }
    }

    phase_start = MPI_Wtime();
    MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
    *comm_time += MPI_Wtime() - phase_start;
    if (dims[0] < 0) return EXIT_FAILURE;

    if ((size_t)dims[0] * (size_t)dims[1] > INT32_MAX / 3) {
        if (rank == 0) fprintf(stderr, "Imagem grande demais para MPI_Scatterv.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    total_pixels = dims[0] * dims[1];
    counts = malloc((size_t)processes * sizeof(*counts));
    displs = malloc((size_t)processes * sizeof(*displs));
    if (!counts || !displs) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    distribution_opt(total_pixels, processes, counts, displs);
    local_bytes = counts[rank];
    local_pixels = local_bytes / 3;
    local_pixels_data = malloc((size_t)(local_bytes ? local_bytes : 1));
    if (!local_pixels_data) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    phase_start = MPI_Wtime();
    MPI_Scatterv(rank == 0 ? im->pix : NULL, counts, displs, MPI_BYTE,
                 local_pixels_data, local_bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    *comm_time += MPI_Wtime() - phase_start;

    if (rank == 0) {
        phase_start = MPI_Wtime();
        lut = build_lut_opt(im, colors, &palette_count);
        *comp_time += MPI_Wtime() - phase_start;
    }

    phase_start = MPI_Wtime();
    MPI_Bcast(&palette_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    *comm_time += MPI_Wtime() - phase_start;
    if (rank != 0) {
        lut = malloc(LUT_BYTES);
        if (!lut) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    phase_start = MPI_Wtime();
    MPI_Bcast(lut, (int)LUT_BYTES, MPI_BYTE, 0, MPI_COMM_WORLD);
    *comm_time += MPI_Wtime() - phase_start;

    phase_start = MPI_Wtime();
    replace_pixels_opt(local_pixels_data, local_pixels, lut);
    *comp_time += MPI_Wtime() - phase_start;

    phase_start = MPI_Wtime();
    MPI_Gatherv(local_pixels_data, local_bytes, MPI_BYTE,
                rank == 0 ? im->pix : NULL, counts, displs, MPI_BYTE,
                0, MPI_COMM_WORLD);
    *comm_time += MPI_Wtime() - phase_start;

    if (rank == 0) {
        phase_start = MPI_Wtime();
        if (!write_ppm(im, (char *)output)) {
            fprintf(stderr, "Erro ao escrever a imagem: %s\n", output);
            status = EXIT_FAILURE;
        } else {
            printf("%s -> %s: %d cores, %d processos\n",
                   input, output, palette_count, processes);
            fflush(stdout);
        }
        *io_time += MPI_Wtime() - phase_start;
    }

    phase_start = MPI_Wtime();
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    *comm_time += MPI_Wtime() - phase_start;

    free(lut);
    free(local_pixels_data);
    free(displs);
    free(counts);
    if (rank == 0) free(im);
    return status;
}

int main(int argc, char **argv)
{
    int rank, processes, colors, requested, result = EXIT_SUCCESS;
    double computation_time = 0.0, communication_time = 0.0, io_time = 0.0;
    double start, phase_start, local_total;
    double local_timings[4], *all_timings = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    if (argc < 5 || (argc - 3) % 2 != 0 ||
        !positive_int_opt(argv[1], &colors) ||
        !positive_int_opt(argv[2], &requested)) {
        if (rank == 0) usage_opt(argv[0]);
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    if (requested != processes) {
        if (rank == 0)
            fprintf(stderr, "Solicitados %d processos, mas MPI iniciou %d.\n",
                    requested, processes);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    for (int argument = 3; argument < argc; argument += 2)
        if (process_image_opt(argv[argument], argv[argument + 1], colors,
                              rank, processes, &computation_time,
                              &communication_time, &io_time))
            result = EXIT_FAILURE;

    phase_start = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    communication_time += MPI_Wtime() - phase_start;
    local_total = MPI_Wtime() - start;

    printf("RANK_TIMING rank=%d total=%.9f comp=%.9f comm=%.9f io=%.9f\n",
           rank, local_total, computation_time, communication_time, io_time);
    fflush(stdout);

    local_timings[0] = local_total;
    local_timings[1] = computation_time;
    local_timings[2] = communication_time;
    local_timings[3] = io_time;
    if (rank == 0)
        all_timings = malloc((size_t)processes * 4 * sizeof(*all_timings));
    MPI_Gather(local_timings, 4, MPI_DOUBLE, all_timings, 4, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        double *values = malloc((size_t)processes * sizeof(*values));
        if (!all_timings || !values) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        for (int metric = 0; metric < 4; metric++) {
            for (int process = 0; process < processes; process++)
                values[process] = all_timings[4 * process + metric];
            print_stats_opt(metric == 0 ? "total" :
                            metric == 1 ? "comp" :
                            metric == 2 ? "comm" : "io",
                            values, processes);
        }
        free(values);
    }

    free(all_timings);
    MPI_Finalize();
    return result;
}


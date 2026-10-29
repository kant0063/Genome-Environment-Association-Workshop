################################################################################
# candidate_genes_generic.R
# Species-agnostic candidate gene search around significant association markers
# (GWAS / EGWAS / GEA results), using a VCF, a gene annotation (GFF3 or GTF) and,
# optionally, functional annotation tables, protein sequences and LD estimates.
#
# Usage:  Rscript candidate_genes_generic.R [input_dir] [output_dir]
#   input_dir  : folder searched (recursively) for inputs; default "."
#   output_dir : default "candidate_genes_out"
# Inputs are found by file-name pattern (gzipped or not); any can be set explicitly in FM.
#
#   required  VCF                    *.vcf / *.vcf.gz
#   required  association results    *EGWAS* / *GWAS* / *GEA* / *assoc*  (.csv/.tsv/.txt)
#   required  gene annotation        *.gff3 / *.gff / *.gtf  (gz ok; files named *exon* skipped)
#   optional  functional annotation  *annotation_info* / *defline* / *annot* / *description*:
#                                    any table with a column of gene or transcript IDs
#                                    (the ID column is detected automatically)
#   optional  proteins (BLAST only)  *protein*primary* / *pep*.fa / *protein*.fa
#   optional  LD decay               *LD_decay*.csv with columns set, half_decay_bp;
#                                    set = chromosome name, subgenome label or "genome"
#   optional  chromosome map         *chr_map*.tsv: 2 columns, no header: VCF name, annotation name
#
# Association table: one row per significant marker x trait (x set). Recognised column
# names (case-insensitive): marker SNP/marker/rs/ID; chromosome CHROM/CHR/chromosome;
# position POS/BP/position; trait trait/variable/phenotype; set set/group/population;
# q q/qvalue/FDR/padj; p p/pvalue/p_wald; PVE PVE/R2. Override with FM$col_* if needed.
# Only a marker column (matching VCF IDs) or chromosome + position is required.
#
# Priority per marker: 2 if significant in any set other than FM$full_set (within-group
# signal), +1 if significant for more than one trait. Markers with priority >=
# FM$min_priority are searched; if none qualify, the FM$top_n best by q (or p) are used.
#
# Reported per gene: distance; SNP position (CDS, 5'UTR, 3'UTR, intron, noncoding exon,
# upstream, downstream; "5'/3'_untranslated" when UTRs and exons are not annotated);
# nearest-gene and report-window flags; description from the annotation file; joined
# functional annotation; optional BLASTP best hits. Per marker: MAF, observed
# heterozygosity, low-MAF and excess-heterozygosity flags (any ploidy).
#
# Species examples for the FM block:
#   finger millet (Phytozome Ecoracana v1.1; TASSEL VCF numbering 1-19): chr_map file
#     1->1A, 2->1B, ..., 18->9B; subgenome_regex = "([AB])$";
#     ld_window = c(A = 88000, B = 168000, genome = 106000); report_window = 25000;
#     inbreeding_F = 0.5
#   wheat (IWGSC chr1A..chr7D): subgenome_regex = "([ABD])$"
#   diploids with Ensembl/NCBI/Phytozome GFF3 or GTF: defaults; set window to the LD decay
#     distance for your panel; inbreeding_F = 0 for outcrossers, ~0.5-0.9 for selfers
################################################################################
suppressPackageStartupMessages(library(data.table))

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "."
out_dir <- if (length(args) >= 2) args[2] else "candidate_genes_out"

## ---- Settings ---------------------------------------------------------------------
FM <- list(
  # input paths: NULL = auto-detect in in_dir
  vcf = NULL, assoc = NULL, gff = NULL, annotation = NULL, proteins_fa = NULL,
  ld_summary = NULL, chr_map = NULL,
  # association table columns: NULL = auto-detect
  col_snp = NULL, col_chr = NULL, col_pos = NULL, col_trait = NULL, col_set = NULL,
  col_q = NULL, col_p = NULL, col_pve = NULL,
  full_set     = "full",   # set name meaning "all individuals"
  # genome structure
  subgenome_regex = NULL,  # capture group giving the subgenome from the chromosome name, e.g. "([AB])$"
  gene_types   = c("gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"),
  tx_types     = c("mRNA", "transcript", "lnc_RNA", "ncRNA", "pseudogenic_transcript"),
  # search windows (+/- bp)
  window       = 50000,    # used when no LD value applies
  ld_window    = NULL,     # named by chromosome, subgenome label or "genome", e.g. c(A = 88000, B = 168000)
  min_window   = 0,
  max_window   = 500000,
  report_window = 25000,   # genes within this distance are flagged (in_report_window)
  min_priority = 2,
  top_n        = 20,
  # marker quality flags
  low_maf      = 0.05,
  inbreeding_F = 0,        # expected F: 0 outcrossing, higher for selfers; het > 2pq(1-F) is flagged
  # optional BLAST: named makeblastdb protein database prefixes
  blast_dbs    = c(),      # e.g. c(TAIR10 = "blastdb/TAIR10_pep", rice = "blastdb/rice_pep")
  threads      = 4)

## ---- helpers --------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tab <- function(f) file.path(out_dir, f)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
cat_cmd <- function(f) paste(if (grepl("gz$", f)) "zcat" else "cat", shQuote(f))
norm_chr <- function(x) {                        # "Chr01" / "chr1" / "1" -> "1"
  x <- sub("^(chromosome|chrom|chr)[_.-]?", "", as.character(x), ignore.case = TRUE)
  sub("^0+(?=[0-9])", "", toupper(x), perl = TRUE) }
strip_pref <- function(x) sub("^(gene|transcript|mrna|cds|protein)[:-]", "", x, ignore.case = TRUE)
url_dec <- function(x) { for (p in list(c("%2C", ","), c("%3B", ";"), c("%3D", "="), c("%20", " "), c("%25", "%")))
                           x <- gsub(p[1], p[2], x, fixed = TRUE); x }
db_exists <- function(db) any(file.exists(paste0(db, c(".pin", ".pal", ".psq", ".pdb", ".00.pin"))))
tbl_str <- function(x) { t <- table(x, useNA = "ifany"); paste(names(t), t, sep = "=", collapse = " ") }

all_files <- list.files(in_dir, recursive = TRUE, full.names = TRUE)
find_in <- function(key, patterns, exclude = NULL, required = FALSE, multi = FALSE) {
  if (!is.null(FM[[key]])) {
    miss <- FM[[key]][!file.exists(FM[[key]])]
    if (length(miss)) stop(key, " not found: ", paste(miss, collapse = ", "))
    return(FM[[key]]) }
  hits <- character(0)
  for (p in patterns) {
    f <- all_files[grepl(p, basename(all_files), ignore.case = TRUE)]
    if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f), ignore.case = TRUE)]
    if (length(f)) { hits <- c(hits, f[order(nchar(basename(f)), f)]); if (!multi) break }
  }
  if (!length(hits)) {
    if (required) stop("No ", key, " file in ", normalizePath(in_dir), " (patterns: ", paste(patterns, collapse = " | "), ")")
    return(NULL) }
  hits <- unique(hits); hits <- hits[!duplicated(tools::md5sum(hits))]   # drop identical copies
  if (multi) return(hits)
  if (length(hits) > 1) message("  ", key, ": ", length(hits), " matches, using ", hits[1])
  hits[1] }

FM$vcf   <- find_in("vcf", "\\.vcf(\\.gz)?$", required = TRUE)
FM$assoc <- find_in("assoc", c("EGWAS.*\\.(csv|tsv|txt)(\\.gz)?$", "GWAS.*\\.(csv|tsv|txt)(\\.gz)?$",
                               "GEA.*\\.(csv|tsv|txt)(\\.gz)?$", "assoc.*\\.(csv|tsv|txt)(\\.gz)?$"), required = TRUE)
FM$gff   <- find_in("gff", c("gene[._]gff3", "\\.gff3?([._]gz)?$", "_gff3?([._]gz)?$", "\\.gtf([._]gz)?$", "_gtf([._]gz)?$"),
                    exclude = "exon", required = TRUE)
FM$annotation  <- find_in("annotation", c("annotation_info", "defline", "annot", "description", "functional"),
                          exclude = "gff|gtf|\\.fa|fasta|_fa", multi = TRUE)
FM$proteins_fa <- find_in("proteins_fa", c("protein.*primary", "pep.*\\.fa", "protein.*\\.fa", "protein.*fa"))
FM$ld_summary  <- find_in("ld_summary", "LD_decay.*\\.(csv|tsv|txt)$")
FM$chr_map     <- find_in("chr_map", "chr_map.*\\.(tsv|txt)$")
fk <- c("vcf", "assoc", "gff", "annotation", "proteins_fa", "ld_summary", "chr_map")
msg("Inputs:", paste(fk, sapply(FM[fk], function(x) if (is.null(x)) "-" else paste(basename(x), collapse = ",")),
                     sep = "=", collapse = "; "))

## ---- 1. chromosome harmonisation ----------------------------------------------------------
chr_lookup <- character(0)
if (!is.null(FM$chr_map)) {
  cm <- fread(FM$chr_map, header = FALSE, colClasses = "character")
  chr_lookup <- setNames(norm_chr(cm[[2]]), norm_chr(cm[[1]]))
  msg("Chromosome map:", nrow(cm), "entries")
}
to_key <- function(x) { k <- norm_chr(x); ifelse(k %in% names(chr_lookup), chr_lookup[k], k) }
subgenome <- function(key) {
  if (is.null(FM$subgenome_regex)) return(rep(NA_character_, length(key)))
  m <- regmatches(key, regexec(FM$subgenome_regex, key))
  vapply(m, function(z) if (length(z) >= 2) z[2] else NA_character_, "") }

## ---- 2. gene annotation (GFF3 or GTF) --------------------------------------------------
gff <- fread(cmd = paste(cat_cmd(FM$gff), "| grep -v '^#'"), sep = "\t", header = FALSE, quote = "",
             col.names = c("seqid", "src", "type", "start", "end", "score", "strand", "phase", "attr"),
             colClasses = list(character = c(1:3, 6:9), integer = 4:5))
is_gtf <- !any(grepl("=", head(gff$attr, 500)))
get_att <- function(a, key) {
  pat <- if (is_gtf) paste0("(^|;)\\s*", key, " \"[^\"]*\"") else paste0("(^|;)", key, "=[^;]+")
  m <- regexpr(pat, a, ignore.case = TRUE); out <- rep(NA_character_, length(a))
  v <- regmatches(a, m)
  out[m > 0] <- if (is_gtf) sub(".*\"([^\"]*)\"", "\\1", v) else url_dec(sub("^;?[^=]+=", "", v))
  out }
gff[, key := to_key(seqid)]
gff[, ftype := fcase(type == "CDS", "CDS", type == "exon", "exon",
                     grepl("five_prime_utr|5'?utr", type, ignore.case = TRUE), "UTR5",
                     grepl("three_prime_utr|3'?utr", type, ignore.case = TRUE), "UTR3", default = NA_character_)]
span_by <- function(d, by) d[, .(seqid = seqid[1], key = key[1], start = min(start), end = max(end),
                                 strand = strand[1], attr = attr[1]), by = by]
if (is_gtf) {
  gff[, `:=`(gene_id = get_att(attr, "gene_id"), tx_id = get_att(attr, "transcript_id"))]
  tx <- gff[type %in% FM$tx_types & !is.na(tx_id)]
  if (!nrow(tx)) tx <- merge(span_by(gff[!is.na(tx_id)], "tx_id"),
                             unique(gff[!is.na(tx_id), .(tx_id, gene_id)]), by = "tx_id")
  tx[, tx_name := fcoalesce(get_att(attr, "transcript_name"), tx_id)]
  genes <- gff[type %in% FM$gene_types & !is.na(gene_id)]
  if (!nrow(genes)) genes <- span_by(tx, "gene_id")
  genes[, gene_name := fcoalesce(get_att(attr, "gene_name"), gene_id)]
  feats <- gff[!is.na(ftype) & !is.na(tx_id), .(tx_id, ftype, start, end)]
} else {
  gff[, `:=`(ID = strip_pref(get_att(attr, "ID")), Parent = get_att(attr, "Parent"))]
  tx <- gff[type %in% FM$tx_types]
  tx[, `:=`(tx_id = ID, gene_id = strip_pref(sub(",.*", "", Parent)), tx_name = fcoalesce(get_att(attr, "Name"), ID))]
  genes <- gff[type %in% FM$gene_types][, gene_id := ID]
  if (!nrow(genes) && nrow(tx)) genes <- span_by(tx, "gene_id")
  if (!nrow(tx)) tx <- genes[, .(seqid, key, start, end, strand, attr, tx_id = gene_id, gene_id, tx_name = gene_id)]
  genes[, gene_name := fcoalesce(get_att(attr, "Name"), get_att(attr, "gene_name"), gene_id)]
  ff <- gff[!is.na(ftype) & !is.na(Parent), .(Parent, ftype, start, end)]
  pl <- strsplit(ff$Parent, ",", fixed = TRUE)
  feats <- data.table(tx_id = strip_pref(unlist(pl)), ff[rep(seq_len(nrow(ff)), lengths(pl)), .(ftype, start, end)])
}
tx[, canon := grepl("canonical|MANE_Select", attr, ignore.case = TRUE) | get_att(attr, "longest") %in% "1"]
genes[, description := fcoalesce(get_att(attr, "description"), get_att(attr, "Note"),
                                 get_att(attr, "product"), get_att(attr, "gene_description"))]
genes[, description := sub(" \\[Source:.*\\]$", "", description)]
# primary transcript: flagged canonical/longest, else longest CDS, else longest span
cds_len <- feats[ftype == "CDS", .(cds_len = sum(end - start + 1)), by = tx_id]
tx <- merge(tx, cds_len, by = "tx_id", all.x = TRUE)[is.na(cds_len), cds_len := 0L]
tx[, span := end - start + 1L]
setorder(tx, gene_id, -canon, -cds_len, -span)
genes <- merge(genes, tx[!duplicated(gene_id), .(gene_id, tx_id, tx_name)], by = "gene_id", all.x = TRUE)
genes <- genes[!duplicated(gene_id)]
genes[, subgenome := subgenome(key)]
setkey(feats, tx_id)
msg("Annotation:", basename(FM$gff), if (is_gtf) "(GTF)" else "(GFF3)", "|", nrow(genes), "genes on",
    uniqueN(genes$seqid), "sequences |", sum(!is.na(genes$tx_id)), "with a primary transcript |",
    sum(!is.na(genes$description)), "with a description | features:", tbl_str(feats$ftype))

## ---- 3. VCF positions and chromosome check -----------------------------------------------
vcf <- fread(cmd = paste(cat_cmd(FM$vcf), "| grep -v '^##' | cut -f1-3"), sep = "\t",
             col.names = c("CHROM", "POS", "ID"), colClasses = c("character", "integer", "character"))
vcf[, key := to_key(CHROM)]
ext <- genes[, .(ann_max = max(end), n_genes = .N), by = key]
chk <- vcf[, .(n_variants = .N, vcf_max = max(POS)), by = .(CHROM, key)]
chk[ext, on = "key", `:=`(annotation_max = i.ann_max, n_genes = i.n_genes)]
chk[, subgenome := subgenome(key)]
setnames(chk, "key", "annotation_chr")
fwrite(chk[order(suppressWarnings(as.numeric(CHROM)), CHROM)], tab("chr_name_check.tsv"), sep = "\t")
unm <- chk[is.na(annotation_max)]
msg("VCF:", nrow(vcf), "variants on", nrow(chk), "chromosomes;", nrow(chk) - nrow(unm), "matched to the annotation",
    if (nrow(unm)) paste0("(unmatched, not searched: ", paste(head(unm$CHROM, 10), collapse = ","),
                          if (nrow(unm) > 10) ",...", ")") else "")
if (nrow(unm) == nrow(chk)) stop("No VCF chromosome matches the annotation (VCF e.g. ", chk$CHROM[1],
                                 "; annotation e.g. ", genes$seqid[1], "). Supply a chr_map file.")
bad <- chk[!is.na(annotation_max) & vcf_max > annotation_max * 1.05 + 1e6]
if (nrow(bad)) warning("VCF positions run well past the last annotated gene on: ", paste(bad$CHROM, collapse = ", "),
                       " - check the chromosome map / reference version")

## ---- 4. association table ------------------------------------------------------------
as_raw <- fread(FM$assoc)
pick_col <- function(key, aliases) {
  if (!is.null(FM[[key]])) return(FM[[key]])
  nn <- tolower(gsub("[^A-Za-z0-9]", "", names(as_raw)))
  i <- match(tolower(gsub("[^A-Za-z0-9]", "", aliases)), nn); i <- i[!is.na(i)]
  if (length(i)) names(as_raw)[i[1]] else NULL }
cols <- list(
  snp   = pick_col("col_snp", c("SNP", "marker", "markername", "snp_id", "rs", "rsid", "variant", "ID")),
  chr   = pick_col("col_chr", c("CHROM", "CHR", "chromosome", "chrom", "seqid", "contig")),
  pos   = pick_col("col_pos", c("POS", "BP", "position", "pos_bp", "ps", "basepair")),
  trait = pick_col("col_trait", c("trait", "variable", "phenotype", "env", "envvar", "covariate")),
  set   = pick_col("col_set", c("set", "group", "population", "subset", "panel", "cluster")),
  q     = pick_col("col_q", c("q", "qvalue", "q_value", "FDR", "padj", "p_adj", "adj_p")),
  p     = pick_col("col_p", c("p", "pvalue", "p_value", "p.value", "p_wald", "p_lrt")),
  pve   = pick_col("col_pve", c("PVE", "R2", "rsq", "var_explained")))
msg("Association columns:", paste(names(cols), sapply(cols, function(x) if (is.null(x)) "-" else x), sep = "=", collapse = " "))
getc <- function(k, f, d) if (!is.null(cols[[k]])) f(as_raw[[cols[[k]]]]) else rep(d, nrow(as_raw))
A <- data.table(SNP = getc("snp", as.character, NA_character_), CHROM = getc("chr", as.character, NA_character_),
                POS = getc("pos", as.integer, NA_integer_), trait = getc("trait", as.character, "trait"),
                set = getc("set", as.character, FM$full_set), q = getc("q", as.numeric, NA_real_),
                p = getc("p", as.numeric, NA_real_), PVE = getc("pve", as.numeric, NA_real_))
if (anyNA(A$CHROM) || anyNA(A$POS)) {            # fill positions from VCF IDs
  if (is.null(cols$snp)) stop("Association table needs a marker column or chromosome + position columns")
  i <- match(A$SNP, vcf$ID)
  A[is.na(CHROM), CHROM := vcf$CHROM[i[is.na(A$CHROM)]]]; A[is.na(POS), POS := vcf$POS[i[is.na(A$POS)]]]
  if (anyNA(A$POS)) msg("Markers without a position (dropped):", uniqueN(A$SNP[is.na(A$POS)]))
  A <- A[!is.na(POS)]
}
A[is.na(SNP), SNP := paste(CHROM, POS, sep = ":")]
mnum <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
mk <- A[, .(n_traits = uniqueN(trait), within_group = any(set != FM$full_set),
            sets = paste(sort(unique(set)), collapse = ";"), traits = paste(sort(unique(trait)), collapse = ";"),
            max_PVE = mnum(PVE, max), min_q = mnum(q, min), min_p = mnum(p, min)), by = .(SNP, CHROM, POS)]
mk[, `:=`(annotation_chr = to_key(CHROM), priority = within_group * 2L + as.integer(n_traits > 1))]
mk[, `:=`(subgenome = subgenome(annotation_chr), in_annotation = annotation_chr %in% ext$key)]

## ---- 5. marker MAF and heterozygosity (only association markers are read) -------------
vk <- vcf[mk, on = .(key = annotation_chr, POS), nomatch = NULL, .(CHROM = x.CHROM, POS)]
if (nrow(vk)) {
  kf <- tempfile(fileext = ".tsv"); fwrite(unique(vk), kf, sep = "\t", col.names = FALSE)
  awk <- "awk -F'\\t' 'NR==FNR{k[$1\"\\t\"$2]=1;next} /^#/{next} (($1\"\\t\"$2) in k)'"
  gv <- fread(cmd = paste(cat_cmd(FM$vcf), "|", awk, shQuote(kf), "-"), sep = "\t", header = FALSE, colClasses = "character")
  g  <- sub(":.*", "", as.matrix(gv[, -(1:9)]))                 # GT field, any ploidy, phased or not
  miss <- grepl(".", g, fixed = TRUE) | g == ""
  n_ref <- nchar(gsub("[^0]", "", g)); n_alt <- nchar(gsub("[^1-9]", "", g))
  n_ref[miss] <- NA; n_alt[miss] <- NA
  af <- rowSums(n_alt, na.rm = TRUE) / rowSums(n_ref + n_alt, na.rm = TRUE)
  gs <- data.table(CHROM = gv[[1]], POS = as.integer(gv[[2]]), n_called = rowSums(!miss),
                   MAF = round(pmin(af, 1 - af), 4), het_obs = round(rowSums(n_ref > 0 & n_alt > 0, na.rm = TRUE) / rowSums(!miss), 4))
  gs[, het_excess := het_obs > 2 * MAF * (1 - MAF) * (1 - FM$inbreeding_F)]
  gs[, key := to_key(CHROM)]
  mk[gs, on = .(annotation_chr = key, POS), `:=`(MAF = i.MAF, het_obs = i.het_obs, het_excess = i.het_excess)]
  mk[, low_MAF := MAF < FM$low_maf]
  msg("Genotypes:", ncol(g), "samples;", nrow(gs), "association markers; low MAF:", sum(mk$low_MAF, na.rm = TRUE),
      "; excess heterozygosity:", sum(mk$het_excess, na.rm = TRUE))
} else mk[, `:=`(MAF = NA_real_, het_obs = NA_real_, het_excess = NA, low_MAF = NA)]

## ---- 6. search windows -------------------------------------------------------------------
ldw <- FM$ld_window
if (!is.null(FM$ld_summary)) {
  ld <- fread(FM$ld_summary)
  ldw <- c(ldw[setdiff(names(ldw), ld$set)], setNames(as.numeric(ld$half_decay_bp), ld$set))
  msg("LD windows from", basename(FM$ld_summary))
}
win_for <- function(chr, sg) {                   # chromosome value > subgenome value > genome value > FM$window
  w <- rep(if ("genome" %in% names(ldw)) ldw[["genome"]] else FM$window, length(chr))
  if (length(ldw)) {
    s <- !is.na(sg) & sg %in% names(ldw); w[s] <- ldw[sg[s]]
    ldc <- setNames(ldw, to_key(names(ldw))); c2 <- chr %in% names(ldc); w[c2] <- ldc[chr[c2]] }
  pmin(pmax(unname(w), FM$min_window), FM$max_window) }
mk[, window_bp := win_for(annotation_chr, subgenome)]
msg("Search windows (+/- bp):", tbl_str(mk$window_bp))

ord <- if (any(!is.na(mk$min_q))) mk$min_q else mk$min_p
mk <- mk[order(-priority, ord)]
setcolorder(mk, c("SNP", "CHROM", "POS", "annotation_chr", "subgenome", "in_annotation", "window_bp", "priority",
                  "within_group", "n_traits", "traits", "sets", "max_PVE", "min_q", "min_p", "MAF", "het_obs",
                  "low_MAF", "het_excess"))
fwrite(mk, tab("candidate_markers_prioritised.csv"))
msg(nrow(mk), "markers;", sum(!mk$in_annotation), "on chromosomes absent from the annotation;",
    "priority:", tbl_str(mk$priority))

## ---- 7. functional annotation tables --------------------------------------------------
ann_cols <- character(0)
norm_id <- function(x, level) {                  # level 1: drop ".v1.1"-style suffix; level 2: also ".1" isoform
  x <- strip_pref(sub("\\.v[0-9.]+$", "", x)); if (level == 2) x <- sub("\\.[0-9]+$", "", x); x }
id_fields <- c(tx_name = 0.5, tx_id = 0.5, gene_name = 0, gene_id = 0)   # small preference for transcript IDs
for (af in FM$annotation) {
  sep <- if (grepl("csv", basename(af), ignore.case = TRUE)) "," else "\t"
  an <- tryCatch(fread(cmd = cat_cmd(af), sep = sep, quote = "", colClasses = "character", fill = TRUE, header = FALSE),
                 error = function(e) NULL)
  if (is.null(an) || nrow(an) < 2) { msg("Could not read", basename(af)); next }
  first <- unlist(an[1]); all_ids <- unlist(genes[, .(gene_id, gene_name, tx_id, tx_name)])
  if (any(norm_id(first, 2) %in% norm_id(all_ids, 2))) {        # no header row
    stem <- sub("([._](txt|tsv|csv|gz|[0-9]))+$", "", basename(af))
    setnames(an, paste0(sub(".*[._]", "", stem), "_", seq_along(an)))
  } else { setnames(an, make.unique(sub("^#", "", ifelse(first == "", paste0("V", seq_along(first)), first)))); an <- an[-1] }
  setnames(an, gsub("[^A-Za-z0-9_.]", "_", names(an)))
  smp <- an[1:min(5000, nrow(an))]
  best <- list(score = 0); id_like <- character(0)
  min_hit <- max(1, min(10, 0.2 * nrow(genes)))   # matches needed to accept an ID column
  for (lev in 1:2) {
    for (cc in names(an)) {
      v <- unique(norm_id(smp[[cc]], lev))
      for (fld in names(id_fields)) {
        n <- sum(v %in% norm_id(genes[[fld]][!is.na(genes[[fld]])], lev))
        if (n >= min_hit) id_like <- union(id_like, cc)
        if (n + id_fields[[fld]] > best$score) best <- list(score = n + id_fields[[fld]], n = n, col = cc, fld = fld, lev = lev) } }
    if (best$score >= min_hit) break }
  if (best$score < min_hit) { msg("No gene/transcript ID column found in", basename(af), "- skipped"); next }
  an[, .k := norm_id(get(best$col), best$lev)]
  an <- an[!duplicated(.k)]
  keep <- setdiff(names(an), c(id_like, ".k"))
  keep <- keep[vapply(keep, function(k) uniqueN(an[[k]]) > 1, TRUE)]     # drop constant columns
  if (!length(keep)) next
  for (k in keep) an[, (k) := sub("^\\([0-9]+ of [0-9]+\\) ", "", get(k))]   # Phytozome "(1 of 2)" prefixes
  newn <- make.unique(c(names(genes), keep))[-(seq_along(names(genes)))]
  setnames(an, keep, newn)
  idx <- match(norm_id(genes[[best$fld]], best$lev), an$.k)
  for (k in newn) genes[, (k) := an[[k]][idx]]
  ann_cols <- c(ann_cols, newn)
  msg("Annotation table", basename(af), ": ID column", best$col, "->", best$fld, ";", sum(!is.na(idx)), "genes matched;",
      length(newn), "columns added")
}

## ---- 8. genes within the window of priority markers ------------------------------------
snp_context <- function(p, txid, strand, s, e) {
  if (p < s) return(if (strand == "-") "downstream" else "upstream")
  if (p > e) return(if (strand == "-") "upstream" else "downstream")
  f <- if (is.na(txid)) feats[0] else feats[.(txid), nomatch = NULL]
  if (!nrow(f)) return("genic")
  hit <- f[start <= p & end >= p]$ftype
  if ("CDS" %in% hit) return("CDS")
  if ("UTR5" %in% hit) return("5'UTR")
  if ("UTR3" %in% hit) return("3'UTR")
  cds <- f[ftype == "CDS"]; has_exon <- any(f$ftype == "exon")
  if (nrow(cds) && (p < min(cds$start) || p > max(cds$end))) {
    five <- (p < min(cds$start)) == (strand != "-")
    if ("exon" %in% hit) return(if (five) "5'UTR" else "3'UTR")
    if (!has_exon) return(if (five) "5'_untranslated" else "3'_untranslated")
  }
  if ("exon" %in% hit) return(if (nrow(cds)) "exon" else "noncoding_exon")
  "intron" }

top <- mk[priority >= FM$min_priority & in_annotation == TRUE]
if (!nrow(top)) { top <- head(mk[in_annotation == TRUE], FM$top_n); msg("No marker reaches min_priority; using top", nrow(top)) }
setkey(genes, key)
cg <- rbindlist(lapply(seq_len(nrow(top)), function(i) {
  p <- top$POS[i]; w <- top$window_bp[i]
  g <- genes[.(top$annotation_chr[i]), nomatch = NULL][end >= p - w & start <= p + w]
  if (!nrow(g)) return(NULL)
  g <- copy(g)
  g[, `:=`(SNP = top$SNP[i], CHROM_vcf = top$CHROM[i], POS = p, traits = top$traits[i], sets = top$sets[i],
           priority = top$priority[i], MAF = top$MAF[i], low_MAF = top$low_MAF[i], het_excess = top$het_excess[i],
           window_bp = w, dist_bp = pmax(0L, start - p, p - end))]
  g[, snp_position := mapply(snp_context, p, tx_id, strand, start, end)]
  g[, `:=`(nearest = dist_bp == min(dist_bp), in_report_window = dist_bp <= FM$report_window)]
  g[, c("SNP", "CHROM_vcf", "POS", "seqid", "subgenome", "traits", "sets", "priority", "MAF", "low_MAF", "het_excess",
        "gene_id", "gene_name", "tx_id", "start", "end", "strand", "window_bp", "dist_bp", "snp_position",
        "nearest", "in_report_window", "description", ann_cols), with = FALSE] }), fill = TRUE)
if (nrow(cg)) setnames(cg, "seqid", "annotation_chr")
msg("Markers searched:", nrow(top), "; candidate genes:", nrow(cg), "(", sum(cg$in_report_window), "within +/-",
    FM$report_window, "bp ); markers with no gene in window:", nrow(top) - uniqueN(cg$SNP))

## ---- 9. optional BLASTP of primary-transcript proteins -----------------------------------
if (nrow(cg) && length(FM$blast_dbs) && !is.null(FM$proteins_fa) && nzchar(Sys.which("blastp"))) {
  fa  <- readLines(if (grepl("gz$", FM$proteins_fa)) gzfile(FM$proteins_fa) else FM$proteins_fa)
  hdr <- grepl("^>", fa); rec <- cumsum(hdr); h <- fa[hdr]
  ids <- lapply(h, function(x) {
    tok <- strsplit(sub("^>", "", x), "\\s+")[[1]]
    kv  <- sub("^(transcript|locus|gene|parent|ID)[:=]", "", tok[grepl("^(transcript|locus|gene|parent|ID)[:=]", tok, ignore.case = TRUE)], ignore.case = TRUE)
    strip_pref(unique(c(tok[1], kv, sub("\\.[0-9]+$", "", tok[1])))) })
  want <- unique(na.omit(cg$tx_id)); want_g <- unique(cg$gene_id)
  qid <- vapply(ids, function(v) { m <- intersect(v, want); if (length(m)) m[1] else NA_character_ }, "")
  if (all(is.na(qid))) qid <- vapply(ids, function(v) { m <- intersect(v, want_g); if (length(m)) m[1] else NA_character_ }, "")
  keep <- !is.na(qid) & !duplicated(qid)
  sel <- fa[keep[rec]]; is_h <- grepl("^>", sel)
  sel[is_h] <- paste0(">", qid[keep]); sel[!is_h] <- gsub("\\*", "", sel[!is_h])
  qfa <- tab("candidate_proteins.fa"); writeLines(sel, qfa)
  msg("Proteins extracted:", sum(keep), "of", uniqueN(cg$gene_id), "candidate genes")
  for (db in names(FM$blast_dbs)) {
    if (!db_exists(FM$blast_dbs[[db]])) { msg("BLAST db not found, skipped:", FM$blast_dbs[[db]]); next }
    out <- tab(paste0("blast_", db, ".tsv"))
    system(sprintf("blastp -query %s -db %s -evalue 1e-10 -max_target_seqs 3 -num_threads %d -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' > %s",
                   shQuote(qfa), shQuote(FM$blast_dbs[[db]]), FM$threads, shQuote(out)))
    if (file.exists(out) && file.size(out) > 0) {
      b <- fread(out, sep = "\t", header = FALSE, quote = "",
                 col.names = c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"))
      b <- b[order(qseqid, -bitscore)][!duplicated(qseqid)]
      nm <- paste0("blast_", db, c("_hit", "_pident", "_evalue", "_desc"))
      setnames(b, c("sseqid", "pident", "evalue", "stitle"), nm)
      key_col <- if (any(b$qseqid %in% cg$tx_id)) "tx_id" else "gene_id"
      cg <- merge(cg, b[, c("qseqid", nm), with = FALSE], by.x = key_col, by.y = "qseqid", all.x = TRUE, sort = FALSE)
      msg("BLAST", db, ":", nrow(b), "queries with hits")
    }
  }
}

if (nrow(cg)) setorder(cg, -priority, SNP, dist_bp)
fwrite(cg, tab("candidate_genes.csv"))
if (nrow(cg)) msg("Candidate genes by SNP position:", tbl_str(cg$snp_position),
                  if (!is.null(FM$subgenome_regex)) paste("| by subgenome:", tbl_str(cg$subgenome)) else "")
msg("Done ->", normalizePath(out_dir))

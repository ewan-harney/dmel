doMatch_4 <- function (rnaMeta, chipMeta, region, method, annotation, 
                     fragLength = 180, promoter.length = 5000, conditions=2) 
{
    checkInputParam <- function(rnaMeta, chipMeta, region, method, 
                                annotation) {
        if (missing(rnaMeta) || is.null(rnaMeta) || missing(chipMeta) || 
            is.null(chipMeta) || missing(region) || is.null(region) || 
            missing(method) || is.null(method) || missing(annotation) || is.null(annotation)) {
            stop("arguments rnaMeta, chipMeta, region,\n           method, annotation must be specified")
        }
        if (class(rnaMeta) != "data.frame") {
            stop("rnaMeta must be a data.frame")
        }
        if (!c("condition") %in% colnames(rnaMeta)) {
            stop("rnaMeta has to contain \"condition\" column that indicates\n           the condition of experiment.")
        }
        if (!c("files") %in% colnames(rnaMeta)) {
            stop("rnaMeta has to contain \"files\" column indicates the paths of\n           cprresponing abundance.tsv file calculated from Kallisto")
        }
        if (nrow(rnaMeta) <= 1) {
            stop("rnaMeta has to contain at least 2 rows\n           (one for each biological condition)")
        }
        if (class(chipMeta) != "data.frame") {
            stop("chipMeta must be a data.frame")
        }
        if (!c("mark") %in% colnames(chipMeta)) {
            stop("chipMeta has to contain \"mark\" column indicates the markers of\n           histone modifications")
        }
        if (!c("condition") %in% colnames(chipMeta)) {
            stop("chipMeta has to contain \"condition\" column that indicates the\n           condition of experiment")
        }
        if (!c("files") %in% colnames(chipMeta)) {
            stop("chipMeta has to contain \"files\" column indicates the paths of\n           aligned bam files")
        }
        if (!is.character(chipMeta$files)) {
            stop("the class of the paths of the aligned bam files should be\n           character")
        }
        if (nrow(chipMeta) <= 1) {
            stop("chipMeta has to contain at least 2 rows\n           (one for each biological condition)")
        }
        if (length(conditions) == 3 ) {
            if (length(unique(chipMeta[chipMeta$condition%in%conditions,]$condition)) != 3 | length(unique(rnaMeta[rnaMeta$condition%in%conditions,]$condition)) != 
                3) {
                stop("three distinct biological conditions\n           expected in chipMeta and rnaMeta")
            }
            if (length(unique(c(chipMeta[chipMeta$condition%in%conditions,]$condition, rnaMeta[rnaMeta$condition%in%conditions,]$condition))) != 
                3) {
                stop("identical biological conditions\n           expected in chipMeta and rnaMeta")
            }
        } else if (length(conditions) == 2 ) {
            if (length(unique(chipMeta[chipMeta$condition%in%conditions,]$condition)) != 2 | length(unique(rnaMeta[rnaMeta$condition%in%conditions,]$condition)) != 
                2) {
                stop("two distinct biological conditions\n           expected in chipMeta and rnaMeta")
            }
            if (length(unique(c(chipMeta[chipMeta$condition%in%conditions,]$condition, rnaMeta[rnaMeta$condition%in%conditions,]$condition))) != 
                2) {
                stop("identical biological conditions\n           expected in chipMeta and rnaMeta")
            }  
        }
        
        if (length(region) != 1) {
            stop("region has to be specified as \"promoter\" or \"genebody\"")
        }
        if (!is.character(region)) {
            stop("region has to be specified as \"promoter\" or \"genebody\"")
        }
        if (!region %in% c("promoter", "genebody")) {
            stop("region has to be specified as \"promoter\" or \"genebody\"")
        }
        if (region == "promoter") {
            if (length(method) != 1) {
                stop("method has to be specified as \"weighted.mean\" or\n             \"highest\" if region is set as \"promoter\"")
            }
            if (!is.character(method)) {
                stop("method has to be specified as \"weighted.mean\" or\n             \"highest\" if region is set as \"promoter\"")
            }
            if (!method %in% c("weighted.mean", "highest")) {
                stop("method has to be specified as \"weighted.mean\" or\n             \"highest\" if region is set as \"promoter\"")
            }
        }
    }
    checkInputParam(rnaMeta = rnaMeta, chipMeta = chipMeta, region = region, 
                    method = method, annotation = annotation)
    message("get RNA counts")
    est_counts <- function(x) {
        ex <- utils::read.table(x, sep = "\t", header = TRUE, 
                                row.names = "target_id")
        ex[, "est_counts", drop = FALSE]
    }
    files <- NULL
    for (i in conditions) {
        a <- rnaMeta[rnaMeta$condition == i, ]$files
        names(a) <- paste(i, "_REP", 1:length(a), sep = "")
        files <- c(a, files)
    }
    fileNames <- names(files)
    exprsList <- list()
    for (f in fileNames) {
        x <- est_counts(files[f])
        names(x) <- f
        exprsList[[f]] <- x
    }
    exprs <- do.call("cbind", exprsList)
    #ensembl <- biomaRt::useMart("ensembl", dataset = annotation, 
    #    host = host)
    N <- as.data.frame(annotation)
    #N <- data.frame(biomaRt::getBM(attributes = c("ensembl_transcript_id", 
    #    "ensembl_gene_id", "ensembl_gene_id", "chromosome_name", 
    #    "start_position", "end_position", "transcript_start", 
    #    "transcript_end", "transcription_start_site", "transcript_length", 
    #    "strand"), filters = c("ensembl_transcript_id"), 
    #    values = rownames(exprs), mart = ensembl))
    exprsRNA <- merge(N, exprs, by.x = "ensembl_transcript_id", 
                      by.y = 0)
    exprsRNA <- exprsRNA[exprsRNA$chromosome_name %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chrX"), ]
    chipMeta <- chipMeta[chipMeta$condition%in%conditions,]
    aln <- GenomicAlignments::readGAlignments(file = chipMeta$files[1])
    message("get chip counts")
    if (region == "promoter") {
        if ("chr2L" %in% as.character(GenomicRanges::seqnames(aln))) {
            tss <- GenomicRanges::GRanges(seqnames = paste0(exprsRNA$chromosome_name), 
                                          ranges = IRanges::IRanges(start = exprsRNA$transcription_start_site, 
                                                                    width = 1))
        }
        promoter <- suppressMessages(GenomicRanges::resize(tss,
                                                           fix = "center", width = promoter.length))
        names(promoter) <- exprsRNA$ensembl_transcript_id
        files <- NULL
        for (i in unique(paste(chipMeta$mark, "_HM_", chipMeta$condition, 
                               sep = ""))) {
            a <- chipMeta[paste(chipMeta$mark, "_HM_", chipMeta$condition, 
                                sep = "") == i, ]$files
            names(a) <- paste(i, "_REP", 1:length(a), sep = "")
            files <- c(a, files)
        }
        fileNames <- names(files)
        exprsChIP.list <- NULL
        for (f in fileNames) {
            a <- bam2counts(bamFile = files[f], region = promoter, 
                            fragLength = fragLength)
            exprsChIP.list[[f]] <- a
        }
        exprsChIP <- as.data.frame(exprsChIP.list)
        match.rc <- merge(exprsRNA, exprsChIP, by.x = "ensembl_transcript_id", 
                          by.y = 0)
        chip <- match.rc[, colnames(match.rc) %in% c("ensembl_gene_id", 
                                                     colnames(exprsChIP))]
        rna <- match.rc[, colnames(match.rc) %in% c("ensembl_gene_id", 
                                                    colnames(exprs))]
        res.rna <- stats::aggregate(. ~ ensembl_gene_id, data = rna, 
                                    FUN = sum)
        if (method == "weighted.mean") {
            res.chip <- stats::aggregate(. ~ ensembl_gene_id, 
                                         data = chip, FUN = function(x) stats::weighted.mean(x + 
                                                                                                 0.01, x + 0.01))
            matched.wtmean <- merge(res.rna, res.chip, by = "ensembl_gene_id")
            res.rna <- stats::aggregate(. ~ ensembl_gene_id, 
                                        data = res.rna, FUN = as.integer)
            res.chip <- stats::aggregate(. ~ ensembl_gene_id, 
                                         data = res.chip, FUN = as.integer)
            matched.wtmean <- stats::aggregate(. ~ ensembl_gene_id, 
                                               data = matched.wtmean, FUN = as.integer)
            return(list(res.rna = res.rna, res.chip = res.chip, 
                        matched.data = matched.wtmean))
        }
        else {
            res.chip <- stats::aggregate(. ~ ensembl_gene_id, 
                                         data = chip, FUN = max)
            matched.highest <- merge(res.rna, res.chip, by = "ensembl_gene_id")
            res.rna <- stats::aggregate(. ~ ensembl_gene_id, 
                                        data = res.rna, FUN = as.integer)
            res.chip <- stats::aggregate(. ~ ensembl_gene_id, 
                                         data = res.chip, FUN = as.integer)
            matched.highest <- stats::aggregate(. ~ ensembl_gene_id, 
                                                data = matched.highest, FUN = as.integer)
            return(list(res.rna = res.rna, res.chip = res.chip, 
                        matched.data = matched.highest))
        }
    }
    else {
        rna <- exprsRNA[, colnames(exprsRNA) %in% c("ensembl_gene_id", 
                                                    colnames(exprs))]
        res.rna <- stats::aggregate(. ~ ensembl_gene_id, data = rna, 
                                    FUN = sum)
        genebody.region <- exprsRNA[colnames(exprsRNA) %in% c("ensembl_gene_id", 
                                                              "chromosome_name", "start_position", "end_position")]
        genebody.region <- genebody.region[!duplicated(genebody.region$ensembl_gene_id), 
        ]
        if ("chr2L" %in% as.character(GenomicRanges::seqnames(aln))) {
            genebody <- GenomicRanges::GRanges(seqnames = paste0(genebody.region$chromosome_name),ranges = IRanges::IRanges(start = genebody.region$start_position, 
                                                                                                                             end = genebody.region$end_position))
        }
        else {
            genebody <- GenomicRanges::GRanges(seqnames = genebody.region$chromosome_name, 
                                               ranges = IRanges::IRanges(start = genebody.region$start_position, 
                                                                         end = genebody.region$end_position))
        }
        names(genebody) <- genebody.region$ensembl_gene_id
        files <- NULL
        for (i in unique(paste(chipMeta$mark, "_HM_", chipMeta$condition, 
                               sep = ""))) {
            a <- chipMeta[paste(chipMeta$mark, "_HM_", chipMeta$condition, 
                                sep = "") == i, ]$files
            names(a) <- paste(i, "_REP", 1:length(a), sep = "")
            files <- c(a, files)
        }
        fileNames <- names(files)
        exprsChIP.list <- NULL
        for (f in fileNames) {
            a <- bam2counts(bamFile = files[f], region = genebody, 
                            fragLength = fragLength)
            exprsChIP.list[[f]] <- a
        }
        res.chip <- as.data.frame(exprsChIP.list)
        res.chip$ensembl_gene_id <- rownames(res.chip)
        matched.genebody <- merge(res.rna, res.chip, by = "ensembl_gene_id")
        res.rna <- stats::aggregate(. ~ ensembl_gene_id, data = res.rna, 
                                    FUN = as.integer)
        res.chip <- stats::aggregate(. ~ ensembl_gene_id, 
                                     data = res.chip, FUN = as.integer)
        matched.genebody <- stats::aggregate(. ~ ensembl_gene_id, 
                                             data = matched.genebody, FUN = as.integer)
        return(list(res.rna = res.rna, res.chip = res.chip, matched.data = matched.genebody))
    }
}

#!/bin/env Rscript
# ---------------------------------------------------------------------------- #
# ncbi.r: Download genomic files from the NCBI genome archive
# version: 1.1.0
# https://github.com/J-Moravec/ncbi.r
# ---------------------------------------------------------------------------- #
basename_sans_ext = function(x){
    tools::file_path_sans_ext(basename(x))
    }


replace_ext = function(x, ext){
    paste0(tools::file_path_sans_ext(x), ".", ext)
    }


download = function(x, dir = ".", types = c("GENOME_FASTA", "GENOME_GTF"), overwrite = FALSE,
                    quiet = TRUE, timeout = 3600){

    dest = file.path(dir, paste0(x, ".zip"))

    if(file.exists(dest) && !overwrite)
        return(invisible(dest))

    formats = c("GENOME_FASTA", "GENOME_GFF", "GENOME_GTF", "RNA_FASTA",
                "CDS_FASTA", "PROT_FASTA", "SEQUENCE_REPORT")
    types = match.arg(types, formats, several.ok = TRUE)

    base = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
    types = paste0(types, collapse=",")

    url = paste0(
        base,
        x, "/download?include_annotation_type=",
        types
        )

    # Remove file if something stopped the download
    ok = FALSE
    on.exit(
        if(!ok){
            file.remove(dest)
            stop("Download of '", x, "' from NCBI failed.", call. = FALSE)
            },
        add = TRUE
        )

    op = options("timeout" = timeout)
    err = download.file(url, dest, quiet = quiet)
    ok = TRUE
    options(op)

    if(err != 0){
        file.remove(dest)
        stop(paste0("Download of '", x, "' from NCBI failed."))
        }

    # If invalid ID is provided, the download progresses normally but the file is malformed.
    # The usual size of such file seems to be 855 bytes.
    if(file.size(dest) < 1000){
        file.remove(dest)
        stop(paste0("The requested ID '", x, "' doesn't exist."))
        }

    invisible(dest)
    }


extract = function(x, dir = ".", keep = FALSE){
    ext = tools::file_ext(x)
    prefix = basename_sans_ext(x)

    if(ext != "zip")
        stop(paste0("Unrecognised extension '", ext, "'. Only zip archives are supported.")) 

    workdir = file.path(dir, prefix)
    unzip(x, exdir = workdir, unzip = "unzip")

    if(!keep)
        file.remove(x)

    invisible(workdir)
    }


get_assembly_name = function(x){
    y = readLines(x)
    regmatches(y, regexpr("assemblyName\":\"([^\"]*)", y)) |> substring(16)
    }


format_mapping_table = function(id, assembly_name, formats = NULL){
    genome_fasta = paste(id, assembly_name, "genomic.fna", sep = "_")

    id_ext = function(ext){
        paste0(id, ".", ext)
        }

    mapping = list(
        c("GENOME_FASTA", genome_fasta, id_ext("fna")),
        c("GENOME_GFF", "genomic.gff", id_ext("gff")),
        c("GENOME_GTF", "genomic.gtf", id_ext("gtf")),
        c("RNA_FASTA", "rna.fna", id_ext("rna.fna")),
        c("CDS_FASTA", "cds_from_genomic.fna", id_ext("cds.fna")),
        c("PROT_FASTA", "protein.faa", id_ext("faa")),
        c("SEQUENCE_REPORT", "sequence_report.jsonl", id_ext("jsonl"))
        ) |> do.call(what = rbind.data.frame) |> setNames(c("formats", "old", "new"))
    mapping[["old"]] = file.path("ncbi_dataset", "data", id, mapping[["old"]])

    if(!is.null(formats))
        mapping = mapping[match(formats, mapping[["formats"]], nomatch = 0),]

    mapping
    }


"%nin%" = Negate("%in%")

flatten = function(
    x, dir = ".",
    types = c("GENOME_FASTA", "GENOME_GFF", "GENOME_GTF", "RNA_FASTA",
                "CDS_FASTA", "PROT_FASTA", "SEQUENCE_REPORT"),
    keep = FALSE, overwrite = FALSE
    ){
    formats = match.arg(types, several.ok = TRUE)
    prefix = basename_sans_ext(x)

    # early exit
    new_files = file.path(dir, format_mapping_table(prefix, "phony", formats)$new)
    if(all(file.exists(new_files)) && !overwrite)
        return(invisible(new_files))


    y = extract(x, dir = dir, keep = keep)
    assembly_name = get_assembly_name(
        file.path(y, "ncbi_dataset", "data", "assembly_data_report.jsonl")
        )
    mapping = format_mapping_table(prefix, assembly_name, formats)

    # safety check
    files = list.files(y, recursive = TRUE)
    if(!all(mapping$old %in% files)){
        missing = mapping$formats[mapping$old %nin% files]
        stop("Some of the requested formats are not in the NCBI archive:\n", toString(missing))
        }

    # move and remove
    old_files = file.path(y, mapping$old)
    file.rename(old_files, new_files) # no copying involved, should be faster

    # always remove the folder
    if(TRUE)
        unlink(y, recursive = TRUE, force = TRUE)

    invisible(new_files)
    }


gzip = function(x, out = NULL, overwrite = FALSE, keep = FALSE){
    if(is.null(out))
        out = paste0(x, ".gz")

    if(x == out)
        stop("Input and output paths are identical.")

    if(file.exists(out) && !overwrite)
        return(invisible())

    if(file.exists(out))
        file.remove(out)

    completed = FALSE
    input_connection = file(x, open = "rb")
    on.exit(expr = {
        close(input_connection)
        if(completed && !keep) file.remove(x)
        }, add = TRUE)

    output_connection = gzfile(out, open = "wb")
    on.exit(expr = {
        close(output_connection)
        if(!completed) file.remove(out)
        }, add = TRUE)

    repeat {
        b = readBin(input_connection, what = raw(0L), size = 1L, n = 1e7)
        if (length(b) == 0L) break
        writeBin(b, con = output_connection, size = 1L)
        }

    completed = TRUE
    invisible()
    }


# ---------------------------------------------------------------------------- #
# rargs: Copy-pastable argument parser
# version: 1.0.0
# https://github.com/J-Moravec/rargs
# ---------------------------------------------------------------------------- #


#' Get a name of a script
#'
#' Get the name of the script's filename when run through Rscript
#'
#' For instance, for a script `script.r` in the `folder` folder,
#' it could be run as `Rscript folder/script.r`. In that case,
#' the `get_scriptname` returns the `script.r`.
get_scriptname = function(){
    args = commandArgs(FALSE)

    file_arg = grep("--file=", args, fixed=TRUE, value=TRUE)[1]

    # not run throught script
    if(length(file_arg) == 0)
        return(NULL)

    sub("^--file=", "", file_arg)
    }


#' Define an option
#'
#' This is a helper function to define an option
opt = function(name, short = NULL, default = NULL, flag = FALSE){
    # if name starts with --, it is a long form
    # otherwise it is a positional argument
    long = NULL
    if(startsWith(name, "--")){
        long = name
        name = substring(name, 3)
        }

    if(flag && !is.null(long))
        default = FALSE

    list("name" = name, "short" = short, "long" = long, "default" = default, "flag" = flag)
    }


#' Parse arguments
#'
#' A POSIX-compatible argument parser.
#'
#' @param args **optional** a character vector of arguments, if not provided the commnad-line
#' arguments are obtained using the `commandArgs()` function.
#' @param options a list of one or more options obtained from the `opt` helper function.
#' @return a list of parsed arguments
parse_args = function(args = NULL, options){

    split_short_form = function(x){
        f = \(y){ if(grepl("^-[^-]", y)) paste0("-", strsplit(y, "")[[1]][-1]) else y}
        lapply(x, f) |> unlist()
        }

    short_to_long = function(args, options){
        id = match(args, options$short, nomatch = 0)
        args[id != 0] = unlist(options$long[id])
        args
        }

    check_unknown = function(args, options){
        args = grep("^-", args, value = TRUE)
        opts = unlist(c(options$long, options$short))
        unknown = args[!args %in% opts]

        if(length(unknown) != 0)
            stop("Unknown arguments: ", paste0(unknown, collapse = ","))
        }

    if(is.null(args))
        args = commandArgs(TRUE)

    # add --help to options
    options = c(options, opt("--help", "-h", flag = TRUE) |> list())
    names(options) = sapply(options, getElement, "name")
    options = as.data.frame(do.call(rbind, options))
    options[] = lapply(options, setNames, rownames(options)) # fix missing names

    # remove everything after --
    dashdash = c()
    if("--" %in% args){
        dashdash_id = which(args == "--")
        dashdash = args[-seq_len(dashdash_id)]
        args = args[seq_len(dashdash_id - 1)]
        }

    args = split_short_form(args)
    args = short_to_long(args, options)
    check_unknown(args, options)

    positional = c()
    pars = options$default

    # if arguments contain help, stop parsing
    if("--help" %in% args){
        pars$help = TRUE
        return(pars)
        }

    while(length(args) > 0){
        id = match(args[1], options$long)

        if(is.na(id)){
            positional = c(positional, args[1])
            args = args[-1]
            next
            }

        if(options$flag[[id]]){
            pars[id] = TRUE
            args = args[-1]
            next
            }

        if(length(args) < 2 || args[2] %in% options$long)
            stop("Not enough arguments for ", args[1], call. = FALSE)

        pars[id] = args[2]
        args = args[-c(1:2)]
        }

    # assign positionals to named args and the rest into pars$positional
    positional = c(positional, dashdash)
    pos = lengths(options$long) == 0
    n_pos = sum(pos)
    pars[pos] = as.list(positional)[seq_len(n_pos)]

    if(n_pos < length(positional))
        pars$positional = positional[-seq_len(n_pos)]

    pars
    }


usage = function(){
    prog = get_scriptname()
    blnk = strrep(" ", nchar(prog))
    cat(paste0(
"Usage: ", prog, " [options] ID\n",
"Download and flatten NCBI genome archive. By default, the archive is\n",
"downloaded, the requested files are extracted and gzipped. If any of\n",
"the options download, flatten, and gzip are specified, this behaviour\n",
"is supressed.\n\n",
"  -d,  --download   download genomic archive from NCBI\n",
"  -f,  --flatten    flatten the archive\n",
"  -g,  --gzip       gzips output files to save space\n",
"  -o,  --overwrite  overwrite existing files\n",
"  -k,  --keep       keep the genomic archive\n",
"       --dir        directory where the files will be downloaded\n",
"  -t,  --type=TYPES comma separated list of types, see details\n",
"  -h,  --help       display this help and exit\n",
"\n",
"ID is the NCBI id of genome that will be downloaded from the NCBI genome\n",
"website, see: https://www.ncbi.nlm.nih.gov/datasets/genome/\n",
"TYPES are comma separated list of types, allowed types are:\n",
"GENOME_FASTA, GENOME_GFF, GENOME_GTF, RNA_FASTA, CDS_FASTA,\n",
"PROT_FASTA, and SEQUENCE_REPORT. Spaces are not allowed.\n",
"By default, GENOME_FASTA and GENOME_GTF are downloaded.\n\n",
"Examples:\n",
"  Rscript ", prog, " GCF_000091225.2\n",
"          ", blnk, " downloads the reference genome fasta and gtf file for\n",
"          ", blnk, " Eencephalitozoon cuniculi GB-M1\n",
"  Rscript ", prog, " -dfgt GENOME_FASTA,GENOME_GTF GCF_000091225.2\n",
"  Rscript ", blnk, " equivalent to the above\n",
"  Rscript ", prog, " -fgk GCF_000091225.2.zip\n",
"  Rscript ", blnk, " extract and preserve already existing archive\n",
"  Rscript ", prog, " -dkt PROT_FASTA --dir PATH GCF_000091225.2\n",
"  Rscript ", blnk, " download only the protein sequences into path\n\n"
        ))
    }


opts = list(
    opt("--download", "-d", flag = TRUE),
    opt("--flatten", "-f", flag = TRUE),
    opt("--gzip", "-g", flag = TRUE),
    opt("--overwrite", "-o", flag = TRUE),
    opt("--keep", "-k", flag = TRUE),
    opt("--dir", default = "."),
    opt("--type", "-t", default = c("GENOME_FASTA,GENOME_GTF")),
    opt("id")
    )


main = function(){
    args = parse_args(options = opts)

    if(args$help){
        usage()
        return(invisible())
        }

    if(is.null(args$id)){
        usage()
        stop("Not enough argument.", call. = FALSE)
        }


    types = c("GENOME_FASTA", "GENOME_GFF", "GENOME_GTF", "RNA_FASTA",
             "CDS_FASTA", "PROT_FASTA", "SEQUENCE_REPORT")
    if(identical(args$type, "all")){
        args$type = types
        } else {
        args$type = strsplit(args$type, ",", fixed = TRUE)[[1]]
        }

    if(any(args$type %nin% types))
        stop("Unknown types:", toString(args$type[args$type %nin% type]))

    # default behaviour
    default = c("download", "flatten", "gzip")
    if(all(unlist(args[default]) == FALSE)) args[default] = TRUE

    # targets
    targets = file.path(args$dir, format_mapping_table(args$id, "phony", args$type)$new)

    if(args$gzip)
        targets = paste0(targets, ".gz")

    redo = args$overwrite || !all(file.exists(targets))

    # Nothing needs to be done
    if(!redo){
        return(invisible(targets))
        }

    var = args$id
    if(args$download)
        var = download(var, dir = args$dir, types = args$type, overwrite = args$overwrite)

    if(args$flatten)
        var = flatten(var, dir = args$dir, types = args$type, keep = args$keep,
                overwrite = args$overwrite)

    if(args$gzip)
        sapply(var, gzip, overwrite = args$overwrite)

    invisible()
    }

if(sys.nframe() == 0){
    main()
    }

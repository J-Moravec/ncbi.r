#!/bin/env Rscript
#' ncbi.r
#'
#' Download, extract, and flatten NCBI genome archive.
#'
#' See:
#'
#' Rscript ncbi.r --help
#'
#' for usage information


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


extract = function(x, local = FALSE){
    ext = tools::file_ext(x)
    prefix = basename_sans_ext(x)

    if(ext != "zip")
        stop(paste0("Unrecognised extension '", ext, "'. Only zip archives are supported.")) 

    workdir = file.path(if(local) dirname(x) else tempdir(), prefix)
    unzip(x, exdir = workdir, unzip = "unzip")
    }


flatten = function(x, local = FALSE, clean = TRUE){
    prefix = basename_sans_ext(x)

    workdir = file.path(if(local) dirname(x) else tempdir(), prefix)
    files = file.path(workdir, unzip(x, list = TRUE)$Name)

    lapply(files, move, dir = dirname(x), prefix = prefix)

    if(clean)
        unlink(workdir, recursive = TRUE, force = TRUE)

    invisible()
    }


"%nin%" = Negate("%in%")


move = function(x, dir, prefix){
    ext = tools::file_ext(x)

    if(ext %nin% c("fna", "gtf"))
        return()

    new_path = file.path(dir, paste0(prefix, ".", ext))
    file.copy(x, new_path)
    }


gzip = function(x, out = NULL, overwrite=FALSE, keep = FALSE){
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
        "Usage: ", prog, " [options] file.zip\n",
        "Download, extract, and flatten NCBI genome archive.\n\n",
        "  -d,  --download  download genomic archive from NCBI\n",
        "  -o,  --overwrite overwrite archive if it exists\n",
        "  -e,  --extract   extract genomic archive\n",
        "  -l,  --locally   extract the archive locally\n",
        "  -f,  --flatten   flatten the extracted archive\n",
        "  -g,  --gzip      gzips output files to save space\n",
        "  -h,  --help      display this help and exit\n",
        "\n",
        "The option can be specified in the POSIX-like format. Both long and\n",
        "short options are allowed. By default, ", prog, " downloads, extracts,\n",
        "and flattens the specified archive. When any argument is specified,\n",
        "this behaviour is suppressed. This is useful when the archive is already\n",
        "downloaded.\n\n",
        "The file must be in the format [NCBI ID].zip, where [NCBI ID] is the ID\n",
        "of the genome that will be downloaded from NCBI genome website, see:\n",
        "https://www.ncbi.nlm.nih.gov/datasets/genome/ for more information.\n\n",
        "If the '-el' options are specified, the archive is extracted into the\n",
        "local directory. To flatten this extracted archive, '-fl' needs\n",
        "to be specified.\n\n",
        "Examples:\n",
        "  Rscript ", prog, " file.zip       downloads, extracts, and flattens\n",
        "  Rscript ", prog, " -ref file.zip  equivalent to above\n",
        "  Rscript ", prog, " -ef file.zip   only extract and flattens\n",
        "  Rscript ", prog, " -el file.zip   archive is extracted to a local path\n",
        "  Rscript ", prog, " -e  file.zip   archive is extracted to a temp path,\n",
        "          ", blnk, "                this path is deleted when R session\n",
        "          ", blnk, "                ends, so no output is produced.\n",
        "  Rscript ", prog, " -g file.zip    gzip the fna and gtf files in the zip\n",
        "          ", blnk, "                producing file.fna.gz and file.gtf.gz\n\n"
        ))
    }


opts = list(
    opt("--download", "-d", flag = TRUE),
    opt("--extract", "-e", flag = TRUE),
    opt("--flatten", "-f", flag = TRUE),
    opt("--local", "-l", flag = TRUE),
    opt("--overwrite", "-o", flag = TRUE),
    opt("--gzip", "-g", flag = TRUE),
    opt("file")
    )


main = function(){
    args = parse_args(options = opts)

    if(args$help){
        usage()
        return(invisible())
        }

    if(is.null(args$file))
        stop("Not enough argument.", call. = FALSE)

    if(tools::file_ext(args$file) != "zip")
        stop(paste0("Unrecognised extension '", ext, "'. Only zip archives are supported.")) 

    # implement requested default behaviour
    default = c("download", "extract", "flatten", "gzip")
    if(all(unlist(args[default]) == FALSE)) args[default] = TRUE

    targets = replace_ext(args$file, c("fna", "gtf"))
    if(args$gzip)
        targets = paste0(targets, ".gz")
    redo = args$overwrite || !all(file.exists(targets))

    if(args$download)
        download(args$file, overwrite = args$overwrite)

    if(args$extract && args$flatten && redo){
        extract(args$file, args$local)
        flatten(args$file, args$local)
        }

    if(args$extract && !args$flatten)
        extract(args$file, args$local)

    if(!args$extract && args$flatten)
        flatten(args$file, args$local)

    if(args$gzip){
        gzip(replace_ext(args$file, "fna"), overwrite = args$overwrite)
        gzip(replace_ext(args$file, "gtf"), overwrite = args$overwrite)
        }
    }


if(sys.nframe() == 0){
    main()
    }

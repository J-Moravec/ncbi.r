# Encephalitozoon cuniculi GB-M1
# The smallest genome with all available file types
# Smaller genomes (viruses) do not have transcripts.
workdir = "test-workdir"

id = "GCF_000091225.2"
zip = "test-data/GCF_000091225.2.zip"
zip_copy = file.path(workdir, basename(zip))

if(!dir.exists(workdir))
    dir.create(workdir)

TEST_SET("Can downloads file from NCBI", {
    # download all types
    types = c("GENOME_FASTA", "GENOME_GFF", "GENOME_GTF", "RNA_FASTA", "CDS_FASTA",
              "PROT_FASTA", "SEQUENCE_REPORT")

    file = download(id, dir = workdir, types = types)
    TEST(file.exists(file))
    TEST(file.size(file) == 3359895)
    file.copy(file, zip)
    })


TEST_SET("Can extract the NCBI zip file", {
    dir = extract(zip, workdir, keep = TRUE)

    TEST(file.exists(zip)) # original file was not removed
    TEST(dir.exists(dir))
    TEST(dir.size(dir) == 16049816)
    dir.remove(dir)

    file.copy(zip, zip_copy)
    TEST(file.exists(zip))
    TEST(file.exists(zip_copy))

    dir = extract(zip_copy, workdir, keep = FALSE)
    TEST(file.exists(zip))
    TEST(!file.exists(zip_copy))
    TEST(dir.exists(dir))
    dir.remove(dir)
    })


TEST_SET("Can flatten the NCBI zip file", {
    file.copy(zip, zip_copy)
    sizes = readRDS(file.path("test-data", paste0(id, ".sizes.rds")))

    files = flatten(zip_copy, dir = workdir, keep = TRUE)
    TEST(all(file.size(files) == sizes[basename(files)]))

    # even if files exist, early exit ignores the keep argument
    files = flatten(zip_copy, dir = workdir, keep = FALSE)
    TEST(all(file.size(files) == sizes[basename(files)]))
    TEST(file.exists(zip_copy))
    file.remove(files)

    types = c("GENOME_FASTA", "GENOME_GTF", "CDS_FASTA")
    matching_files = paste0(id, c(".fna", ".gtf", ".cds.fna"))

    files = flatten(zip_copy, dir = workdir, types = types, keep = TRUE)
    TEST(all(basename(files) == matching_files))
    file.remove(files)

    # Order corresponds to requested types
    files = flatten(zip_copy, dir = workdir, types = rev(types), keep = TRUE)
    TEST(all(basename(files) == rev(matching_files)))
    file.remove(files)
    })

file.copy(zip, zip_copy) # restore the file to prevent redownload


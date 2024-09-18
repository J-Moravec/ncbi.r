# To limit redownloading, the genome and other intensive operations
# files are not deleted at the end of the test

testdir = "test-data"
if(!dir.exists(testdir))
    dir.create(testdir)

TEST_SET("Can downloads file from NCBI", {
    # download all types
    types = c("GENOME_FASTA", "GENOME_GFF", "GENOME_GTF", "RNA_FASTA", "CDS_FASTA",
              "PROT_FASTA", "SEQUENCE_REPORT")

    file = download("GCF_034140825.1", dir = testdir, types = types)
    TEST(file.exists(file))
    TEST(file.size(file) == 195458653)
    })

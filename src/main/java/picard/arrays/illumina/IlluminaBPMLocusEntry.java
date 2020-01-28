package picard.arrays.illumina;

/**
 * A simple class to represent a locus entry in an Illumina BPM file.
 */
public class IlluminaBPMLocusEntry {
    String ilmnId;      // IlmnID (probe identifier) of locus
    String name;        // Name (variant identifier) of locus
    int index;          // Index of this entry.
    String ilmnStrand;  // Illumina Strand value
    String snp;         // SNP value for locus (e.g., [A/C])
    String chrom;       // Chromosome for the locus (e.g., XY)
    String ploidy;      // Ploidy
    String species;
    int mapInfo;        // Mapping location of locus
    String customerStrand;   // ? Customer Strand
    int addressA;       // AddressA ID of locus
    int addressB;       // AddressB ID of locus (0 if none)

    String genomeBuild;
    String source;
    String sourceVersion;
    String sourceStrand;

    int expClusters;
    int intensityOnly;
    int assayType;      // Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)

    float fracA;
    float fracC;
    float fracT;
    float fracG;

    String refStrand;

    // Not part of the locusEntry record in the BPM, added here for convenience
    int normalizationId;

    public IlluminaBPMLocusEntry() {
        ilmnId = "";
        name = "";
        index = -1;
        ilmnStrand = "";
        snp = "";
        chrom = "";
        ploidy = "";
        species = "";
        mapInfo = -1;
        customerStrand = "";

        addressA = -1;
        addressB = -1;

        genomeBuild = "";
        source = "";
        sourceVersion = "";
        sourceStrand = "";

        sourceStrand = "";

        expClusters = -1;
        intensityOnly = -1;
        assayType = -1;

        fracA = 0.0f;
        fracC = 0.0f;
        fracT = 0.0f;
        fracG = 0.0f;

        refStrand = "";

        normalizationId = -1;
    }

}

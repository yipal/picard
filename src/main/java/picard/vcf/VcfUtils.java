package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import java.io.File;

/**
 * Created by farjoun on 4/1/17.
 */
public class VcfUtils {

    public static String UNCOMPRESSED_VCF_ENDING = IOUtil.VCF_FILE_EXTENSION;
    public static String COMPRESSED_VCF_ENDING = IOUtil.COMPRESSED_VCF_FILE_EXTENSION;
    public static String BCF_ENDING = IOUtil.BCF_FILE_EXTENSION;

    /**
     * Checks if the suffix is one of those that are allowed for the various
     * formats that contain variants (currently vcf and bcf)
     */
    static public boolean isVariantFile(final File file){
        final String name = file.getName();

        return name.endsWith(UNCOMPRESSED_VCF_ENDING) ||
                name.endsWith(COMPRESSED_VCF_ENDING) ||
                name.endsWith(BCF_ENDING);
    }
}

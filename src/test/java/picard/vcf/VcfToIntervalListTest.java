package picard.vcf;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by lichtens on 7/11/17.
 */
public class VcfToIntervalListTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("testdata/picard/vcf");
    public static final String smallM2VcfMore = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2_more_variants.vcf";

    public String getCommandLineProgramName() {
        return VcfToIntervalList.class.getSimpleName();
    }

    @Test
    public void testExcludingFiltered() throws IOException {
        final File outputFile = File.createTempFile("vcftointervallist_", ".interval_list");
        outputFile.deleteOnExit();
        final List<String> arguments = new ArrayList<>();
        arguments.add("I=" + smallM2VcfMore);
        arguments.add("O=" + outputFile.getAbsolutePath());
        runPicardCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<Interval> intervals = IntervalList.fromFile(outputFile).getIntervals();

        // 11 total, 10 defined unique intervals (two variants are adjacent), one is filtered in INFO and
        // one is filtered in FORMAT FT, but only INFO counts
        Assert.assertEquals(intervals.size(), 11 - 1 - 1);
    }

    @Test
    public void testIncludingFiltered() throws IOException {
        final File outputFile = File.createTempFile("vcftointervallist_incl_", ".interval_list");
        outputFile.deleteOnExit();
        final List<String> arguments = new ArrayList<>();
        arguments.add("I=" + smallM2VcfMore);
        arguments.add("O=" + outputFile.getAbsolutePath());
        arguments.add(VcfToIntervalList.INCLUDE_FILTERED_SHORT_NAME + "=true");

        runPicardCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<Interval> intervals = IntervalList.fromFile(outputFile).getIntervals();

        // 11 total, 10 defined unique intervals (two variants are adjacent)
        Assert.assertEquals(intervals.size(), 11 - 1 );
    }
}

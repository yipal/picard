package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_CIGAR;
import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_START;

/**
 * Tests related to code in AbstractAlignmentMerger
 */
public class AbstractAlignmentMergerTest {
    final Log log = Log.getInstance(AbstractAlignmentMergerTest.class);

    @Test public void tesOverlappedReadClippingWithNonOverlappedReads() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(110);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 200, false, false, "110M", "110M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "110M");
        Assert.assertEquals(r2.getAlignmentStart(), 200);
        Assert.assertEquals(r2.getCigarString(), "110M");
    }

    @Test public void testBasicOverlappedReadClipping() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(110);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 90, false, false, "110M", "110M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "100M10S");
        Assert.assertEquals(r2.getAlignmentStart(), 100);
        Assert.assertEquals(r2.getCigarString(), "10S100M");
    }

    @Test public void testOverlappedReadClippingWithExistingSoftClipping() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(120);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 95, false, false, "110M10S", "15S105M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "100M20S");
        Assert.assertEquals(r2.getAlignmentStart(), 100);
        Assert.assertEquals(r2.getCigarString(), "20S100M");
    }

    @Test public void testOverlappedReadClippingWithExistingSoftClippingAndHardClipping() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(120);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 95, false, false, "110M10S5H", "5H15S105M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "100M20S"); // Should ideally be 100M20S5H
        Assert.assertEquals(r2.getAlignmentStart(), 100);
        Assert.assertEquals(r2.getCigarString(), "20S100M"); // Should ideally be 5H20S100M
    }


    @DataProvider(name="transferAlignmentData")
    Object [][] transferAlignmentData(){
        List<Object[]> tests = new ArrayList<>();

        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(120);
        final SAMRecord alignedRec = set.addFrag( "q1", 0, 100, false, false, "5H15S105M", null,30);

        final SAMRecord unalignedRec = set.addFrag("q1", 0, 100, true, true, NO_ALIGNMENT_CIGAR, null, 30);
        unalignedRec.setAttribute("MC", "MATE_CIGAR");
        unalignedRec.setMateReferenceIndex(-1);

        final SAMRecord expectedRec = set.addFrag( "q1", 0, 100, false, false, "5H15S105M", null, 30);
        expectedRec.setAttribute("MC", "MATE_CIGAR");

        tests.add(new Object[] {expectedRec, alignedRec, unalignedRec, false, Collections.singletonList("MC"), AbstractAlignmentMerger.UnmappingReadStrategy.COPY_TO_TAG});

        return tests.toArray(new Object[][]{});

    }


    @Test(dataProvider = "transferAlignmentData")
    public void testTransferAlignmentInfoToFragment(SAMRecord expected, SAMRecord aligned, SAMRecord unaligned, final boolean isContaminant, List<String> attributesToRetain, AbstractAlignmentMerger.UnmappingReadStrategy unmappingReadStrategy){

        AbstractAlignmentMerger.actuallyTransferAlignmentInfoToFragment(unaligned, aligned, isContaminant, false, attributesToRetain, Collections.emptyList(), Collections.emptySet(), Collections.emptySet(), false, 0, 0, unmappingReadStrategy, log);

        Assert.assertEquals(unaligned, expected);

    }
}

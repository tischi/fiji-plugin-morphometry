package tests;

import de.embl.cba.morphometry.commands.SpindleMorphometryCommand;
import de.embl.cba.morphometry.spindle.SpindleMeasurements;
import loci.common.DebugTools;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class TestSpindleMorphometry
{
	@Test
	public < R extends RealType< R > > void testSmallSpindle( )
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthSmall.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		assertEquals( 7.5, spindleWidth, 1.0 );
	}

	@Test
	public < R extends RealType< R > > void testLargeSpindle()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthLarge.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		assertEquals( 8.79, spindleWidth, 1.0 );
	}


	@Test
	public < R extends RealType< R > > void testSpindleWithBrightOtherDNA()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/BrightOtherDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 2; // This is the other way around than in the others!
		command.dnaChannelIndexOneBased = 1;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		assertEquals( spindleWidth, 9.25, 1.0 );
	}


	@Test
	public < R extends RealType< R > > void testDimDNA()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/DimDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.minimalDynamicRange = 20;
		command.showIntermediateResults = false;
		command.saveResults = true;
		command.outputDirectory = new File("/Users/tischer/" +
				"Documents/fiji-plugin-morphometry/src" +
				"/test/resources/test-data/spindle/output" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final String comment = ( String) measurements.get( 0 ).get(
				SpindleMeasurements.COMMENT );

		assertEquals( comment, SpindleMeasurements.ANALYSIS_INTERRUPTED_LOW_DYNAMIC_DNA );
	}


	public static void main( String[] args )
	{
//		new TestSpindleMorphometry().testDimDNA();
//		new TestSpindleMorphometry().testSpindleWithBrightOtherDNA();
		new TestSpindleMorphometry().testSmallSpindle();
//		new TestSpindleMorphometry().testLargeSpindle();
	}

}
